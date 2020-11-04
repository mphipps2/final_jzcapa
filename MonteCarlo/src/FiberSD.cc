//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Author: Michael Phipps

#include "FiberSD.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FiberSD::FiberSD(G4String sdName, G4int modNum, G4bool optical)
  :G4VSensitiveDetector(sdName),
  m_modNum(modNum),
  OPTICAL(optical),
  REDUCED_TREE(false),
  ZDC(false),
  RPD(false),
  m_cherenkovVec(0) {
  collectionName.insert(sdName);
  HCID = -1;

  if( sdName.contains("R") ) RPD = true;
  if( sdName.contains("Z") ) ZDC = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FiberSD::~FiberSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::HistInitialize(){
  std::string name = GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::Initialize(G4HCofThisEvent* HCE){

  fiberCollection = new FiberHitsCollection(SensitiveDetectorName,
					      m_modNum);

  std::string name = collectionName[0];

  if(REDUCED_TREE){
    // Grab the vector from AnalysisManager if reduced tree mode is active
    AnalysisManager* analysisManager = AnalysisManager::getInstance();
    m_cherenkovVec = analysisManager->GetFiberVector(ZDC,RPD,m_modNum);
    m_timeVec      = analysisManager->GetTimeVector(ZDC,RPD,m_modNum);

    // Clear it and resize in case this isn't the first event
    m_cherenkovVec->clear();
    m_cherenkovVec->resize( m_nFibers, 0 );
    m_timeVec->clear();
    m_timeVec->resize( 128, 0 );
  }

  m_nCherenkovs = 0;
  m_nHits = 0;

  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID( name );}

  HCE->AddHitsCollection( HCID, fiberCollection );
  G4cout << " HCID " << HCID << " name " << name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FiberSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

  //Get the number of Cherenkov photons created in this step
  int capturedPhotons = 0;
  const std::vector<const G4Track*>* secVec = aStep->GetSecondaryInCurrentStep();
  for(uint i = 0; i < secVec->size(); i++){
    if( secVec->at(i)->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
      capturedPhotons++;
    }//end if photon
  }//end secondary track loop
  m_nCherenkovs += capturedPhotons; // Record the total in case OPTICAL is true

  //Don't record hits that didn't produce cherenkov photons
  // if(capturedPhotons == 0) return true;

  G4int rodNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);

  G4ThreeVector pos = aStep->GetTrack()->GetPosition();
  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();

  // If OPTICAL is true, determine if the photon has reached the top of the topOfVolume
  // and add the hit to the collection if it has
  if(OPTICAL){
    if( aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() && pos.y() >= m_topOfVolume - 0.1*mm){
      if(REDUCED_TREE){
        m_cherenkovVec->at(rodNum)++;

        // This vector ranges from 0 to 64ns in 0.5ns bins. The 128th bin is overflow
        double time = aStep->GetTrack()->GetGlobalTime();
        int bin = (time <= 64*ns) ? time/0.5 : 127;
        m_timeVec->at(bin)++;

        m_nHits++;
        return true;
      }
      FiberHit* newHit = new FiberHit();
      newHit->setPos      ( pos );
      newHit->setOrigin   ( aStep->GetTrack()->GetVertexPosition() );
      newHit->setMomentum ( aStep->GetPreStepPoint()->GetMomentum() );
      newHit->setEnergy   ( aStep->GetPreStepPoint()->GetTotalEnergy() );
      newHit->setTime     ( aStep->GetTrack()->GetGlobalTime() );
      newHit->setRodNb    ( rodNum );

      fiberCollection->insert ( newHit );
      m_nHits++;

      aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
      return true;
    }
  }else{ // Otherwise record all hits
    if(REDUCED_TREE){
      m_cherenkovVec->at(rodNum) += capturedPhotons;

      // This vector ranges from 0 to 64ns in 0.5ns bins. The 128th bin is overflow
      double time = aStep->GetTrack()->GetGlobalTime();
      int bin = (time <= 64*ns) ? time/0.5 : 127;
      m_timeVec->at(bin) += capturedPhotons;

      m_nHits++;
      return true;
    }
    FiberHit* newHit = new FiberHit();
    newHit->setCharge      ( aStep->GetPreStepPoint()->GetCharge() );
    newHit->setTrackID     ( aStep->GetTrack()->GetTrackID() );
    newHit->setModNb       ( m_modNum );
    newHit->setRodNb       ( rodNum );
    newHit->setEdep        ( aStep->GetTotalEnergyDeposit() );
    newHit->setOrigin      ( aStep->GetTrack()->GetVertexPosition() );
    newHit->setPos         ( pos );
    newHit->setParticle    ( particle );
    newHit->setEnergy      ( aStep->GetPreStepPoint()->GetTotalEnergy() );
    newHit->setMomentum    ( aStep->GetPreStepPoint()->GetMomentum() );
    newHit->setTime        ( aStep->GetTrack()->GetGlobalTime() );
    newHit->setNCherenkovs ( capturedPhotons );

    fiberCollection->insert ( newHit );
    m_nHits++;
    return true;
  }

  return false; //Something failed if we got here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::EndOfEvent(G4HCofThisEvent*)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
