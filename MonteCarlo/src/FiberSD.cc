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

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

#include <string>
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FiberSD::FiberSD(G4String sdName, G4int modNum, G4bool optical)
  :G4VSensitiveDetector(sdName), m_modNum(modNum), OPTICAL(optical) {
  collectionName.insert(sdName);
  HCID = -1;

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

  m_nCherenkovs = 0;

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

  //      Figure out if this is necessary
  //
  G4int    totalRodNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
  G4String radNum_s = "7" +  aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() + "0";
  G4int    rodNum;
  G4int    radNum;

	rodNum = totalRodNum;
	radNum = (std::stoi (radNum_s)*100)+rodNum;
  // ^^^^ Figure out if this is necessary ^^^^

  G4ThreeVector pos = aStep->GetTrack()->GetPosition();
  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();

  // If OPTICAL is true, determine if the photon has reached the top of the topOfVolume
  // and add the hit to the collection if it has
  if(OPTICAL){
    if( aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() && pos.y() >= m_topOfVolume){
      FiberHit* newHit = new FiberHit();
      newHit->setPos      ( pos );
      newHit->setOrigin   ( aStep->GetTrack()->GetVertexPosition() );
      newHit->setMomentum ( aStep->GetPreStepPoint()->GetMomentum() );
      newHit->setEnergy   ( aStep->GetPreStepPoint()->GetTotalEnergy() );
      newHit->setRodNb    ( rodNum );

      fiberCollection->insert ( newHit );

      aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
      return true;
    }
  }else{ // Otherwise record all hits
    FiberHit* newHit = new FiberHit();
    newHit->setCharge      ( aStep->GetPreStepPoint()->GetCharge() );
    newHit->setTrackID     ( aStep->GetTrack()->GetTrackID() );
    newHit->setModNb       ( m_modNum );
    newHit->setRadNb       ( radNum );
    newHit->setRodNb       ( rodNum );
    newHit->setEdep        ( aStep->GetTotalEnergyDeposit() );
    newHit->setOrigin      ( aStep->GetTrack()->GetVertexPosition() );
    newHit->setPos         ( pos );
    newHit->setParticle    ( particle );
    newHit->setEnergy      ( aStep->GetPreStepPoint()->GetTotalEnergy() );
    newHit->setMomentum    ( aStep->GetPreStepPoint()->GetMomentum() );
    newHit->setNCherenkovs ( capturedPhotons );

    fiberCollection->insert ( newHit );
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::EndOfEvent(G4HCofThisEvent*)
{

  //  G4int NbHits = fiberCollection->entries();

  /*
  if(verboseLevel>0) {
      std::cout << " if verbose loop" << std::endl;
      std::cout << "\n-------->Hits Collection: in this event they are " << NbHits
		<< " hits in the calorimeter cells: " << std::endl;
      for (G4int i=0;i<NbHits;i++) {
	if (i %100 == 0) std::cout << " i " << i << std::endl;
	(*fiberCollection)[i]->Print();
      }
  }
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
