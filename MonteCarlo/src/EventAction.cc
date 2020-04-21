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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
/// \author Chad Lantz
/// \date 16 April 2020

#include "EventAction.hh"
#include "FiberHit.hh"
#include "FiberSD.hh"
#include "AnalysisManager.hh"

#include "G4UnitsTable.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction( )
: G4UserEventAction(){
  hitsCollID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fEventNo = evt->GetEventID();
  m_analysisManager = AnalysisManager::getInstance();
  hitsCollID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt){
  G4PrimaryVertex* pVert = evt->GetPrimaryVertex();
  m_RPDdblVec = m_analysisManager->GetRPDdoubleVectors( );
  m_RPDintVec = m_analysisManager->GetRPDintVectors   ( );
  m_ZDCdblVec = m_analysisManager->GetZDCdoubleVectors( );
  m_ZDCintVec = m_analysisManager->GetZDCintVectors   ( );
  CLUSTER = m_analysisManager->GetClusterFlag();

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  FiberHitsCollection* HC = 0;
  G4int nCollections =  HCE->GetNumberOfCollections();
  if(HCE) {
    while (hitsCollID < nCollections) {
      HC = (FiberHitsCollection*)(HCE->GetHC(hitsCollID));
      G4String name = HC->GetSDname();

      FiberSD* sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
      if( sd->OpticalIsOn() ) ProcessOpticalHitCollection( HC );
      else                    ProcessHitCollection( HC );
      hitsCollID++;

    }// end while < nCollections

    //Use the base class to fill individual data points
    //There's probably a better way to do this, but I just want it
    //to work right now
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    // fill ntuples  //
    for(int i = 0; i < analysisManager->GetNofNtuples(); i++){
      analysisManager->FillNtupleDColumn(i,1, pVert->GetX0() );
      analysisManager->FillNtupleDColumn(i,2, pVert->GetY0() );
      analysisManager->FillNtupleDColumn(i,3, pVert->GetZ0() );

      analysisManager->FillNtupleIColumn(i,4, fEventNo );
    }

    //Use our custom class to finish the job
    m_analysisManager->FillNtuples();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ProcessHitCollection( FiberHitsCollection* HC ){
  G4int prevTrackId = 0;
  G4int prevRadiatorNo = 0;
  G4int nCherenkovsSum = 0;
  G4double eDepSum = 0.0;
  G4String name = HC->GetSDname();


  int n_hit = HC->entries();
  G4cout << name << " nHits = " << n_hit << G4endl;
  for ( G4int i = 0 ; i < n_hit; i++){

    G4double      eDep          = (*HC)[i]->getEdep();

    G4int         radiatorNo    = (*HC)[i]->getRadNb();
    G4int         nCherenkovs   = (*HC)[i]->getNCherenkovs(); // This is the number of cherenkovs in a single step within the SD
    G4int         trackID       = (*HC)[i]->getTrackID();

    //Sum energy from all steps a particle takes in a single scoring volume
    if (trackID == prevTrackId || radiatorNo == prevRadiatorNo || eDep != 0) {
      nCherenkovsSum += nCherenkovs;
      eDepSum += eDep;

      prevTrackId = trackID;
      prevRadiatorNo = radiatorNo;
      continue;
    }// end summation


    G4ThreeVector position      = (*HC)[i]->getPos();
    G4ThreeVector origin        = (*HC)[i]->getOrigin();
    G4ThreeVector momentum      = (*HC)[i]->getMomentum();
    G4double      energy        = (*HC)[i]->getEnergy();
    G4double      velocity      = (*HC)[i]->getVelocity();
    G4double      beta          = (*HC)[i]->getBeta();
    G4double      charge        = (*HC)[i]->getCharge();

    G4int         rodNo         = (*HC)[i]->getRodNb();
    G4int         modNb         = (*HC)[i]->getModNb();
    G4int         pid           = (*HC)[i]->getParticle()->GetPDGEncoding();

    if(name.compare(0,3,"RPD") == 0){ //RPD hits
      int rpdNo = atoi( name.substr(3,1).c_str() );
      if(!CLUSTER){
        //doubles
        m_RPDdblVec->at(rpdNo-1).at(0). push_back( position.x() );
        m_RPDdblVec->at(rpdNo-1).at(1). push_back( position.y() );
        m_RPDdblVec->at(rpdNo-1).at(2). push_back( position.z() );
        m_RPDdblVec->at(rpdNo-1).at(3). push_back( momentum.x() );
        m_RPDdblVec->at(rpdNo-1).at(4). push_back( momentum.y() );
        m_RPDdblVec->at(rpdNo-1).at(5). push_back( momentum.z() );
        m_RPDdblVec->at(rpdNo-1).at(6). push_back( energy       );
        m_RPDdblVec->at(rpdNo-1).at(7). push_back( velocity     );
        m_RPDdblVec->at(rpdNo-1).at(8). push_back( beta         );
        m_RPDdblVec->at(rpdNo-1).at(9). push_back( eDepSum      );
        m_RPDdblVec->at(rpdNo-1).at(10).push_back( charge       );

        //ints
        m_RPDintVec->at(rpdNo-1).at(0).push_back( modNb          );
        m_RPDintVec->at(rpdNo-1).at(1).push_back( radiatorNo     );
        m_RPDintVec->at(rpdNo-1).at(2).push_back( rodNo          );
        m_RPDintVec->at(rpdNo-1).at(3).push_back( nCherenkovsSum );
        m_RPDintVec->at(rpdNo-1).at(4).push_back( trackID        );
        m_RPDintVec->at(rpdNo-1).at(5).push_back( pid            );

      } else{
        //doubles
        m_RPDdblVec->at(rpdNo-1).at(0).push_back( position.x() );
        m_RPDdblVec->at(rpdNo-1).at(1).push_back( position.y() );
        m_RPDdblVec->at(rpdNo-1).at(2).push_back( position.z() );

        //ints
        m_RPDintVec->at(rpdNo-1).at(0).push_back( nCherenkovsSum );
        m_RPDintVec->at(rpdNo-1).at(1).push_back( rodNo );
      }// end if !CLUSTER
    // end if RPD
    } else{
      if( name.compare(0,3,"ZDC") == 0 ){//ZDC hitsCollID, check to be sure/symmetric
        int zdcNo = atoi( name.substr(3,1).c_str() );
        if(!CLUSTER){
          //doubles
          m_ZDCdblVec->at(zdcNo-1).at(0). push_back( position.x() );
          m_ZDCdblVec->at(zdcNo-1).at(1). push_back( position.y() );
          m_ZDCdblVec->at(zdcNo-1).at(2). push_back( position.z() );
          m_ZDCdblVec->at(zdcNo-1).at(3). push_back( momentum.x() );
          m_ZDCdblVec->at(zdcNo-1).at(4). push_back( momentum.y() );
          m_ZDCdblVec->at(zdcNo-1).at(5). push_back( momentum.z() );
          m_ZDCdblVec->at(zdcNo-1).at(6). push_back( energy       );
          m_ZDCdblVec->at(zdcNo-1).at(7). push_back( velocity     );
          m_ZDCdblVec->at(zdcNo-1).at(8). push_back( beta         );
          m_ZDCdblVec->at(zdcNo-1).at(9). push_back( eDepSum      );
          m_ZDCdblVec->at(zdcNo-1).at(10).push_back( charge       );

          //ints
          m_ZDCintVec->at(zdcNo-1).at(0).push_back( modNb          );
          m_ZDCintVec->at(zdcNo-1).at(1).push_back( radiatorNo     );
          m_ZDCintVec->at(zdcNo-1).at(2).push_back( rodNo          );
          m_ZDCintVec->at(zdcNo-1).at(3).push_back( nCherenkovsSum );
          m_ZDCintVec->at(zdcNo-1).at(4).push_back( trackID        );
          m_ZDCintVec->at(zdcNo-1).at(5).push_back( pid            );

        } else{
          //doubles
          m_ZDCdblVec->at(zdcNo-1).at(0). push_back( position.x() );
          m_ZDCdblVec->at(zdcNo-1).at(1). push_back( position.y() );
          m_ZDCdblVec->at(zdcNo-1).at(2). push_back( position.z() );

          //ints
          m_ZDCintVec->at(zdcNo-1).at(0).push_back( nCherenkovsSum );
          m_ZDCintVec->at(zdcNo-1).at(1).push_back( radiatorNo  );
        }// end if !CLUSTER
      }// end if ZDC
    }// end else (!RPD)
    nCherenkovsSum = 0;
    eDepSum = 0.0;
  }//end of hit loop
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ProcessOpticalHitCollection ( FiberHitsCollection* HC ){
  G4int prevTrackId = 0;
  G4int nCherenkovsSum = 0;
  G4String name = HC->GetSDname();
  FiberSD* sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
  G4double topOfVolume = sd->GetTopOfVolume();

  int n_hit = HC->entries();
  G4cout << name << " nOpticalHits = " << n_hit << G4endl;
  for ( G4int i = 0 ; i < n_hit; i++){


    G4int         nCherenkovs   = (*HC)[i]->getNCherenkovs();

    //Grab nCherenkovs and skip anything that isn't a photon
    if( (*HC)[i]->getParticle() != G4OpticalPhoton::OpticalPhotonDefinition() ){
      nCherenkovsSum += nCherenkovs;
      continue;
    }

    //Skip any photons that we've already recorded
    G4ThreeVector position = (*HC)[i]->getPos();
    G4int         trackID  = (*HC)[i]->getTrackID();
    if( position.y() < topOfVolume - 0.2*mm || trackID == prevTrackId ){
      continue;
    }
    trackID = prevTrackId;

    G4ThreeVector origin   = (*HC)[i]->getOrigin();
    G4ThreeVector momentum = (*HC)[i]->getMomentum();
    G4int         rodNo    = (*HC)[i]->getRodNb();

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if( name.compare(0,3,"ZDC") == 0 ){//ZDC hitsCollID, check to be sure/symmetric
      int zdcNo = atoi( name.substr(3,1).c_str() );
      m_ZDCdblVec->at(zdcNo-1).at(0). push_back( origin.x()   );
      m_ZDCdblVec->at(zdcNo-1).at(1). push_back( origin.y()   );
      m_ZDCdblVec->at(zdcNo-1).at(2). push_back( origin.z()   );
      m_ZDCdblVec->at(zdcNo-1).at(3). push_back( momentum.x() );
      m_ZDCdblVec->at(zdcNo-1).at(4). push_back( momentum.y() );
      m_ZDCdblVec->at(zdcNo-1).at(5). push_back( momentum.z() );

      analysisManager->FillNtupleIColumn( zdcNo - 1, 5, nCherenkovsSum );
      m_ZDCintVec->at(zdcNo-1).at(0).push_back( rodNo          );
    }//end fill ZDC vectors
    if( name.compare(0,3,"RPD") == 0 ){//RPD hitsCollID, check to be sure/symmetric
      int rpdNo = atoi( name.substr(3,1).c_str() );
      m_RPDdblVec->at(rpdNo-1).at(0). push_back( origin.x()   );
      m_RPDdblVec->at(rpdNo-1).at(1). push_back( origin.y()   );
      m_RPDdblVec->at(rpdNo-1).at(2). push_back( origin.z()   );
      m_RPDdblVec->at(rpdNo-1).at(3). push_back( momentum.x() );
      m_RPDdblVec->at(rpdNo-1).at(4). push_back( momentum.y() );
      m_RPDdblVec->at(rpdNo-1).at(5). push_back( momentum.z() );

      analysisManager->FillNtupleIColumn( rpdNo - 1 + m_ZDCdblVec->size() , 5, nCherenkovsSum );
      m_RPDintVec->at(rpdNo-1).at(0). push_back( rodNo );
    }//end fill RPD vectors
  }// end hit loop
}
