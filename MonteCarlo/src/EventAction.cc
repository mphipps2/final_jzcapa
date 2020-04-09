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
// Author: Chad Lantz
//
// \file RunAction.cc
// \brief Implementation of the RunAction class

#include "EventAction.hh"
#include "FiberHit.hh"
#include "Analysis.hh"

#include "G4UnitsTable.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PrimaryVertex.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction( )
: G4UserEventAction(){
  hitsCollID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::SetnZDCs(G4int nZDCs){
  //Make a new vector of vectors of pointers to vectors of doubles
  fZDCdblVec = new std::vector< std::vector< std::vector<double>* > >(nZDCs);
  fZDCintVec = new std::vector< std::vector< std::vector< int  >* > >(nZDCs);

  //If cluster flag is false, we have more vectors to store
  G4int nDblVecs = (!CLUSTER) ? 11 : 0;
  G4int nIntVecs = (!CLUSTER) ?  7 : 2;

  for(G4int i = 0; i < nZDCs; i++){

    fZDCdblVec->at(i).resize(nDblVecs);
    for(G4int j = 0; j < nDblVecs; j++){
      fZDCdblVec->at(i).at(j) = new std::vector<double>;
    }//end double vector loop

    fZDCintVec->at(i).resize(nIntVecs);
    for(G4int j = 0; j < nIntVecs; j++){
      fZDCintVec->at(i).at(j) = new std::vector< int  >;
    }//end int vector loop
  }//end ZDC loop
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void EventAction::SetnRPDs(G4int nRPDs){
  //Make a new vector of vectors of pointers to vectors of doubles
  fRPDdblVec = new std::vector< std::vector< std::vector<double>* > >(nRPDs);
  fRPDintVec = new std::vector< std::vector< std::vector< int  >* > >(nRPDs);

  //If cluster flag is false, we have more vectors to store
  G4int nDblVecs = (!CLUSTER) ? 11 : 3;
  G4int nIntVecs = (!CLUSTER) ?  7 : 2;

  for(G4int i = 0; i < nRPDs; i++){

    fRPDdblVec->at(i).resize(nDblVecs);
    for(G4int j = 0; j < nDblVecs; j++){
      fRPDdblVec->at(i).at(j) = new std::vector<double>;
    }//end double vector loop

    fRPDintVec->at(i).resize(nIntVecs);
    for(G4int j = 0; j < nIntVecs; j++){
      fRPDintVec->at(i).at(j) = new std::vector< int  >;
    }//end int vector loop
  }//end RPD loop
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fEventNo = evt->GetEventID();
  hitsCollID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt){
  G4PrimaryVertex* pVert = evt->GetPrimaryVertex();
  gunPosX = pVert->GetX0();
  gunPosY = pVert->GetY0();
  gunPosZ = pVert->GetZ0();

  // Last step in volume?
  G4double LastStepInVolume = pVert->GetPosition().z();

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  FiberHitsCollection* HC = 0;
  G4int nCollections =  HCE->GetNumberOfCollections();
  int IDholder = 0;
  if(HCE) {
    while (hitsCollID < nCollections) {
      HC = (FiberHitsCollection*)(HCE->GetHC(hitsCollID));
      G4String name = HC->GetSDname();

      G4int prevTrackId = 0;
      G4int prevRadiatorNo = 0;
      G4int TotalnCherenkovs = 0;
      G4double TotaleDep = 0.0;

      int n_hit = HC->entries();
      std::cout << printf("%s hits collection, nHits = %d", name.c_str(),n_hit) << std::endl;
      for ( G4int i = 0 ; i < n_hit; i++){

        G4ThreeVector position      = (*HC)[i]->getPos();
        G4ThreeVector momentum      = (*HC)[i]->getMomentum();
        G4double      energy        = (*HC)[i]->getEnergy();
        G4double      velocity      = (*HC)[i]->getVelocity();
        G4double      beta          = (*HC)[i]->getBeta();
        G4double      eDep          = (*HC)[i]->getEdep();
        G4double      charge        = (*HC)[i]->getCharge();

        G4int         radiatorNo    = (*HC)[i]->getRadNb();
        G4int         rodNo         = (*HC)[i]->getRodNb();
        G4int         modNb         = (*HC)[i]->getModNb();
        G4int         trackID       = (*HC)[i]->getTrackID();
        G4int         pid           = (*HC)[i]->getParticle()->GetPDGEncoding();
        G4int         nCherenkovs   = (*HC)[i]->getNCherenkovs(); // This is the number of cherenkovs in a single step within the SD

        //Sum energy from all steps a particle takes in a single scoring volume
        if (trackID == prevTrackId || radiatorNo == prevRadiatorNo || eDep != 0) {
          TotalnCherenkovs += nCherenkovs;
          TotaleDep += eDep;

          prevTrackId = trackID;
          prevRadiatorNo = radiatorNo;
          continue;
        }// end summation

        if(name.compare(0,3,"RPD")){ //RPD hits
          int rpdNo = atoi( name.substr(3,1).c_str() );
          if(!CLUSTER){
            //doubles
            fRPDdblVec->at(rpdNo).at(0)-> push_back( position.x() );
            fRPDdblVec->at(rpdNo).at(1)-> push_back( position.y() );
            fRPDdblVec->at(rpdNo).at(2)-> push_back( position.z() );
            fRPDdblVec->at(rpdNo).at(3)-> push_back( momentum.x() );
            fRPDdblVec->at(rpdNo).at(4)-> push_back( momentum.y() );
            fRPDdblVec->at(rpdNo).at(5)-> push_back( momentum.z() );
            fRPDdblVec->at(rpdNo).at(6)-> push_back( energy       );
            fRPDdblVec->at(rpdNo).at(7)-> push_back( velocity     );
            fRPDdblVec->at(rpdNo).at(8)-> push_back( beta         );
            fRPDdblVec->at(rpdNo).at(9)-> push_back( TotaleDep    );
            fRPDdblVec->at(rpdNo).at(10)->push_back( charge       );

            //ints
            fRPDintVec->at(rpdNo).at(0)->push_back( fEventNo     );
            fRPDintVec->at(rpdNo).at(1)->push_back( modNb        );
            fRPDintVec->at(rpdNo).at(2)->push_back( radiatorNo   );
            fRPDintVec->at(rpdNo).at(3)->push_back( rodNo        );
            fRPDintVec->at(rpdNo).at(4)->push_back( TotalnCherenkovs  );
            fRPDintVec->at(rpdNo).at(5)->push_back( trackID      );
            fRPDintVec->at(rpdNo).at(6)->push_back( pid          );

          } else{
            //doubles
            fRPDdblVec->at(rpdNo).at(0)->push_back( position.x() );
            fRPDdblVec->at(rpdNo).at(1)->push_back( position.y() );
            fRPDdblVec->at(rpdNo).at(2)->push_back( position.z() );

            //ints
            fRPDintVec->at(rpdNo).at(0)->push_back( rodNo );
            fRPDintVec->at(rpdNo).at(1)->push_back( TotalnCherenkovs );
          }// end if !CLUSTER
        // end if RPD
        } else{
          if( name.compare(0,3,"ZDC") ){//ZDC hitsCollID, check to be sure/symmetric
            int zdcNo = atoi( name.substr(3,1).c_str() );
            if(!CLUSTER){
              //doubles
              fZDCdblVec->at(zdcNo).at(0)-> push_back( position.x() );
              fZDCdblVec->at(zdcNo).at(1)-> push_back( position.y() );
              fZDCdblVec->at(zdcNo).at(2)-> push_back( position.z() );
              fZDCdblVec->at(zdcNo).at(3)-> push_back( momentum.x() );
              fZDCdblVec->at(zdcNo).at(4)-> push_back( momentum.y() );
              fZDCdblVec->at(zdcNo).at(5)-> push_back( momentum.z() );
              fZDCdblVec->at(zdcNo).at(6)-> push_back( energy       );
              fZDCdblVec->at(zdcNo).at(7)-> push_back( velocity     );
              fZDCdblVec->at(zdcNo).at(8)-> push_back( beta         );
              fZDCdblVec->at(zdcNo).at(9)-> push_back( TotaleDep    );
              fZDCdblVec->at(zdcNo).at(10)->push_back( charge       );

              //ints
              fZDCintVec->at(zdcNo).at(0)->push_back( fEventNo     );
              fZDCintVec->at(zdcNo).at(1)->push_back( modNb        );
              fZDCintVec->at(zdcNo).at(2)->push_back( radiatorNo   );
              fZDCintVec->at(zdcNo).at(3)->push_back( rodNo        );
              fZDCintVec->at(zdcNo).at(4)->push_back( TotalnCherenkovs  );
              fZDCintVec->at(zdcNo).at(5)->push_back( trackID      );
              fZDCintVec->at(zdcNo).at(6)->push_back( pid          );

            } else{
              //ints
              fZDCintVec->at(zdcNo).at(0)->push_back( radiatorNo  );
              fZDCintVec->at(zdcNo).at(1)->push_back( TotalnCherenkovs );
            }// end if !CLUSTER
          }// end if ZDC
        }// end else (!RPD)
        TotalnCherenkovs = 0;
        TotaleDep = 0.0;
        hitsCollID++;
      }//end of hit loop

      // fill ntuples  //
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      for(int i = 0; i < analysisManager->GetNofNtuples(); i++){
        analysisManager->FillNtupleDColumn(i,1,LastStepInVolume);
        analysisManager->AddNtupleRow(i);
      }


      // Clear RPD vectors //
      for(uint i = 0; i < fRPDdblVec->size(); i++){
        //Clear double vectors
        for(uint j = 0; j < fRPDdblVec->at(0).size(); j++){
          fRPDdblVec->at(i).at(j)->clear();
        }
        //Clear int vectors
        for(uint j = 0; j < fRPDintVec->at(0).size(); j++){
          fRPDintVec->at(i).at(j)->clear();
        }
      }
      // Clear ZDC vectors //
      for(uint i = 0; i < fZDCdblVec->size(); i++){
        //Clear double vectors
        for(uint j = 0; j < fZDCdblVec->at(0).size(); j++){
          fZDCdblVec->at(i).at(j)->clear();
        }
        //Clear int vectors
        for(uint j = 0; j < fZDCintVec->at(0).size(); j++){
          fZDCintVec->at(i).at(j)->clear();
        }
      }


    }// end while < nCollections
  }// end if HCE
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
