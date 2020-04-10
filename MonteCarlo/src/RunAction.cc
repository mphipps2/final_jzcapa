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
// $Id: RunAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4HCtable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction()
{
  // Get number of ZDCs and RPDs from DetectorConstruction
  const DetectorConstruction* constDetectorConstruction
    = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  DetectorConstruction* detectorConstruction
    = const_cast<DetectorConstruction*>(constDetectorConstruction);

  int    nZDCs   = detectorConstruction->GetnZDCs();
  int    nRPDs   = detectorConstruction->GetnRPDs();
  G4bool CLUSTER = detectorConstruction->GetClusterFlag();

  // Get event action
  const EventAction* constEventAction
    = static_cast<const EventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
  EventAction* eventAction
    = const_cast<EventAction*>(constEventAction);

  eventAction->SetClusterFlag(CLUSTER);
  eventAction->SetnZDCs(nZDCs);
  eventAction->SetnRPDs(nRPDs);

  // The sizes of these vectors will tell us how many of each detector we have
  // and what branches we want on the trees. We also need the addresses so we can
  // write the contents to the trees.
  std::vector< std::vector< std::vector<double>* > >*  RPDdblVec
    = eventAction->GetRPDdoubleVectors( );
  std::vector< std::vector< std::vector< int  >* > >*  RPDintVec
    = eventAction->GetRPDintVectors( );
  std::vector< std::vector< std::vector<double>* > >*  ZDCdblVec
    = eventAction->GetZDCdoubleVectors( );
  std::vector< std::vector< std::vector< int  >* > >*  ZDCintVec
    = eventAction->GetZDCintVectors( );


  //Create the analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;
  char name[20];
  //Create ZDC trees and branches
  for(int zdcNo = 0; zdcNo < nZDCs; zdcNo++){
    sprintf(name,"ZDC%dtree",zdcNo);
    analysisManager->CreateNtuple( name, "ZDC data");
    if(!CLUSTER){
      //Make double branches
      analysisManager->CreateNtupleDColumn( zdcNo, "lastStep"                                 );
      analysisManager->CreateNtupleDColumn( zdcNo, "lastSteptest"                             );
      analysisManager->CreateNtupleDColumn( zdcNo, "gunPosX"                                  );
      analysisManager->CreateNtupleDColumn( zdcNo, "gunPosY"                                  );
      analysisManager->CreateNtupleDColumn( zdcNo, "gunPosZ"                                  );
      analysisManager->CreateNtupleDColumn( zdcNo, "x",           *ZDCdblVec->at(zdcNo).at(0) );
      analysisManager->CreateNtupleDColumn( zdcNo, "y",           *ZDCdblVec->at(zdcNo).at(1) );
      analysisManager->CreateNtupleDColumn( zdcNo, "z",           *ZDCdblVec->at(zdcNo).at(2) );
      analysisManager->CreateNtupleDColumn( zdcNo, "Px",          *ZDCdblVec->at(zdcNo).at(3) );
      analysisManager->CreateNtupleDColumn( zdcNo, "Py",          *ZDCdblVec->at(zdcNo).at(4) );
      analysisManager->CreateNtupleDColumn( zdcNo, "Pz",          *ZDCdblVec->at(zdcNo).at(5) );
      analysisManager->CreateNtupleDColumn( zdcNo, "energy",      *ZDCdblVec->at(zdcNo).at(6) );
      analysisManager->CreateNtupleDColumn( zdcNo, "velocity",    *ZDCdblVec->at(zdcNo).at(7) );
      analysisManager->CreateNtupleDColumn( zdcNo, "beta",        *ZDCdblVec->at(zdcNo).at(8) );
      analysisManager->CreateNtupleDColumn( zdcNo, "eDep",        *ZDCdblVec->at(zdcNo).at(9) );

      //Make int branches
      analysisManager->CreateNtupleIColumn( zdcNo, "EventNo",     *ZDCintVec->at(zdcNo).at(0) );
      analysisManager->CreateNtupleIColumn( zdcNo, "modNo",       *ZDCintVec->at(zdcNo).at(1) );
      analysisManager->CreateNtupleIColumn( zdcNo, "radNo",       *ZDCintVec->at(zdcNo).at(2) );
      analysisManager->CreateNtupleIColumn( zdcNo, "rodNo",       *ZDCintVec->at(zdcNo).at(3) );
      analysisManager->CreateNtupleIColumn( zdcNo, "nCherenkovs", *ZDCintVec->at(zdcNo).at(4) );
      analysisManager->CreateNtupleIColumn( zdcNo, "trackID",     *ZDCintVec->at(zdcNo).at(5) );
      analysisManager->CreateNtupleIColumn( zdcNo, "pid",         *ZDCintVec->at(zdcNo).at(6) );
      analysisManager->CreateNtupleIColumn( zdcNo, "charge",      *ZDCintVec->at(zdcNo).at(7) );

    } else { // There's only two branches to save space on cluster jobs
      analysisManager->CreateNtupleDColumn( zdcNo, "lastStep"                                 );
      analysisManager->CreateNtupleDColumn( zdcNo, "gunPosX"                                  );
      analysisManager->CreateNtupleDColumn( zdcNo, "gunPosY"                                  );
      analysisManager->CreateNtupleDColumn( zdcNo, "gunPosZ"                                  );
      analysisManager->CreateNtupleIColumn( zdcNo, "radNo",       *ZDCintVec->at(zdcNo).at(0) );
      analysisManager->CreateNtupleIColumn( zdcNo, "nCherenkovs", *ZDCintVec->at(zdcNo).at(1) );
    }//end if !CLUSTER
  }//end ZDC loop

  //Create RPD trees and branches
  for(int rpdNo = 0; rpdNo < nRPDs; rpdNo++){
    sprintf(name,"RPD%dtree",rpdNo);
    analysisManager->CreateNtuple( name, "RPD data");
    int nTuple = nZDCs + rpdNo;
    if(!CLUSTER){
      //Make double branches
      analysisManager->CreateNtupleDColumn( nTuple, "lastStep"                                 );
      analysisManager->CreateNtupleDColumn( nTuple, "lastSteptest"                             );
      analysisManager->CreateNtupleDColumn( nTuple, "gunPosX"                                  );
      analysisManager->CreateNtupleDColumn( nTuple, "gunPosY"                                  );
      analysisManager->CreateNtupleDColumn( nTuple, "gunPosZ"                                  );
      analysisManager->CreateNtupleDColumn( nTuple, "x",           *RPDdblVec->at(rpdNo).at(0) );
      analysisManager->CreateNtupleDColumn( nTuple, "y",           *RPDdblVec->at(rpdNo).at(1) );
      analysisManager->CreateNtupleDColumn( nTuple, "z",           *RPDdblVec->at(rpdNo).at(2) );
      analysisManager->CreateNtupleDColumn( nTuple, "Px",          *RPDdblVec->at(rpdNo).at(3) );
      analysisManager->CreateNtupleDColumn( nTuple, "Py",          *RPDdblVec->at(rpdNo).at(4) );
      analysisManager->CreateNtupleDColumn( nTuple, "Pz",          *RPDdblVec->at(rpdNo).at(5) );
      analysisManager->CreateNtupleDColumn( nTuple, "energy",      *RPDdblVec->at(rpdNo).at(6) );
      analysisManager->CreateNtupleDColumn( nTuple, "velocity",    *RPDdblVec->at(rpdNo).at(7) );
      analysisManager->CreateNtupleDColumn( nTuple, "beta",        *RPDdblVec->at(rpdNo).at(8) );
      analysisManager->CreateNtupleDColumn( nTuple, "eDep",        *RPDdblVec->at(rpdNo).at(9) );

      //Make int branches
      analysisManager->CreateNtupleIColumn( nTuple, "EventNo",     *RPDintVec->at(rpdNo).at(0) );
      analysisManager->CreateNtupleIColumn( nTuple, "modNo",       *RPDintVec->at(rpdNo).at(1) );
      analysisManager->CreateNtupleIColumn( nTuple, "radNo",       *RPDintVec->at(rpdNo).at(2) );
      analysisManager->CreateNtupleIColumn( nTuple, "rodNo",       *RPDintVec->at(rpdNo).at(3) );
      analysisManager->CreateNtupleIColumn( nTuple, "nCherenkovs", *RPDintVec->at(rpdNo).at(4) );
      analysisManager->CreateNtupleIColumn( nTuple, "trackID",     *RPDintVec->at(rpdNo).at(5) );
      analysisManager->CreateNtupleIColumn( nTuple, "pid",         *RPDintVec->at(rpdNo).at(6) );
      analysisManager->CreateNtupleIColumn( nTuple, "charge",      *RPDintVec->at(rpdNo).at(7) );

    } else { // There's only two branches to save space on cluster jobs
      //Make double branches
      analysisManager->CreateNtupleDColumn( nTuple, "lastStep"                                 );
      analysisManager->CreateNtupleDColumn( nTuple, "gunPosX"                                  );
      analysisManager->CreateNtupleDColumn( nTuple, "gunPosY"                                  );
      analysisManager->CreateNtupleDColumn( nTuple, "gunPosZ"                                  );
      analysisManager->CreateNtupleDColumn( nTuple, "x",           *RPDdblVec->at(rpdNo).at(0) );
      analysisManager->CreateNtupleDColumn( nTuple, "y",           *RPDdblVec->at(rpdNo).at(1) );
      analysisManager->CreateNtupleDColumn( nTuple, "z",           *RPDdblVec->at(rpdNo).at(2) );

      //Make int branches
      analysisManager->CreateNtupleIColumn( nTuple, "rodNo",       *RPDintVec->at(rpdNo).at(0) );
      analysisManager->CreateNtupleIColumn( nTuple, "nCherenkovs", *RPDintVec->at(rpdNo).at(1) );
    }//end if !CLUSTER
  }//end RPD loop
}//end constructor

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(__attribute__((unused)) const G4Run* run)
{
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //if(analysisManager->GetFileName() == "") 
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{

	G4int nofEvents = run->GetNumberOfEvent();
	if (nofEvents == 0) return;

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
