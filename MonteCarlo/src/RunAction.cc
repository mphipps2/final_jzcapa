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
#include "FiberSD.hh"

#include "G4SDManager.hh"
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
  //Create the analysis manager
  m_analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << m_analysisManager->GetType() << " type Analysis Manager" << G4endl;
  char name[20];


  // Get number of ZDCs and RPDs from DetectorConstruction
  const DetectorConstruction* constDetectorConstruction
    = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  DetectorConstruction* detectorConstruction
    = const_cast<DetectorConstruction*>(constDetectorConstruction);

  int    nZDCs   = detectorConstruction->GetnZDCs();
  int    nRPDs   = detectorConstruction->GetnRPDs();
  CLUSTER = detectorConstruction->GetClusterFlag();

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
  m_ZDCdblVec = eventAction->GetZDCdoubleVectors( );
  m_ZDCintVec = eventAction->GetZDCintVectors( );
  m_RPDdblVec = eventAction->GetRPDdoubleVectors( );
  m_RPDintVec = eventAction->GetRPDintVectors( );

  FiberSD* sd;
  //Create ZDC trees and branches
  for(int zdcNo = 0; zdcNo < nZDCs; zdcNo++){
    sprintf(name,"ZDC%d_SD",zdcNo);
    //Find out from the SD if optical is on for this detector
    sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
      if( sd->OpticalIsOn() ) MakeZDCOpticalTree( zdcNo, zdcNo );
    else MakeZDCTree( zdcNo, zdcNo );
    m_analysisManager->FinishNtuple();
  }//end ZDC loop

  //Create RPD trees and branches
  for(int rpdNo = 0; rpdNo < nRPDs; rpdNo++){
    sprintf(name,"RPD%d_SD",rpdNo);
    int nTuple = nZDCs + rpdNo;
    //Find out from the SD if optical is on for this detector
    sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
      if( sd->OpticalIsOn() ) MakeRPDOpticalTree( nTuple, rpdNo );
    else MakeRPDTree( nTuple, rpdNo );
    m_analysisManager->FinishNtuple();
  }//end RPD loop




}//end constructor

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::MakeZDCTree( G4int nTupleNo, G4int zdcNo ){
  char name[20];
  sprintf(name,"ZDC%dtree",zdcNo);
  m_analysisManager->CreateNtuple( name, "ZDC data");
  if(!CLUSTER){
    //Make double branches
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastSteptest"                              );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_ZDCdblVec->at(zdcNo).at(3) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_ZDCdblVec->at(zdcNo).at(4) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_ZDCdblVec->at(zdcNo).at(5) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "energy",      m_ZDCdblVec->at(zdcNo).at(6) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "velocity",    m_ZDCdblVec->at(zdcNo).at(7) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "beta",        m_ZDCdblVec->at(zdcNo).at(8) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "eDep",        m_ZDCdblVec->at(zdcNo).at(9) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "charge",      m_ZDCdblVec->at(zdcNo).at(10) );

    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_ZDCintVec->at(zdcNo).at(0) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "radNo",       m_ZDCintVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(2) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_ZDCintVec->at(zdcNo).at(3) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_ZDCintVec->at(zdcNo).at(4) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_ZDCintVec->at(zdcNo).at(5) );

  } else { //Fewer vector branches to save storage space
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );

    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_ZDCintVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "radNo",       m_ZDCintVec->at(zdcNo).at(0) );
  }//end if !CLUSTER
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::MakeZDCOpticalTree( G4int nTupleNo, G4int zdcNo ){
  char name[20];
  sprintf(name,"ZDC%dtree",zdcNo);
  m_analysisManager->CreateNtuple( name, "ZDC data");

  //Make double branches
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_ZDCdblVec->at(zdcNo).at(6) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_ZDCdblVec->at(zdcNo).at(7) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_ZDCdblVec->at(zdcNo).at(8) );

  //Make int branches
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(zdcNo).at(0) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::MakeRPDTree       ( G4int nTupleNo, G4int rpdNo ){
  char name[20];
  sprintf(name,"RPD%dtree",rpdNo);
  m_analysisManager->CreateNtuple( name, "RPD data");
  if(!CLUSTER){
    //Make double branches
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_RPDdblVec->at(rpdNo).at(3) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_RPDdblVec->at(rpdNo).at(4) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_RPDdblVec->at(rpdNo).at(5) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "energy",      m_RPDdblVec->at(rpdNo).at(6) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "velocity",    m_RPDdblVec->at(rpdNo).at(7) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "beta",        m_RPDdblVec->at(rpdNo).at(8) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "eDep",        m_RPDdblVec->at(rpdNo).at(9) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "charge",      m_RPDdblVec->at(rpdNo).at(10));

    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_RPDintVec->at(rpdNo).at(0) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "radNo",       m_RPDintVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(2) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_RPDintVec->at(rpdNo).at(3) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_RPDintVec->at(rpdNo).at(4) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_RPDintVec->at(rpdNo).at(5) );

  } else { //Fewer vector branches to save storage space
    //Make double branches
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );

    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_RPDintVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );
  }//end if !CLUSTER
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::MakeRPDOpticalTree( G4int nTupleNo, G4int rpdNo ){
  char name[20];
  sprintf(name,"RPD%dtree",rpdNo);
  m_analysisManager->CreateNtuple( name, "RPD data");
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );

  //Make int branches
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(__attribute__((unused)) const G4Run* run)
{

  // Open an output file
  //if(analysisManager->GetFileName() == "")
  m_analysisManager->OpenFile("output.root");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{

	G4int nofEvents = run->GetNumberOfEvent();
	if (nofEvents == 0) return;

  m_analysisManager->Write();
  m_analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
