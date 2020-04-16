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
// $Id: RunAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "FiberSD.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager* AnalysisManager::analysisManager = NULL;

AnalysisManager* AnalysisManager::getInstance(void)
{
    if (analysisManager == NULL) {
        analysisManager = new AnalysisManager();
    }
    return analysisManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager::AnalysisManager()
 : m_FactoryOn(false)
{
  // Get a pointer to DetectorConstruction
  // const DetectorConstruction* constDetectorConstruction
  //   = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // m_detectorConstruction
  //   = const_cast<DetectorConstruction*>(constDetectorConstruction);

  m_analysisManager = G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisManager::~AnalysisManager()
{
  delete m_ZDCdblVec;
  delete m_RPDdblVec;
  delete m_ZDCintVec;
  delete m_RPDintVec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in AnalysisManager.hh
  G4cout << "Using " << m_analysisManager->GetType() << " type Analysis Manager" << G4endl;
  m_analysisManager->SetVerboseLevel(1);
  m_analysisManager->SetNtupleMerging(true);

  // Open an output file
  //
  G4bool fileOpen;
  if(m_analysisManager->GetFileName() == ""){
    fileOpen = m_analysisManager->OpenFile("output");
  } else {
    fileOpen = m_analysisManager->OpenFile();
  }
  if (! fileOpen) {
    G4cerr << "\n---> AnalysisManager::Book(): cannot open "
           << m_analysisManager->GetFileName() << G4endl;
    return;
  }

  DetectorConstruction* detectorConstruction = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  int    nZDCs   = detectorConstruction->GetnZDCs();
  int    nRPDs   = detectorConstruction->GetnRPDs();
  CLUSTER = detectorConstruction->GetClusterFlag();

  m_ZDCdblVec = new std::vector< std::vector< std::vector<double> > >(nZDCs);
  m_ZDCintVec = new std::vector< std::vector< std::vector< int  > > >(nZDCs);
  m_RPDdblVec = new std::vector< std::vector< std::vector<double> > >(nRPDs);
  m_RPDintVec = new std::vector< std::vector< std::vector< int  > > >(nRPDs);

  //Create ZDC trees and branches
  G4cout << "Number of ZDCs " << nZDCs << G4endl;
  G4cout << "Number of RPDs " << nRPDs << G4endl;
  FiberSD* sd;
  char name[20];
  for(int zdcNo = 0; zdcNo < nZDCs; zdcNo++){
    sprintf(name,"ZDC%d_SD",zdcNo+1);
    //Find out from the SD if optical is on for this detector
    sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
      if( sd->OpticalIsOn() ) MakeZDCOpticalTree( zdcNo, zdcNo );
    else MakeZDCTree( zdcNo, zdcNo );
  }//end ZDC loop
  //Create RPD trees and branches
  for(int rpdNo = 0; rpdNo < nRPDs; rpdNo++){
    sprintf(name,"RPD%d_SD",rpdNo+1);
    int nTuple = nZDCs + rpdNo;
    //Find out from the SD if optical is on for this detector
    sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );
      if( sd->OpticalIsOn() ) MakeRPDOpticalTree( nTuple, rpdNo );
    else MakeRPDTree( nTuple, rpdNo );
  }//end RPD loop

  m_FactoryOn = true;

  G4cout << "\n----> Output file is open in "
         << m_analysisManager->GetFileName() << "."
         << m_analysisManager->GetFileType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::Save()
{
  if (! m_FactoryOn) return;

  m_analysisManager->Write();
  m_analysisManager->CloseFile();

  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;

  delete m_analysisManager;
  m_FactoryOn = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillNtuples(){

  // fill ntuples  //
  for(int i = 0; i < m_analysisManager->GetNofNtuples(); i++){
    G4cout << "Writing nTuple " << i << G4endl;
    m_analysisManager->AddNtupleRow(i);
  }

  // Clear ZDC vectors //
  for(uint i = 0; i < m_ZDCdblVec->size(); i++){
    std::cout << "Clearing ZDC" << i << std::endl;
      m_ZDCdblVec->at(i).clear();
      m_ZDCintVec->at(i).clear();
  }
  // Clear RPD vectors //
  for(uint i = 0; i < m_RPDdblVec->size(); i++){
    std::cout << "Clearing RPD" << i << std::endl;
    m_RPDdblVec->at(i).clear();
    m_RPDintVec->at(i).clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeZDCTree( G4int nTupleNo, G4int zdcNo ){
  char name[20];
  sprintf(name,"ZDC%dtree",zdcNo+1);
  m_analysisManager->CreateNtuple( name, "ZDC data");
  if(!CLUSTER){
    G4cout << "making non cluster tree" << G4endl;
    //Resize the vector for the number of branches storing double vectors
    m_ZDCdblVec->at(zdcNo).resize(11);
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

    //Resize the vector for the number of branches storing int vectors
    m_ZDCintVec->at(zdcNo).resize(6);
    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_ZDCintVec->at(zdcNo).at(0) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "radNo",       m_ZDCintVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(2) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_ZDCintVec->at(zdcNo).at(3) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_ZDCintVec->at(zdcNo).at(4) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_ZDCintVec->at(zdcNo).at(5) );

  } else { //Fewer vector branches to save storage space
    std::cout << "making cluster tree" << std::endl;
    //Resize the vector for the number of branches storing double vectors
    m_ZDCdblVec->at(zdcNo).resize(6);
    //Make double branches
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );

    //Resize the vector for the number of branches storing int vectors
    m_ZDCintVec->at(zdcNo).resize(1);
    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_ZDCintVec->at(zdcNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "radNo",       m_ZDCintVec->at(zdcNo).at(0) );
  }//end if !CLUSTER
  m_analysisManager->FinishNtuple( );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeZDCOpticalTree( G4int nTupleNo, G4int zdcNo ){
  char name[20];
  sprintf(name,"ZDC%dtree",zdcNo+1);
  m_analysisManager->CreateNtuple( name, "ZDC data");
  G4cout << "making optical tree" << G4endl;

  //Resize the vector for the number of branches storing double vectors
  m_ZDCdblVec->at(zdcNo).resize(6);
  //Make double branches
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_ZDCdblVec->at(zdcNo).at(3) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_ZDCdblVec->at(zdcNo).at(4) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_ZDCdblVec->at(zdcNo).at(5) );

  //Resize the vector for the number of branches storing int vectors
  m_ZDCintVec->at(zdcNo).resize(1);
  //Make int branches
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(0) );

  m_analysisManager->FinishNtuple( nTupleNo );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeRPDTree( G4int nTupleNo, G4int rpdNo ){
  char name[20];
  sprintf(name,"RPD%dtree",rpdNo+1);
  m_analysisManager->CreateNtuple( name, "RPD data");
  if(!CLUSTER){
    //Resize the vector for the number of branches storing double vectors
    m_RPDdblVec->at(rpdNo).resize(11);
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

    //Resize the vector for the number of branches storing int vectors
    m_RPDintVec->at(rpdNo).resize(6);
    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_RPDintVec->at(rpdNo).at(0) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "radNo",       m_RPDintVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(2) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_RPDintVec->at(rpdNo).at(3) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_RPDintVec->at(rpdNo).at(4) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_RPDintVec->at(rpdNo).at(5) );

  } else { //Fewer vector branches to save storage space
    //Resize the vector for the number of branches storing double vectors
    m_RPDdblVec->at(rpdNo).resize(3);
    //Make double branches
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );

    //Resize the vector for the number of branches storing int vectors
    m_RPDintVec->at(rpdNo).resize(2);
    //Make int branches
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_RPDintVec->at(rpdNo).at(1) );
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );
  }//end if !CLUSTER

  m_analysisManager->FinishNtuple( nTupleNo );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeRPDOpticalTree( G4int nTupleNo, G4int rpdNo ){
  char name[20];
  sprintf(name,"RPD%dtree",rpdNo+1);
  //Resize the vector for the number of branches storing double vectors
  m_RPDdblVec->at(rpdNo).resize(3);
  //Make double branches
  m_analysisManager->CreateNtuple( name, "RPD data");
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStep"                                  );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "lastStepTest"                              );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosX"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosY"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "gunPosZ"                                   );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
  m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );

  //Resize the vector for the number of branches storing int vectors
  //Make int branches
  m_RPDintVec->at(rpdNo).resize(1);
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "EventNo"                                   );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );
  m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );

  m_analysisManager->FinishNtuple( nTupleNo );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::CreateVectors(G4int nZDCs, G4int nRPDs){
  //Make a new vector of vectors of pointers to vectors of doubles
  m_ZDCdblVec = new std::vector< std::vector< std::vector<double> > >(nZDCs);
  m_ZDCintVec = new std::vector< std::vector< std::vector< int  > > >(nZDCs);
  m_RPDdblVec = new std::vector< std::vector< std::vector<double> > >(nRPDs);
  m_RPDintVec = new std::vector< std::vector< std::vector< int  > > >(nRPDs);

  // Resize to the largest number of branches we will fill from vectors
  for(G4int i = 0; i < nZDCs; i++){
    m_ZDCdblVec->at(i).resize(10);
    m_ZDCintVec->at(i).resize(8);
  }//end ZDC loop
  for(G4int i = 0; i < nRPDs; i++){
    m_RPDdblVec->at(i).resize(10);
    m_RPDintVec->at(i).resize(8);

  }//end RPD loop
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
