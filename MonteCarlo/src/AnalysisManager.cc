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
//
/// \file AnalysisManager.hh
/// \brief Definition of the AnalysisManager class
/// \author Chad Lantz
/// \date 16 April 2020

#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "SteppingAction.hh"
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
  m_analysisManager = G4AnalysisManager::Instance();
  m_lastStepVec = new std::vector< G4ThreeVector >;
  m_gunPos = new G4ThreeVector;

  m_Pi0Vert = new std::vector< G4ThreeVector >;
  m_Pi0Mom = new std::vector< G4ThreeVector >;
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

void AnalysisManager::Book( G4String fileName )
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
  if(fileName == ""){
    fileOpen = m_analysisManager->OpenFile("output");
  } else {
    fileOpen = m_analysisManager->OpenFile(fileName);
  }
  if (! fileOpen) {
    G4cerr << "\n---> AnalysisManager::Book(): cannot open "
           << m_analysisManager->GetFileName() << G4endl;
    return;
  }

  //Get information about the detector configuration from DetectorConstruction
  DetectorConstruction* detectorConstruction = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  int nZDCs = detectorConstruction->GetnZDCs();
  int nRPDs = detectorConstruction->GetnRPDs();
  OPTICAL   = detectorConstruction->GetOpticalFlag();
  PI0   = detectorConstruction->GetPI0Flag();

  //Make vectors for the detectors we have
  //Indecies are [module#][dataType][dataPoint]
  m_ZDCdblVec = new std::vector< std::vector< std::vector<double> > >(nZDCs);
  m_ZDCintVec = new std::vector< std::vector< std::vector< int  > > >(nZDCs);
  m_RPDdblVec = new std::vector< std::vector< std::vector<double> > >(nRPDs);
  m_RPDintVec = new std::vector< std::vector< std::vector< int  > > >(nRPDs);

  SteppingAction* steppingAction = (SteppingAction*)G4RunManager::GetRunManager()->GetUserSteppingAction();
  steppingAction->SetLastStepVec( m_lastStepVec, &m_lastStepPidVec );
  steppingAction->SetPi0Mom( m_Pi0Mom );
  steppingAction->SetPi0Vertex( m_Pi0Vert );


  //Create the event data tree
  MakeEventDataTree();

  //Create ZDC trees and branches
  FiberSD* sd;
  char name[20];
  for(int zdcNo = 0; zdcNo < nZDCs; zdcNo++){
    int nTuple = zdcNo + 1;

    //Find out from the SD if optical is on for this detector
    sprintf(name,"ZDC%d_SD",zdcNo+1);
    sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );

    MakeZDCTree( nTuple, zdcNo, sd->GetFiberVec(), sd->OpticalIsOn() );

  }//end ZDC loop

  //Create RPD trees
  for(int rpdNo = 0; rpdNo < nRPDs; rpdNo++){
    int nTuple = nZDCs + 1 + rpdNo;

    //Find out from the SD if optical is on for this detector
    sprintf(name,"RPD%d_SD",rpdNo+1);
    sd = (FiberSD*)G4SDManager::GetSDMpointer()->FindSensitiveDetector( name );

    MakeRPDTree( nTuple, rpdNo, sd->GetFiberVec(), sd->OpticalIsOn() );

  }//end RPD loop

  //Turn on the data factory I guess
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

  // Fill the last step vectors

  //G4cout << "Before lastStep loop" << G4endl;
  for(uint i = 0; i < m_lastStepVec->size(); i++){
    m_lastStepXVec.push_back( m_lastStepVec->at(i).x() );
    m_lastStepYVec.push_back( m_lastStepVec->at(i).y() );
    m_lastStepZVec.push_back( m_lastStepVec->at(i).z() );
  }
  for(uint i = 0; i < m_Pi0Mom->size(); i++){
    m_Pi0MomX.push_back( m_Pi0Mom->at(i).x() );
    m_Pi0MomY.push_back( m_Pi0Mom->at(i).y() );
    m_Pi0MomZ.push_back( m_Pi0Mom->at(i).z() );
  }

  for(uint i = 0; i < m_Pi0Vert->size(); i++){
    m_Pi0VertX.push_back( m_Pi0Vert->at(i).x() );
    m_Pi0VertY.push_back( m_Pi0Vert->at(i).y() );
    m_Pi0VertZ.push_back( m_Pi0Vert->at(i).z() );
  }



  m_analysisManager->FillNtupleIColumn( 0, 0, m_eventNo );

  m_analysisManager->FillNtupleDColumn( 0, 1, m_gunPos->x() );
  m_analysisManager->FillNtupleDColumn( 0, 2, m_gunPos->y() );
  m_analysisManager->FillNtupleDColumn( 0, 3, m_gunPos->z() );


  // fill ntuples  //
  for(int i = 0; i < m_analysisManager->GetNofNtuples(); i++){
    m_analysisManager->AddNtupleRow(i);
  }

  // Clear ZDC vectors //
  for(uint i = 0; i < m_ZDCdblVec->size(); i++){
    for(uint j = 0; j < m_ZDCdblVec->at(i).size(); j++){
      m_ZDCdblVec->at(i).at(j).clear();
    }
    for(uint j = 0; j < m_ZDCintVec->at(i).size(); j++){
      m_ZDCintVec->at(i).at(j).clear();
    }
  }
  // Clear RPD vectors //
  for(uint i = 0; i < m_RPDdblVec->size(); i++){
    for(uint j = 0; j < m_RPDdblVec->at(i).size(); j++){
      m_RPDdblVec->at(i).at(j).clear();
    }
    for(uint j = 0; j < m_RPDintVec->at(i).size(); j++){
      m_RPDintVec->at(i).at(j).clear();
    }
  }
  m_lastStepVec->clear();
  m_lastStepXVec.clear();
  m_lastStepYVec.clear();
  m_lastStepZVec.clear();
  m_lastStepPidVec.clear();

  m_Pi0Mom->clear();
  m_Pi0MomX.clear();
  m_Pi0MomY.clear();
  m_Pi0MomZ.clear();

  m_Pi0Vert->clear();
  m_Pi0VertX.clear();
  m_Pi0VertY.clear();
  m_Pi0VertZ.clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillZDCnCherenkovs( int zdcNo, int nCherenkovs ){
  // nTuple 0 is the EventData tree, and ZDC numbering starts at 1,
  // so the ZDC number is also the nTuple number
  m_analysisManager->FillNtupleIColumn( zdcNo, 0, nCherenkovs );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillRPDnCherenkovs( int rpdNo, int nCherenkovs ){
  // RPD trees are made last, so their nTuple number is rpdNo + nZDCs
  m_analysisManager->FillNtupleIColumn( rpdNo + m_ZDCintVec->size(), 0, nCherenkovs );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeEventDataTree( ){
  m_analysisManager->CreateNtuple( "EventData", "Event Data");

  //Integer
  m_analysisManager->CreateNtupleIColumn( 0, "EventNo"  );

  //Doubles
  m_analysisManager->CreateNtupleDColumn( 0, "gunPosX"  );
  m_analysisManager->CreateNtupleDColumn( 0, "gunPosY"  );
  m_analysisManager->CreateNtupleDColumn( 0, "gunPosZ"  );

  //std::vector< double >
  m_analysisManager->CreateNtupleDColumn( 0, "lastStepX",   m_lastStepXVec   );
  m_analysisManager->CreateNtupleDColumn( 0, "lastStepY",   m_lastStepYVec   );
  m_analysisManager->CreateNtupleDColumn( 0, "lastStepZ",   m_lastStepZVec   );


  //std::vector< int >
  m_analysisManager->CreateNtupleIColumn( 0, "lastStepPID", m_lastStepPidVec );

  if(PI0){
      //std::vector< double >
      m_analysisManager->CreateNtupleDColumn( 0, "pi0Px",   m_Pi0MomX );
      m_analysisManager->CreateNtupleDColumn( 0, "pi0Py",   m_Pi0MomY );
      m_analysisManager->CreateNtupleDColumn( 0, "pi0Pz",   m_Pi0MomZ );

      //std::vector< double >
      m_analysisManager->CreateNtupleDColumn( 0, "pi0Vx",   m_Pi0VertX );
      m_analysisManager->CreateNtupleDColumn( 0, "pi0Vy",   m_Pi0VertY );
      m_analysisManager->CreateNtupleDColumn( 0, "pi0Vz",   m_Pi0VertZ );
  }

  m_analysisManager->FinishNtuple( );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeEventGenTree( std::vector< std::vector<int>* > &intVec , std::vector< std::vector<double>* > &dblVec ){

  m_eventGenNtupleNo = m_analysisManager->GetNofNtuples();

  m_analysisManager->CreateNtuple( "eventGen", "Event Generator Data");

  //Integer
  m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "nPart"  );
  m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "nSpec"  );
  m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "model"  );

  //Doubles
  m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "impactParameter"  );

  ////std::vector< double >
  m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "px", *dblVec[0] );
  m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "py", *dblVec[1] );
  m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo, "pz", *dblVec[2] );
  m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo,  "E", *dblVec[3] );
  m_analysisManager->CreateNtupleDColumn( m_eventGenNtupleNo,  "m", *dblVec[4] );


  //vector< int > branches
  m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo,      "pdgid", *intVec[0] );
  m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "CRMCstatus", *intVec[1] );
  m_analysisManager->CreateNtupleIColumn( m_eventGenNtupleNo, "keptStatus", *intVec[2] );

  m_analysisManager->FinishNtuple( );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::FillEventGenTree( int nPart, int nSpec, int model, double impactParam ){

  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 0, nPart );
  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 1, nSpec );
  m_analysisManager->FillNtupleIColumn( m_eventGenNtupleNo, 2, model );

  //Doubles
  m_analysisManager->FillNtupleDColumn( m_eventGenNtupleNo, 3, impactParam );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeZDCTree( G4int nTupleNo, G4int zdcNo, std::vector< int >* nCherenkovVec, G4bool thisIsOptical ){
  char name[20];
  sprintf(name,"ZDC%dtree",zdcNo+1);
  m_analysisManager->CreateNtuple( name, "ZDC data");

  // If nCherenkovVec is defined for this SD the user has set the ReducedTree flag
  // Otherwise generate the full tree
  if( nCherenkovVec ){

    // This branch is a vector nFibers long that contains the number of cherenkovs
    // produced in each fiber during an event
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *nCherenkovVec );

  } else { //Full tree

    if(thisIsOptical){
      // This branch is a vector nFibers long that contains the number of cherenkovs
      // produced in each fiber during an event. Make this first so it's columnId is 0
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );

      //vector< int > branch
      m_ZDCintVec->at(zdcNo).resize(1);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(0) );

      //Resize the vector for the number of optical branches storing double vectors
      m_ZDCdblVec->at(zdcNo).resize(7);
      //vector< double > branches
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_ZDCdblVec->at(zdcNo).at(0) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_ZDCdblVec->at(zdcNo).at(1) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_ZDCdblVec->at(zdcNo).at(2) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_ZDCdblVec->at(zdcNo).at(3) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_ZDCdblVec->at(zdcNo).at(4) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_ZDCdblVec->at(zdcNo).at(5) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_ZDCdblVec->at(zdcNo).at(6) );

    //The next branches are only created if the optical flag is off
    } else {
      m_ZDCdblVec->at(zdcNo).resize(12);
      //Non-optical vector< double > branches
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
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_ZDCdblVec->at(zdcNo).at(11) );

      //Non-optical vector< int > branches
      m_ZDCintVec->at(zdcNo).resize(5);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_ZDCintVec->at(zdcNo).at(0) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_ZDCintVec->at(zdcNo).at(1) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_ZDCintVec->at(zdcNo).at(2) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_ZDCintVec->at(zdcNo).at(3) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_ZDCintVec->at(zdcNo).at(4) );
    }
  }//end if nCherenkovVec
  m_analysisManager->FinishNtuple( );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::MakeRPDTree( G4int nTupleNo, G4int rpdNo, std::vector< int >* nCherenkovVec, G4bool thisIsOptical ){
  char name[20];
  sprintf(name,"RPD%dtree",rpdNo+1);
  m_analysisManager->CreateNtuple( name, "RPD data");

  // If nCherenkovVec is defined for this SD the user has set the ReducedTree flag
  // Otherwise generate the full tree
  if( nCherenkovVec ){

    // This branch is a vector nFibers long that contains the number of cherenkovs
    // produced in each fiber during an event
    m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", *nCherenkovVec );

  } else { //Full tree


    //This is done first so it can have an columnId of 0 if optical is on
    if(thisIsOptical){
      // This branch is a vector nFibers long that contains the number of cherenkovs
      // produced in each fiber during an event. Make this first so it's columnId is 0
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs"                               );

      //vector< int > branch
      m_RPDintVec->at(rpdNo).resize(1);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );

      //Resize the vector for the number of optical branches storing double vectors
      m_RPDdblVec->at(rpdNo).resize(7);
      //vector< double > branches
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "x",           m_RPDdblVec->at(rpdNo).at(0) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "y",           m_RPDdblVec->at(rpdNo).at(1) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "z",           m_RPDdblVec->at(rpdNo).at(2) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Px",          m_RPDdblVec->at(rpdNo).at(3) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Py",          m_RPDdblVec->at(rpdNo).at(4) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "Pz",          m_RPDdblVec->at(rpdNo).at(5) );
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_RPDdblVec->at(rpdNo).at(6) );

    //The next branches are only created if the optical flag is off
    } else {
      m_RPDdblVec->at(rpdNo).resize(12); //Increase the size to store more branches
      //Non-optical vector< double > branches
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
      m_analysisManager->CreateNtupleDColumn( nTupleNo, "time",        m_RPDdblVec->at(rpdNo).at(11));

      //Non-optical vector< int > branches
      m_RPDintVec->at(rpdNo).resize(5);
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "rodNo",       m_RPDintVec->at(rpdNo).at(0) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "modNo",       m_RPDintVec->at(rpdNo).at(1) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "nCherenkovs", m_RPDintVec->at(rpdNo).at(2) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "trackID",     m_RPDintVec->at(rpdNo).at(3) );
      m_analysisManager->CreateNtupleIColumn( nTupleNo, "pid",         m_RPDintVec->at(rpdNo).at(4) );
    }
  }//end if nCherenkovVec
  m_analysisManager->FinishNtuple( nTupleNo );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
