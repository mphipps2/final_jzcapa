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
// $Id: example.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file example.cc
/// \brief Main program of the  example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "SharedData.hh"

#include "PhysicsList.hh"
/*
#ifdef G4MULTITHREADED
#include "MyMTRunManager.hh"
#else
#include "MyRunManager.hh"
#endif
*/
#include "MyRunManager.hh"

#include "G4UImanager.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include <TString.h>
#include <TEnv.h>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  // Get some arguments for RunManager
  TString cfgName = "config/config.cfg";
  
  if( argc == 4 ){
    TString arg;
    arg = TString( argv[3] );
    if( arg.Contains("config") || arg.Contains("cfg") )
      cfgName = arg;
  }
  std::string outputName;
  outputName = "analysis/temp.root";
  if(argc == 3) {
    TString arg;
    arg = TString( argv[2]);
    if( arg.Contains(".root")  ) {
      outputName = arg;
    }
  }

  // Create SharedData and Run Manager
  SharedData* sharedData = new SharedData( outputName, cfgName.Data() );
  sharedData->Initialize();

  // Construct the default run manager
  //
  /*
    Does not work for now. Need to create multiple root files
    for various nodes. Otherwise there is a crash.
    (At least I think that is the reason)
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new MyMTRunManager( sharedData );
  runManager->SetNumberOfThreads(7);
#else
  G4RunManager* runManager = new MyRunManager( sharedData );
#endif
  */
  G4RunManager* runManager = new MyRunManager( sharedData );
  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new DetectorConstruction( sharedData ) );
  
  // Physics list
  TEnv* config = sharedData->GetConfig();
  std::string physicsListName = config->GetValue("physicsList","FTFP_BERT");

  // for now passing constructors directly
  // if need to modify something (SetSomething)
  // then take it out and add via this pointer.
  //  G4VModularPhysicsList* physicsList = NULL;
  runManager->SetUserInitialization( new PhysicsList(physicsListName, sharedData));
    
  // User action initialization
  runManager->SetUserInitialization(new ActionInitialization(sharedData));
  
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !
  sharedData->Finalize();
  std::cerr << "\nAbout to delete sharedData - " << sharedData << std::endl;  
  // WHY THIS NOT WORKING???
  // delete sharedData; 
  
  delete visManager;
  delete runManager;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
