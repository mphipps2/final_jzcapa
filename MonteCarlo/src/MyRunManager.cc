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
// $Id: MyRunManager.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file MyRunManager.cc
/// \brief Implementation of the MyRunManager class

#include "MyRunManager.hh"
#include "SharedData.hh"

#include "G4Timer.hh"
 
#include "G4RunManagerKernel.hh"
 
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VPersistencyManager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <sstream>

#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyRunManager::MyRunManager()
  : m_sd( NULL )
{}

MyRunManager::MyRunManager( SharedData* sd )
  : m_sd( sd )
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


MyRunManager::~MyRunManager()
{}

void MyRunManager::TerminateOneEvent()
{
  StackPreviousEvent(currentEvent);
  currentEvent = 0;

  // This is new, call shared data's end of event
  m_sd->EndOfEvent();
  
  numberOfEventProcessed++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

