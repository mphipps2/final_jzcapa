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
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4HCtable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Version.hh"

#if G4VERSION_NUMBER > 999
    #include "G4AccumulableManager.hh"
#else
    #include "G4ParameterManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(SharedData *sd)
: G4UserRunAction()
{
    fSharedData = sd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //  G4int nbEventInRun = run->GetNumberOfEventToBeProcessed();

  fSharedData->AddOutputToTree("ID",&fTrackID_v);
  fSharedData->AddOutputToTree("ModNb",&fModNb_v);
  fSharedData->AddOutputToTree("RadNb",&fRadNb_v);
  fSharedData->AddOutputToTree("RodNb",&fRodNb_v);
  fSharedData->AddOutputToTree("EDep",&fEdep_v);
  fSharedData->AddOutputToTree("Pid",&fPid_v);
  fSharedData->AddOutputToTree("X",&fX_v);
  fSharedData->AddOutputToTree("Y",&fY_v);
  fSharedData->AddOutputToTree("Z",&fZ_v);
  fSharedData->AddOutputToTree("Px",&fPx_v);
  fSharedData->AddOutputToTree("Py",&fPy_v);
  fSharedData->AddOutputToTree("Pz",&fPz_v);
  fSharedData->AddOutputToTree("EventNo",&fEventNo_v);
  fSharedData->AddOutputToTree("Energy",&fEnergy_v);
  fSharedData->AddOutputToTree("Charge",&fCharge_v);
  fSharedData->AddOutputToTree("Velocity",&fVelocity_v);
  fSharedData->AddOutputToTree("NCherenkovs",&fNCherenkovs_v);
  fSharedData->AddOutputToTree("Beta",&fBeta_v);
  
  // reset parameters to their initial values
#if G4VERSION_NUMBER > 999
    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
#else
    G4ParameterManager* parameterManager = G4ParameterManager::Instance();
#endif
  parameterManager->Reset();   

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge parameters 
#if G4VERSION_NUMBER > 999
    G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
#else
    G4ParameterManager* parameterManager = G4ParameterManager::Instance();
#endif
    parameterManager->Merge();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

