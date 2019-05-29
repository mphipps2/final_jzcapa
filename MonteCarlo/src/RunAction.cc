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

void RunAction::BeginOfRunAction(__attribute__((unused)) const G4Run* run)
{ 
  // inform the runManager to save random number seed
  
  long seeds[2];
  long systime = time(NULL);
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds); 
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //G4int nbEventInRun = run->GetNumberOfEventToBeProcessed();
	
  fSharedData->AddOutputToZDCTree("ID",&fTrackID_v);
  fSharedData->AddOutputToZDCTree("ModNb",&fModNb_v);
  fSharedData->AddOutputToZDCTree("RadNb",&fRadNb_v);
  fSharedData->AddOutputToZDCTree("RodNb",&fRodNb_v);
  fSharedData->AddOutputToZDCTree("EDep",&fEdep_v);
  fSharedData->AddOutputToZDCTree("Pid",&fPid_v);
  fSharedData->AddOutputToZDCTree("X",&fX_v);
  fSharedData->AddOutputToZDCTree("Y",&fY_v);
  fSharedData->AddOutputToZDCTree("Z",&fZ_v);
  fSharedData->AddOutputToZDCTree("Px",&fPx_v);
  fSharedData->AddOutputToZDCTree("Py",&fPy_v);
  fSharedData->AddOutputToZDCTree("Pz",&fPz_v);
  fSharedData->AddOutputToZDCTree("EventNo",&fEventNo_v);
  fSharedData->AddOutputToZDCTree("Energy",&fEnergy_v);
  fSharedData->AddOutputToZDCTree("Charge",&fCharge_v);
  fSharedData->AddOutputToZDCTree("Velocity",&fVelocity_v);
  fSharedData->AddOutputToZDCTree("NCherenkovs",&fNCherenkovs_v);
  fSharedData->AddOutputToZDCTree("Beta",&fBeta_v);
  
  fSharedData->AddOutputToRPDTree("ID",&fTrackID_v2);
  fSharedData->AddOutputToRPDTree("ModNb",&fModNb_v2);
  fSharedData->AddOutputToRPDTree("RadNb",&fRadNb_v2);
  fSharedData->AddOutputToRPDTree("RodNb",&fRodNb_v2);
  fSharedData->AddOutputToRPDTree("EDep",&fEdep_v2);
  fSharedData->AddOutputToRPDTree("Pid",&fPid_v2);
  fSharedData->AddOutputToRPDTree("X",&fX_v2);
  fSharedData->AddOutputToRPDTree("Y",&fY_v2);
  fSharedData->AddOutputToRPDTree("Z",&fZ_v2);
  fSharedData->AddOutputToRPDTree("Px",&fPx_v2);
  fSharedData->AddOutputToRPDTree("Py",&fPy_v2);
  fSharedData->AddOutputToRPDTree("Pz",&fPz_v2);
  fSharedData->AddOutputToRPDTree("EventNo",&fEventNo_v2);
  fSharedData->AddOutputToRPDTree("Energy",&fEnergy_v2);
  fSharedData->AddOutputToRPDTree("Charge",&fCharge_v2);
  fSharedData->AddOutputToRPDTree("Velocity",&fVelocity_v2);
  fSharedData->AddOutputToRPDTree("NCherenkovs",&fNCherenkovs_v2);
  fSharedData->AddOutputToRPDTree("Beta",&fBeta_v2);
  
  fSharedData->AddOutputToFiberTree("ID",&fTrackID_v3);
  fSharedData->AddOutputToFiberTree("ModNb",&fModNb_v3);
  fSharedData->AddOutputToFiberTree("RadNb",&fRadNb_v3);
  fSharedData->AddOutputToFiberTree("RodNb",&fRodNb_v3);
  fSharedData->AddOutputToFiberTree("EDep",&fEdep_v3);
  fSharedData->AddOutputToFiberTree("Pid",&fPid_v3);
  fSharedData->AddOutputToFiberTree("X",&fX_v3);
  fSharedData->AddOutputToFiberTree("Y",&fY_v3);
  fSharedData->AddOutputToFiberTree("Z",&fZ_v3);
  fSharedData->AddOutputToFiberTree("Px",&fPx_v3);
  fSharedData->AddOutputToFiberTree("Py",&fPy_v3);
  fSharedData->AddOutputToFiberTree("Pz",&fPz_v3);
  fSharedData->AddOutputToFiberTree("EventNo",&fEventNo_v3);
  fSharedData->AddOutputToFiberTree("Energy",&fEnergy_v3);
  fSharedData->AddOutputToFiberTree("Charge",&fCharge_v3);
  fSharedData->AddOutputToFiberTree("Velocity",&fVelocity_v3);
  fSharedData->AddOutputToFiberTree("NCherenkovs",&fNCherenkovs_v3);
  fSharedData->AddOutputToFiberTree("Beta",&fBeta_v3);
  
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

