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
// $Id: SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "Analysis.hh"
#include "DetectorConstruction.hh"
#include "FiberSD.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction( )
: G4UserSteppingAction(){
	// Get number of OPTICAL flag from DetectorConstruction
	const DetectorConstruction* constDetectorConstruction
		= static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	DetectorConstruction* detectorConstruction
		= const_cast<DetectorConstruction*>(constDetectorConstruction);

	OPTICAL = detectorConstruction->GetOpticalFlag();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{delete G4AnalysisManager::Instance();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(__attribute__((unused)) const G4Step* theStep  )
{
	G4Track* theTrack = theStep->GetTrack();

	if(theTrack->GetParentID()==0 && theStep->IsLastStepInVolume()){
		lastStep = theTrack->GetPosition().getZ();
    auto analysisManager = G4AnalysisManager::Instance();
    for(int i = 0; i < analysisManager->GetNofNtuples(); i++){
      analysisManager->FillNtupleDColumn(i,0,lastStep);
    }
	}

	// If we are tracking a photon
	if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ){
		FiberSD* sd = (FiberSD*)theTrack->GetVolume()->GetLogicalVolume()->GetSensitiveDetector();
		// If World OPTICAL is on
		if(OPTICAL){
			// Kill photons only in SD volumes with OPTICAL off
			if(sd != 0 && !sd->OpticalIsOn() ){
				theTrack->SetTrackStatus( fStopAndKill );
			}
		} else { // World OPTICAL is off
			// Kill all photons except those in SD volumes with OPTICAL on
			if(sd == 0 || !sd->OpticalIsOn() ){
				theTrack->SetTrackStatus( fStopAndKill );
			}
		}
	}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
