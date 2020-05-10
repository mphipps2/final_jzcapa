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
/// \file include/PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
//
//

#include "PrimaryGeneratorMessenger.hh"

#include <sstream>
#include <iostream>

#include "G4OpticalSurface.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction * Gen)
:G4UImessenger(),fGenerator(Gen)
{
  fGeneratorDir = new G4UIdirectory("/beam/");
  fGeneratorDir->SetGuidance("");

  fLHCDir = new G4UIdirectory("/beam/LHC/");
  fLHCDir->SetGuidance("");

  fSPSDir = new G4UIdirectory("/beam/SPS/");
  fSPSDir->SetGuidance("");

  fFNALDir = new G4UIdirectory("/beam/FNAL/");
  fFNALDir->SetGuidance("");

  fBeamTypeCmd = new G4UIcmdWithAString("/beam/type", this);
  fBeamTypeCmd->SetGuidance("Set the beam type. LHC, SPS, or FNAL");
  fBeamTypeCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fBeamTypeCmd->SetToBeBroadcasted(false);

  fVerticalCrossingCmd = new G4UIcmdWithADoubleAndUnit("/beam/LHC/verticalCrossingAngle", this);
  fVerticalCrossingCmd->SetGuidance("Set vertical crossing angle for LHC beam");
  fVerticalCrossingCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fVerticalCrossingCmd->SetToBeBroadcasted(false);
  fVerticalCrossingCmd->SetDefaultUnit("mrad");

  fHorizontalCrossingCmd = new G4UIcmdWithADoubleAndUnit("/beam/LHC/horizontalCrossingAngle", this);
  fHorizontalCrossingCmd->SetGuidance("Set horizontal crossing angle for LHC beam");
  fHorizontalCrossingCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fHorizontalCrossingCmd->SetToBeBroadcasted(false);
  fHorizontalCrossingCmd->SetDefaultUnit("mrad");

  fProjectBeamCmd = new G4UIcmdWithADoubleAndUnit("/beam/projectBeam", this);
  fProjectBeamCmd->SetGuidance("Set the plane in Z which the beam will be projected to");
  fProjectBeamCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fProjectBeamCmd->SetToBeBroadcasted(false);
  fProjectBeamCmd->SetDefaultUnit("mm");

  fBeamPosCmd = new G4UIcmdWith3VectorAndUnit("/beam/pos", this);
  fBeamPosCmd->SetGuidance("Set the origin position of the beam");
  fBeamPosCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fBeamPosCmd->SetToBeBroadcasted(false);
  fBeamPosCmd->SetParameterName("X","Y","Z",true);
  fBeamPosCmd->SetDefaultUnit("mm");



}

/*
 *
 */
PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger(){
  delete fBeamTypeCmd;
  delete fVerticalCrossingCmd;
  delete fHorizontalCrossingCmd;
  delete fProjectBeamCmd;
  delete fBeamPosCmd;
}

/*
 *
 */
void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == fBeamTypeCmd){
    newValue.toLower();
    fGenerator->SetBeamType(newValue);
  }
  else if(command == fVerticalCrossingCmd){
    fGenerator->SetVerticalCrossingAngle( fVerticalCrossingCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fHorizontalCrossingCmd){
    fGenerator->SetHorizontalCrossingAngle( fHorizontalCrossingCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fProjectBeamCmd){
    fGenerator->SetProjectionPlane( fProjectBeamCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fBeamPosCmd){
    fGenerator->SetBeamPos( new G4ThreeVector( fBeamPosCmd->GetNew3VectorValue(newValue) ) );
  }
}
