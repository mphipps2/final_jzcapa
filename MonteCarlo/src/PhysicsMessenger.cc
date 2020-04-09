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
/// \file src/PhysicsMessenger.cc
/// \brief Implementation of the PhysicsMessenger class
//
//

#include "PhysicsMessenger.hh"
#include "PhysicsList.hh"

#include "G4OpticalSurface.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


PhysicsMessenger::PhysicsMessenger(PhysicsList * pList)
:G4UImessenger(),fPhysicsList(pList)
{

  fPhysicsDir = new G4UIdirectory("/Physics/");
  fPhysicsDir->SetGuidance("PhysicsList options");

  fSelectListCmd = new G4UIcmdWithAString("/Physics/SelectList", this);
  fSelectListCmd->SetGuidance("Choose a physics list to use");
  fSelectListCmd->SetParameterName("List Name",false);
  fSelectListCmd->SetDefaultValue("FTFP_BERT");
}

/*
 *
 */
PhysicsMessenger::~PhysicsMessenger(){
  delete fSelectListCmd;
}

/*
 *
 */
void PhysicsMessenger::SetNewValue(G4UIcommand* command,G4String newValue){
  if (command == fSelectListCmd) {
    fPhysicsList->SetList(newValue);
  }
}
