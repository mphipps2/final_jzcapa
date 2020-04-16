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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
/// \author Chad Lantz
/// \date 16 April 2020

#include "DetectorMessenger.hh"

#include "G4OpticalSurface.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det)
{
  fDetDir = new G4UIdirectory("/Detector/");
  fDetDir->SetGuidance("Geometry creation and modification");

  fRPDDir = new G4UIdirectory("/Detector/ZDC/");
  fRPDDir->SetGuidance("ZDC creation and modification");

  fZDCDir = new G4UIdirectory("/Detector/RPD/");
  fZDCDir->SetGuidance("RPD creation and modification");

  fClusterCmd = new G4UIcmdWithABool("/Detector/Cluster", this);
  fClusterCmd->SetGuidance("Set cluster flag");
  fClusterCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fClusterCmd->SetDefaultValue(true);
  fClusterCmd->SetToBeBroadcasted(false);

  fOpticalCmd = new G4UIcmdWithABool("/Detector/Optical", this);
  fOpticalCmd->SetGuidance("Set optical flag");
  fOpticalCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fOpticalCmd->SetDefaultValue(true);
  fOpticalCmd->SetToBeBroadcasted(false);

  fOverlapsCmd = new G4UIcmdWithABool("/Detector/Overlaps", this);
  fOverlapsCmd->SetGuidance("Set check overlaps flag for current ZDC");
  fOverlapsCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fOverlapsCmd->SetDefaultValue(true);
  fOverlapsCmd->SetToBeBroadcasted(false);

  fForcePositionCmd = new G4UIcmdWithABool("/Detector/ForcePosition", this);
  fForcePositionCmd->SetGuidance("Force detector positions");
  fForcePositionCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fForcePositionCmd->SetDefaultValue(true);
  fForcePositionCmd->SetToBeBroadcasted(false);

  fSetRunNumberCmd = new G4UIcmdWithAnInteger("/Detector/RunNumber", this);
  fSetRunNumberCmd->SetGuidance("Set run number if loading configuration from file");
  fSetRunNumberCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fSetRunNumberCmd->SetToBeBroadcasted(false);

  fConfigFileCmd = new G4UIcmdWithAString("/Detector/ConfigFile", this);
  fConfigFileCmd->SetGuidance("Set configuration file name");
  fConfigFileCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fConfigFileCmd->SetToBeBroadcasted(false);

  fPrintDebugCmd = new G4UIcmdWithAString("/Detector/PrintDebugStatement", this);
  fPrintDebugCmd->SetGuidance("Print a string to the console to help debug");
  fPrintDebugCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fPrintDebugCmd->SetDefaultValue("DEBUG");
  fPrintDebugCmd->SetToBeBroadcasted(false);

  fSetWorldDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/Detector/SetWorldDimensions", this);
  fSetWorldDimensionsCmd->SetGuidance("Set world volume dimensions");
  fSetWorldDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fSetWorldDimensionsCmd->SetToBeBroadcasted(false);
  fSetWorldDimensionsCmd->SetParameterName("X","Y","Z",true);
  fSetWorldDimensionsCmd->SetDefaultValue(G4ThreeVector(2.,2.,4.));
  fSetWorldDimensionsCmd->SetDefaultUnit("m");


  //ZDC commands
  fZDCAddCmd = new G4UIcmdWithoutParameter("/Detector/ZDC/Add", this);
  fZDCAddCmd->SetGuidance("Add a ZDC to the world");
  fZDCAddCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCAddCmd->SetToBeBroadcasted(false);

  fZDCPositionCmd = new G4UIcmdWith3VectorAndUnit("/Detector/ZDC/Position", this);
  fZDCPositionCmd->SetGuidance("Set current ZDC position");
  fZDCPositionCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCPositionCmd->SetToBeBroadcasted(false);
  fZDCPositionCmd->SetParameterName("X","Y","Z",true);
  fZDCPositionCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  fZDCPositionCmd->SetDefaultUnit("mm");

  fZDCFiberDiametersCmd = new G4UIcmdWith3VectorAndUnit("/Detector/ZDC/FiberDiameters", this);
  fZDCFiberDiametersCmd->SetGuidance("Set current ZDC fiber core, cladding, and buffer diameters");
  fZDCFiberDiametersCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCFiberDiametersCmd->SetToBeBroadcasted(false);
  fZDCFiberDiametersCmd->SetParameterName("Core","Cladding","Buffer",true);
  fZDCFiberDiametersCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  fZDCFiberDiametersCmd->SetDefaultUnit("mm");

  fZDCAbsorberDimensionsCmd = new G4UIcmdWith3VectorAndUnit("/Detector/ZDC/AbsorberDimensions", this);
  fZDCAbsorberDimensionsCmd->SetGuidance("Set current ZDC absorber dimensions");
  fZDCAbsorberDimensionsCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCAbsorberDimensionsCmd->SetToBeBroadcasted(false);
  fZDCAbsorberDimensionsCmd->SetParameterName("X","Y","Z",true);
  fZDCAbsorberDimensionsCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  fZDCAbsorberDimensionsCmd->SetDefaultUnit("mm");

  fZDCnAbsorbersCmd = new G4UIcmdWithAnInteger("/Detector/ZDC/nAbsorbers", this);
  fZDCnAbsorbersCmd->SetGuidance("Set number of absorbers in the current ZDC");
  fZDCnAbsorbersCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCnAbsorbersCmd->SetDefaultValue(11);
  fZDCnAbsorbersCmd->SetToBeBroadcasted(false);

  fZDCHousingThicknessCmd =
    new G4UIcmdWithADoubleAndUnit("/Detector/ZDC/HousingThickness",this);
  fZDCHousingThicknessCmd->SetGuidance("Set housing thickness of the current ZDC");
  fZDCHousingThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCHousingThicknessCmd->SetToBeBroadcasted(false);
  fZDCHousingThicknessCmd->SetDefaultUnit("mm");

  fZDCGapThicknessCmd =
    new G4UIcmdWithADoubleAndUnit("/Detector/ZDC/GapThickness",this);
  fZDCGapThicknessCmd->SetGuidance("Set radiator gap thickness of the current ZDC");
  fZDCGapThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCGapThicknessCmd->SetToBeBroadcasted(false);
  fZDCGapThicknessCmd->SetDefaultUnit("mm");

  fZDCOpticalFlagCmd = new G4UIcmdWithABool("/Detector/ZDC/Optical", this);
  fZDCOpticalFlagCmd->SetGuidance("Set optical flag for current ZDC");
  fZDCOpticalFlagCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCOpticalFlagCmd->SetDefaultValue(true);
  fZDCOpticalFlagCmd->SetToBeBroadcasted(false);

  fZDCOverlapsFlagCmd = new G4UIcmdWithABool("/Detector/ZDC/CheckOverlaps", this);
  fZDCOverlapsFlagCmd->SetGuidance("Set check overlaps flag for current ZDC");
  fZDCOverlapsFlagCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCOverlapsFlagCmd->SetDefaultValue(true);
  fZDCOverlapsFlagCmd->SetToBeBroadcasted(false);

  fZDCHousingMaterialCmd = new G4UIcmdWithAString("/Detector/ZDC/HousingMaterial", this);
  fZDCHousingMaterialCmd->SetGuidance("Set current ZDC housing material");
  fZDCHousingMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCHousingMaterialCmd->SetToBeBroadcasted(false);

  fZDCAbsorberMaterialCmd = new G4UIcmdWithAString("/Detector/ZDC/AbsorberMaterial", this);
  fZDCAbsorberMaterialCmd->SetGuidance("Set current ZDC absorber material");
  fZDCAbsorberMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCAbsorberMaterialCmd->SetToBeBroadcasted(false);

  fZDCSetCurrentCmd = new G4UIcmdWithAnInteger("/Detector/ZDC/SetCurrent", this);
  fZDCSetCurrentCmd->SetGuidance("Select ZDC to be modified");
  fZDCSetCurrentCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCSetCurrentCmd->SetToBeBroadcasted(false);

  fZDCDuplicateCmd = new G4UIcmdWithAnInteger("/Detector/ZDC/Duplicate", this);
  fZDCDuplicateCmd->SetGuidance("Select duplicate the ZDC with the given index");
  fZDCDuplicateCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fZDCDuplicateCmd->SetDefaultValue(0);
  fZDCDuplicateCmd->SetToBeBroadcasted(false);


  //RPD commands
  fRPDAddCmd = new G4UIcmdWithoutParameter("/Detector/RPD/Add", this);
  fRPDAddCmd->SetGuidance("Add an RPD to the world");
  fRPDAddCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDAddCmd->SetToBeBroadcasted(false);

  fRPDPositionCmd = new G4UIcmdWith3VectorAndUnit("/Detector/RPD/Position", this);
  fRPDPositionCmd->SetGuidance("Set current RPD position");
  fRPDPositionCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDPositionCmd->SetToBeBroadcasted(false);
  fRPDPositionCmd->SetParameterName("X","Y","Z",true);
  fRPDPositionCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  fRPDPositionCmd->SetDefaultUnit("mm");

  fRPDFiberDiametersCmd = new G4UIcmdWith3VectorAndUnit("/Detector/RPD/FiberDiameters", this);
  fRPDFiberDiametersCmd->SetGuidance("Set current RPD fiber core, cladding, and buffer diameters");
  fRPDFiberDiametersCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDFiberDiametersCmd->SetToBeBroadcasted(false);
  fRPDFiberDiametersCmd->SetParameterName("Core","Cladding","Buffer",true);
  fRPDFiberDiametersCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  fRPDFiberDiametersCmd->SetDefaultUnit("mm");

  fRPDHousingThicknessCmd =
    new G4UIcmdWithADoubleAndUnit("/Detector/RPD/HousingThickness",this);
  fRPDHousingThicknessCmd->SetGuidance("Set housing thickness of the current RPD");
  fRPDHousingThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDHousingThicknessCmd->SetToBeBroadcasted(false);
  fRPDHousingThicknessCmd->SetDefaultUnit("mm");

  fRPDSetFiberPitchCmd =
    new G4UIcmdWithADoubleAndUnit("/Detector/RPD/FiberPitch",this);
  fRPDSetFiberPitchCmd->SetGuidance("Set fiber pitch of the current RPD");
  fRPDSetFiberPitchCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDSetFiberPitchCmd->SetToBeBroadcasted(false);
  fRPDSetFiberPitchCmd->SetDefaultUnit("mm");

  fRPDSetTileSizeCmd =
    new G4UIcmdWithADoubleAndUnit("/Detector/RPD/TileSize",this);
  fRPDSetTileSizeCmd->SetGuidance("Set tile size of the current RPD");
  fRPDSetTileSizeCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDSetTileSizeCmd->SetToBeBroadcasted(false);
  fRPDSetTileSizeCmd->SetDefaultUnit("mm");

  fRPDMinWallThicknessCmd =
    new G4UIcmdWithADoubleAndUnit("/Detector/RPD/MinWallThickness",this);
  fRPDMinWallThicknessCmd->SetGuidance("Set minimum distance between fibers of the current RPD");
  fRPDMinWallThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDMinWallThicknessCmd->SetToBeBroadcasted(false);
  fRPDMinWallThicknessCmd->SetDefaultUnit("mm");

  fRPDTypeCmd = new G4UIcmdWithAString("/Detector/RPD/RPDtype", this);
  fRPDTypeCmd->SetGuidance("Set design type for the current RPD (cms or panflute)");
  fRPDTypeCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDTypeCmd->SetToBeBroadcasted(false);

  fRPDOpticalFlagCmd = new G4UIcmdWithABool("/Detector/RPD/Optical", this);
  fRPDOpticalFlagCmd->SetGuidance("Set optical flag for current RPD");
  fRPDOpticalFlagCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDOpticalFlagCmd->SetDefaultValue(true);
  fRPDOpticalFlagCmd->SetToBeBroadcasted(false);

  fRPDOverlapsFlagCmd = new G4UIcmdWithABool("/Detector/RPD/CheckOverlaps", this);
  fRPDOverlapsFlagCmd->SetGuidance("Set check overlaps flag for current RPD");
  fRPDOverlapsFlagCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDOverlapsFlagCmd->SetDefaultValue(true);
  fRPDOverlapsFlagCmd->SetToBeBroadcasted(false);

  fRPDSetCurrentCmd = new G4UIcmdWithAnInteger("/Detector/RPD/SetCurrent", this);
  fRPDSetCurrentCmd->SetGuidance("Select RPD to be modified");
  fRPDSetCurrentCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDSetCurrentCmd->SetToBeBroadcasted(false);

  fRPDDuplicateCmd = new G4UIcmdWithAnInteger("/Detector/RPD/Duplicate", this);
  fRPDDuplicateCmd->SetGuidance("Select duplicate the RPD with the given index");
  fRPDDuplicateCmd->AvailableForStates(G4State_PreInit, G4State_Init, G4State_Idle);
  fRPDDuplicateCmd->SetDefaultValue(0);
  fRPDDuplicateCmd->SetToBeBroadcasted(false);

}

/*
 *
 */
DetectorMessenger::~DetectorMessenger(){

  delete fClusterCmd;
  delete fOpticalCmd;
  delete fOverlapsCmd;
  delete fForcePositionCmd;
  delete fSetRunNumberCmd;
  delete fConfigFileCmd;
  delete fPrintDebugCmd;
  delete fSetWorldDimensionsCmd;


  //ZDC Commands
  delete fZDCAddCmd;
  delete fZDCPositionCmd;
  delete fZDCFiberDiametersCmd;
  delete fZDCAbsorberDimensionsCmd;
  delete fZDCnAbsorbersCmd;
  delete fZDCHousingThicknessCmd;
  delete fZDCGapThicknessCmd;
  delete fZDCOpticalFlagCmd;
  delete fZDCOverlapsFlagCmd;
  delete fZDCHousingMaterialCmd;
  delete fZDCAbsorberMaterialCmd;
  delete fZDCSetCurrentCmd;
  delete fZDCDuplicateCmd;


  //RPD commands
  delete fRPDAddCmd;
  delete fRPDPositionCmd;
  delete fRPDFiberDiametersCmd;
  delete fRPDHousingThicknessCmd;
  delete fRPDSetFiberPitchCmd;
  delete fRPDSetTileSizeCmd;
  delete fRPDMinWallThicknessCmd;
  delete fRPDTypeCmd;
  delete fRPDOpticalFlagCmd;
  delete fRPDOverlapsFlagCmd;
  delete fRPDSetCurrentCmd;
  delete fRPDDuplicateCmd;


}

/*
 *
 */
void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fClusterCmd) {
    fDetector->SetClusterFlag( fClusterCmd->GetNewBoolValue(newValue) );
  }
  if (command == fOpticalCmd) {
    fDetector->SetClusterFlag( fOpticalCmd->GetNewBoolValue(newValue) );
  }
  if (command == fOverlapsCmd) {
    fDetector->SetClusterFlag( fOverlapsCmd->GetNewBoolValue(newValue) );
  }
  else if(command == fForcePositionCmd){
    fDetector->ForcePosition( fForcePositionCmd->GetNewBoolValue(newValue) );
    std::cout << "Forcing detector position" << std::endl;
  }
  else if(command == fConfigFileCmd){
    fDetector->SetConfigFileName( newValue );
  }
  else if(command == fPrintDebugCmd){
    std::cout << newValue << std::endl;
  }
  else if(command == fSetRunNumberCmd){
    fDetector->SetRunNumber( fSetRunNumberCmd->GetNewIntValue(newValue) );
  }
  else if(command == fSetWorldDimensionsCmd){
    fDetector->SetWorldDimensions( new G4ThreeVector( fSetWorldDimensionsCmd->GetNew3VectorValue(newValue) ) );
  }


  //ZDC Commands
  else if(command == fZDCAddCmd){
    fDetector->AddZDC();
  }
  else if(command == fZDCPositionCmd){
    fDetector->SetZDCPosition( new G4ThreeVector( fZDCPositionCmd->GetNew3VectorValue(newValue) ) );
  }
  else if(command == fZDCFiberDiametersCmd){
    fDetector->SetZDCFiberDiameters( new G4ThreeVector( fZDCFiberDiametersCmd->GetNew3VectorValue(newValue) ) );
  }
  else if(command == fZDCAbsorberDimensionsCmd){
    fDetector->SetZDCAbsorberDimensions( new G4ThreeVector( fZDCAbsorberDimensionsCmd->GetNew3VectorValue(newValue) ) );
  }
  else if(command == fZDCnAbsorbersCmd){
    fDetector->SetZDCnAbsorbers( fZDCnAbsorbersCmd->GetNewIntValue(newValue) );
  }
  else if(command == fZDCHousingThicknessCmd){
    fDetector->SetZDCHousingThickness( fZDCHousingThicknessCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fZDCGapThicknessCmd){
    fDetector->SetZDCGapThickness( fZDCGapThicknessCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fZDCOpticalFlagCmd){
    fDetector->SetZDCOpticalFlag( fZDCOpticalFlagCmd->GetNewBoolValue(newValue) );
  }
  else if(command == fZDCOverlapsFlagCmd){
    fDetector->SetZDCOverlapsFlag( fZDCOverlapsFlagCmd->GetNewBoolValue(newValue) );
  }
  else if(command == fZDCHousingMaterialCmd){
    fDetector->SetZDCHousingMaterial( newValue );
  }
  else if(command == fZDCAbsorberMaterialCmd){
    fDetector->SetZDCAbsorberMaterial( newValue );
  }
  else if(command == fZDCSetCurrentCmd){
    fDetector->SetCurrentZDC( fZDCSetCurrentCmd->GetNewIntValue(newValue) );
  }
  else if(command == fZDCDuplicateCmd){
    fDetector->DuplicateZDC( fZDCDuplicateCmd->GetNewIntValue(newValue) );
  }


  //RPD commands
  else if(command == fRPDAddCmd){
    fDetector->AddRPD();
  }
  else if(command == fRPDPositionCmd){
    fDetector->SetRPDPosition( new G4ThreeVector( fRPDPositionCmd->GetNew3VectorValue(newValue) ) );
  }
  else if(command == fRPDFiberDiametersCmd){
    fDetector->SetRPDFiberDimensions( new G4ThreeVector( fRPDFiberDiametersCmd->GetNew3VectorValue(newValue) ) );
  }
  else if(command == fRPDHousingThicknessCmd){
    fDetector->SetRPDHousingThickness( fRPDHousingThicknessCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fRPDSetFiberPitchCmd){
    fDetector->SetRPDFiberPitch( fRPDSetFiberPitchCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fRPDSetTileSizeCmd){
    fDetector->SetRPDTileSize( fRPDSetTileSizeCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fRPDMinWallThicknessCmd){
    fDetector->SetRPDMinWallThickness( fRPDMinWallThicknessCmd->GetNewDoubleValue(newValue) );
  }
  else if(command == fRPDTypeCmd){
    fDetector->SetRPDDetectorType( newValue );
  }
  else if(command == fRPDOpticalFlagCmd){
    fDetector->SetRPDOpticalFlag( fRPDOpticalFlagCmd->GetNewBoolValue(newValue) );
  }
  else if(command == fRPDOverlapsFlagCmd){
    fDetector->SetRPDOverlapsFlag( fRPDOverlapsFlagCmd->GetNewBoolValue(newValue) );
  }
  else if(command == fRPDSetCurrentCmd){
    fDetector->SetCurrentRPD( fRPDSetCurrentCmd->GetNewIntValue(newValue) );
  }
  else if(command == fRPDDuplicateCmd){
    fDetector->DuplicateRPD( fRPDDuplicateCmd->GetNewIntValue(newValue) );
  }
}
