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
// Author: Michael Phipps

#include "DetectorConstruction.hh"
#include "SharedData.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <iostream>
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), m_sd(NULL),
    m_solidWorld(NULL), m_logicWorld(NULL), m_physWorld(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction( SharedData* sd )
  : G4VUserDetectorConstruction(), m_sd( sd ), 
    m_solidWorld(NULL), m_logicWorld(NULL), m_physWorld(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
  if ( m_physWorld ) {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
  }

  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction :: DefineBorderProperties()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //----------------------------------------------     
  // Set Some Values
  //----------------------------------------------
  G4double modSizeX[5]; G4double modSizeY[5]; G4double modSizeZ[5]; G4double modCasingThickness[5]; G4double modWidth[5]; G4String modAbsorberMat[5]; G4double modAbsorberThickness[5]; G4double modAbsorberHeight[5]; G4double modAbsorberWidth[5]; G4double modCoreDiameter[5]; G4double modCladdingThickness[5]; G4int modNRadiators[5]; G4int modNAbsorbers[5]; G4int modType[5]; bool cladding[5]; G4double modRadiatorGapLength[5]; G4double stripPitch[5]; G4int modNStripsPerGap[5];
  
  G4int    nModules            = config->GetValue( "nModules",4);
  modType[0]            = config->GetValue( "mod1Type",1);
  modType[1]            = config->GetValue( "mod2Type",2);
  modType[2]            = config->GetValue( "mod3Type",3);
  modType[3]            = config->GetValue( "mod4Type",3);
  modType[4]            = config->GetValue( "mod5Type",3);    
  //  std::cout << " nMods " << nModules << " modtype1 " << modType[0] << " mod2 " << modType[1] << " mod3 " << modType[2] << " mod4 " << modType[3] << std::endl;
  // Calculate modSizeX, Y and Z for each mod. Mod size includes casings
  // In custom case, this is calculated by adding constituent lengths
  // For types 1-3, this is hard-coded directly
  for (int i = 0; i < nModules; ++i) {
    if (modType[i] == 4) {
      char variable[256];
      std::string modCladding;
      sprintf(variable,"mod%dCasingThickness",i+1);
      modCasingThickness[i] = config->GetValue(variable,.605);
      sprintf(variable,"mod%dAbsorberThickness",i+1);
      modAbsorberThickness[i] = config->GetValue(variable,10.);
      sprintf(variable,"mod%dAbsorberHeight",i+1);
      modAbsorberHeight[i] = config->GetValue(variable,180.);
      sprintf(variable,"mod%dAbsorberWidth",i+1);
      modAbsorberWidth[i] = config->GetValue(variable,100.);
      sprintf(variable,"mod%dRadiatorGapLength",i+1);
      modRadiatorGapLength[i] = config->GetValue(variable,2.);
      sprintf(variable,"mod%dCoreDiameter",i+1);
      modCoreDiameter[i] = config->GetValue(variable,1.5);
      sprintf(variable,"mod%dCladding",i+1);
      modCladding = config->GetValue(variable,"true");
      sprintf(variable,"mod%dCladdingThickness",i+1);
      modCladdingThickness[i] = config->GetValue(variable,0.605);
      sprintf(variable,"mod%dNRadiators",i+1);
      modNRadiators[i] = config->GetValue(variable,12);
      sprintf(variable,"mod%dNStripsPerGap",i+1);
      modNStripsPerGap[i] = config->GetValue(variable,52);
      modNAbsorbers[i] = modNRadiators[i] - 1;
      stripPitch[i] = modCoreDiameter[i] + 2*modCladdingThickness[i];
      if (modNRadiators[i] != 0) modWidth[i] = modNStripsPerGap[i]*stripPitch[i];
      else modWidth[i] = modAbsorberWidth[i];
      std::transform(modCladding.begin(), modCladding.end(), modCladding.begin(), ::tolower);
      cladding[i] = modCladding == "true" ?  true : false;
      if (cladding[i]) modCladdingThickness[i] = 0.;
      if (modNRadiators[i] == 0) {
	modNAbsorbers[i] = 1; // the case where you are defining a solid absorber block with no active channels
      }
      modSizeX[i] = modWidth[i] + 2*modCasingThickness[i];
      modSizeY[i] = 2*modCasingThickness[i]+modAbsorberHeight[i];
      modSizeZ[i] = 2*modCasingThickness[i]+modNRadiators[i]*modRadiatorGapLength[i] + modNAbsorbers[i]*modAbsorberThickness[i];
      
    }
    else {
      modSizeX[i] = 90.78;
      modSizeY[i] = 200.;
      modSizeZ[i] = 150.;
    }
  }

  
  // Option to switch on/off checking of volumes overlaps
  //
  bool checkOverlaps = true;
  
  //----------------------------------------------     
  // World
  //----------------------------------------------
  G4double maxModSizeX = 0; 
  G4double maxModSizeY = 0;
  G4double totalModSizeZ = 0;
  for (int i = 0; i < nModules; ++i) {
    if (modSizeX[i] > maxModSizeX) maxModSizeX = modSizeX[i];
    if (modSizeY[i] > maxModSizeY) maxModSizeY = modSizeY[i];
    totalModSizeZ += modSizeZ[i];
    //    std::cout << " modSizeZ " << modSizeZ[i] << std::endl;
  }

  G4double worldSizeX       = 1.1 * maxModSizeX * mm;    // mm
  G4double worldSizeY       = 1.1 * maxModSizeY * mm;    // mm
  G4double worldSizeZ       = 1.2 * totalModSizeZ * mm;
  //  std::cout << " totalmodsizeZ " << totalModSizeZ << " worldSizeZ " << worldSizeZ << std::endl;
  
  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");
  
  printf( "Building world with x %5.1f y %5.1f z %5.1f\n",
	  worldSizeX, worldSizeY, worldSizeZ );
  
  m_solidWorld =    
    new G4Box("World",              //its name
	      0.5*worldSizeX,       //its size
	      0.5*worldSizeY,
	      0.5*worldSizeZ );   
  
  m_logicWorld =                         
    new G4LogicalVolume(m_solidWorld,     //its solid
                        g4Air,            //its material
                        "World");         //its name
                                   
  m_physWorld = 
    new G4PVPlacement(0,                  //no rotation
                      G4ThreeVector(),    //at (0,0,0)
                      m_logicWorld,       //its logical volume
                      "World",            //its name
                      0,                  //its mother  volume
                      false,              //no boolean operation
                      0,                  //copy number
                      checkOverlaps);     //overlaps checking
  
  //----------------------------------------------     
  // Build EMCal Module
  //----------------------------------------------
    // air spacing between modules
  G4double dZ = 10 * mm;

  G4double zCoord = (totalModSizeZ / 2) + ((dZ*(nModules-1))/2) - modSizeZ[0]/2;
  //  std::cout << " zCoord " << zCoord << " dZ " << dZ << std::endl;  
  for(int i=0; i < nModules; ++i) {
    G4ThreeVector globalPos = G4ThreeVector(0, 0, zCoord ); // 1 cm b/w modules
    if (modType[i] == 1) {
      ModType1 *mod = new ModType1(i,globalPos,m_logicWorld,m_sd);
      mod->Construct();
    }
    else if (modType[i] == 2) {
      ModType2 *mod = new ModType2(i,globalPos,m_logicWorld,m_sd);
      mod->Construct();
    }
    else if (modType[i] == 3) {
      ModType3 *mod = new ModType3(i,globalPos,m_logicWorld,m_sd);
      mod->Construct();
    }
    else if (modType[i] == 4) {
      ModTypeCustom *mod = new ModTypeCustom(i,globalPos,m_logicWorld,m_sd);
      mod->Construct();
    }
    if (i != nModules - 1) zCoord -= ( modSizeZ[i]/2 + modSizeZ[i+1]/2 + dZ );
    //    std::cout << " z coord " << zCoord << " incremented by " << modSizeZ[i]/2 + modSizeZ[i+1]/2 + dZ  << std::endl;
  }
  
  return m_physWorld;
}
