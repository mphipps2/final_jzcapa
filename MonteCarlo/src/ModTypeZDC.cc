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
// Michael Phipps
// For an explanation of the hierarchy scheme see: https://twiki.cern.ch/twiki/bin/view/Atlas/ZdcSimulation#Geometry_Implementation_Develope

#include "ModTypeZDC.hh"
#include "QuartzSD.hh"
#include "CherenkovSD.hh"
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
#include "G4SDManager.hh"

#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
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

ModTypeZDC::ModTypeZDC(const int cn,const G4ThreeVector& pos,
	      G4LogicalVolume* mother, SharedData* sd)
  : m_modNum( cn ),  m_pos( pos ), m_logicMother( mother ),
    m_sd( sd ),
    m_matQuartz(0),
    m_matAbsorber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeZDC::ModTypeZDC()
  : m_modNum( 0 ), m_pos(G4ThreeVector()), m_logicMother(NULL),
    m_sd(NULL),
    m_matQuartz(0), m_matAbsorber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeZDC::~ModTypeZDC()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeZDC::Construct(){
  DefineMaterials();
  ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeZDC::DefineMaterials()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();
  std::string simCherenkov_string = config->GetValue("simCherenkov","false");
  std::transform(simCherenkov_string.begin(), simCherenkov_string.end(), simCherenkov_string.begin(), ::tolower);
  m_simCherenkov = simCherenkov_string == "true" ? true : false;
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //----------------------------------------------
  // Define Materials
  //----------------------------------------------
  m_matQuartz    = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  // Define materials not in NIST
  G4double fractionMass;
  G4int ncomponents;

  G4Element*  O   = new G4Element ("O"  ,"O"  ,  8.0 ,  16.0*g/mole);
  //  G4Element*  Si  = new G4Element ("Si" ,"Si" , 14.0 ,  28.085*g/mole);
  G4Element*  W   = new G4Element ("W"  ,"W"  , 74.0 ,  183.84*g/mole);
  G4Element*  Fe  = new G4Element ("Fe" ,"Fe" , 26.0 ,  55.845 *g/mole);
  G4Element*  C   = new G4Element ("C"  ,"C"  ,  6.0 ,  12.0107*g/mole);
  G4Element*  Ni  = new G4Element ("Ni" ,"Ni" , 28.0 ,  58.6934*g/mole);
	G4Element*  N  = new G4Element ("N" ,"N" , 14.0 ,  28.014*g/mole);


  // Absorber composition:  savannah.cern.ch/task/download.php?file_id=22925
  m_matAbsorber = new G4Material("Tungsten",18.155*g/cm3,ncomponents=3);
  m_matAbsorber->AddElement(W,  fractionMass=0.948);
  m_matAbsorber->AddElement(Ni, fractionMass=0.037);
  m_matAbsorber->AddElement(Fe, fractionMass=0.015);

  m_matHousing  = new G4Material("Steel", 7.9*gram/cm3,ncomponents=2);
  m_matHousing->AddElement(Fe  , fractionMass=0.98);
  m_matHousing->AddElement(C   , fractionMass=0.02);

  m_matAir = new G4Material("Air",1.290*mg/cm3,ncomponents=2);
  m_matAir->AddElement(N, fractionMass=.70);
  m_matAir->AddElement(O, fractionMass=.30);
  //----------------------------------------------
  // Define Material Properties
  //----------------------------------------------
  const G4int NUMENTRIES = 2;

  G4double ephoton         [NUMENTRIES] = {2.00*eV,4.80*eV};

  G4double rindexCore[NUMENTRIES] = {1.46,1.46};

  //Fill in the Marterial properties table for each material.
  //Guide for undestanding Optical processes at http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Photo
  // Only needed if direct Cherenkov production turned on
  G4MaterialPropertiesTable *quartzMPT = new G4MaterialPropertiesTable();
  if (m_simCherenkov) {
    quartzMPT->AddProperty("RINDEX",ephoton,rindexCore,NUMENTRIES);
    //  quartzMPT->AddProperty("ABSLENGTH",ephoton,absorptionQuartz,NUMENTRIES);
    m_matQuartz->SetMaterialPropertiesTable(quartzMPT);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeZDC :: DefineBorderProperties(){
  const G4int NUMENTRIES = 2;
  G4double ephoton[NUMENTRIES] = {2.00*eV,4.80*eV};
  // Get Config
  //  TEnv* config = m_sd->GetConfig();

  //----------------------------------------------
  // Housing Skin
  //----------------------------------------------
  G4double housingReflectivity      [NUMENTRIES] = {0.4, 0.4};
  G4double housingEfficiency        [NUMENTRIES] = {0.10, 0.10};

  G4MaterialPropertiesTable* housingMPT
    = new G4MaterialPropertiesTable();
  G4OpticalSurface* housingOS =
    new G4OpticalSurface("HousingOpSurface",unified, polished, dielectric_metal);
  if (m_simCherenkov) {
    housingMPT->AddProperty("REFLECTIVITY", ephoton, housingReflectivity, NUMENTRIES);
    housingMPT->AddProperty("EFFICIENCY"  , ephoton, housingEfficiency  , NUMENTRIES);
    housingOS->SetMaterialPropertiesTable( housingMPT );
  }

  new G4LogicalSkinSurface("housingSkinSurface", m_SteelLogical, housingOS );

  //----------------------------------------------
  // Absorber Skin
  //----------------------------------------------
  G4double absorberReflectivity      [NUMENTRIES] = {0.0, 0.0};
  G4MaterialPropertiesTable* absorberMPT
    = new G4MaterialPropertiesTable();
  G4OpticalSurface* absorberOS =
    new G4OpticalSurface("AbsorberOpSurface",glisur, polished, dielectric_metal);
  if (m_simCherenkov) {
    absorberMPT->AddProperty("REFLECTIVITY", ephoton, absorberReflectivity, NUMENTRIES);
    absorberOS->SetMaterialPropertiesTable( absorberMPT );
  }

  new G4LogicalSkinSurface("absorberSkinSurface", m_WLogical, absorberOS );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeZDC::ConstructDetector()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();

  //----------------------------------------------
  // Set Some Values
  //----------------------------------------------
  float modCasingThickness[5]; float modAbsorberThickness[5]; float modAbsorberHeight[5]; float modAbsorberWidth[5]; G4String modAbsorberMat[5]; float modRadiatorGapLength[5]; float modCoreDiameter[5];  float modCladdingThickness[5]; int modNRadiators[5]; int modNAbsorbers[5];  bool cladding[5];  int modType[5]; int modNStripsPerGap[5];

  int nModules = config->GetValue("nModules",2);
    // Option to switch on/off checking of volumes overlaps
  bool checkOverlaps = config->GetValue("checkOverlaps",false);

  for (int i = 0; i < nModules; ++i) {
    char variable[256];
    sprintf(variable,"mod%dType",i+1);
    modType[i] = config->GetValue(variable,5);
    if (modType[i] == 5){
      std::string modCladding;

      sprintf(variable,"mod%dCasingThickness",5);
      modCasingThickness[i] = config->GetValue(variable,7.94);
      sprintf(variable,"mod%dAbsorberThickness",5);
      modAbsorberThickness[i] = config->GetValue(variable,10.);
      sprintf(variable,"mod%dAbsorberHeight",5);
      modAbsorberHeight[i] = config->GetValue(variable,180.);
      sprintf(variable,"mod%dAbsorberWidth",5);
      modAbsorberWidth[i] = config->GetValue(variable,100.);
      sprintf(variable,"mod%dAbsorberMat",5);
      modAbsorberMat[i] = config->GetValue(variable,"W");
      sprintf(variable,"mod%dRadiatorGapLength",5);
      modRadiatorGapLength[i] = config->GetValue(variable,2.);
      sprintf(variable,"mod%dCoreDiameter",5);
      modCoreDiameter[i] = config->GetValue(variable,1.5);
      sprintf(variable,"mod%dCladdingThickness",5);
      modCladdingThickness[i] = config->GetValue(variable,0.1);
      sprintf(variable,"mod%dNStripsPerGap",5);
      modNStripsPerGap[i] = config->GetValue(variable,52);
      sprintf(variable,"mod%dCladding",5);
      modCladding = config->GetValue(variable,"true");
      std::transform(modCladding.begin(), modCladding.end(), modCladding.begin(), ::tolower);
      cladding[i] = modCladding == "true" ? true : false;
      if (!cladding[i]) modCladdingThickness[i] = 0.;
      sprintf(variable,"mod%dNRadiators",i+1);
      modNRadiators[i] = config->GetValue(variable,12);
      modNAbsorbers[i] = modNRadiators[i] - 1;
      if (modNRadiators[i] == 0) modNAbsorbers[i] = 1; // the case where you are defining a solid absorber block with no active channels
    }
  }

  // geometric constants
  float zPitch;
  float xStartStrip; // middle of left most strip -- note this strip doesn't actually exist since the sets on the edge have 5 strips instead of 6
  float stripPitch; // distance between center of each rod ie) the diameter of one strip
  float zStartW; // position where first tungsten plate gets placed
  float zStartRad; // position where first radiator gap gets placed

  float modLengthZ = modAbsorberThickness[m_modNum]*modNAbsorbers[m_modNum] + modRadiatorGapLength[m_modNum]*modNRadiators[m_modNum];
  zPitch = modAbsorberThickness[m_modNum]+modRadiatorGapLength[m_modNum];
  zStartW = -1*modLengthZ/2 + modRadiatorGapLength[m_modNum] + modAbsorberThickness[m_modNum]/2;
  zStartRad = -1*modLengthZ/2 + (modRadiatorGapLength[m_modNum])/2;
  stripPitch = modCoreDiameter[m_modNum] + 2*modCladdingThickness[m_modNum];
  xStartStrip = -1*modNStripsPerGap[m_modNum]*stripPitch/2 + stripPitch/2;
  //  std::cout << " xStartStrip " << xStartStrip << std::endl;
  float modWidthX = modNStripsPerGap[m_modNum]*stripPitch;
  if (modWidthX == 0) modWidthX = modAbsorberWidth[m_modNum]; // the case where you are defining a solid absorber block with no active channels
  float modHeightY = modAbsorberHeight[m_modNum];
  float boxLengthZ = modLengthZ + modCasingThickness[m_modNum]*2;

  //  std::cout << " modcustom boxLengthZ " << boxLengthZ << std::endl;
  float boxWidthX = modWidthX + modCasingThickness[m_modNum]*2;
  float boxHeightY = modHeightY + modCasingThickness[m_modNum]*2;

  //----------------------------------------------
  // Housing
  //----------------------------------------------

  m_SteelBox         =  new G4Box("SteelCasing", boxWidthX*mm/2.0 ,boxHeightY*mm/2.0  , boxLengthZ*mm/2.0);
  m_ModuleBox        =  new G4Box("ModuleCasing", modWidthX*mm/2.0 ,modHeightY*mm/2.0  , modLengthZ*mm/2.0);
  m_SteelLogical        = new G4LogicalVolume(m_SteelBox           ,m_matHousing, "Steel_Logical");
  m_ModuleLogical       = new G4LogicalVolume(m_ModuleBox          ,m_matAir, "Module_Logical");

  char name[256];
  int cn = m_modNum;
  sprintf(name,"Steel_Box_Physical_i %d", cn);
  G4RotationMatrix* nullRotation = new G4RotationMatrix();
  G4ThreeVector  pos;
  pos = G4ThreeVector(0,0,0);
  m_SteelBoxPhysical = new G4PVPlacement(nullRotation,m_pos,m_SteelLogical,name,m_logicMother,false,cn,checkOverlaps);

  sprintf(name,"Module_Physical_i %d", cn);
  m_ModulePhysical = new G4PVPlacement(nullRotation,pos,m_ModuleLogical,name,m_SteelLogical,false,cn,checkOverlaps);

  G4VisAttributes* moduleColor = new G4VisAttributes( G4Colour::Gray() );
  m_ModuleLogical->SetVisAttributes( moduleColor );


  //----------------------------------------------
  // Quartz
  //----------------------------------------------

  m_StripTube 		= new G4Tubs( "Strip_Tube",
							0.0*mm,
							(modCoreDiameter[m_modNum]/2.0-0.005)*mm,
							modHeightY*mm/2.0 ,
							0.0*deg,
							360.0*deg);
  m_StripLogical 	= new G4LogicalVolume(m_StripTube
							,m_matQuartz,
							"Strip_Logical");

  if (cladding[m_modNum]) {

    m_CladdingTube  = new G4Tubs( "Cladding_Tube",
							modCoreDiameter[m_modNum]/2.0*mm - (0.005)*mm,
							(modCoreDiameter[m_modNum]/2.0*mm-0.0005*mm+modCladdingThickness[m_modNum])*mm, modHeightY*mm/2.0, 0.0*deg,
							360.0*deg);
    m_CladdingLogical = new G4LogicalVolume(m_CladdingTube, m_matQuartz, "Cladding_Logical");
  }

  G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Cyan() );
  quartzColor->SetForceSolid(true);
  m_StripLogical->SetVisAttributes( quartzColor );
  if (cladding[m_modNum]) m_CladdingLogical->SetVisAttributes( quartzColor );


  //----------------------------------------------
  // Plates
  //----------------------------------------------


  m_W                           = new G4Box ("W", modWidthX*mm/2.0 ,     modAbsorberHeight[m_modNum]*mm/2.0  , modAbsorberThickness[m_modNum]*mm/2.0);
  m_WLogical                    = new G4LogicalVolume(m_W          ,m_matAbsorber, "W_Logical");
  G4VisAttributes* airColor  = new G4VisAttributes( G4Colour::White() );
   airColor->SetForceSolid(true);
  m_ModuleLogical->SetVisAttributes( airColor );

  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );
   absorberColor->SetForceSolid(true);
  m_WLogical->SetVisAttributes( absorberColor );

  int modOffset = 100000; //  just a convention to organize copy numbers for different modules

  // Quartz strip dimensions: radiator gap, group, strip
  // populate all 8 modules with quartz strips
  G4RotationMatrix* stripRotation = new G4RotationMatrix();
  stripRotation->rotateX(90.*deg);
  cn = 0;
  for(int K=0;K<modNRadiators[m_modNum];K++) {
    for(int M=0;M<modNStripsPerGap[m_modNum];M++) {

	  sprintf(name,"Strip_a %d", cn);
      m_StripPhysical[K][M] = new G4PVPlacement(
									stripRotation,
									G4ThreeVector((xStartStrip+(M*stripPitch))*mm,0,(zStartRad+K*zPitch)*mm),
									m_StripLogical,
									name,
									m_ModuleLogical,
									false,
									cn,
									checkOverlaps);
      if (cladding[m_modNum]) {

	sprintf(name,"Cladding_a %d", cn);
	m_CladdingPhysical[K][M] = new G4PVPlacement(
									stripRotation,
									G4ThreeVector((xStartStrip+(M*stripPitch))*mm,0,(zStartRad+K*zPitch)*mm),m_CladdingLogical,
									name,
									m_ModuleLogical,
									false,
									cn,
									checkOverlaps);
      }
      ++cn;
    }
  }
  cn = 0;
  //  W plates (non pixel modules): Physical plates that span length of each absorber gap
  for(int K=0;K<modNAbsorbers[m_modNum];K++) {    // 11 layers of plates
    char volName[256];
    sprintf(volName,"W_i %d",K + m_modNum*modOffset);
    m_WPhysical[K] = new G4PVPlacement(nullRotation,G4ThreeVector(0,0,zStartW*mm+K*zPitch*mm),m_WLogical,name,m_ModuleLogical,false,cn,checkOverlaps);
    ++cn;
  }

  //----------------------------------------------
  // Define Surface/Border Properties
  //----------------------------------------------
  DefineBorderProperties();

  //----------------------------------------------
  // SD and Scoring Volumes
  //----------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Note one SD object for each module
  char quartzSDname[256];
  sprintf( quartzSDname, "QuartzSD%d", m_modNum+1);
  QuartzSD* aQuartzSD = new QuartzSD( quartzSDname, m_sd, m_modNum );
  aQuartzSD->HistInitialize();
  SDman->AddNewDetector( aQuartzSD );
  m_StripLogical->SetSensitiveDetector( aQuartzSD );



  /*
  if (m_simCherenkov) {
    //Note one SD object for each module
    char CherenkovSDname[256];
    sprintf( CherenkovSDname, "CherenkovSD%d", m_modNum);
    CherenkovSD* aCherenkovSD = new CherenkovSD( CherenkovSDname, m_sd, m_modNum);
    aCherenkovSD->HistInitialize();
    SDman->AddNewDetector( aCherenkovSD );
    //    m_AirScoringLogical->SetSensitiveDetector( aCherenkovSD );
  }
  */
  std::cout << "ModTypeZDC construction finished: SD name " << quartzSDname << std::endl;

}
