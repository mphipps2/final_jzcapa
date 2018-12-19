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
// For an explanation of the hierarchy scheme see: https://twiki.cern.ch/twiki/bin/view/Atlas/ZdcSimulation#Geometry_Implementation_Develope

#include "ModType3.hh"
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

ModType3::ModType3(const int cn,const G4ThreeVector& pos,
	      G4LogicalVolume* mother, SharedData* sd)
  : m_modNum( cn ),  m_pos( pos ), m_logicMother( mother ),
    m_sd( sd ),
    m_matQuartz(0),
    m_matAbsorber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModType3::ModType3()
  : m_modNum( 0 ), m_pos(G4ThreeVector()), m_logicMother(NULL),
    m_sd(NULL), 
    m_matQuartz(0), m_matAbsorber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModType3::~ModType3()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType3::Construct(){
  DefineMaterials();
  ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType3::DefineMaterials()
{
  // Get Config
  //  TEnv* config = m_sd->GetConfig();
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  //----------------------------------------------     
  // Define Materials
  //----------------------------------------------
  //  m_matHousing   = nist->FindOrBuildMaterial("G4_Al");
  m_matQuartz    = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  //  m_matAbsorber  = nist->FindOrBuildMaterial("G4_W");

  // Define materials not in NIST
  G4double fractionMass;
  G4int ncomponents;
  
  G4Element*  O   = new G4Element ("O"  ,"O"  ,  8.0 ,  16.0*g/mole);
  //  G4Element*  Si  = new G4Element ("Si" ,"Si" , 14.0 ,  28.085*g/mole);
  G4Element*  W   = new G4Element ("W"  ,"W"  , 74.0 ,  183.84*g/mole);
  G4Element*  Fe  = new G4Element ("Fe" ,"Fe" , 26.0 ,  55.845 *g/mole);
  G4Element*  C   = new G4Element ("C"  ,"C"  ,  6.0 ,  12.0107*g/mole);
  G4Element*  Ni  = new G4Element ("Ni" ,"Ni" , 28.0 ,  58.6934*g/mole);
  
  // Absorber composition:  savannah.cern.ch/task/download.php?file_id=22925
  m_matAbsorber = new G4Material("Tungsten",18.155*g/cm3,ncomponents=3);
  m_matAbsorber->AddElement(W,  fractionMass=0.948);
  m_matAbsorber->AddElement(Ni, fractionMass=0.037);
  m_matAbsorber->AddElement(Fe, fractionMass=0.015);

  m_matHousing  = new G4Material("Steel", 7.9*gram/cm3,ncomponents=2);
  m_matHousing->AddElement(Fe  , fractionMass=0.98);
  m_matHousing->AddElement(C   , fractionMass=0.02);

  m_matAir = new G4Material("Air",1.290*mg/cm3,ncomponents=2);
  m_matAir->AddElement(Ni, fractionMass=.70);
  m_matAir->AddElement(O, fractionMass=.30);
  //----------------------------------------------     
  // Define Material Properties
  //----------------------------------------------
  const G4int NUMENTRIES = 2;
  
  G4double ephoton         [NUMENTRIES] = {2.00*eV,4.80*eV};

  G4double rindexQuartz    [NUMENTRIES] = {1.46,1.46};
  G4double absorptionQuartz[NUMENTRIES] = {46*m,46*m};

  //Fill in the Marterial properties table for each material.
  //Guide for undestanding Optical processes at http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Photo
  G4MaterialPropertiesTable *quartzMPT = new G4MaterialPropertiesTable();
  quartzMPT->AddProperty("RINDEX",ephoton,rindexQuartz,NUMENTRIES);
  quartzMPT->AddProperty("ABSLENGTH",ephoton,absorptionQuartz,NUMENTRIES);   
  m_matQuartz->SetMaterialPropertiesTable(quartzMPT);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType3 :: DefineBorderProperties(){
  const G4int NUMENTRIES = 2;
  G4double ephoton[NUMENTRIES] = {2.00*eV,4.80*eV};

  //----------------------------------------------     
  // Housing Skin
  //----------------------------------------------
  G4double housingReflectivity      [NUMENTRIES] = {0.4, 0.4};
  G4double housingEfficiency        [NUMENTRIES] = {0.10, 0.10};

  G4MaterialPropertiesTable* housingMPT
    = new G4MaterialPropertiesTable();
  housingMPT->AddProperty("REFLECTIVITY", ephoton, housingReflectivity, NUMENTRIES);
  housingMPT->AddProperty("EFFICIENCY"  , ephoton, housingEfficiency  , NUMENTRIES);
  G4OpticalSurface* housingOS =
    new G4OpticalSurface("HousingOpSurface",unified, polished, dielectric_metal);
  housingOS->SetMaterialPropertiesTable( housingMPT );

  new G4LogicalSkinSurface("housingSkinSurface", m_SteelLogical, housingOS );
  
 
  //----------------------------------------------     
  // Absorber Skin
  //----------------------------------------------
  G4double absorberReflectivity      [NUMENTRIES] = {0.0, 0.0};
  G4MaterialPropertiesTable* absorberMPT
    = new G4MaterialPropertiesTable();
  absorberMPT->AddProperty("REFLECTIVITY", ephoton, absorberReflectivity, NUMENTRIES);
  G4OpticalSurface* absorberOS =
    new G4OpticalSurface("AbsorberOpSurface",glisur, polished, dielectric_metal);
  //      new G4OpticalSurface("AbsorberOpSurface",unified, polished, dielectric_metal);
  absorberOS->SetMaterialPropertiesTable( absorberMPT );

  new G4LogicalSkinSurface("absorberSkinSurface", m_WLogical, absorberOS );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType3::ConstructDetector()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();
  
  //----------------------------------------------     
  // Set Some Values
  //----------------------------------------------  
  //  G4double gapThicknessZ       = 2.;
  //  G4double modSizeX; G4double modSizeY; G4double modSizeZ;   G4double modAbsorberThickness; G4double modCoreDiameter; G4double modCoreIndexRefraction;
  G4String modAbsorberMat;

    // Option to switch on/off checking of volumes overlaps
  bool checkOverlaps = config->GetValue("checkOverlaps",true);
  
  // geometric constants
  const float zPitch = 1.20;
  const int nRadGaps = 12; // radiator gaps
  const int nAbsGaps = 11; // absorber gaps
  const int nStripSets = 9; // # of strip sets in each radiator gap
  const int nStrips = 6; // # of strips/set in each radiator gap  
  const float zStartAir  = -6.6  ; // (13.4/2) + 1.0 + 0.1
  const float xStartStrip = -4.5535; // middle of left most strip -- note this strip doesn't actually exist since the sets on the edge have 5 strips instead of 6  
  const float largePitchStrip = 1.044625; // distance between each set
  const float smallPitchStrip = .15; // distance between center of each in same set ie) the diameter of one strip
  const float zStartW = -6.0; // (13.4/2) + 0.5
  //  const float bufferSpace = 0.000001; // very small buffer space between strips so there are no floating point overlaps between volumes
  //----------------------------------------------     
  // variables ending in prime are
  // for construction of G4Para
  // the sizes there have to be adjusted for
  // depending on theta
  //----------------------------------------------
  
  //----------------------------------------------     
  // Housing
  //----------------------------------------------
  
  m_SteelBox         =  new G4Box("SteelCasing", 9.078*cm/2.0 ,19.0*cm/2.0  , 15.0*cm/2.0);
  //  m_AirScoringBox    =  new G4Box("AirScoring", 9.078*cm/2.0 ,2.0*cm/2.0, 15.0*cm/2.0);
  m_ModuleBox        =  new G4Box("ModuleCasing", 8.957*cm/2.0 ,18.0*cm/2.0  , 13.4*cm/2.0);
  
  m_SteelLogical        = new G4LogicalVolume(m_SteelBox           ,m_matHousing, "Steel_Logical");
  //  m_AirScoringLogical   = new G4LogicalVolume(m_AirScoringBox      ,m_matAir, "Air_Logical"); 
  m_ModuleLogical       = new G4LogicalVolume(m_ModuleBox          ,m_matAir, "Module_Logical");  

  char name[256];
  int cn = m_modNum;    
  sprintf(name,"Steel_Box_Physical_i %d", cn);
  G4RotationMatrix* nullRotation = new G4RotationMatrix();
  G4ThreeVector  pos;
  pos = G4ThreeVector(0,0,0);  
  m_SteelBoxPhysical = new G4PVPlacement(nullRotation,m_pos,m_SteelLogical,name,m_logicMother,false,cn,checkOverlaps);
  // This was added as a dummy scoring volume for run mode where cherenkov light is simulated directly
  sprintf(name,"Air_Scoring_Physical_i %d", cn);

  /*  
  pos = G4ThreeVector(0,9.0*cm,0);
  m_AirScoringPhysical  = new G4PVPlacement(nullRotation, pos, m_AirScoringLogical, name, m_SteelLogical, false, cn, checkOverlaps);
  */
  sprintf(name,"Module_Physical_i %d", cn);
  m_ModulePhysical = new G4PVPlacement(nullRotation,pos,m_ModuleLogical,name,m_SteelLogical,false,cn,checkOverlaps);	        
  //  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");  
  
  G4VisAttributes* moduleColor = new G4VisAttributes( G4Colour::Gray() );
  m_ModuleLogical->SetVisAttributes( moduleColor );

	
  //----------------------------------------------     
  // Quartz
  //----------------------------------------------

  //  m_StripTube                     = new G4Tubs( "Strip_Tube", 0.0*mm     , 1.5*mm/2.0  , 18.0*cm/2.0, 0.0*deg, 360.0*deg);
    m_StripTube                     = new G4Tubs( "Strip_Tube", 0.0*mm     , 1.4995*mm/2.0  , 18.0*cm/2.0, 0.0*deg, 360.0*deg);
  m_StripLogical                         = new G4LogicalVolume(m_StripTube          ,m_matQuartz, "Strip_Logical");

  G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Cyan() );
  quartzColor->SetForceSolid(true);  
  m_StripLogical->SetVisAttributes( quartzColor );


  //----------------------------------------------     
  // Plates
  //----------------------------------------------
  // Variables used as input to G4Para
  //  G4double absorberSizeX  = 89.57; // fix mod x
  //  G4double absorberSizeY  = moduleSizeY;
  //  G4double absorberSizeZ  = 10.0;
    

  m_W                           = new G4Box ("W", 8.957*cm/2.0 ,     18.0*cm/2.0  , 1.0*cm/2.0);
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
  for(int K=0;K<nRadGaps;K++) {    // 12 layers of strips
    for(int l=0;l<nStripSets;l++) {   // each layer has 9 sets of strips
      for(int M=0;M<nStrips;M++) { // each non-edge set has 6 strips while sets along edge have 5 (ie one compartment goes: 5,6,6,6,6,6,6,6,5)
	if ( (l == 8 && M == 5) || (l == 0 && M == 0)) continue; // edge strip sets only have 5 strips (not 6)
	//	cn = M + l*nStrips + K*nStrips*nStripSets + m_modNum*modOffset;
	sprintf(name,"Strip_a %d", cn);
	m_StripPhysical[K][l][M] = new G4PVPlacement(stripRotation,G4ThreeVector((xStartStrip+(l*largePitchStrip)+(M*smallPitchStrip))*cm,0,(zStartAir+K*zPitch)*cm),m_StripLogical,name,m_ModuleLogical,false,cn,checkOverlaps);
	++cn;
      }
    }
  }
  cn = 0;    
  //  W plates (non pixel modules): Physical plates that span length of each absorber gap
  for(int K=0;K<nAbsGaps;K++) {    // 11 layers of plates	
    char volName[256];
    sprintf(volName,"W_i %d",K + m_modNum*modOffset);
    m_WPhysical[K] = new G4PVPlacement(nullRotation,G4ThreeVector(0,0,zStartW*cm+K*zPitch*cm),m_WLogical,name,m_ModuleLogical,false,cn,checkOverlaps);
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
  if (config->GetValue("cherenkovYield",false)) {
    std::cout <<" adding CherenkovYield! " << std::endl;
    //Note one SD object for each module
    char CherenkovSDname[256];
    sprintf( CherenkovSDname, "CherenkovSD%d", m_modNum);
    CherenkovSD* aCherenkovSD = new CherenkovSD( CherenkovSDname, m_sd, m_modNum);
    aCherenkovSD->HistInitialize();
    SDman->AddNewDetector( aCherenkovSD );
    //    m_AirScoringLogical->SetSensitiveDetector( aCherenkovSD );    
  }
  */
  std::cout << "ModType3 construction finished: SD name " << quartzSDname << std::endl;

}
