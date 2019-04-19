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

#include "ModType2.hh"
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

ModType2::ModType2(const int cn,const G4ThreeVector& pos,
	      G4LogicalVolume* mother, SharedData* sd)
  : m_modNum( cn ),  m_pos( pos ), m_logicMother( mother ),
    m_sd( sd ),
    m_matQuartz(0),
    m_matAbsorber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModType2::ModType2()
  : m_modNum( 0 ), m_pos(G4ThreeVector()), m_logicMother(NULL),
    m_sd(NULL), 
    m_matQuartz(0), m_matAbsorber(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModType2::~ModType2()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType2::Construct(){
  DefineMaterials();
  ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType2::DefineMaterials()
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

void ModType2 :: DefineBorderProperties(){
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

  new G4LogicalSkinSurface("absorberSkinSurface", m_WLatEdgeLogical, absorberOS );
  new G4LogicalSkinSurface("absorberSkinSurface", m_WLatCenterLogical, absorberOS );
  new G4LogicalSkinSurface("absorberSkinSurface", m_WLongLogical, absorberOS );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModType2::ConstructDetector()
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
  const G4double pitch = 1.044625 * cm; // pitch in both x and y
  const G4double zPitch = 1.20 * cm;
  //  const int nSides = 2; // # of arms
  //  const int nMods = 4; // modules/side
  const int nRadGaps = 12; // radiator gaps
  //  const int nRadGaps = 2; // radiator gaps  
  const int nAbsGaps = 11; // absorber gaps
   //  const int nAbsGaps = 1; // absorber gaps    
  const int nLatPlates = 9; // nHolesX (machined holes in tungsten) + 1
  const int nHolesX = 8; // nHolesX (machined holes in tungsten)
  const int nHolesY = 16; // number of machined holes in Y
  //const int nHolesY = 2; // number of machined holes in Y  
  const int nStripSets = 9; // # of strip sets in each radiator gap
  //  const int nWGroups = 10; // # of strip sets + 1
  const int nStrips = 6; // # of strips/set in each radiator gap
  const G4double xStartPix = -3.6561875 * cm; // middle of first pixel
  const G4double yStartPix = -7.725 * cm;
  const G4double zStartWLat = -6.0 * cm  ; // (13.4/2) + 0.5
  //  const float zStartPixW = -6.0  ; // (13.4/2) + 0.5
  const G4double zStartAir  = -6.6 * cm ; // (13.4/2) + 1.0 + 0.1
  const G4double xStartWLatEdge = -4.1035 * cm;
  const G4double xStartWLatCenter = -4.1785 * cm;
  const G4double xStartStrip = -4.5535 * cm; // middle of left most strip -- note this strip doesn't actually exist since the sets on the edge have 5 strips instead of 6  
  const G4double largePitchStrip = 1.044625 * cm; // distance between each set
  const G4double smallPitchStrip = .15 * cm; // distance between center of each in same set ie) the diameter of one strip
  
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
  m_AirScoringBox    =  new G4Box("AirScoring", 9.078*cm/2.0 ,2.0*cm/2.0, 15.0*cm/2.0);
  m_ModuleBox        =  new G4Box("ModuleCasing", 8.957*cm/2.0 ,18.0*cm/2.0  , 13.4*cm/2.0);
  
  m_SteelLogical        = new G4LogicalVolume(m_SteelBox           ,m_matHousing, "Steel_Logical");
  m_AirScoringLogical   = new G4LogicalVolume(m_AirScoringBox      ,m_matAir, "Air_Logical"); 
  m_ModuleLogical       = new G4LogicalVolume(m_ModuleBox          ,m_matAir, "Module_Logical");  

  //  G4double moduleSizeX = 92.0;
  //  G4double moduleSizeY = 180.;
  //  G4double moduleSizeZ = 140.;
  //  G4double housingThickness   = 2.;    
  char name[256];
  int cn = m_modNum;    
  sprintf(name,"Steel_Box_Physical_i %d", cn);
  G4RotationMatrix* nullRotation = new G4RotationMatrix();
  m_SteelBoxPhysical = new G4PVPlacement(nullRotation,m_pos,m_SteelLogical,name,m_logicMother,false,cn,checkOverlaps);
  // This was added as a dummy scoring volume for run mode where cherenkov light is simulated directly
  sprintf(name,"Air_Scoring_Physical_i %d", cn);
  G4ThreeVector  pos;
  pos = G4ThreeVector(0,9.0*cm,0);
  //  m_AirScoringPhysical  = new G4PVPlacement(nullRotation, pos, m_AirScoringLogical, name, m_SteelLogical, false, cn, checkOverlaps);
  sprintf(name,"Module_Physical_i %d", cn);
  pos = G4ThreeVector(0,0,0);
  m_ModulePhysical = new G4PVPlacement(nullRotation,pos,m_ModuleLogical,name,m_SteelLogical,false,cn,checkOverlaps);	        
  //  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");  
  
  G4VisAttributes* moduleColor = new G4VisAttributes( G4Colour::Gray() );
  m_ModuleLogical->SetVisAttributes( moduleColor );

	
  //----------------------------------------------     
  // Quartz
  //----------------------------------------------

  m_PixelTubeW                   = new G4Tubs( "Pixel_Tube_W", 0.0*mm     , 1.0*mm/2.0  , 1.0*cm/2.0, 0.0*deg, 360.0*deg);
  m_PixelTubeRad                 = new G4Tubs( "Pixel_Tube_Rad", 0.0*mm     , 1.0*mm/2.0  , 0.2*cm/2.0, 0.0*deg, 360.0*deg);
  m_AirTubeW                     = new G4Tubs( "Air_Tube_W", 0.0*mm     , 1.0*mm/2.0  , 1.0*cm/2.0, 0.0*deg, 360.0*deg);
  m_AirTubeRad                   = new G4Tubs( "Air_Tube_Rad", 0.0*mm     , 1.0*mm/2.0  , 0.2*cm/2.0, 0.0*deg, 360.0*deg);
  m_StripTube                     = new G4Tubs( "Strip_Tube", 0.0*mm     , 1.5*mm/2.0  , 18.0*cm/2.0, 0.0*deg, 360.0*deg);
  m_StripLogical                         = new G4LogicalVolume(m_StripTube          ,m_matQuartz, "Strip_Logical");
  m_AirGapLogical                       = new G4LogicalVolume(m_AirGap             ,m_matAir, "Air_Gap_Logical");
  m_PixelWLogical                       = new G4LogicalVolume(m_PixelTubeW        ,m_matQuartz, "Pixel_W_Logical");
  m_PixelRadLogical                     = new G4LogicalVolume(m_PixelTubeRad      ,m_matQuartz, "Pixel_Rad_Logical");
  m_PixelAirWLogical                   = new G4LogicalVolume(m_AirTubeW          ,m_matAir, "Pixel_Air_W_Logical");
  m_PixelAirRadLogical                 = new G4LogicalVolume(m_AirTubeRad        ,m_matAir, "Pixel_Air_Rad_Logical");


  G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Cyan() );
  quartzColor->SetForceSolid(true);  
  m_StripLogical->SetVisAttributes( quartzColor );
  m_PixelWLogical->SetVisAttributes( quartzColor );
  m_PixelRadLogical->SetVisAttributes( quartzColor );



  //----------------------------------------------     
  // Plates
  //----------------------------------------------
  // Variables used as input to G4Para
  //  G4double absorberSizeX  = 89.57; // fix mod x
  //  G4double absorberSizeY  = moduleSizeY;
  //  G4double absorberSizeZ  = 10.0;
    
  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );
  absorberColor->SetForceSolid(true);
  
  // first 3 W volumes used for pixel modules
  m_WLatCenter                   = new G4Box ("W_Lat_Center", 0.9*cm/2.0 ,     18.0*cm/2.0  , 1.0*cm/2.0);
  m_WLatEdge                     = new G4Box ("W_Lat_Edge", 0.75*cm/2.0 ,    18.0*cm/2.0  , 1.0*cm/2.0);
  m_WLong                         = new G4Box ("W_Long", 0.144625*cm/2.0, 18.0*cm/2.0  , 13.4*cm/2.0);
  m_AirGap                        = new G4Box ("Air_Gap", 0.144625*cm/2.0 ,18.0*cm/2.0  , 0.2*cm/2.0);
  m_WLatCenterLogical                  = new G4LogicalVolume(m_WLatCenter        ,m_matAbsorber, "W_Lat_Center_Logical");
  m_WLatEdgeLogical                    = new G4LogicalVolume(m_WLatEdge          ,m_matAbsorber, "W_Lat_Edge_Logical" );
  m_WLongLogical                        = new G4LogicalVolume(m_WLong              ,m_matAbsorber, "W_Long_Logical");
  m_AirGapLogical                 = new G4LogicalVolume(m_AirGap, m_matAir, "Air_Gap_Logical");
      
  m_WLatCenterLogical->SetVisAttributes( absorberColor );
  m_WLatEdgeLogical->SetVisAttributes( absorberColor );
  m_WLongLogical->SetVisAttributes( absorberColor );  

  G4VisAttributes* airColor2  = new G4VisAttributes( G4Colour::White() );
  airColor2->SetForceSolid(true);
  m_PixelAirWLogical->SetVisAttributes( airColor2 );
  m_PixelAirRadLogical->SetVisAttributes( airColor2 );

  G4VisAttributes* airColor  = new G4VisAttributes( G4Colour::Blue() );
  airColor->SetForceSolid(true);
  m_AirGapLogical->SetVisAttributes( airColor );
  
  __attribute__((unused)) int modOffset = 100000; //  just a convention to organize copy numbers for different modules
  cn = 0;
  // Create Pixel and air cylinder that fill radiator gap in pixel modules
  for(int K=0;K<nRadGaps;K++) {    // 12 radiator gaps; 11 absorber gaps
    for(int l=0;l<nHolesY;l++) {    // 16 pixel holes/Longitudinal W plates in Y dimension
      for (int M=0; M<nHolesX;M++) {  // 8 pixels/Longitudinal W plates in X dimension
	if (l==0) {
	  //	  cn = M + nHolesX*K + modOffset*m_modNum;
	  sprintf(name,"Air_Gap_i %d", cn);
	  G4double zSpacing = zStartAir+K*zPitch;	  
	  m_AirPhysical[K][M] = new G4PVPlacement(nullRotation, G4ThreeVector(0,0,zSpacing ), m_AirGapLogical, name, m_WLongLogical, false, cn, checkOverlaps);
	}
	if ( K==0 && l==0) {
	  //	  cn = M + nHolesX*K + modOffset*m_modNum;
	  sprintf(name,"W_Long_i %d", cn);
	  G4double xSpacing = xStartPix+M*pitch;	  
	  m_WLongPhysical[M] = new G4PVPlacement(nullRotation,G4ThreeVector(xSpacing,0,0), m_WLongLogical, name, m_ModuleLogical, false, cn, checkOverlaps);
	}
	// for Hadronic Module (ie Type2)
	// Note the convention for volume names. "_i" == inactive volumes, "_a" == active volumes, "_d" == dead volumes. Beyond that, each volume should be labeled with a component name andended with a volume number. The volume number is set s.t. the 10000 digit indicates side, the 1000 digit denotes module. This convention is used later in the SD classes for ID purposes.
	//	std::cout << " K " << K << " l " << l << " M " << M << std::endl;
	if (l < 10) {
	  //	  cn = M + l*nHolesX + K*nHolesY*nHolesX + m_modNum*modOffset;
	  sprintf(name,"Pixel_Rad_i %d", cn);
	  G4double ySpacing = yStartPix+l*pitch;	  
	  m_PixelRadPhysical[K][l][M] = new G4PVPlacement(nullRotation,G4ThreeVector(0,ySpacing,0),m_PixelRadLogical,name,m_AirGapLogical,false,cn,checkOverlaps);	    
	}
	else {
	  //	  cn = M + l*nHolesX + K*nHolesY*nHolesX + m_modNum*modOffset;
	  sprintf(name,"Air_Rad_i %d", cn);
	  G4double ySpacing = yStartPix+l*pitch;
	  m_PixelRadPhysical[K][l][M] = new G4PVPlacement(nullRotation,G4ThreeVector(0,ySpacing,0),m_PixelAirRadLogical,name,m_AirGapLogical,false,cn,checkOverlaps);	    
	}
	++cn;
      }
    }
  }
  
  cn = 0;
  // Create Pixel and air cylinder and air gaps that fill absorber gap in pixel modules
  for(int K=0;K<nAbsGaps;K++) {    // 12 radiator gaps; 11 absorber gaps
    for(int l=0;l<nHolesY;l++) {    // 16 pixel holes/Longitudinal W plates in Y dimension
      for (int M=0; M<nHolesX;M++) {  // 8 machined holes in x
	// for EM Module
	if (l > 0 && l < 9) {
	  //	  cn = M + l*nHolesX + K*nHolesY*nHolesX + m_modNum*modOffset;
	  sprintf(name,"Pixel_W_i %d", cn);
	  G4double ySpacing = yStartPix+l*pitch;
	  G4double zSpacing = zStartWLat+K*zPitch;	  
	  m_PixelWPhysical[K][l][M] = new G4PVPlacement(nullRotation,G4ThreeVector(0,ySpacing,zSpacing),m_PixelWLogical,name,m_WLongLogical,false,cn,checkOverlaps);	    
	} 
	else {
	  //	  cn = M + l*nHolesX + K*nHolesY*nHolesX + m_modNum*modOffset;
	  sprintf(name,"Air_W_i %d", cn);
	  G4double ySpacing = yStartPix+l*pitch;
	  G4double zSpacing = zStartWLat+K*zPitch ;	  
	  m_PixelWPhysical[K][l][M] = new G4PVPlacement(nullRotation,G4ThreeVector(0,ySpacing,zSpacing),m_PixelAirWLogical,name,m_WLongLogical,false,cn,checkOverlaps);	    
	}
        ++cn;
      }
    }
  }

  cn =0;
  // Quartz strip dimensions: side, module, radiator gap, group, strip
  // populate all 8 modules with quartz strips
  G4RotationMatrix* stripRotation = new G4RotationMatrix();
  stripRotation->rotateX(90.*deg);
  for(int K=0;K<nRadGaps;K++) {    // 12 layers of strips
    for(int l=0;l<nStripSets;l++) {   // each layer has 9 sets of strips
      for(int M=0;M<nStrips;M++) { // each non-edge set has 6 strips while sets along edge have 5 (ie one compartment goes: 5,6,6,6,6,6,6,6,5)
	if ( (l == 8 && M == 5) || (l == 0 && M == 0)) continue; // edge strip sets only have 5 strips (not 6)
	//	cn = M + l*nStrips + K*nStrips*nStripSets + m_modNum*modOffset;
	sprintf(name,"Strip_a %d", cn);
	G4double xSpacing = xStartStrip+(l*largePitchStrip)+(M*smallPitchStrip);
	G4double zSpacing = zStartAir+K*zPitch;	
	m_StripPhysical[K][l][M] = new G4PVPlacement(stripRotation,G4ThreeVector(xSpacing,0,zSpacing),m_StripLogical,name,m_ModuleLogical,false,cn,checkOverlaps);
	++cn;
      }
    }
  }
    

  // Here we populate physical volumes in pixel modules (EM_C, H0_C and H0_A) ... note: EM_A doesn't have pixels --- physical volumes include longitudinal W plates, air gaps and pixel segments in both absorber and radiator
  //  Note: Longitudinal W plates (only in pixel modules) used as tool to carve away pixels through inheritance, rather than boolean subtraction (which is much more computationally taxing)
  cn = 0;
  for(int K=0;K<nAbsGaps;K++) {    // 12 radiator gaps; 11 absorber gaps
    for(int l=0;l<nHolesY;l++) {    // 8 pixels/Longitudinal W plates in X dimension; 9 lateral W plates
      for (int M=0; M<nLatPlates;M++) {  // 9: ie # of pixels in x direction
	// for lateral plates
	if (l==0) {
	  //	  cn = M + K*nLatPlates + m_modNum*modOffset;
	  sprintf(name,"W_Lat_i %d", cn);
	  if (M == 0 || M == 8) {
	    G4double xSpacing = xStartWLatEdge;
	    G4double zSpacing = zStartWLat+K*zPitch;
	    if (M == 0) m_WLatPhysical[K][M] = new G4PVPlacement(nullRotation,G4ThreeVector(xSpacing,0,zSpacing),m_WLatEdgeLogical,name,m_ModuleLogical,false,cn,checkOverlaps);
	    else m_WLatPhysical[K][M] = new G4PVPlacement(nullRotation,G4ThreeVector(-1*xSpacing,0,zSpacing),m_WLatEdgeLogical,name,m_ModuleLogical,false,cn,checkOverlaps);	    
	  }  
	  else {
	    G4double xSpacing = xStartWLatCenter+M*pitch;
	    G4double zSpacing = zStartWLat+K*zPitch;
	    m_WLatPhysical[K][M] = new G4PVPlacement(nullRotation,G4ThreeVector(xSpacing,0,zSpacing),m_WLatCenterLogical,name,m_ModuleLogical,false,cn,checkOverlaps);	    
	  }
	  ++cn;
	}
      }
    }
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
    //Note one SD object for each module
    char CherenkovSDname[256];
    sprintf( CherenkovSDname, "CherenkovSD%d", m_modNum);
    CherenkovSD* aCherenkovSD = new CherenkovSD( CherenkovSDname, m_sd, m_modNum);
    aCherenkovSD->HistInitialize();
    SDman->AddNewDetector( aCherenkovSD );
    m_AirScoringLogical->SetSensitiveDetector( aCherenkovSD );    
  }
  */
  std::cout << "ModType2 construction finished: SD name " << quartzSDname << std::endl;

}
