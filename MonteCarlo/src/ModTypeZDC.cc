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
#include "FiberSD.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4MaterialTable.hh"

#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeZDC::ModTypeZDC(const int cn,
	      G4LogicalVolume* mother, G4ThreeVector* pos)
  : m_modNum( cn ), m_pos( pos ), m_fiberDiam (new G4ThreeVector(1.5,0.,0.)),
	m_absDim (new G4ThreeVector(90.,180.,11.)), m_logicMother( mother )
{
	m_materials = Materials::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	ModTypeZDC::ModTypeZDC(const int cn, ModTypeZDC* right)
		: m_modNum( cn )
	{
		m_nAbsorbers 			 = right->m_nAbsorbers;
		m_pos				 			 = new G4ThreeVector(*right->m_pos);
		m_fiberDiam 	 		 = new G4ThreeVector(*right->m_fiberDiam);
		m_absDim 		 			 = new G4ThreeVector(*right->m_absDim);
		m_HousingThickness = right->m_HousingThickness;
		m_GapThickness 		 = right->m_GapThickness;
		OPTICAL 					 = right->OPTICAL;
		CHECK_OVERLAPS 		 = right->CHECK_OVERLAPS;
		m_matAbsorber 		 = right->m_matAbsorber;
		m_matHousing			 = right->m_matHousing;
		m_logicMother			 = right->m_logicMother;
		m_materials				 = right->m_materials;
	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeZDC::ModTypeZDC(const int cn, G4LogicalVolume* mother)
  : m_modNum( cn ), m_pos(new G4ThreeVector(0.,0.,0.)), m_fiberDiam (new G4ThreeVector(1.5,0.,0.)),
	m_absDim (new G4ThreeVector(90.,180.,11.)), m_logicMother( mother )
{
	m_materials = Materials::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeZDC::ModTypeZDC()
  : m_modNum( 0 ), m_pos(new G4ThreeVector(0.,0.,0.)), m_fiberDiam (new G4ThreeVector(1.5,0.,0.)),
	m_absDim (new G4ThreeVector(90.,180.,11.)), m_logicMother(NULL)
{
	m_materials = Materials::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeZDC::~ModTypeZDC()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeZDC::Construct(){
  ConstructDetector();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeZDC::ConstructDetector()
{
	bool BUFFER, CLADDING;
	G4double fiberMaxDia;
	//FiberDimension x=Core, y=Cladding, z=Buffer diameters
	BUFFER   = ( m_fiberDiam->z() == 0.0 ) ? false : true;
	CLADDING = ( m_fiberDiam->y() == 0.0 ) ? false : true;
	if 		 (BUFFER)  fiberMaxDia = m_fiberDiam->z();
	else if(CLADDING)fiberMaxDia = m_fiberDiam->y();
	else						 fiberMaxDia = m_fiberDiam->x();

  // geometric constants
	float zPitch;
  float xStartStrip; // middle of left most strip -- note this strip doesn't actually exist since the sets on the edge have 5 strips instead of 6
  float stripPitch;  // distance between center of each rod ie) the diameter of one strip
  float zStartW; 		 // position where first tungsten plate gets placed
  float zStartRad; 	 // position where first radiator gap gets placed

  float modLengthZ = m_absDim->z()*m_nAbsorbers + m_GapThickness*m_nAbsorbers;
  zPitch = m_absDim->z()+ m_GapThickness;
  zStartW = -1*modLengthZ/2 + m_GapThickness + m_absDim->z()/2;
  zStartRad = -1*modLengthZ/2 + m_GapThickness/2;
  stripPitch = fiberMaxDia;
  xStartStrip = -1*fmod(m_absDim->x(),fiberMaxDia)*stripPitch/2 + stripPitch/2;
  float modWidthX = fmod(m_absDim->x(),fiberMaxDia)*stripPitch;
  if (modWidthX == 0) modWidthX = m_absDim->x(); // the case where you are defining a solid absorber block with no active channels
  float modHeightY = m_absDim->y();
	float boxWidthX  = modWidthX  + m_HousingThickness*2;
	float boxHeightY = modHeightY + m_HousingThickness*2;
  float boxLengthZ = modLengthZ + m_HousingThickness*2;


  //----------------------------------------------
  // Housing
  //----------------------------------------------

  m_HousingBox     = new G4Box("SteelCasing",  boxWidthX*mm/2.0 ,boxHeightY*mm/2.0  , boxLengthZ*mm/2.0);
  m_ModuleBox      = new G4Box("ModuleCasing", modWidthX*mm/2.0 ,modHeightY*mm/2.0  , modLengthZ*mm/2.0);
  m_HousingLogical = new G4LogicalVolume(m_HousingBox           ,m_matHousing, 		 "Housing_Logical");
  m_ModuleLogical  = new G4LogicalVolume(m_ModuleBox            ,m_materials->Air, "Module_Logical");

  char name[256];
  int cn = m_modNum;
  sprintf(name,"ZDC%d_Case_Physical", m_modNum);
  G4RotationMatrix* nullRotation = new G4RotationMatrix();
  G4ThreeVector  pos;
  pos = G4ThreeVector(0,0,0);
  m_HousingPhysical = new G4PVPlacement(nullRotation,*m_pos,m_HousingLogical,name,m_logicMother,false,cn,CHECK_OVERLAPS);

  sprintf(name,"ZDC%d_Air_Physical", cn);
  m_ModulePhysical = new G4PVPlacement(nullRotation,pos,m_ModuleLogical,name,m_HousingLogical,false,cn,CHECK_OVERLAPS);

  G4VisAttributes* moduleColor = new G4VisAttributes( G4Colour::Gray() );
  m_ModuleLogical->SetVisAttributes( moduleColor );


  //----------------------------------------------
  // Quartz
  //----------------------------------------------

  m_FiberCoreTube =
		new G4Tubs( "Fiber_Core_Tube",
							  0.0*mm,
								(m_fiberDiam->x()/2.0-0.005)*mm,
								modHeightY*mm/2.0 ,
								0.0*deg,
								360.0*deg);
  m_FiberCoreLogical =
		new G4LogicalVolume(m_FiberCoreTube,
												m_materials->pQuartz,
												"FiberCore_Logical");

	G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Cyan() );
  quartzColor->SetForceSolid(true);
  m_FiberCoreLogical->SetVisAttributes( quartzColor );

  if ( CLADDING ) {

    m_CladdingTube  =
			new G4Tubs( "Fiber_Cladding_Tube",
									m_fiberDiam->x()/2.0*mm - (0.005)*mm,
									(m_fiberDiam->x()/2.0*mm-0.0005*mm+m_fiberDiam->y())*mm, modHeightY*mm/2.0, 0.0*deg,
									360.0*deg);
    m_CladdingLogical =
			new G4LogicalVolume(m_CladdingTube,
													m_materials->pQuartz,
													"Fiber_Cladding_Logical");
		m_CladdingLogical->SetVisAttributes( quartzColor );
  }

  //----------------------------------------------
  // Plates
  //----------------------------------------------

  m_W = new G4Box("W",
									modWidthX*mm/2.0 ,
									m_absDim->y()/2.0,
									m_absDim->y()/2.0);
  m_WLogical =
		new G4LogicalVolume(m_W,
												m_matAbsorber,
												"W_Logical");

  G4VisAttributes* airColor = new G4VisAttributes( G4Colour::White() );
  airColor->SetForceSolid(true);
  m_ModuleLogical->SetVisAttributes( airColor );

  G4VisAttributes* absorberColor = new G4VisAttributes( G4Colour::Red() );
  absorberColor->SetForceSolid(true);
  m_WLogical->SetVisAttributes( absorberColor );

  // Quartz strip dimensions: radiator gap, group, strip
  // populate all 8 modules with quartz strips
  G4RotationMatrix* stripRotation = new G4RotationMatrix();
  stripRotation->rotateX(90.*deg);
  cn = 0;
	m_FiberCorePhysical.resize(m_nAbsorbers + 1);
	m_FiberCladdingPhysical.resize(m_nAbsorbers + 1);
  for(int K = 0; K < m_nAbsorbers + 1; K++ ){
    for(int M = 0; M < fmod(m_absDim->x(),fiberMaxDia); M++ ){

	  	sprintf(name,"ZDC%d_Core%d", m_modNum, cn);
      m_FiberCorePhysical.at(K).push_back( new G4PVPlacement(
															  stripRotation,
																G4ThreeVector((xStartStrip+(M*stripPitch))*mm,0,(zStartRad+K*zPitch)*mm),
																m_FiberCoreLogical,
																name,
																m_ModuleLogical,
																false,
																cn,
																CHECK_OVERLAPS) );
      if ( CLADDING ) {

				sprintf(name,"ZDC%d_Cladding%d", m_modNum, cn);
				m_FiberCladdingPhysical.at(K).push_back(new G4PVPlacement(
																					stripRotation,
																					G4ThreeVector((xStartStrip+(M*stripPitch))*mm,0,(zStartRad+K*zPitch)*mm),
																					m_CladdingLogical,
																					name,
																					m_ModuleLogical,
																					false,
																					cn,
																					CHECK_OVERLAPS) );
      }
      ++cn;
    }
  }
  cn = 0;
  //  W plates (non pixel modules): Physical plates that span length of each absorber gap
  for(int K=0;K<m_nAbsorbers;K++) {    // 11 layers of plates
    char volName[256];
    sprintf(volName,"ZDC%d_Absorber%d",m_modNum, K);
    m_WPhysical.push_back( new G4PVPlacement(nullRotation,
													 G4ThreeVector(0,0,zStartW*mm+K*zPitch*mm),
													 m_WLogical,name,
													 m_ModuleLogical,
													 false,
													 cn,
													 CHECK_OVERLAPS) );
    ++cn;
  }

  //----------------------------------------------
  // SD and Scoring Volumes
  //----------------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Note one SD object for each module
  char fiberSDname[256];
  sprintf( fiberSDname, "ZDC%d_SD", m_modNum);
  FiberSD* aFiberSD = new FiberSD( fiberSDname, m_modNum, OPTICAL );
  aFiberSD->HistInitialize();
  SDman->AddNewDetector( aFiberSD );
  m_FiberCoreLogical->SetSensitiveDetector( aFiberSD );


  std::cout << "ModTypeZDC construction finished: SD name " << fiberSDname << std::endl;

}

void ModTypeZDC::SetHousingMaterial(G4String material)
{
	material.toLower();
			 if( material == "steel"   ) m_matHousing = m_materials->Steel;
	else if( material == "aluminum") m_matHousing = m_materials->Al;
	else G4cout << "Invalid housing material selection, defaulting to steel" << G4endl;
}

void ModTypeZDC::SetAbsorberMaterial(G4String material)
{
	material.toLower();
			 if( material == "pure"  		 ) m_matHousing = m_materials->pureW;
	else if( material == "composite" ) m_matHousing = m_materials->NiW;
	else G4cout << "Invalid absorber material selection, defaulting to composite" << G4endl;
}
