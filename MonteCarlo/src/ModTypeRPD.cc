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

#include "ModTypeRPD.hh"
#include "FiberSD.hh"


#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4UserLimits.hh"


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

ModTypeRPD::ModTypeRPD(const int cn,const G4ThreeVector& pos,
	      G4LogicalVolume* mother)
  : m_modNum( cn ),  m_pos( pos ), m_logicMother( mother ),
    m_matQuartz(0)
{
	materials = Materials::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD(const int cn, ModTypeRPD* right)
  : m_modNum( cn )
{
	m_modNum 					 = right->m_modNum;
	m_pos 						 = new G4ThreeVector(right->m_pos);
	m_fiberDiam 			 = new G4ThreeVector(right->m_fiberDiam);
	m_HousingThickness = right->m_HousingThickness;
	m_fiberPitch 			 = right->m_fiberPitch;
	m_tileSize 				 = right->m_tileSize;
	m_logicMother 		 = right->m_logicMother;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD()
  : m_modNum( 0 ), m_pos(G4ThreeVector()), m_logicMother(NULL),
    m_matQuartz(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::~ModTypeRPD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::Construct(){
  DefineMaterials();
  ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::DefineMaterials()
{
	materials->UseOpticalMaterials(true);
	materials->DefineOpticalProperties();

  //----------------------------------------------
  // Define Materials
  //----------------------------------------------

	//Quartz
	m_matQuartz = materials->pQuartz;

	// m_quartzOS = materials->quartzOS;
	// m_opSurface = materials->opSurface;

	//Aluminum
	m_Al = materials->Al;

	//Air
	m_Air = materials->Air;

	//Polyethylene/Clad
	m_Poly = materials->Polyethylene;

	//Optical Grease
	m_Grease = materials->Grease;

	//Wavelength Shifter/Fiber
	m_PMMA = materials->PMMA;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::ConstructDetector()
{
	bool pan_flute = false;
	bool test_tile = false;

	if(pan_flute || test_tile){
		for(int i =0;i<8;i++){
		rpd_comp[i]=false;
		}
	}
	else for(int i =0;i<8;i++){
		rpd_comp[i]=true;
		}

	// RPD components
	bool tile_flag = rpd_comp[0];				// quartz tiles
	bool fiber_flag = rpd_comp[1];			// readout fibers
	bool cladding_flag = rpd_comp[2];		// fiber cladding
	bool grease_flag = rpd_comp[3];			// optical grease
	bool foil_flag = rpd_comp[4];				// aluminum foil
	bool fr_bck_flag = rpd_comp[5]; 		// front and back aluminum casing
	bool vertical_flag = rpd_comp[6]; 	// veritcal aluminum dividers
	bool air_detec_flag = rpd_comp[7];	// photo detectors

  // Get Config
  TEnv* config = m_sd->GetConfig();

	//retrieve RPD parameters
  float tileX 		= config->GetValue("tileXsize",20);
  float tileY 		= config->GetValue("tileYsize",20);
  float tileZ 		= config->GetValue("tileZsize",10);
  float halfX_gap 	= ((config->GetValue("modGapX",1.58))/2)*mm;
	float fiber_diam 	= config->GetValue("mod6CoreDiameter",1);
	float foil_thickness 	= config->GetValue("foil_thickness",0.016);


	//manually create some RPD parameters
	float halfY_gap 	= (2.5*foil_thickness/2)*mm;
	float hole_center_offset = 0.05;
	float case_thickness = 0.9*2.0*halfX_gap;
	float core_diam;
	float grease_offset;
	if(grease_flag) grease_offset = 0.05;
	else grease_offset = 0.0;

	if(cladding_flag) core_diam = 0.970*(fiber_diam)*mm;
	else 							core_diam = (fiber_diam)*mm;

	float fiberHeightY[4];
	float foilHeightY[4];

  char name[256];

	//create some rotation matrices
	G4RotationMatrix* stripRotation = new G4RotationMatrix();
	stripRotation->rotateX(90.*deg);
	G4RotationMatrix* nullRotation = new G4RotationMatrix();

  // Option to switch on/off checking of volumes overlaps
  bool checkOverlaps = config->GetValue("checkOverlaps",false);

	//create tile that has 4 holes to be used repeatedly later
  m_tile_no_fiber_hole = new G4Box("tile_no_hole",(tileX/2)*mm, (tileY/2)*mm, (tileZ/2)*mm);
	m_fiber_subtract = new G4Tubs( "fiber_sub",
						0.0*mm,
						((fiber_diam/2.0)+ grease_offset )*mm, //1mm diameter
						tileY*mm ,
						0.0*deg,
						360.0*deg);

	m_tile = new G4SubtractionSolid("tile",m_tile_no_fiber_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 1*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset)*mm));
  m_tile = new G4SubtractionSolid("tile",m_tile,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 2*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
	m_tile = new G4SubtractionSolid("tile",m_tile,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 3*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
	m_tile = new G4SubtractionSolid("tile",m_tile,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 4*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));


  m_tileLogical	 = new G4LogicalVolume(m_tile, m_matQuartz, "tile_Logical");

  G4VisAttributes* quartzColor  = new G4VisAttributes(  G4Colour::Cyan());//
	//quartzColor->SetForceSolid(true);
	m_tileLogical->SetVisAttributes( quartzColor );

	int cn = 0, cn_fiber = 0;

	float RPD_centerX = m_pos.getX();
	float RPD_centerY = m_pos.getY();
	float RPD_centerZ = m_pos.getZ();

	float RPD_startX = RPD_centerX + 3*halfX_gap + 1.5*tileX;
	float RPD_startY = RPD_centerY + 3*halfY_gap + 1.5*tileY;


//positioning variables
///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float factor1 = 1.5;
float hole_shift		= grease_offset + hole_center_offset;

float tileZcenter[4];
for(int i=0; i < 4; i++){
	tileZcenter[i] = RPD_centerZ + (i * (factor1*fiber_diam + foil_thickness) );
}
float fiberZcenter[4];
for(int i=0; i < 4; i++){
	fiberZcenter[i] = tileZcenter[i] + tileZ/2.0 - fiber_diam/2.0 - hole_shift;
}

float foil_gap = (fiberZcenter[1]-(fiber_diam/2.0)) - (fiberZcenter[0] + (fiber_diam/2.0));

float foilZcenter[4];
for(int i=0; i < 4; i++){
	foilZcenter[i] = fiberZcenter[i] + fiber_diam/2.0 + foil_gap/2.0;
}

float final_foil_pos = foilZcenter[3]+foil_thickness/2.0;
float assembly_midZ = RPD_centerZ + (final_foil_pos - (RPD_centerZ + (tileZ/2.0)))/2.0;
float half_Z_length = final_foil_pos-assembly_midZ;


int fiber_cnt=0;

//create fibers/cladding/grease with correct lengths
for(int k=0;k<4;k++) {
	fiberHeightY[k]=((k+1)*tileY)+(k*2*halfY_gap);

	sprintf(name,"fiber_a %d", k);
	m_fiber[k] 		= new G4Tubs( name,
						0.0*mm,
						(core_diam/2.0)*mm, //1mm diameter
						fiberHeightY[k]*mm/2.0 ,
						0.0*deg,
						360.0*deg);
if(grease_flag){
	sprintf(name,"fibergrease_a %d", k);
  m_fibergrease[k] 		= new G4Tubs( name,
						(fiber_diam/2.0)*mm,
						((fiber_diam+grease_offset)/2.0)*mm,
						fiberHeightY[k]*mm/2.0 ,
						0.0*deg,
						360.0*deg);
					}

if(cladding_flag){
	sprintf(name,"fiberclad_a %d", k);
	m_fiberclad[k] 		= new G4Tubs( name,
							(core_diam/2.0)*mm,
							((fiber_diam)/2.0)*mm,
							fiberHeightY[k]*mm/2.0 ,
							0.0*deg,
							360.0*deg);
						}

if(grease_flag){
	sprintf(name,"fibergreaseLogical_a %d", k);

	m_fibergreaseLogical[k] 	= new G4LogicalVolume(m_fibergrease[k]
						,m_Grease,//m_Poly
						name);
	m_fibergreaseLogical[k]->SetVisAttributes( G4Colour(1,1,1,0.3) );//G4VisAttributes::Invisible
	}

if(cladding_flag){
	sprintf(name,"fibercladLogical_a %d", k);

	m_fibercladLogical[k] 	= new G4LogicalVolume(m_fiberclad[k]
						,m_Poly,//m_Poly
						name);
	m_fibercladLogical[k]->SetVisAttributes( G4Colour(1,0,0,0.2));//G4VisAttributes::Invisible
	}
}

for(int i=0;i<4;i++) {
	for(int j=0;j<4;j++) {
		for(int k=0;k<4;k++) {
				sprintf(name,"fiberLogical_%d_%d_%d",k,j,i);

						m_fiberLogical[fiber_cnt] 	= new G4LogicalVolume(m_fiber[k]
											,m_PMMA,
											name);
						m_fiberLogical[fiber_cnt]->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.2) );

						m_fiberLogical[fiber_cnt]->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

							fiber_cnt++;
					}
				}
			}

// create alum containment

//create foil sections
for(int k=0;k<4;k++) {
foilHeightY[k]=(2*(k+1)-2) * halfY_gap + (k+1)*tileY;
sprintf(name,"foil_V %d", k);

m_foilV[k] 		= new G4Box(name,(tileX/2.0)*mm, (foilHeightY[k]/2.0)*mm, (0.75*foil_thickness/2.0)*mm);//

sprintf(name,"foilVlog %d", k);

m_foilVLogical[k] 	= new G4LogicalVolume(m_foilV[k]
					,m_Al,
					name);
m_foilVLogical[k]->SetVisAttributes(G4Colour(0.9,0.0,0.0,0.4));//G4VisAttributes::Invisible
}


m_foilVfront 		= new G4Box("m_foilVfront",(tileX/2.0)*mm, (tileY/2.0)*mm, (0.75*foil_thickness/2.0)*mm);//
m_foilVfrontLogical 	= new G4LogicalVolume(m_foilVfront
				,m_Al,
				"foilVfrntlog");
m_foilVfrontLogical->SetVisAttributes( G4Colour(0.9,0.0,0.0,0.4) );//


m_foilH 			= new G4Box("foil_H",(tileX/2.0)*mm, (0.75*foil_thickness/2.0)*mm, ((tileZ+(foil_gap)-hole_center_offset)/2.0) *mm) ;// foil/2.0
m_foilHLogical 	= new G4LogicalVolume(m_foilH
				,m_Al,
				"foil_H_log");
m_foilHLogical->SetVisAttributes( G4Colour(0.9,0.0,0.0,0.4));//


m_foilHtop 			= new G4Box("foil_H_top",(tileX/2.0)*mm, (0.75*foil_thickness/2.0)*mm, ((tileZ)/2.0) *mm) ;//+(foil_gap)-hole_center_offset
m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 1*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset)*mm));
m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 2*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 3*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));
m_foilHtop_hole = new G4SubtractionSolid("foil_H_top_hole",m_foilHtop_hole,m_fiber_subtract,stripRotation, G4ThreeVector(( tileX/2.0 - 4*(tileX/5)  )*mm,0,(tileZ/2.0 - fiber_diam/2.0 - grease_offset - hole_center_offset )*mm));

m_foilHtopLogical 	= new G4LogicalVolume(m_foilHtop_hole
				,m_Al,
				"foil_Htop_log");
m_foilHtopLogical->SetVisAttributes(  G4Colour(0.9,0.0,0.0,0.4) );//

//create alum case pieces

//vertical dividers
m_AlcaseV 			= new G4Box("case_V",(case_thickness/2.0)*mm, ((fiberHeightY[3]+case_thickness)/2.0)*mm, ((final_foil_pos)-(RPD_centerZ-(tileZ/2.0))+case_thickness)/2.0 *mm );
m_AlcaseVLogical 	= new G4LogicalVolume(m_AlcaseV
					,m_Al,
					"case_V_log");
m_AlcaseVLogical->SetVisAttributes( G4Colour(0.0,0.0,0.9,0.4));

//front and back plate
m_Alcase 			= new G4Box("case_fr_bck",((4*tileX+8*halfX_gap+case_thickness)/2.0)*mm, ((fiberHeightY[3]+case_thickness)/2.0)*mm, ((0.99*case_thickness)/2.0)*mm );
m_AlcaseLogical 	= new G4LogicalVolume(m_Alcase
					,m_Al,
					"case_fr_bck_log");
m_AlcaseLogical->SetVisAttributes( G4Colour(0.0,0.0,0.9,0.4));//


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// if Optical is turned on we need to build air detector above each fibers

if (config->GetValue("OPTICAL_ON",false) == 1){

float air_detect_thickness = 0.1;
float correction = 0.0;

m_air_detect = new G4Tubs( "air_detect",
					0.0*mm,
					((core_diam/2.0)-correction)*mm,
					((air_detect_thickness/2.0))*mm ,
					0.0*deg,
					360.0*deg);


int air_detecT_cnt = 0;
for(int i=0;i<4;i++) {
	for(int j=0;j<4;j++) {
		for(int k=0;k<4;k++) {
			sprintf(name,"photo_detect_log_%d_%d_%d",i, j, k);
			m_air_detect_Logical[air_detecT_cnt]	 = new G4LogicalVolume(m_air_detect, m_PMMA, name);
			m_air_detect_Logical[air_detecT_cnt]->SetVisAttributes( G4Colour(0.0,1.0,0.0,5.0) );

			//sprintf(name,"photo_detect %d", k);
			sprintf(name,"photo_detect_phys_%d_%d_%d",i, j, k);//comment

			if(air_detec_flag){

			m_air_detectPhysical[air_detecT_cnt]		  = new G4PVPlacement(
					nullRotation,
					G4ThreeVector(( 0.0 ) *mm ,
									 ( 0.0 ) 	 *mm ,
									 (fiberHeightY[k]/2 - (air_detect_thickness/2.0) - correction)*mm),
					m_air_detect_Logical[air_detecT_cnt],
					name,
					m_fiberLogical[air_detecT_cnt],
					false,
					air_detecT_cnt,
					checkOverlaps);



					air_detecT_cnt++;

				}
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

//place vertical dividers
for(int k=0;k<5;k++) {

	sprintf(name,"case_V %d", k);

	if(vertical_flag){
		m_AlcaseVPhysical[k] = new G4PVPlacement(
						nullRotation,
						G4ThreeVector((( 	(RPD_startX + tileX/2.0) + halfX_gap  ) - k*(2*halfX_gap+tileX)  ) 	*mm ,
													 (RPD_centerY )   	*mm ,
													 (assembly_midZ )			*mm),
						m_AlcaseVLogical,
						name,
						m_logicMother,
						false,
						k,
						checkOverlaps);
						}
}
//place front/back plate
for(int k=0;k<2;k++) {

	sprintf(name,"case_fr_bck %d", k);

	if(fr_bck_flag){
		m_AlcasePhysical[k] = new G4PVPlacement(
						nullRotation,
						G4ThreeVector(  (RPD_centerX) 	*mm ,
														(RPD_centerY)   	*mm ,
													 (assembly_midZ - half_Z_length - case_thickness + 2.0*k*(half_Z_length+case_thickness) )			*mm),
						m_AlcaseLogical,
						name,
						m_logicMother,
						false,
						k,
						checkOverlaps);
					}
}
//place tiles, fibers, foil
  for(int j=0;j<4;j++) {
    for(int i=0;i<4;i++) {



		if(tile_flag){
			sprintf(name,"tilephys_%d_%d", i,j);
    	m_tilePhysical[i][j] = new G4PVPlacement(
							nullRotation,
							G4ThreeVector( ( 	RPD_startX - (i*( tileX+(2*halfX_gap) ) ) )  	*mm ,
										   			( 	RPD_startY - (j*( tileY+(2*halfY_gap) ) ) )	*mm ,
											 					tileZcenter[j] *mm),
							m_tileLogical,
							name,
							m_logicMother,
							false,
							cn,
							checkOverlaps);
						}



			if(foil_flag){
				sprintf(name,"foil_V_%d_%d",  i,j);
				m_foilVPhysical[i][j] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector( ( 	RPD_startX - (j*( tileX+(2*halfX_gap) ) ) )  	*mm ,
											   		   ( RPD_startY + tileY/2 - (fiberHeightY[i]/2)  ) 		*mm ,
												 			 ( foilZcenter[i]) *mm),
								m_foilVLogical[i],
								name,
								m_logicMother,
								false,
								cn,
								checkOverlaps);

			sprintf(name,"foil_V_front_%d_%d",  i,j);

				m_foilVfrontPhysical[i][j] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector( ( 	RPD_startX - (i*( tileX+(2*halfX_gap) ) ) )  	*mm ,
											   		  ( 	RPD_startY - (j*( tileY+(2*halfY_gap) ) ) ) 		*mm ,
												 			 ( tileZcenter[j] - tileZ/2.0 -foil_gap/2.0) *mm),//fiberZcenter[i] + fiber_diam/2.0 + foil_gap/2.0;
								m_foilVfrontLogical,
								name,
								m_logicMother,
								false,
								cn,
								checkOverlaps);


			sprintf(name,"foil_H_%d_%d", i,j);

				m_foilHPhysical[i][j] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector( ( 	RPD_startX - (j*( tileX+(2*halfX_gap) ) ) )  	*mm ,
											   		   ( RPD_startY + tileY/2 + 1.3*halfY_gap - (i+1) * (tileY + 2 * halfY_gap ))   		*mm ,
												 			 (tileZcenter[i]  - hole_center_offset/2.0)			*mm),
								m_foilHLogical,
								name,
								m_logicMother,
								false,
								i,
								checkOverlaps);

			if(i==0){

			sprintf(name,"foil_H_top_%d", j);

				m_foilHtopPhysical[j] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector( ( 	RPD_startX - (j*( tileX+(2*halfX_gap) ) ) )  	*mm ,
											   		   (RPD_startY + tileY/2 + 1.5*foil_thickness)   		*mm ,
												 			 (tileZcenter[0] *mm)),
								m_foilHtopLogical,
								name,
								m_logicMother,
								false,
								j,
								checkOverlaps);
							}
				}

		for(int k=0;k<4;k++) {
			if(fiber_flag){
				sprintf(name,"fiberphys_%d_%d_%d",i, j, k);
				m_fiberPhysical[cn_fiber]		  = new G4PVPlacement(
							stripRotation,
							G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,
										   ( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
											 (fiberZcenter[k]) *mm),
							m_fiberLogical[cn_fiber],
							name,
							m_logicMother,
							false,
							cn_fiber,
							checkOverlaps);
						}
			if(cladding_flag){
				sprintf(name,"cladphys_%d_%d_%d",i, j, k);
				m_fibercladPhysical[cn_fiber]		  = new G4PVPlacement(
							stripRotation,
							G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,
												( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
												(fiberZcenter[k]) *mm),
							m_fibercladLogical[k],
							name,
							m_logicMother,
							false,
							cn_fiber,
							checkOverlaps);
						}
			if(grease_flag){
				sprintf(name,"greasephys_%d_%d_%d",i, j, k);
				m_fibergreasePhysical[cn_fiber]		  = new G4PVPlacement(
							stripRotation,
							G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,
													( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
													(fiberZcenter[k]) *mm),
							m_fibergreaseLogical[k],
							name,
							m_logicMother,
							false,
							cn_fiber,
							checkOverlaps);
									}

							cn_fiber++;
		}

			//
      // std::cout  << std:: endl << "tilex = "
			// 	<<  RPD_startX - (i*( tileX+(2*halfX_gap) ) )
			// 	<< ", tiley = "
			// 	<<  RPD_startY - (j*( tileY+(2*halfY_gap) ) )
			// 	<< ", tilez = "
			// 	<< tileZcenter[j]
			// 	<< ", ("
			// 	<< i
			// 	<< ","
			// 	<< j
			// 	<< ")"
			// 	<< std::endl << std::endl;


      cn++;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// TEST SETUP
if(test_tile){
		float testtileY = 100;
		float testtileZ = 50;
		float testgreaseY = 30;
		float testgreaseZ = 2.5;
		float testPDY = 2.5;
		float testcoreZ = 5;
		float testblckZ 	= 15;
		float offsetZ 	= testtileZ/2+testblckZ/2;
		float offsetY 	= 5;

	  m_test_tile 		= new G4Box("test_tile",(100/2)*mm, (testtileY/2)*mm, (testtileZ/2)*mm);
		m_test_alum 		= new G4Box("test_alum",(100/2)*mm, (testtileY/2)*mm, (testblckZ/2)*mm);
		m_test_wls 			= new G4Box("test_wls",(100/2)*mm, (testgreaseY/2)*mm, (testcoreZ/2)*mm);
		m_test_PD 			= new G4Box("test_PD",(100/2)*mm, (testPDY/2)*mm, (testcoreZ/2)*mm);
		m_test_clad 		= new G4Box("test_clad",(100/2)*mm, (testgreaseY/2)*mm, (testgreaseZ/2)*mm);
		m_test_grease 	= new G4Box("test_grease",(100/2)*mm, (testgreaseY/2)*mm, (testgreaseZ/2)*mm);
		m_test_block 		= new G4Box("test_blck1",(100/2)*mm, ((testgreaseY)/2)*mm, (testblckZ/2)*mm);


		m_test_tileLogical	 = new G4LogicalVolume(m_test_tile, m_matQuartz, "testtile_Logical");
		m_test_alumLogical	 = new G4LogicalVolume(m_test_alum, m_Al, "testalum_Logical");
		m_test_wlsLogical	 = new G4LogicalVolume(m_test_wls, m_PMMA, "testwls_Logical");
		m_test_PDLogical	 = new G4LogicalVolume(m_test_PD, m_PMMA, "testPD_Logical");
		m_test_cladLogical	 = new G4LogicalVolume(m_test_clad, m_Poly, "testclad_Logical");
		m_test_greaseLogical	 = new G4LogicalVolume(m_test_grease, m_Grease, "testgrease_Logical");
		m_test_blockLogical	 = new G4LogicalVolume(m_test_block, m_Al, "testblck_Logical");

		m_test_tileLogical->SetVisAttributes( G4Colour(0.0,1.0,0.0,0.3) );
		m_test_alumLogical->SetVisAttributes( G4Colour(0.8,0.0,0.0,0.3) );
		m_test_wlsLogical->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.3) );
		m_test_PDLogical->SetVisAttributes( G4Colour(0.0,0.0,1.0,0.3) );
		m_test_cladLogical->SetVisAttributes( G4Colour(0.3,0.5,0.0,0.3) );
		m_test_greaseLogical->SetVisAttributes( G4Colour(0.0,0.6,0.8,0.3) );
		m_test_blockLogical->SetVisAttributes( G4Colour(0.8,0.0,0.0,0.3) );

		m_test_tilePhysical = new G4PVPlacement(
						nullRotation,
						G4ThreeVector(( 0.0  ) 	*mm ,
													 (0.0 )   	*mm ,
													 (assembly_midZ)			*mm),
						m_test_tileLogical,
						"testtilephys",
						m_logicMother,
						false,
						1,
						checkOverlaps);

		m_test_greasePhysical[0] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector(( 0.0  ) 	*mm ,
															 (testtileY/2 + testgreaseY/2 - offsetY)   	*mm ,
															 (assembly_midZ - testcoreZ/2 - 1.5*testgreaseZ +offsetZ)			*mm),
								m_test_greaseLogical,
								"testgrease1phys",
								m_logicMother,
								false,
								0,
								checkOverlaps);

		m_test_greasePhysical[1] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector(( 0.0  ) 	*mm ,
															 (testtileY/2 + testgreaseY/2- offsetY )   	*mm ,
															 (assembly_midZ + testcoreZ/2 + 1.5*testgreaseZ+offsetZ )			*mm),
								m_test_greaseLogical,
								"testgrease2phys",
								m_logicMother,
								false,
								1,
								checkOverlaps);

		m_test_cladPhysical[0] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector(( 0.0  ) 	*mm ,
															 (testtileY/2 + testgreaseY/2- offsetY )   	*mm ,
															 (assembly_midZ - testcoreZ/2 - testgreaseZ/2+offsetZ )			*mm),
								m_test_cladLogical,
								"testclad1phys",
								m_logicMother,
								false,
								0,
								checkOverlaps);

		m_test_cladPhysical[1] = new G4PVPlacement(
								nullRotation,
								G4ThreeVector(( 0.0  ) 	*mm ,
															 (testtileY/2 + testgreaseY/2- offsetY )   	*mm ,
															 (assembly_midZ+ testcoreZ/2 + testgreaseZ/2+offsetZ )			*mm),
								m_test_cladLogical,
								"testclad2phys",
								m_logicMother,
								false,
								1,
								checkOverlaps);

		m_test_wlsPhysical = new G4PVPlacement(
							nullRotation,
							G4ThreeVector(( 0.0  ) 	*mm ,
														 (testtileY/2 + testgreaseY/2- offsetY )   	*mm,
														 (assembly_midZ + offsetZ)			*mm),
							m_test_wlsLogical,
							"testtilewls",
							m_logicMother,
							false,
							1,
							checkOverlaps);

		m_test_PDPhysical = new G4PVPlacement(
						nullRotation,
						G4ThreeVector(( 0.0  ) 	*mm ,
													 (testgreaseY/2 - testPDY/2  )   	*mm ,
													 (0.0)			*mm),
						m_test_PDLogical,
						"photo_detect",
						m_test_wlsLogical,
						false,
						1,
						checkOverlaps);

		m_test_alumPhysical = new G4PVPlacement(
								nullRotation,
								G4ThreeVector(( 0.0  ) 	*mm ,
															 (0.0 - offsetY )   	*mm ,
															 (assembly_midZ + testblckZ/2 + testtileZ/2 )			*mm),
								m_test_alumLogical,
								"testalumphys",
								m_logicMother,
								false,
								1,
								checkOverlaps);

	m_test_blockPhysical[0] = new G4PVPlacement(
							nullRotation,
							G4ThreeVector(( 0.0  ) 	*mm ,
														 (testtileY/2 + testgreaseY/2 )   	*mm ,
														(assembly_midZ + testblckZ/2 + 2.0* testgreaseZ + testcoreZ/2 +offsetZ)			*mm),
							m_test_blockLogical,
							"testblck1phys",
							m_logicMother,
							false,
							1,
							checkOverlaps);

	m_test_blockPhysical[1] = new G4PVPlacement(
							nullRotation,
							G4ThreeVector(( 0.0  ) 	*mm ,
														 (testtileY/2 + testgreaseY/2   )   	*mm ,
														 (assembly_midZ - testblckZ/2 - 2.0* testgreaseZ - testcoreZ/2 +offsetZ)			*mm),
							m_test_blockLogical,
							"testblck2phys",
							m_logicMother,
							false,
							1,
							checkOverlaps);


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PAN FLUTE RPD
if(pan_flute){

float foil_spacing = 1;
int m_PFrpd_cnt=0;
int m_PFfoil_cnt=0;
float PFfoil_thickness = 0.1;
float rod_diam = tileX/4.0;
float readout_Y = 0.5;
float pan_flute_start = RPD_centerX + 7.5*rod_diam + 7.5*foil_spacing;

for(int k=0;k<4;k++) {
	fiberHeightY[k]=((k+1)*tileY)+(k*2*halfY_gap) + readout_Y;

	sprintf(name,"m_PFrpd_%d", k);
	m_PFrpd[k] 		= new G4Tubs( name,
						0.0*mm,
						(rod_diam/2.0)*mm,
						fiberHeightY[k]*mm/2.0 ,
						0.0*deg,
						360.0*deg);
	//cylindrical foil dividers
  sprintf(name,"m_PFrpd_%d", k);
	m_test_foil[k] 		= new G4Tubs( name,
						(rod_diam/2.0 + PFfoil_thickness/2)*mm,
						(rod_diam/2.0 + PFfoil_thickness)*mm,
						fiberHeightY[k]*mm/2.0 ,
						0.0*deg,
						360.0*deg);
}

	m_PFdetec 		= new G4Tubs( "m_PFdetec",
						0.0*mm,
						(rod_diam/2.0)*mm,
						(readout_Y/2.0)*mm ,
						0.0*deg,
						360.0*deg);

	for(int j=0;j<4;j++) {
		for(int i=0;i<4;i++) {
			for(int k=0;k<4;k++) {
			sprintf(name,"m_PFrpd_log_%d_%d_%d",j,i,k);

					m_PFrpdLogical[m_PFrpd_cnt] 	= new G4LogicalVolume(m_PFrpd[(k+i)%4]
										,m_matQuartz,
										name);
										if(k==0) m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Cyan() 	 );
					else 			if(k==1) m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Red() );
					else 			if(k==2) m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Green() );
					else 			if(k==3) m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Magenta() );

					m_PFrpdLogical[m_PFrpd_cnt]->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

			sprintf(name,"m_foil_log_%d_%d_%d",j,i,k);

					m_test_foilLogical[m_PFrpd_cnt] 	= new G4LogicalVolume(m_test_foil[(k+i)%4]
										,m_Al,
										name);
					m_test_foilLogical[m_PFrpd_cnt]->SetVisAttributes(  G4VisAttributes::Invisible );//G4Colour(1.0,0.0,0.0,0.3)

			sprintf(name,"m_PFdetec_log_%d_%d_%d",j,i,k);

					m_PFdetecLogical[m_PFrpd_cnt] 	= new G4LogicalVolume(m_PFdetec
										,m_matQuartz,
										name);
					m_PFdetecLogical[m_PFrpd_cnt]->SetVisAttributes(   G4VisAttributes::Invisible);//G4Colour(1.0,0.0,0.0,0.3)

			sprintf(name,"m_PFrpd_pyhs_%d_%d_%d",j,i,k);
					m_PFrpdPhysical[m_PFrpd_cnt]		  = new G4PVPlacement(
								stripRotation,
								G4ThreeVector( ( pan_flute_start - (m_PFfoil_cnt*( rod_diam + foil_spacing ) ) )  *mm ,
												 			( RPD_startY + tileY/2 - (fiberHeightY[(k+i)%4]/2) - readout_Y/2   )	 *mm ,
												 			(tileZcenter[0]-k*( rod_diam + foil_spacing)) *mm),
								m_PFrpdLogical[m_PFrpd_cnt],
								name,
								m_logicMother,
								false,
								m_PFrpd_cnt,
								checkOverlaps);

			sprintf(name,"m_PFfoil_pyhs_%d_%d_%d",j,i,k);
					m_test_foilPhysical[m_PFrpd_cnt]		  = new G4PVPlacement(
								stripRotation,
								G4ThreeVector( ( pan_flute_start - (m_PFfoil_cnt*( rod_diam + foil_spacing ) ) )  *mm ,
												 			( RPD_startY + tileY/2 - (fiberHeightY[(k+i)%4]/2) - readout_Y/2   )	 *mm ,
												 			(tileZcenter[0]-k*( rod_diam + foil_spacing)) *mm),
								m_test_foilLogical[m_PFrpd_cnt],
								name,
								m_logicMother,
								false,
								m_PFrpd_cnt,
								checkOverlaps);

			sprintf(name,"photo_det_phys_%d_%d_%d",j,i,k);
					m_PFdetecPhysical[m_PFrpd_cnt]		  = new G4PVPlacement(
								nullRotation,
								G4ThreeVector( 0.0   *mm ,
												 			(0.0) *mm,
															((fiberHeightY[(k+i)%4]/2) - readout_Y/2   )	 *mm) ,
								m_PFdetecLogical[m_PFrpd_cnt],
								name,
								m_PFrpdLogical[m_PFrpd_cnt],
								false,
								m_PFrpd_cnt,
								checkOverlaps);




								m_PFrpd_cnt++;
							}
						m_PFfoil_cnt++;

					}
				}
}


  //----------------------------------------------
  // SD and Scoring Volumes
  //----------------------------------------------


 G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Note one SD object for each module
  char quartzSDname[256];
  sprintf( quartzSDname, "RPD_SD%d", m_modNum+1);

if (config->GetValue("OPTICAL_ON",false) == 0){
 	RpdSD* aRpdSD = new RpdSD( quartzSDname, m_sd, m_modNum );
  aRpdSD->HistInitialize();
  SDman->AddNewDetector( aRpdSD );
  m_tileLogical->SetSensitiveDetector( aRpdSD );
}

   char fiberSDname[256];
  sprintf( fiberSDname, "RPD%d_SD", m_modNum+1);


  FiberSD* aFiberSD = new FiberSD( fiberSDname, m_sd, m_modNum+3 );
  aFiberSD->HistInitialize();
  SDman->AddNewDetector( aFiberSD );


if(pan_flute){
	for(int i=0;i<64;i++){
	 m_PFdetecLogical[i]->SetSensitiveDetector( aFiberSD );
 	}
}
else if(test_tile) m_test_PDLogical->SetSensitiveDetector( aFiberSD );
else{
	if (config->GetValue("OPTICAL_ON",false) == 1){
		for(int i=0;i<64;i++){
				m_air_detect_Logical[i]->SetSensitiveDetector( aFiberSD );}}
	else{
		for(int k=0;k<64;k++) {
				m_fiberLogical[k]->SetSensitiveDetector( aFiberSD );
  		}
		}
	}

  std::cout << "ModTypeRPD construction finished: SD name " << quartzSDname << std::endl;

  std::cout << "Fiber      construction finished: SD name " << fiberSDname << std::endl;

}
