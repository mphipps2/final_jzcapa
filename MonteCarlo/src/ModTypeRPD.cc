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

ModTypeRPD::ModTypeRPD(const int cn, G4LogicalVolume* mother, G4ThreeVector* pos )
  : m_modNum( cn ),  m_pos( pos ), m_logicMother( mother )
{
	materials = Materials::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD(const int cn, ModTypeRPD* right)
  : m_modNum( cn )
{
	m_pos 						 = new G4ThreeVector(*right->m_pos);
	m_fiberDiam 			 = new G4ThreeVector(*right->m_fiberDiam);
	m_HousingThickness = right->m_HousingThickness;
	m_fiberPitch 			 = right->m_fiberPitch;
	m_tileSize 				 = right->m_tileSize;
	m_minWallThickness = right->m_minWallThickness;
  m_detType 				 = right->m_detType;
  OPTICAL 					 = right->OPTICAL;
  CHECK_OVERLAPS 		 = right->CHECK_OVERLAPS;
	materials					 = right->materials;
	m_logicMother 		 = right->m_logicMother;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD()
  : m_modNum( 0 ), m_pos(new G4ThreeVector(0.,0.,0.)), m_logicMother(NULL)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::~ModTypeRPD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::Construct(){
  DefineMaterials();

	if(m_detType == "cms"){
		ConstructCMSDetector();
	}else {
		ConstructPrototypeDetector();
	}
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

void ModTypeRPD::ConstructPrototypeDetector()
{
	//retrieve RPD parameters
	float fiber_diam 	= m_fiberDiam->x(); // Just core for now

	char name[256];

	float RPD_centerX = m_pos->getX();
	float RPD_centerY = m_pos->getY();
	float RPD_centerZ = m_pos->getZ();

	int n_rows = 4;
	int n_columns = 4;
	int n_cycles_per_tile = (m_tileSize/fiber_diam)/n_rows;   //Divide to round down to a whole number
	int n_fibers_per_tile = n_cycles_per_tile*2*n_rows;  //The pattern will be repeated twice per tile

	//If you asked for a
	if(.707*m_minWallThickness < fiber_diam){
		m_minWallThickness = 1.414*fiber_diam;
		std::cout << "Minimum wall thickness entered is too small!!!" << std::endl;
		std::cout << "Using " << m_minWallThickness << "mm instead" << std::endl;
	}
	// If the remaining space in x leaves less than 1mm for each wall, remove a cycle and calculate
	// the new wall thickness
	calculate_wall_thickness:
	float wall_thickness = ((m_tileSize - n_fibers_per_tile*fiber_diam)/n_fibers_per_tile)*mm;
	if( wall_thickness < m_minWallThickness*mm ){
		--n_cycles_per_tile;
		n_fibers_per_tile = n_cycles_per_tile*n_rows;
		goto calculate_wall_thickness;
	}

	// Distance on center of fibers
	float pitch = wall_thickness + fiber_diam;
	// Distance in X and Z from one fiber to another diagonally from it
	float offset = pitch/2;
	// Distance between fiber and aluminum sheath
	float gap = 0.01*mm;
	// Will be filled with fiber height for each row
	float fiber_height;
	// Distance from the top of the RPD area of interest to the readout
	float distance_to_readout = 0*mm;
	// Thickness of the readout puck
	float readout_thickness = 0.1*mm;
	// Positions of the current fiber for pattern1 and pattern2
	float posx1, posz1, posy, posx2, posz2;
	// Count the number of fibers as we go along
	m_PFrpd_cnt = 0;

	//create some rotation matrices
	G4RotationMatrix* stripRotation = new G4RotationMatrix();
	stripRotation->rotateX(90.*deg);
	G4RotationMatrix* nullRotation = new G4RotationMatrix();

	// The SD caps that will go on the fibers. Each fiber will get a separate SD,
	// so we wait to create the logical volume in the loop
	m_PFdetec = new G4Tubs( "m_PFdetec",
													0.0*mm,
													(fiber_diam/2.0)*mm,
													(readout_thickness/2.0)*mm ,
													0.0*deg,
													360.0*deg);

	for(int row = 0; row < n_rows; row++){
		fiber_height = m_tileSize*(row+1) + distance_to_readout;

		// The fiber solid. One length for each row
		sprintf(name,"m_PFrpd_%d", row);
		m_PFrpd[row] = new G4Tubs(name,
														0.0*mm,
														(fiber_diam/2.0)*mm,
														fiber_height*mm/2.0 ,
														0.0*deg,
														360.0*deg);

		//Cylindrical foil dividers
		sprintf(name,"m_PFfoil_%d", row);
		m_PFrpd_foil[row] = new G4Tubs(name,
																(fiber_diam/2.0 + gap/2.0)*mm,
																(fiber_diam/2.0 + gap/2.0 + wall_thickness/4)*mm,
																fiber_height*mm/2.0 ,
																0.0*deg,
																360.0*deg);

		sprintf(name,"m_PFfoil_log_%d",row);
		m_PFrpd_foilLogical[row] 	= new G4LogicalVolume(m_PFrpd_foil[row]
							,m_Al,
							name);
		m_PFrpd_foilLogical[row]->SetVisAttributes( G4Colour(1.0,0.0,0.0,0.3) );// G4Colour(1.0,0.0,0.0,0.3) or G4VisAttributes::Invisible

		for(int col = 0; col < n_columns; col++){
			//Now we're in the realm of working on a single tile
			//we have to cycle through the two patterns until the tile
			//is filled in X
			for(int cycle = 0; cycle < n_cycles_per_tile; cycle++){

				for(int fiber = 0; fiber < n_columns; fiber++){
					// !!!!!!!Position calculations assume an even number of rows and columns!!!!!!! //////

					//Start at RPD center + tile width* number of tiles + cycle number * cycle width + stack number in cycle * pitch
					posx1 = RPD_centerX + m_tileSize*((n_columns/2) - col ) - (cycle*n_columns*pitch) - fiber*pitch - pitch/4;  //ARIC ADDED - pitch/4
					posx2 = posx1 - offset;

					//Start at Z center - half the stack depth + moving front to back and looping back to front
					posz1 = RPD_centerZ - pitch*n_columns/2 + pitch*((row + fiber	  )%4);
					posz2 = RPD_centerZ - pitch*n_columns/2 + pitch*((row + fiber + 2)%4) + offset; //Pattern is offset by 2

					//Start at RPDY center + distance to bottom of top tile + half the fiber height
					posy = RPD_centerY + ((n_rows/2 - 1) - row)*m_tileSize + fiber_height/2;


						std::cout << "Rod# = " << m_PFrpd_cnt << ", posx1 = " << posx1 << ", posy = " << posy  << "fheight = " << fiber_height << std::endl;
						std::cout << "Rod# = " << m_PFrpd_cnt + 1 << ", posx2 = " << posx2 << ", posy = " << posy  << std::endl;


					//----------------------- Place the rods -----------------------//

					sprintf(name,"m_PFrpd_log_%d_%d_%d_%d_0",row,col,cycle,fiber);
					m_PFrpdLogical[m_PFrpd_cnt] =
						new G4LogicalVolume(m_PFrpd[row],
																m_matQuartz,
																name);
					m_PFrpdLogical[m_PFrpd_cnt]->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

					sprintf(name,"m_PFrpd_pyhs_%d_%d_%d_%d_0",row,col,cycle,fiber);
					m_PFrpdPhysical[m_PFrpd_cnt] =
						new G4PVPlacement(stripRotation,
															G4ThreeVector(posx1*mm, posy*mm, posz1*mm),
															m_PFrpdLogical[m_PFrpd_cnt],
															name,
															m_logicMother,
															false,
															m_PFrpd_cnt,
															CHECK_OVERLAPS);


					sprintf(name,"m_PFrpd_log_%d_%d_%d_%d_1",row,col,cycle,fiber);
					m_PFrpdLogical[m_PFrpd_cnt+1] =
						new G4LogicalVolume(m_PFrpd[row],
																m_matQuartz,
																name);
					m_PFrpdLogical[m_PFrpd_cnt+1]->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,10*ms));

					sprintf(name,"m_PFrpd_pyhs_%d_%d_%d_%d_1",row,col,cycle,fiber);
					m_PFrpdPhysical[m_PFrpd_cnt+1] =
						new G4PVPlacement(stripRotation,
															G4ThreeVector(posx2*mm, posy*mm, posz2*mm),
															m_PFrpdLogical[m_PFrpd_cnt+1],
															name,
															m_logicMother,
															false,
															m_PFrpd_cnt+1,
															CHECK_OVERLAPS);

				 if(row==0){
								 m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Cyan()		);
								 m_PFrpdLogical[m_PFrpd_cnt+1]->SetVisAttributes( G4Colour::Cyan()		);}
				 if(row==1){
								 m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Red()			);
								 m_PFrpdLogical[m_PFrpd_cnt+1]->SetVisAttributes( G4Colour::Red()			);}
				 if(row==2){
								 m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Green()		);
								 m_PFrpdLogical[m_PFrpd_cnt+1]->SetVisAttributes( G4Colour::Green()		);}
				 if(row==3){
								 m_PFrpdLogical[m_PFrpd_cnt]->SetVisAttributes( G4Colour::Magenta()	);
								 m_PFrpdLogical[m_PFrpd_cnt+1]->SetVisAttributes( G4Colour::Magenta()	);}

					//----------------------- Place the walls/gap filler -----------------------//
					sprintf(name,"m_PFfoil_pyhs_%d_%d_%d_%d_0",row,col,cycle,fiber);
					m_PFrpd_foilPhysical[m_PFrpd_cnt] =
						new G4PVPlacement(stripRotation,
															G4ThreeVector(posx1*mm, posy*mm, posz1*mm),
															m_PFrpd_foilLogical[row],
															name,
															m_logicMother,
															false,
															m_PFrpd_cnt,
															CHECK_OVERLAPS);

					sprintf(name,"m_PFfoil_pyhs_%d_%d_%d_%d_1",row,col,cycle,fiber);
					m_PFrpd_foilPhysical[m_PFrpd_cnt+1] =
						new G4PVPlacement(stripRotation,
															G4ThreeVector(posx2*mm, posy*mm, posz2*mm),
															m_PFrpd_foilLogical[row],
															name,
															m_logicMother,
															false,
															m_PFrpd_cnt+1,
															CHECK_OVERLAPS);

					//----------------------- Make and place the photodector volumes -----------------------//
					sprintf(name,"m_PFdetec_log_%d_%d_%d_%d_0",row,col,cycle,fiber);
					m_PFdetecLogical[m_PFrpd_cnt]	=
						new G4LogicalVolume(m_PFdetec,
																m_matQuartz,
																name);
					m_PFdetecLogical[m_PFrpd_cnt]->SetVisAttributes(G4Colour(0.0,6.0,4.0,0.3));//G4Colour(1.0,0.0,0.0,0.3)

					sprintf(name,"photo_det_phys_%d_%d_%d_%d_0",row,col,cycle,fiber);
					m_PFdetecPhysical[m_PFrpd_cnt] =
						new G4PVPlacement(nullRotation,
															G4ThreeVector(0*mm, 0*mm, ((fiber_height/2-readout_thickness/2))*mm),
															m_PFdetecLogical[m_PFrpd_cnt],
															name,
															m_PFrpdLogical[m_PFrpd_cnt],
															false,
															m_PFrpd_cnt,
															CHECK_OVERLAPS);

					sprintf(name,"m_PFdetec_log_%d_%d_%d_%d_1",row,col,cycle,fiber);
					m_PFdetecLogical[m_PFrpd_cnt+1]	=
						new G4LogicalVolume(m_PFdetec,
																m_matQuartz,
																name);
					m_PFdetecLogical[m_PFrpd_cnt+1]->SetVisAttributes(G4Colour(0.0,6.0,4.0,0.3));//G4Colour(1.0,0.0,0.0,0.3)

					sprintf(name,"photo_det_phys_%d_%d_%d_%d_1",row,col,cycle,fiber);
					m_PFdetecPhysical[m_PFrpd_cnt+1] =
						new G4PVPlacement(nullRotation,
															G4ThreeVector(0*mm, 0*mm, ((fiber_height/2-readout_thickness/2))*mm),
															m_PFdetecLogical[m_PFrpd_cnt+1],
															name,
															m_PFrpdLogical[m_PFrpd_cnt+1],
															false,
															m_PFrpd_cnt+1,
															CHECK_OVERLAPS);

					m_PFrpd_cnt++;
					m_PFrpd_cnt++;

					//std::cout << "FIBER END" << std::endl;
				}//end fiber
				//std::cout << "CYCLE END" << std::endl;
			}//end cycle
			//std::cout << "COL END" << std::endl;
		}//end column
		//std::cout << "ROW END" << std::endl;
	}//end row


	//----------------------------------------------
	// SD and Scoring Volumes
	//----------------------------------------------


	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	//Note one SD object for each module
	char fiberSDname[256];
	sprintf( fiberSDname, "RPD%d_SD", m_modNum+1);

	FiberSD* aFiberSD = new FiberSD( fiberSDname, m_modNum, OPTICAL );
	aFiberSD->HistInitialize();
	SDman->AddNewDetector( aFiberSD );
	m_tileLogical->SetSensitiveDetector( aFiberSD );



	std::cout << "Prototype RPD construction finished: SD name " << fiberSDname << std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::ConstructCMSDetector()
{
	bool test_tile = false;

	if(test_tile){
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



	//retrieve RPD parameters
  float tileX 		= 20;
  float tileY 		= 20;
  float tileZ 		= 20;
  float halfX_gap 	= (1.58/2)*mm;
	float fiber_diam 	= 1*mm;
	float foil_thickness 	= 0.016*mm;


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

	float RPD_centerX = m_pos->getX();
	float RPD_centerY = m_pos->getY();
	float RPD_centerZ = m_pos->getZ();

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

if ( OPTICAL ){

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
					CHECK_OVERLAPS);



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
						CHECK_OVERLAPS);
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
						CHECK_OVERLAPS);
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
							CHECK_OVERLAPS);
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
								CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);


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
								CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);
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
							CHECK_OVERLAPS);
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
							CHECK_OVERLAPS);
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
							CHECK_OVERLAPS);
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
						CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);

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
							CHECK_OVERLAPS);

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
						CHECK_OVERLAPS);

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
								CHECK_OVERLAPS);

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
							CHECK_OVERLAPS);

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
							CHECK_OVERLAPS);


}

  //----------------------------------------------
  // SD and Scoring Volumes
  //----------------------------------------------


 G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Note one SD object for each module
  char fiberSDname[256];
  sprintf( fiberSDname, "RPD%d_SD", m_modNum+1);

  FiberSD* aFiberSD = new FiberSD( fiberSDname, m_modNum+1, OPTICAL );
  aFiberSD->HistInitialize();
  SDman->AddNewDetector( aFiberSD );



if(test_tile) m_test_PDLogical->SetSensitiveDetector( aFiberSD );
else{
	if (OPTICAL){
		for(int i=0;i<64;i++){
				m_air_detect_Logical[i]->SetSensitiveDetector( aFiberSD );}}
	else{
		for(int k=0;k<64;k++) {
				m_fiberLogical[k]->SetSensitiveDetector( aFiberSD );
  		}
		}
	}

  std::cout << "ModTypeRPD construction finished: SD name " << fiberSDname << std::endl;

  std::cout << "Fiber      construction finished: SD name " << fiberSDname << std::endl;

}
