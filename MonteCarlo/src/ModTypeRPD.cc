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
#include "QuartzSD.hh"
#include "RpdSD.hh"
#include "FiberSD.hh"
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

ModTypeRPD::ModTypeRPD(const int cn,const G4ThreeVector& pos,
	      G4LogicalVolume* mother, SharedData* sd)
  : m_modNum( cn ),  m_pos( pos ), m_logicMother( mother ),
    m_sd( sd ),
    m_matQuartz(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModTypeRPD::ModTypeRPD()
  : m_modNum( 0 ), m_pos(G4ThreeVector()), m_logicMother(NULL),
    m_sd(NULL), 
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

  //----------------------------------------------     
  // Define Material Properties
  //----------------------------------------------
  const G4int NUMENTRIES = 2;
  
  G4double ephoton         [NUMENTRIES] = {2.00*eV,4.80*eV};

  G4double rindexCore[NUMENTRIES] = {1.46,1.46};
  //  G4double absorptionCore[NUMENTRIES] = {46*m,46*m};
  //  G4double rindexCladding[NUMENTRIES] = {1.46,1.46};
  //  G4double absorptionCladding[NUMENTRIES] = {46*m,46*m};  

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModTypeRPD::ConstructDetector()
{
  // Get Config
  TEnv* config = m_sd->GetConfig();

  float tileX 		= config->GetValue("tileXsize",20);
  float tileY 		= config->GetValue("tileYsize",20);
  float tileZ 		= config->GetValue("tileZsize",10);
  float halfX_gap 	= ((config->GetValue("modGapX",1.58))/2)*mm;
  float halfY_gap 	= ((config->GetValue("modGapY",0.5))/2)*mm;
  
  char name[256];
  
    // Option to switch on/off checking of volumes overlaps
  bool checkOverlaps = config->GetValue("checkOverlaps",false);
  
  G4RotationMatrix* nullRotation = new G4RotationMatrix();

  m_tile		 = new G4Box("tile",(tileX/2)*mm, (tileY/2)*mm, (tileZ/2)*mm);
  m_tileLogical	 = new G4LogicalVolume(m_tile, m_matQuartz, "tile_Logical");		
 
  G4VisAttributes* quartzColor  = new G4VisAttributes( G4Colour::Cyan() );
	//quartzColor->SetForceSolid(true);  
	m_tileLogical->SetVisAttributes( quartzColor );

	int cn = 0, cn_fiber = 0;
	
	float RPD_centerX = m_pos.getX();
	float RPD_centerY = m_pos.getY();
	float RPD_centerZ = m_pos.getZ();
	
	float RPD_startX = RPD_centerX + 3*halfX_gap + 1.5*tileX;
	float RPD_startY = RPD_centerY + 3*halfY_gap + 1.5*tileY;
	
	// Readout fibers behind RPD
	G4RotationMatrix* stripRotation = new G4RotationMatrix();
	stripRotation->rotateX(90.*deg);
	float fiberHeightY[4]; 
	
	
	for(int k=0;k<4;k++) {
		fiberHeightY[k]=((k+1)*tileY)+(k*2*halfY_gap);
		
		sprintf(name,"fiber_a %d", k);
		
		m_fiber[k] 		= new G4Tubs( name, 
							0.0*mm, 
							1.5/2*mm, 
							fiberHeightY[k]*mm/2.0 , 
							0.0*deg, 
							360.0*deg);
							
		sprintf(name,"fiberLogical_a %d", k);					
		
		m_fiberLogical[k] 	= new G4LogicalVolume(m_fiber[k]          
							,m_matQuartz, 
							name);
		m_fiberLogical[k]->SetVisAttributes( G4Colour::Green() );				
	}
	
	
	
	
  for(int j=0;j<4;j++) {  
    for(int i=0;i<4;i++) { 
      
		sprintf(name,"tile %d", cn);
    m_tilePhysical[i][j] = new G4PVPlacement(
							nullRotation,
							G4ThreeVector( ( RPD_startX - (i*( tileX+(2*halfX_gap) ) ) )  	*mm ,   
										   ( RPD_startY - (j*( tileY+(2*halfY_gap) ) ) )	*mm ,
											 RPD_centerZ 									*mm),
							m_tileLogical,
							name,
							m_logicMother,
							false,
							cn,
							checkOverlaps);  
		
		for(int k=0;k<4;k++) { 
			//sprintf(name,"fiber %d_%d_%d_%d",i, j, k, cn_fiber );
			sprintf(name,"%d%d%d",i, j, k);
			
			m_fiberPhysical[cn_fiber]		  = new G4PVPlacement(   //m_fiberPhysical[cn_fiber][j]
							stripRotation,
							G4ThreeVector( ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) *mm ,   
										   ( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )	 *mm ,
											 RPD_centerZ + tileZ + (k * 2)					 *mm),
							m_fiberLogical[k],
							name,
							m_logicMother,
							false,
							cn_fiber,
							checkOverlaps);
							
							cn_fiber++;
							
				std::cout << cn_fiber-1 
				<< " height = "
				<< fiberHeightY[k]
				<< ", name = " 
				<< name
				<<  ": Fx = " 
				<< ( RPD_startX - (j*( tileX+(2*halfX_gap) ) ) ) + (tileX/2) - ((i+1)*tileX/5) 
				<< ", Fy = " 
				<<  ( RPD_startY + tileY/2 - (fiberHeightY[k]/2)  )
				<< ", Fz = " 
				<< RPD_centerZ + tileZ + (k * 2) 
				<< ", ("
				<< i 
				<< "," 
				<< j 
				<< ","
				<< k
				<< ")" 
				<< std::endl;	
							
		}

							
      std::cout  << std:: endl << "Sx = " 
				<<  RPD_startX - (i*( tileX+(2*halfX_gap) ) )  
				<< ", Sy = " 
				<<  RPD_startY - (j*( tileY+(2*halfY_gap) ) ) 
				<< ", Sz = " 
				<< RPD_centerZ 
				<< ", ("
				<< i 
				<< "," 
				<< j 
				<< ")" 
				<< std::endl << std::endl;	
      ++cn;
    }
  }
		
	
	
	
  



		
  //----------------------------------------------     
  // Define Surface/Border Properties
  //----------------------------------------------  

//  DefineBorderProperties();

  //----------------------------------------------     
  // SD and Scoring Volumes
  //----------------------------------------------
  

 G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Note one SD object for each module
  char quartzSDname[256];
  sprintf( quartzSDname, "RPD_SD%d", m_modNum+1);

 RpdSD* aRpdSD = new RpdSD( quartzSDname, m_sd, m_modNum );
  aRpdSD->HistInitialize();
  SDman->AddNewDetector( aRpdSD );
  m_tileLogical->SetSensitiveDetector( aRpdSD );

   char fiberSDname[256];
  sprintf( fiberSDname, "Fiber_SD%d", m_modNum+2);
  
  
  FiberSD* aFiberSD = new FiberSD( fiberSDname, m_sd, m_modNum+3 );
  aFiberSD->HistInitialize();
  SDman->AddNewDetector( aFiberSD );
  for(int k=0;k<4;k++) {					
		m_fiberLogical[k]->SetSensitiveDetector( aFiberSD );
  }
  
  

  
  
  
  
  
  
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
  std::cout << "ModTypeRPD construction finished: SD name " << quartzSDname << std::endl;

}
