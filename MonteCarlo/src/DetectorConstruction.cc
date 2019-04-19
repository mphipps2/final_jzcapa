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

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "PurgMagTabulatedField3D.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4PVParameterised.hh"
#include "G4ThreeVector.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"


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
	int ZDC_SETUP = config->GetValue("ZDC_SETUP", 0);
	int runNum = config->GetValue( "RunNumber", -1);
	
	m_sd->LoadConfigurationFile(runNum);
	m_sd->LoadAlignmentFile(runNum);

	// Create variables to be used in beamtest 2018 simulation
	
	G4ThreeVector zdc1Pos,zdc2Pos, rpdPos;
	bool 		mag_on = false, targ_in = false, lead_in = false, 
				ZDC1 = false, ZDC2 = false, RPD = false;
	G4double 	mag_zOffset=0, density=0, fractionmass=0, leadblockZ=0, worldZoffset=16000*mm,
				zdc1X=0,zdc1Y=0,zdc1Z=0,zdc2X=0,zdc2Y=0,zdc2Z=0,rpdX=0,rpdY=0,rpdZ=0,
				tableX_shift=0, tableY_shift=0;
	std::string detector[3];
	G4int   	ncomponents;  
	Survey		*srvy_zdc1 = m_sd->GetSurvey("ZDC1"), 
				*srvy_zdc2 = m_sd->GetSurvey("ZDC2"), 
				*srvy_rpd = m_sd->GetSurvey("RPD");;
	Alignment	*align_run 	= m_sd->GetAlignment();

	//################################ SURVEY/ALIGNMENT_SETUP
	if(ZDC_SETUP == 1){
	
	//table(-2250,500) -> rpd/beam(0,0)	where 100=0.1cm in table coordinates
	tableX_shift = (-2250.0 - (align_run->x_table)  )/100*mm ;
	tableY_shift = (500.0   - (align_run->y_table)  )/100*mm ;
		
	rpdX  = (srvy_rpd ->x_pos)   *1000.0			 					*mm;	
	zdc1X = (( (srvy_zdc1->x_pos)*1000.0 ) - rpdX + tableX_shift    )	*mm;
	zdc2X = (( (srvy_zdc2->x_pos)*1000.0 ) - rpdX + tableX_shift	)	*mm;
	rpdX  =  tableX_shift;
		
	rpdY  = (srvy_rpd ->y_pos)   *1000.0								*mm;	
	zdc1Y = (( (srvy_zdc1->y_pos)*1000.0 ) - 320 - rpdY + tableY_shift)	*mm;
	zdc2Y = (( (srvy_zdc2->y_pos)*1000.0 ) - 320 - rpdY + tableY_shift)	*mm;
	rpdY  =  tableY_shift;
	
	rpdZ  = (( (srvy_rpd ->z_pos)*1000.0 ) -(worldZoffset) )*mm;
	zdc1Z = (( (srvy_zdc1->z_pos)*1000.0 ) -(worldZoffset) )*mm;
	zdc2Z = (( (srvy_zdc2->z_pos)*1000.0 ) -(worldZoffset) )*mm;
	
	 
	//-320mm is offset to get from zdc mid to active area mid
	zdc1Pos = G4ThreeVector( zdc1X, zdc1Y, zdc1Z); 
	zdc2Pos = G4ThreeVector( zdc2X, zdc2Y, zdc2Z);
	rpdPos 	= G4ThreeVector( rpdX , rpdY , rpdZ );
  
	
	detector[0]=align_run->upstream_Det;
	detector[1]=align_run->mid_Det;
	detector[2]=align_run->downstream_Det;
	
	for(int i=0; i<3; i++){
		if(detector[i]=="ZDC1") {
			ZDC1 = true;
		}
		if(detector[i]=="ZDC2") {
			ZDC2 = true;
		}
		if(detector[i]=="RPD") {
			RPD = true;
		}
	}
	
	// Assign lead block position in mm
	leadblockZ = ( zdc1Z- 250)*mm; //approximate position
	
	
	targ_in =  	align_run->target_In;
	mag_on 	= 	align_run->magnet_On;
	lead_in = 	align_run->lead_In;
	}
	//################################ SURVEY/ALIGNMENT_END
	
	
  G4NistManager* nist = G4NistManager::Instance();
		
		G4Element* Pb = nist->FindOrBuildElement(82);
		G4Material* Lead = new G4Material("Lead", density= 11.35*g/cm3, ncomponents=1);
		Lead->AddElement(Pb, fractionmass=1.0);  
	
  // Get nist material manager
  // G4NistManager* nist = G4NistManager::Instance();

  //----------------------------------------------     
  // Set Some Values
  //----------------------------------------------
  G4double modSizeX[5]; G4double modSizeY[5]; G4double modSizeZ[5]; G4double modCasingThickness[5]; G4double modWidth[5]; G4String modAbsorberMat[5]; G4double modAbsorberThickness[5]; G4double modAbsorberHeight[5]; G4double modAbsorberWidth[5]; G4double modCoreDiameter[5]; G4double modCladdingThickness[5]; G4int modNRadiators[5]; G4int modNAbsorbers[5]; G4int modType[5]; bool cladding[5]; G4double modRadiatorGapLength[5]; G4double stripPitch[5]; G4int modNStripsPerGap[5];
  
  G4int    nModules            = config->GetValue( "nModules",2);
  modType[0]            = config->GetValue( "mod1Type",5);
  modType[1]            = config->GetValue( "mod2Type",5);
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
	if(ZDC_SETUP==1){
		if (modType[i] == 5) {
		char variable[256];
		std::string modCladding;
		sprintf(variable,"mod%dCasingThickness",5);
		modCasingThickness[i] = config->GetValue(variable,.605);
		sprintf(variable,"mod%dAbsorberThickness",5);
		modAbsorberThickness[i] = config->GetValue(variable,10.);
		sprintf(variable,"mod%dAbsorberHeight",5);
		modAbsorberHeight[i] = config->GetValue(variable,180.);
		sprintf(variable,"mod%dAbsorberWidth",5);
		modAbsorberWidth[i] = config->GetValue(variable,100.);
		sprintf(variable,"mod%dRadiatorGapLength",5);
		modRadiatorGapLength[i] = config->GetValue(variable,2.);
		sprintf(variable,"mod%dCoreDiameter",5);
		modCoreDiameter[i] = config->GetValue(variable,1.5);
		sprintf(variable,"mod%dCladding",5);
		modCladding = config->GetValue(variable,"true");
		sprintf(variable,"mod%dCladdingThickness",5);
		modCladdingThickness[i] = config->GetValue(variable,0.605);
		sprintf(variable,"mod%dNRadiators",5);
		modNRadiators[i] = config->GetValue(variable,12);
		sprintf(variable,"mod%dNStripsPerGap",5);
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
	}
    else {
      modSizeX[i] = 90.78;
      modSizeY[i] = 200.;
      modSizeZ[i] = 150.;
    }
  }

  
  // Option to switch on/off checking of volumes overlaps
  //
  bool checkOverlaps = false;
  
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

  G4double worldSizeX       = 1.1 * maxModSizeX * mm;    	// mm //1.1
  G4double worldSizeY       = 1.1 * maxModSizeY * mm;    	// mm //1.1
  G4double worldSizeZ       = 1.2 * totalModSizeZ * mm; 	//was 1.2
  
   
   if(ZDC_SETUP==1){
		worldSizeZ= 32000*mm;
		//worldZoffset=(worldSizeZ/2)*mm;
		//320mm is offset to get from zdc mid to active area mid
		if( std::abs(zdc1X) < std::abs(zdc2X) ) worldSizeX = 1.1 * 2 * ( std::abs(zdc2X)+ maxModSizeX/2 ) * mm;
		else  															worldSizeX = 1.1 * 2 * ( std::abs(zdc1X) + maxModSizeX/2 ) * mm;
		if( std::abs(zdc1Y) < std::abs(zdc2Y) ) worldSizeY = 1.1 * 2 * ( std::abs(zdc2Y) + maxModSizeY/2  ) * mm; 
		else  															worldSizeY = 1.1 * 2 * ( std::abs(zdc1Y)+ maxModSizeY/2  ) * mm;
   }

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
  
  
  G4VisAttributes* boxVisAtt_world= new G4VisAttributes(G4Colour(0.5,0.5,0.5)); 

		m_logicWorld ->SetVisAttributes(boxVisAtt_world);
  //----------------------------------------------     
  // Build ZDC Setup
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
  
  
  if(ZDC_SETUP==1){

	  std::cout << "(" << zdc1Pos.getX() << ", " << zdc1Pos.getY() << ", " << zdc1Pos.getZ() << ")" << std::endl;
	  std::cout << "(" << zdc2Pos.getX() << ", " << zdc2Pos.getY() << ", " << zdc2Pos.getZ() << ")" << std::endl;
	  std::cout << "(" << rpdPos.getX()  << ", " << rpdPos.getY()  << ", " << rpdPos.getZ()  << ")" << std::endl;
	  
	  ModTypeZDC *mod1 = new ModTypeZDC(0,zdc1Pos,m_logicWorld,m_sd); 
	  ModTypeZDC *mod2 = new ModTypeZDC(1,zdc2Pos,m_logicWorld,m_sd); 
	  ModTypeRPD *mod3 = new ModTypeRPD(0,rpdPos, m_logicWorld,m_sd); 
	  
	 if(ZDC1) mod1->Construct();
	 if(ZDC2) mod2->Construct();
	 if(RPD)  mod3->Construct();
		
	  
		// Setup magnetic field
	if( mag_on )
    {
		mag_zOffset=-9.55*m;
		//Field grid in A9.TABLE. File must be accessible from run directory. 
		G4MagneticField* PurgMagField= new PurgMagTabulatedField3D((std::getenv("JZCaPA") + std::string("/Utils/PurgMag3D.TABLE")).c_str(), mag_zOffset+(worldZoffset/1000.0));
		fField.Put(PurgMagField);
      
		//This is thread-local
		G4FieldManager* pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
           
		//G4cout<< "DeltaStep "<<pFieldMgr->GetDeltaOneStep()/mm <<"mm" <<endl;
		//G4ChordFinder *pChordFinder = new G4ChordFinder(PurgMagField);

		pFieldMgr->SetDetectorField(fField.Get());
		pFieldMgr->CreateChordFinder(fField.Get());
	}
		// Setup lead target
   if( targ_in ){
		G4Box* leadTarget = new G4Box("Target",5*cm, 5*cm, 1.3*cm);
    
		logic_leadTarget
		= new G4LogicalVolume(leadTarget,
                          Lead,
                          "LeadTarget");
    
		new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, (2600-worldZoffset)*mm),
                      logic_leadTarget,
                      "LeadTarget1",
                      m_logicWorld,
                      false,
                      0
                      );

		// Visualization attributes
		G4VisAttributes* boxVisAtt_lead= new G4VisAttributes(G4Colour(1.0,0.0,1.0)); //magenta
		logic_leadTarget ->SetVisAttributes(boxVisAtt_lead);
		
   }
	if( lead_in ){  
		G4Box* leadBlock = new G4Box("LeadBlock",10*cm, 10*cm, 20*cm);
    
		logic_leadBlock
		= new G4LogicalVolume(leadBlock,
                          Lead,
                          "LeadBlock");
    
		new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, (leadblockZ-worldZoffset)*mm),
                      logic_leadBlock,
                      "LeadBlock1",
                      m_logicWorld,
                      false,
                      0
                      );

		// Visualization attributes
		G4VisAttributes* boxVisAtt_leadblk= new G4VisAttributes(G4Colour(1.0,0.0,1.0)); //magenta
		logic_leadBlock ->SetVisAttributes(boxVisAtt_leadblk);
		cout << "Placed Lead Block at Z = " << leadblockZ*mm << "mm" <<  std::endl;
	}
	  
  }//END_ZDC_SETUP
  
  return m_physWorld;
}
