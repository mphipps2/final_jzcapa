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
  : G4VUserDetectorConstruction(),
    m_solidWorld(NULL), m_logicWorld(NULL), m_physWorld(NULL)
{materials = Materials::getInstance();}

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
  ConstructTestBeam()
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructTestBeam()
{

// Create variables to be used in beamtest 2018 simulation

  G4ThreeVector zdc1Pos,zdc2Pos, rpdPos;
  bool mag_on = targ_in = lead_in = ZDC1 = ZDC2 = RPD = false;
  G4double mag_zOffset=0, density=0, fractionmass=0, leadblockZ=0,
           zdc1X=0, zdc1Y=0, zdc1Z=0, zdc2X=0, zdc2Y=0, zdc2Z=0, rpdX=0,
           rpdY=0, rpdZ=0, tableX_shift=0, tableY_shift=0;
  G4double worldZoffset=16000*mm;
  std::string detector[3];
  G4int  ncomponents;
  LoadConfigurationFile(runNum);
  LoadAlignmentFile(runNum);
  Survey *srvy_zdc1 = GetSurvey("ZDC1"),
         *srvy_zdc2 = GetSurvey("ZDC2"),
         *srvy_rpd  = GetSurvey("RPD");;
  Alignment *m_alignment = GetAlignment();

//################################ SURVEY/ALIGNMENT_SETUP

  if( ForceDetectorPosition ) == 0){

    std::cout << "******************************************" << std::endl
              << "        PLACING DETECTORS MANUALLY        " << std::endl
              << "******************************************" << std::endl;

    rpdX = srvy_rpd->x_pos*mm;
    rpdY = srvy_rpd->y_pos*mm;
    rpdZ = srvy_rpd->z_pos*mm;

    zdc1X = srvy_zdc1->x_pos*mm;
    zdc1Y = srvy_zdc1->y_pos*mm;
    zdc1Z = srvy_zdc1->z_pos*mm;

    zdc2X = srvy_zdc2->x_pos*mm;
    zdc2Y = srvy_zdc2->y_pos*mm;
    zdc2Z = srvy_zdc2->z_pos*mm;
  } else{
    //table(-2250,500) -> rpd/beam(0,0)	where 100=0.1cm in table coordinates
    //-320mm is offset to get from zdc mid to active area mid
    tableX_shift = (-2250.0 - (m_alignment->x_table)  )/100*mm ;//2257 more accurate
    tableY_shift = (500.0   - (m_alignment->y_table)  )/100*mm ;//501  more accurate

    rpdX  = (srvy_rpd ->x_pos)   *1000.0*mm;
    zdc1X = (( (srvy_zdc1->x_pos)*1000.0 ) - rpdX + tableX_shift )*mm;
    zdc2X = (( (srvy_zdc2->x_pos)*1000.0 ) - rpdX + tableX_shift )*mm;
    rpdX  =  tableX_shift;

    rpdY  = (srvy_rpd ->y_pos)   *1000.0*mm;
    zdc1Y = (( (srvy_zdc1->y_pos)*1000.0 ) - 320 - rpdY + tableY_shift )*mm;
    zdc2Y = (( (srvy_zdc2->y_pos)*1000.0 ) - 320 - rpdY + tableY_shift )*mm;
    rpdY  =  tableY_shift;

    rpdZ  = (( (srvy_rpd ->z_pos)*1000.0 ) -(worldZoffset) )*mm;
    zdc1Z = (( (srvy_zdc1->z_pos)*1000.0 ) -(worldZoffset) )*mm;
    zdc2Z = (( (srvy_zdc2->z_pos)*1000.0 ) -(worldZoffset) )*mm;
  }


  zdc1Pos = G4ThreeVector( zdc1X, zdc1Y, zdc1Z);
  zdc2Pos = G4ThreeVector( zdc2X, zdc2Y, zdc2Z);
  rpdPos  = G4ThreeVector( rpdX , rpdY , rpdZ );


  detector[0]=m_alignment->upstream_Det;
  detector[1]=m_alignment->mid_Det;
  detector[2]=m_alignment->downstream_Det;

  for(int i=0; i<3; i++){
    if(detector[i]=="ZDC1") {ZDC1 = true;}
    if(detector[i]=="ZDC2") {ZDC2 = true;}
    if(detector[i]=="RPD" ) {RPD  = true;}
  }

  // Assign lead block position in mm, place approx 1 ft in front of ZDC1
  // half ZDC Z width = 90mm
  // 1ft ~ 300mm
  leadblockZ = ( zdc1Z - 90 - 300)*mm;


  targ_in = m_alignment->target_In;
  mag_on  = m_alignment->magnet_On;
  lead_in = m_alignment->lead_In;
  }
  //################################ SURVEY/ALIGNMENT_END


  // Option to switch on/off checking of volumes overlaps
  //
  bool checkOverlaps = false;

   if(TESTBEAM_SETUP==1){
     worldSizeZ= 32000*mm;
     if( std::abs(zdc1X) < std::abs(zdc2X) ) worldSizeX = 1.1 * 2 * ( std::abs(zdc2X) + maxModSizeX/2 )*mm;
     else                                    worldSizeX = 1.1 * 2 * ( std::abs(zdc1X) + maxModSizeX/2 )*mm;
     if( std::abs(zdc1Y) < std::abs(zdc2Y) ) worldSizeY = 1.1 * 2 * ( std::abs(zdc2Y) + maxModSizeY/2 )*mm;
     else                                    worldSizeY = 1.1 * 2 * ( std::abs(zdc1Y) + maxModSizeY/2 )*mm;
        if(!ZDC1 && !ZDC2) {
          worldSizeX = 1.1 * 2 * ( std::abs(rpdX)+ maxModSizeX/2 ) * mm;
          worldSizeY = 1.1 * 2 * ( std::abs(rpdY) + maxModSizeY/2  ) * mm;
        }
   }

  G4Material* g4Air = nist->FindOrBuildMaterial("G4_AIR");

  //Air
  if (config->GetValue("OPTICAL_ON",false) == 1){
    materials->UseOpticalMaterials(true); //set this an an option later ARIC!
    materials->DefineOpticalProperties();
    g4Air = materials->Air;
  }

  printf( "Building world with x %5.1f y %5.1f z %5.1f\n",
          worldSizeX, worldSizeY, worldSizeZ );

  m_solidWorld =
    new G4Box("World",   //its name
        0.5*worldSizeX,  //its size
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

  ModTypeZDC *mod1 = new ModTypeZDC(0,zdc1Pos,m_logicWorld,m_sd);
  ModTypeZDC *mod2 = new ModTypeZDC(1,zdc2Pos,m_logicWorld,m_sd);
  ModTypeRPD *mod3 = new ModTypeRPD(0,rpdPos, m_logicWorld,m_sd);

  if(ZDC1){ mod1->Construct();
    std::cout << "ZDC1 center = " << "(" << zdc1Pos.getX() << ", " << zdc1Pos.getY() << ", " << zdc1Pos.getZ() << ")" << std::endl;
  }
  if(ZDC2){ mod2->Construct();
    std::cout << "ZDC2 center = " << "(" << zdc2Pos.getX() << ", " << zdc2Pos.getY() << ", " << zdc2Pos.getZ() << ")" << std::endl;
  }
  if(RPD){  mod3->Construct();
    std::cout << "RPD center = " << "(" << rpdPos.getX()  << ", " << rpdPos.getY()  << ", " << rpdPos.getZ()  << ")" << std::endl;
  }


  // Setup magnetic field
  if( mag_on ){
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
     G4Box* leadBlock = new G4Box("LeadBlock",(0.5*worldSizeX)*mm, (0.5*150)*mm, 10*cm);

     logic_leadBlock
     = new G4LogicalVolume(leadBlock,
                           Lead,
                           "LeadBlock");

    new G4PVPlacement(0,
                      G4ThreeVector(0.0, 0.0, (leadblockZ)*mm),
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

  }//END_TESTBEAM_SETUP

  return m_physWorld;
}





//***********************************************************************************
// READ Survey_2018.xml and Alignment_2018.xml
//***********************************************************************************
/**
 * @brief Reads the .xml configuration file and load characteristics for all the alignments, immediately sorted into alignment objects
 * @param _inFile
 */

void DetectorConstruction::LoadConfigurationFile( int m_runNumber, std::string _inFile  ){

    m_XMLparser = new XMLSettingsReader();


    if (!m_XMLparser->parseFile(_inFile)) {
      if(!m_XMLparser->parseFile("Survey_2018.xml")){
        std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
        return;
      }
      std::cout << " Found Survey_2018.xml in executable directory " << std::endl;
    }

    std::cout << "Loading .xml Configuration File..." << std::endl;
    std::cout << "Found " << m_XMLparser->getBaseNodeCount("Survey") << " survey entries " << std::endl;

    int first_run, last_run;

    for (int i = 0; i < m_XMLparser->getBaseNodeCount("Survey"); i++) { //this was unsigned int i = 0

        m_XMLparser->getChildValue("Survey",i,"start_run",first_run);
        m_XMLparser->getChildValue("Survey",i,"end_run",last_run);

        //Discard entries for any alignment that does not apply to our run
        if(m_runNumber < first_run || m_runNumber > last_run) continue;
          m_survey = new Survey();

        //If the entry applies, we store it in the vector
        m_XMLparser->getChildValue("Survey",i,"detector",m_survey->detector);
        m_XMLparser->getChildValue("Survey",i,"x_pos",m_survey->x_pos);
        m_XMLparser->getChildValue("Survey",i,"y_pos",m_survey->y_pos);
        m_XMLparser->getChildValue("Survey",i,"z_pos",m_survey->z_pos);
        m_XMLparser->getChildValue("Survey",i,"cos_x",m_survey->cos_x);
        m_XMLparser->getChildValue("Survey",i,"cos_y",m_survey->cos_y);
        m_XMLparser->getChildValue("Survey",i,"cos_z",m_survey->cos_z);

        surveyEntries.push_back(m_survey);
    }

    if(surveyEntries.size() == 0) std::cout << "WARNING: SURVEY NOT FOUND!!!" << std::endl;

    return;
}

/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void DetectorConstruction::LoadAlignmentFile( int m_runNumber, std::string _inFile ){

  bool debug = false;

    m_XMLparser = new XMLSettingsReader();

    if (!m_XMLparser->parseFile(_inFile)) {
        if(!m_XMLparser->parseFile("Alignment_2018.xml")){
          std::cerr << " Data Reader could not parse file : " << _inFile << std::endl;
          return;
        }
        std::cout << " Found Alignment_2018.xml in executable directory " << std::endl;
    }

    m_alignment = new Alignment();

    if(debug){
      std::cout << "Loading .xml Alignment File..." << std::endl;
      std::cout << "Found " << m_XMLparser->getBaseNodeCount("Alignment") << " alignment entries " << std::endl;
      std::cout << "Retrieving the information for run " << m_runNumber << std::endl;
    }

    int run;
    for ( int i = 0; i < m_XMLparser->getBaseNodeCount("Alignment"); i++) {
      m_XMLparser->getChildValue("Alignment",i,"run",run);
      if(run != m_runNumber) continue;
        if(debug){
          std::cout << "Found Run Entry in Alignment file for run " << m_runNumber << std::endl;
        }
        m_XMLparser->getChildValue("Alignment",i,"x_table",m_alignment->x_table);
        m_XMLparser->getChildValue("Alignment",i,"y_table",m_alignment->y_table);
        m_XMLparser->getChildValue("Alignment",i,"upstream_Det",m_alignment->upstream_Det);
        m_XMLparser->getChildValue("Alignment",i,"mid_Det",m_alignment->mid_Det);
        m_XMLparser->getChildValue("Alignment",i,"downstream_Det",m_alignment->downstream_Det);
        m_XMLparser->getChildValue("Alignment",i,"target_In",m_alignment->target_In);
        m_XMLparser->getChildValue("Alignment",i,"lead_In",m_alignment->lead_In);
        m_XMLparser->getChildValue("Alignment",i,"magnet_On",m_alignment->magnet_On);
    }

    if(m_alignment == NULL) std::cout << "WARNING: ALIGNMENT NOT FOUND!!!" << std::endl;
    return;
}

Survey* DetectorConstruction::GetSurvey(std::string name){
  for(unsigned int i = 0; i < surveyEntries.size(); i++){
    if( name == surveyEntries[i]->detector ){ return surveyEntries[i]; }
  }
  Survey* empty=NULL;
  return empty;
}

Alignment* DetectorConstruction::GetAlignment(){
  return m_alignment;
}
