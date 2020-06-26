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
#include "DetectorMessenger.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4Box.hh"
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
  : G4VUserDetectorConstruction(), m_solidWorld(NULL), m_logicWorld(NULL),
  m_physWorld(NULL),OPTICAL(false), ForceDetectorPosition(false)
{
  new DetectorMessenger(this);
  currentRPD = -1;
  currentZDC = -1;
  m_materials = Materials::getInstance();
  m_materials->UseOpticalMaterials(true);
  m_materials->DefineOpticalProperties();
}

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

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ExecuteMacroFile("geometry.mac");

  if( ForceDetectorPosition ){
    ManualConstruction();
  }else{
    ConstructSPSTestBeam();
    std::cout << "built test beam" << std::endl;
  }
  return m_physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructWorldVolume(G4double x, G4double y, G4double z){
  printf( "Building world with x %5.1f y %5.1f z %5.1f\n", x, y, z );

  m_solidWorld =
    new G4Box("World", 0.5*x, 0.5*y, 0.5*z );

  m_logicWorld =
    new G4LogicalVolume(m_solidWorld,     //its solid
                        m_materials->Air, //its material
                        "World");         //its name

  m_physWorld =
    new G4PVPlacement(0,                  //no rotation
                      G4ThreeVector(),    //at (0,0,0)
                      m_logicWorld,       //its logical volume
                      "World",            //its name
                      0,                  //its mother  volume
                      false,              //no boolean operation
                      0,                  //copy number
                      false);             //overlaps checking


  G4VisAttributes* boxVisAtt_world = new G4VisAttributes(G4Colour(0.5,0.5,0.5));

  m_logicWorld ->SetVisAttributes(boxVisAtt_world);

  return m_physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ManualConstruction(){

  std::cout << "******************************************" << std::endl
            << "        PLACING DETECTORS MANUALLY        " << std::endl
            << "******************************************" << std::endl;

  G4ThreeVector* pos;
  for(ModTypeZDC* zdc : m_ZDCvec){
    pos = zdc->GetPosition();
    printf( "ZDC%d center = (%f,%f,%f)\n", zdc->GetModNum(), pos->x(), pos->y(), pos->z() );
    zdc->Construct();
  }

  for(ModTypeRPD* rpd : m_RPDvec){
    pos = rpd->GetPosition();
    printf( "RPD%d center = (%f,%f,%f)\n", rpd->GetModNum(), pos->x(), pos->y(), pos->z() );
    rpd->Construct();
  }

  return m_physWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructSPSTestBeam(){
  // Create variables to be used in beamtest 2018 simulation
  G4ThreeVector zdc1Pos,zdc2Pos, rpdPos;
  bool ZDC1 = false, ZDC2 = false, RPD = false;
  G4double mag_zOffset=0, firstDetZ, detX, detY, detZ,
           tableX_shift=0, tableY_shift=0;
  G4Material* Lead = m_materials->Pb;
  G4double worldSizeX = 180*mm;
  G4double worldSizeY = 1200*mm;
  G4double worldSizeZ = 32000*mm;
  G4double worldZoffset= worldSizeZ/2.0;
  std::string detector[3];

  //################################ World volume construction

  ConstructWorldVolume( worldSizeX, worldSizeY, worldSizeZ );

  //################################ SURVEY/ALIGNMENT_SETUP

  LoadConfigurationFile();
  LoadAlignmentFile();

  //table(-2250,500) -> rpd/beam(0,0)	where 100=0.1cm in table coordinates
  //-320mm is offset to get from zdc mid to active area mid
  tableX_shift = (-2250.0 - (m_alignment->x_table)  )/100*mm ;//2257 more accurate
  tableY_shift = (500.0   - (m_alignment->y_table)  )/100*mm ;//501  more accurate
  Survey *srvy_rpd = GetSurvey("RPD");

  detector[0]=m_alignment->upstream_Det;
  detector[1]=m_alignment->mid_Det;
  detector[2]=m_alignment->downstream_Det;

  for(int i=0; i<3; i++){
    if(detector[i]=="ZDC1") {ZDC1 = true;}
    if(detector[i]=="ZDC2") {ZDC2 = true;}
    if(detector[i]=="RPD" ) {RPD  = true;}
  }

  firstDetZ = worldSizeZ;
  for(Survey* survey : m_surveyEntries){
    detX = (survey->x_pos*1000.0)*mm;
    detY = (survey->y_pos*1000.0)*mm;
    detZ = (survey->z_pos*1000.0 - worldZoffset )*mm;

    if(detZ < firstDetZ) firstDetZ = detZ;

    if( survey->detector == "ZDC1" && ZDC1 ){
      AddZDC( new G4ThreeVector( detX + ( tableX_shift       - srvy_rpd->x_pos )*mm,
                                 detY + ( tableY_shift - 320 - srvy_rpd->y_pos )*mm,
                                 detZ ) );
    } else if( survey->detector == "ZDC2" && ZDC2 ){
      AddZDC( new G4ThreeVector( detX + ( tableX_shift       - srvy_rpd->x_pos )*mm,
                                 detY + ( tableY_shift - 320 - srvy_rpd->y_pos )*mm,
                                 detZ ) );
    }else if( survey->detector == "RPD" && RPD ){
      AddRPD( new G4ThreeVector( tableX_shift, tableY_shift, detZ) );
    }
  }

  G4ThreeVector* pos;
  G4int modNum;
  for(ModTypeZDC* zdc : m_ZDCvec){
    zdc->SetFiberDiameters    ( new G4ThreeVector(1.5*mm,0.0,0.0) );
    zdc->SetAbsorberDimensions( new G4ThreeVector(90.0*mm, 180.0*mm, 11.0*mm) );
    zdc->SetnAbsorbers        ( 11 );
    zdc->SetHousingThickness  ( 4.5*mm );
    zdc->SetGapThickness      ( 2.5*mm );
    zdc->SetHousingMaterial   ( "aluminum" );
    zdc->SetAbsorberMaterial  ( "pure" );

    pos = zdc->GetPosition();
    modNum = zdc->GetModNum();
    printf( "ZDC%d center = (%f,%f,%f)\n", modNum, pos->x(), pos->y(), pos->z() );
    zdc->Construct();
  }

  for(ModTypeRPD* rpd : m_RPDvec){
    rpd->SetDetectorType  ( "cms" );

    pos = rpd->GetPosition();
    modNum = rpd->GetModNum();
    printf( "RPD%d center = (%f,%f,%f)\n", modNum, pos->x(), pos->y(), pos->z() );
    rpd->Construct();
  }
//################################ Get SURVEY/ALIGNMENT_END


  // Setup magnetic field
  if( m_alignment->magnet_On ){
    mag_zOffset=-9.55*m;
    //Field grid in A9.TABLE. File must be accessible from run directory.
    G4MagneticField* PurgMagField= new PurgMagTabulatedField3D((std::getenv("JZCaPA") + std::string("/Utils/PurgMag3D.TABLE")).c_str(), mag_zOffset+(worldZoffset/1000.0));
    fField.Put(PurgMagField);

    //This is thread-local
    G4FieldManager* pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    pFieldMgr->SetDetectorField(fField.Get());
    pFieldMgr->CreateChordFinder(fField.Get());
  }
    // Setup lead target
    if( m_alignment->target_In ){
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
   if( m_alignment->lead_In ){
     // Assign lead block position in mm, place approx 1 ft in front of ZDC1
     // half ZDC Z width = 90mm
     // 1ft ~ 300mm
     G4double leadblockZ = ( firstDetZ - 90 - 300)*mm;

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
  return m_physWorld;
}





//***********************************************************************************
// READ Survey_2018.xml and Alignment_2018.xml
//***********************************************************************************
/**
 * @brief Reads the .xml configuration file and load characteristics for all the alignments, immediately sorted into alignment objects
 * @param _inFile
 */

void DetectorConstruction::LoadConfigurationFile( G4String _inFile  ){

  if(_inFile = ""){
    _inFile = std::getenv("JZCaPA");
    _inFile += "/Utils/Survey_2018.xml";
  }

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

        m_surveyEntries.push_back(m_survey);
    }

    if(m_surveyEntries.size() == 0) std::cout << "WARNING: SURVEY NOT FOUND!!!" << std::endl;

    return;
}

/**
 * @brief Reads the .xml configuration file and load characteristics for all the channels, immediately sorted into detectors objects
 * @param _inFile
 */
void DetectorConstruction::LoadAlignmentFile( G4String _inFile ){
  bool debug = false;

  if( _inFile = ""){
    _inFile = std::getenv("JZCaPA");
    _inFile += "/Utils/Alignment_2018.xml";
  }

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

/*
 * Retrieve a survey based on name
*/
Survey* DetectorConstruction::GetSurvey(G4String name){
  for(unsigned int i = 0; i < m_surveyEntries.size(); i++){
    if( m_surveyEntries[i]->detector.compare( name.c_str() ) ){ return m_surveyEntries[i]; }
  }
  Survey* empty=NULL;
  return empty;
}

/*
 * Add a ZDC module
*/
void DetectorConstruction::AddZDC(G4ThreeVector* position){
  uint newModNum = m_ZDCvec.size()+1;
  m_ZDCvec.push_back(new ModTypeZDC(newModNum, m_logicWorld, position ));
  currentZDC = newModNum;
  printf("Added ZDC%d\n", newModNum);
}

/*
 * Add an RPD module
*/
void DetectorConstruction::AddRPD(G4ThreeVector* position){
  uint newModNum = m_RPDvec.size()+1;
  m_RPDvec.push_back(new ModTypeRPD(newModNum, m_logicWorld, position ));
  currentRPD = newModNum;
  printf("Added RPD%d\n", newModNum);
}

/*
 * Duplicate a ZDC module
*/
void DetectorConstruction::DuplicateZDC( G4int module ){
  uint newModNum = m_ZDCvec.size()+1;
  if((unsigned)module >= newModNum ){
    printf("\n\n Cannot duplicate. ZDC%d does not exist \n\n",module);
    return;
  }
  m_ZDCvec.push_back( new ModTypeZDC( newModNum, m_ZDCvec.at(module-1) ) );
  currentZDC = newModNum;
  printf("Duplicate ZDC%d built from ZDC%d\n", newModNum, module);
}

/*
 * Duplicate an RPD module
*/
void DetectorConstruction::DuplicateRPD( G4int module ){
  uint newModNum = m_RPDvec.size()+1;
  if((unsigned)module >= newModNum ){
    printf("\n\n Cannot duplicate. RPD%d does not exist \n\n",module);
    return;
  }
  m_RPDvec.push_back( new ModTypeRPD( newModNum, m_RPDvec.at(module-1) ) );
  currentRPD = newModNum;
  printf("Duplicate RPD%d built from RPD%d\n", newModNum, module);
}
