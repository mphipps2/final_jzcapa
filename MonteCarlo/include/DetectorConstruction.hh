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
// $Id: DetectorConstruction.hh 69565 2013-05-08 12:35:31Z gcosmo $
//
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "ModTypeZDC.hh"
#include "ModTypeRPD.hh"
#include "G4Cache.hh"
#include "G4MagneticField.hh"
#include "XMLSettingsReader.hh"
#include <vector>
#include "Materials.hh"

class G4Box;
class G4Para;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

class G4UniformMagField;
class PurgMagTabulatedField3D;


class Survey {

 public :

	/** Type of detector - ZDC or RPD **/
    std::string detector;
    /** x_pos for this survey */
    double x_pos;
	/** y_pos for this survey */
    double y_pos;
	/** z_pos for this survey */
    double z_pos;
	/** cos_x for this survey */
    double cos_x;
	/** cos_y for this survey */
    double cos_y;
	/** cos_z for this survey */
    double cos_z;

};

class Alignment {

 public:

    /** X position of the Desy Table **/
    double x_table;
    /** Y position of the Desy Table **/
    double y_table;
    /** First detector met by the beam **/
    std::string upstream_Det;
    /** Second detector met by the beam **/
    std::string mid_Det;
    /** Third detector met by the beam **/
    std::string downstream_Det;
    /** GOLIATH magnet status **/
    bool magnet_On;
    /** Target in **/
    bool target_In;
    /** Lead absorber in **/
    bool lead_In;

};

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  virtual void               DefineBorderProperties();

  virtual G4VPhysicalVolume* ConstructDetector();
  virtual G4VPhysicalVolume* ConstructTestBeam();
  virtual G4VPhysicalVolume* LoadConfigurationFile();
  virtual Survey*            GetSurvey();
  virtual Alignment*         LoadAlignmentFile();

  inline  G4int              GetnZDCs(){return ZDCvec.size();}
  inline  G4int              GetnRPDs(){return RPDvec.size();}
  inline  G4bool             SetClusterFlag(G4bool arg){CLUSTER = arg;}
  inline  G4bool             GetClusterFlag(){return CLUSTER;}


protected:
  Materials*              materials;

protected:
  /* World objects */
  G4Box*                  m_solidWorld;
  G4LogicalVolume*        m_logicWorld;
  G4VPhysicalVolume*      m_physWorld;

  /* Lead target objects */
  G4LogicalVolume*        logic_leadTarget;
  G4LogicalVolume*        logic_leadBlock;

  /* Configuration parser objects and storage */
  XMLSettingsReader*      m_XMLparser;
  Alignment* 	            m_alignment;
  Survey*                 m_survey;
  std::vector < Survey* > surveyEntries;

  /* Number of each detector */
  std::vector< ModTypeZDC* > ZDCvec;
  std::vector< ModTypeRPD* > RPDvec;

  /* Cluster Flag for forwarding to RunAction */
  G4bool CLUSTER;



private:
	G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
