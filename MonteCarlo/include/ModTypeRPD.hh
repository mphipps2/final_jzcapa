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
/// \file ModTypeRPD.hh
/// \brief Implementation of the ModTypeRPD class
/// \author Aric Tate
/// \date February 2019

#ifndef ModTypeRPD_h
#define ModTypeRPD_h 1


#include "globals.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "Materials.hh"

#include <vector>


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;

/// Detector construction class to define materials and geometry.

class ModTypeRPD
{
public:
  ModTypeRPD(const int, ModTypeRPD*);
  ModTypeRPD(const int, G4LogicalVolume*, G4ThreeVector*);
  ~ModTypeRPD();



  virtual void  Construct();
  virtual void  ConstructCMSDetector();
  virtual void  ConstructPanFluteDetector();

  virtual void  DefineMaterials();


  inline  void  SetPosition          ( G4ThreeVector* vec ){ delete m_pos;       m_pos       = vec; }
  inline  void  SetFiberDiameters    ( G4ThreeVector* vec ){ delete m_fiberDiam; m_fiberDiam = vec; }
  inline  void  SetHousingThickness  ( G4double       arg ){ m_HousingThickness = arg; }
  inline  void  SetFiberPitchX       ( G4double       arg ){ m_fiberPitchX      = arg; }
  inline  void  SetFiberPitchZ       ( G4double       arg ){ m_fiberPitchZ      = arg; }
  inline  void  SetTileSize          ( G4double       arg ){ m_tileSize         = arg; }
  inline  void  SetMinWallThickness  ( G4double       arg ){ m_minWallThickness = arg; }
  inline  void  SetReadoutDistance   ( G4double       arg ){ READOUT = true; m_distanceToReadout = arg;}
  inline  void  SetDetectorType      ( G4String       arg ){ m_detType          = arg; }
  inline  void  SetOpticalFlag       ( G4bool         arg ){ OPTICAL            = arg; }
  inline  void  SetOverlapsFlag      ( G4bool         arg ){ CHECK_OVERLAPS     = arg; }


  inline  G4ThreeVector* GetPosition ( ){ return m_pos;    }
  inline  G4int          GetModNum   ( ){ return m_modNum; }

protected:
  const G4int      m_modNum;
  G4ThreeVector*   m_pos;
  G4ThreeVector*   m_fiberDiam;
  G4double         m_HousingThickness;
  G4double         m_fiberPitchX;
  G4double         m_fiberPitchZ;
  G4double         m_tileSize;
  G4double         m_minWallThickness;
  G4double         m_distanceToReadout;
  G4String         m_detType;
  G4bool           OPTICAL;
  G4bool           CHECK_OVERLAPS;
  G4bool           READOUT;
  Materials*       materials;
  G4LogicalVolume* m_logicMother;


protected:

  G4Material*           m_matQuartz;
  G4Material*           m_Al;
  G4Material*           m_Poly;
  G4Material*           m_Air;
  G4Material*           m_PMMA;
  G4Material*           m_Grease;

  G4OpticalSurface*   m_photonDetSurface;

  G4SubtractionSolid* m_tile;



//PAN FLUTE START---------------------------------
  G4VSolid*           m_PFdetec;
  G4LogicalVolume*    m_PFdetecLogical[256];
  G4VPhysicalVolume*  m_PFdetecPhysical[256];


  G4VSolid*           m_PFrpd_housing;
  G4LogicalVolume*    m_PFrpd_housingLogical;
  G4VPhysicalVolume*  m_PFrpd_housingPhysical;

  G4VSolid*           m_PFreadout_air;
  G4LogicalVolume*    m_PFreadout_airLogical;
  G4VPhysicalVolume*  m_PFreadout_airPhysical;


  std::vector< G4VSolid* >                m_PFrpd;
  std::vector< G4LogicalVolume* >         m_PFrpdLogical;
  std::vector< G4VPhysicalVolume* >       m_PFrpdPhysical;

  std::vector< G4VSolid* >                m_PFreadout_fiber;
  std::vector< G4LogicalVolume* >         m_PFreadout_fiberLogical;
  std::vector< G4VPhysicalVolume* >       m_PFreadout_fiberPhysical;

  std::vector< G4VSolid* >                m_PFrpd_channel;
  std::vector< G4SubtractionSolid* >      m_subChannel;
  std::vector< G4LogicalVolume* >         m_PFrpd_channelLogical;
  std::vector< G4VPhysicalVolume* >       m_PFrpd_channelPhysical;

  std::vector< G4AssemblyVolume* >        m_PFrpd_FiberAssy;

  G4int               m_PFrpd_cnt;

//PAN FLUTE STOP----------------------------------



//CMS RPD DESIGN START----------------------------------
//STATIC ARRAYS

  G4VSolid*           m_tile_no_fiber_hole;
  G4VSolid*           m_fiber_subtract;
  G4LogicalVolume*    m_tileLogical;
  G4VPhysicalVolume*  m_tilePhysical[4][4];


  G4VSolid*          m_foilV[4];
  G4LogicalVolume*   m_foilVLogical[4];
  G4VPhysicalVolume* m_foilVPhysical[4][4];

  G4VSolid*          m_foilVfront;
  G4LogicalVolume*   m_foilVfrontLogical;
  G4VPhysicalVolume* m_foilVfrontPhysical[4][4];

  G4VSolid*          m_foilH;
  G4LogicalVolume*   m_foilHLogical;
  G4VPhysicalVolume* m_foilHPhysical[4][4];

  G4SubtractionSolid* m_foilHtop_hole;
  G4VSolid*          m_foilHtop;
  G4LogicalVolume*   m_foilHtopLogical;
  G4VPhysicalVolume* m_foilHtopPhysical[4];

  G4VSolid*          m_AlcaseV;
  G4LogicalVolume*   m_AlcaseVLogical;
  G4VPhysicalVolume* m_AlcaseVPhysical[5];

  G4VSolid*          m_Alcase;
  G4LogicalVolume*   m_AlcaseLogical;
  G4VPhysicalVolume* m_AlcasePhysical[2];

  G4VSolid*          m_fiber[4];
  G4LogicalVolume*   m_fiberLogical[64];
  G4VPhysicalVolume* m_fiberPhysical[64];

  G4VSolid*          m_fiberclad[4];
  G4LogicalVolume*   m_fibercladLogical[4];
  G4VPhysicalVolume* m_fibercladPhysical[64];

  G4VSolid*          m_fibergrease[4];
  G4LogicalVolume*   m_fibergreaseLogical[4];
  G4VPhysicalVolume* m_fibergreasePhysical[64];

  G4VSolid*          m_air_detect;
  G4LogicalVolume*   m_air_detect_Logical[64];
  G4VPhysicalVolume* m_air_detectPhysical[64];

//CMS RPD DESIGN STOP ----------------------------------

//TEST SETUP START -------------------------------------
//STATIC ARRAYS
  G4VSolid*           m_test_tile;
  G4LogicalVolume*    m_test_tileLogical;
  G4VPhysicalVolume*  m_test_tilePhysical;

  G4VSolid*           m_test_alum;
  G4LogicalVolume*    m_test_alumLogical;
  G4VPhysicalVolume*  m_test_alumPhysical;

  G4VSolid*           m_test_wls;
  G4LogicalVolume*    m_test_wlsLogical;
  G4VPhysicalVolume*  m_test_wlsPhysical;

  G4VSolid*           m_test_PD;
  G4LogicalVolume*    m_test_PDLogical;
  G4VPhysicalVolume*  m_test_PDPhysical;

  G4VSolid*           m_test_clad;
  G4LogicalVolume*    m_test_cladLogical;
  G4VPhysicalVolume*  m_test_cladPhysical[2];

  G4VSolid*           m_test_grease;
  G4LogicalVolume*    m_test_greaseLogical;
  G4VPhysicalVolume*  m_test_greasePhysical[2];

  G4VSolid*           m_test_block;
  G4LogicalVolume*    m_test_blockLogical;
  G4VPhysicalVolume*  m_test_blockPhysical[2];
//TEST SETUP STOP ----------------------------------

  //rpd booleans
  bool rpd_comp[8];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
