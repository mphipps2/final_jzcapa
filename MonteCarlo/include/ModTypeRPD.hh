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

#ifndef ModTypeRPD_h
#define ModTypeRPD_h 1


#include "globals.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
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
  ModTypeRPD(const int,
	 const G4ThreeVector&,
	 G4LogicalVolume*);
  ModTypeRPD();
  ~ModTypeRPD();



  virtual void  Construct();

  virtual void  DefineMaterials();
//  virtual void  DefineBorderProperties();

  virtual void  ConstructDetector();

protected:
  const int            m_modNum;
  const G4ThreeVector  m_pos;
  G4LogicalVolume*     m_logicMother;

  Materials*          materials;
protected:

  G4Material*           m_matQuartz;

  G4Material*           m_Al;

  G4Material*           m_Poly;

  G4Material*           m_Air;

  G4Material*           m_PMMA;

  G4Material*           m_Grease;

  G4OpticalSurface* m_photonDetSurface;

  bool               m_simCherenkov;

  G4SubtractionSolid* m_tile;

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
  G4LogicalVolume*   m_fiberLogical[64];//4
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


//test_setup
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

  //pan flute rpd
  G4VSolid*          m_PFrpd[4];
  G4LogicalVolume*   m_PFrpdLogical[64];
  G4VPhysicalVolume* m_PFrpdPhysical[64];

  G4VSolid*           m_test_foil[4];
  G4LogicalVolume*    m_test_foilLogical[64];
  G4VPhysicalVolume*  m_test_foilPhysical[64];

  G4VSolid*           m_PFdetec;
  G4LogicalVolume*    m_PFdetecLogical[64];
  G4VPhysicalVolume*  m_PFdetecPhysical[64];

  //rpd booleans
  bool rpd_comp[8];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
