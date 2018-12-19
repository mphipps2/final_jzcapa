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

#ifndef ModType2_h
#define ModType2_h 1

#include "globals.hh"
#include "G4PVPlacement.hh"

#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;

class SharedData;

/// Detector construction class to define materials and geometry.

class ModType2
{
public:
  ModType2(const int,
	 const G4ThreeVector&,
	 G4LogicalVolume*, SharedData*);
  ModType2();
  ~ModType2();
  
  virtual void  Construct();
  
  virtual void  DefineMaterials();  
  virtual void  DefineBorderProperties();

  virtual void  ConstructDetector();
  
protected:
  const int            m_modNum;  
  const G4ThreeVector  m_pos;
  G4LogicalVolume*     m_logicMother;
  SharedData*          m_sd;

protected:
  G4Material*        m_matHousing;
  G4Material*        m_matQuartz;
  G4Material*        m_matAbsorber;
  G4Material*        m_matAir;
  
  G4VSolid*          m_ModuleBox;
  G4LogicalVolume*   m_ModuleLogical;
  G4VPhysicalVolume* m_ModulePhysical;

  G4VSolid*          m_AirScoringBox;
  G4LogicalVolume*   m_AirScoringLogical;
  G4VPhysicalVolume* m_AirScoringPhysical;
  
  G4VSolid*          m_SteelBox;
  G4LogicalVolume*   m_SteelLogical;
  G4VPhysicalVolume* m_SteelBoxPhysical;

  // Lateral segment of tungsten plate for the seven middle groups (6 quartz strips/group in these sections)
  G4VSolid*          m_WLatCenter;
  G4LogicalVolume*   m_WLatCenterLogical;
  G4VPhysicalVolume* m_WLatPhysical[11][9]; // filled with LatCenter and LatEdge

  // Lateral segment of tungsten plate for the two edge groups (only 5 quartz strips in the edge groups so lateral tungsten layer adjusted accordingly)
  G4VSolid*          m_WLatEdge;
  G4LogicalVolume*   m_WLatEdgeLogical;

  // Longitudinal segments of tungsten plates -- 8/module that run length of entire module in space between the quartz strip groups
  G4VSolid*          m_WLong;
  G4LogicalVolume*   m_WLongLogical;
  G4VPhysicalVolume* m_WLongPhysical[8];

  // Pixel segments carved out of tungsten layers
  G4VSolid*          m_PixelTubeW;
  G4LogicalVolume*   m_PixelWLogical;
  G4PVPlacement*     m_PixelWPhysical[11][16][8];  // filled with PixelWLogical and PixelAirWLogical

  // Air segments carved out of tungsten layers (holes unfilled by pixels)
  G4VSolid*          m_AirTubeW;
  G4LogicalVolume*   m_PixelAirWLogical;

  // Pixel segments carved out of radiator layers (more precisely from "longitudinal W plates"
  G4VSolid*          m_PixelTubeRad;
  G4LogicalVolume*   m_PixelRadLogical;
  G4PVPlacement*     m_PixelRadPhysical[12][16][8];  // filled with PixelRadLogical and PixelAirRadLogical

  // Air segments carved out of radiator layers
  G4VSolid*          m_AirTubeRad;
  G4LogicalVolume*   m_PixelAirRadLogical;

  // Vertical quartz strips (these strips are full -- no partial segments)
  G4VSolid*          m_StripTube;
  G4LogicalVolume*   m_StripLogical;
  G4PVPlacement*     m_StripPhysical[12][9][6];  

  // Rectangular air volumes carved out of longitudinal W plates (entire gap Y -- pixels then get carved out of this)
  G4VSolid*          m_AirGap; // add this to source
  G4LogicalVolume*   m_AirGapLogical;
  G4PVPlacement*     m_AirPhysical[12][8];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

