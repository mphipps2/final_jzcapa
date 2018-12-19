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

#ifndef ModType3_h
#define ModType3_h 1

#include "globals.hh"
#include "G4PVPlacement.hh"

#include <vector>

//class G4Box;
//class G4Para;
//class G4Tubs;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;

class SharedData;

/// Detector construction class to define materials and geometry.

class ModType3
{
public:
  ModType3(const int,
	 const G4ThreeVector&,
	 G4LogicalVolume*, SharedData*);
  ModType3();
  ~ModType3();
  
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

  G4VSolid*          m_W;
  G4LogicalVolume*   m_WLogical;
  G4VPhysicalVolume* m_WPhysical[11];

  // Vertical quartz strips (these strips are full -- no partial segments)
  G4VSolid*          m_StripTube;
  G4LogicalVolume*   m_StripLogical;
  G4VPhysicalVolume* m_StripPhysical[12][9][6];  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

