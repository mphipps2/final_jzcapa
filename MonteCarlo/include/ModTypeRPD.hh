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

#include <vector>


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4Material;

class SharedData;

/// Detector construction class to define materials and geometry.

class ModTypeRPD
{
public:
  ModTypeRPD(const int,
	 const G4ThreeVector&,
	 G4LogicalVolume*, SharedData*);
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
  SharedData*          m_sd;

protected:
  
  G4Material*        m_matQuartz;

  bool               m_simCherenkov;
  
  G4VSolid*          m_tile;
  G4LogicalVolume*   m_tileLogical;
  G4VPhysicalVolume* m_tilePhysical[4][4];
  
  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

