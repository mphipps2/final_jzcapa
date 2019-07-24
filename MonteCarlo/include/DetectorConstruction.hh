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
#include "ModType1.hh"
#include "ModType2.hh"
#include "ModType3.hh"
#include "ModTypeCustom.hh"
#include "ModTypeZDC.hh"
#include "ModTypeRPD.hh"
#include "G4Cache.hh"
#include "G4MagneticField.hh"
#include <vector>
#include "Materials.hh"

class G4Box;
class G4Para;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class SharedData;

class G4UniformMagField;
class PurgMagTabulatedField3D;


/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  DetectorConstruction( SharedData* );
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

  virtual void               DefineBorderProperties();

  virtual G4VPhysicalVolume* ConstructDetector();

protected:
  SharedData*        m_sd;
  Materials*          materials;

protected:
  G4Box*               m_solidWorld;
  G4LogicalVolume*     m_logicWorld;
  G4VPhysicalVolume*   m_physWorld;
  G4LogicalVolume*     logic_leadTarget;
  G4LogicalVolume*     logic_leadBlock;

private:
	G4Cache<G4MagneticField*> fField;  //pointer to the thread-local fields

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
