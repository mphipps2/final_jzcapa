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
/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//
//
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class ModTypeRPD;
class ModTypeZDC;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;


class DetectorMessenger: public G4UImessenger{
  public:

    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    DetectorConstruction*      fDetector;

    G4UIdirectory*             fRPDDir;
    G4UIdirectory*             fZDCDir;

    G4UIcmdWithABool*          fClusterCmd;

    //ZDC Commands
    G4UIcmdWithoutParameter*   fZDCAddCmd;
    G4UIcmdWithAnInteger*      fZDCSetCurrentCmd;
    G4UIcmdWithABool*          fZDCForcePositionCmd;
    G4UIcmdWith3VectorAndUnit* fZDCPositionCmd;
    G4UIcmdWithADouble*        fZDCHousingThicknessCmd;
    G4UIcmdWithAnInteger*      fZDCnAbsorbersCmd;
    G4UIcmdWith3VectorAndUnit* fZDCAbsorberDimensionsCmd;
    G4UIcmdWithADouble*        fZDCFiberCoreDiameterCmd;
    G4UIcmdWithABool*          fZDCUseOpticalMaterialsCmd;


    //RPD commands
    G4UIcmdWithoutParameter*   fRPDAddCmd;
    G4UIcmdWithAnInteger*      fRPDSetCurrentCmd;
    G4UIcmdWith3VectorAndUnit* fRPDPositionCmd;
    G4UIcmdWithADouble*        fRPDSetFiberPitchCmd;
    G4UIcmdWithADouble*        fRPDSetTileSizeCmd;
    G4UIcmdWith3VectorAndUnit* fRPDHousingThicknessCmd;
    G4UIcmdWith3VectorAndUnit* fRPDFiberDimensionsCmd; //Core, Cladding, Buffer
    G4UIcmdWithABool*          fRPDUseOpticalMaterialsCmd;


};
#endif