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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef FiberSD_h
#define FiberSD_h 1

#include "G4VSensitiveDetector.hh"
#include "FiberHit.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class FiberSD : public G4VSensitiveDetector
{
public:
  FiberSD(G4String, G4int,G4bool);
  ~FiberSD();

  void HistInitialize();

  void   Initialize          ( G4HCofThisEvent* );
  G4bool ProcessHits         ( G4Step*, G4TouchableHistory* );
  void   EndOfEvent          ( G4HCofThisEvent* );

  inline G4bool   OpticalIsOn    ( ){ return OPTICAL; }
  inline void     SetTopOfVolume ( G4double _top ){ m_topOfVolume = _top; }
  inline G4double GetTopOfVolume ( ){ return m_topOfVolume; }
  inline int      GetNCherenkovs ( ){ return m_nCherenkovs; }

private:
  int HCID;
  G4double m_modCoreIndexRefraction;
  FiberHitsCollection* fiberCollection;
  G4int m_modNum;
  G4int m_nCherenkovs;
  G4bool OPTICAL;
  G4double m_topOfVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
