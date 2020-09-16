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
// $Id: PrimaryGeneratorAction.hh 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "globals.hh"
#include "AnalysisManager.hh"

#include "TFile.h"
#include "TTree.h"

#include <vector>

class G4GeneralParticleSource;
class G4Event;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);
    virtual void GenerateLHCEvent (G4Event*);
    virtual void GenerateSPSEvent (G4Event*);
    virtual void GenerateFNALEvent(G4Event*);

    virtual void InitializeCRMC(G4String GenModel);
    virtual void OpenInputFile(G4String fileName);
    virtual void ReadEvent();

    inline void SetBeamType               ( G4String arg ){ fBeamType = arg; }
    inline void SetCRMCmodel              ( G4String arg ){ InitializeCRMC( arg ); }
    inline void SetInputFile              ( G4String arg ){ OpenInputFile( arg ); }
    inline void SetVerticalCrossingAngle  ( G4double arg ){ fVertXingAngle = arg; }
    inline void SetHorizontalCrossingAngle( G4double arg ){ fHorizXingAngle = arg; }
    inline void SetProjectionPlane        ( G4double arg ){ PROJECT = true; fProjPlane = arg; }
    inline void SetPseudorapidityCut      ( G4double arg ){ fpsrCut = arg; }
    inline void SetnPrimaries             ( G4int    arg ){ fnPrimaries = arg; }
    inline void SetRandomizeRP            ( G4bool   arg ){ RANDOMIZE_RP = arg; }
    inline void SetBeamPos                ( G4ThreeVector* arg ){ delete fpos; fpos = arg; }


  private:
  	G4GeneralParticleSource*   fParticleGun;
    PrimaryGeneratorMessenger* fGeneratorMessenger;

    G4String       fBeamType;
    G4String       fGenInputFile;
    G4double       fVertXingAngle;
    G4double       fHorizXingAngle;
    G4double       fProjPlane;
    G4double       fpsrCut;
    G4double       fCRMCimpactPar;
    G4int          fnPrimaries;
    G4int          fCRMCnPart;
    G4int          fCurrentEvent;
    G4bool         PROJECT;
    G4bool         INPUT_INITIALIZED;
    G4bool         RANDOMIZE_RP;
    G4ThreeVector* fpos;
    AnalysisManager* m_analysisManager;
    G4RunManager* runManager;
    std::vector< G4PrimaryParticle* > fPrimaryVec;
    std::vector< std::vector<double>* > fdblVec;
    std::vector< std::vector< int  >* > fintVec;
    std::vector<double> *fCRMCpx, *fCRMCpy, *fCRMCpz, *fCRMCenergy, *fCRMCm;
    std::vector< int  > *fCRMCpdgid, *fCRMCstatus, *fCRMCkeptStatus, fCRMCkeptIndex;
    TFile* eventGenFile;
    TTree* eventGenDataTree;
    TTree* eventGenParticleTree;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
