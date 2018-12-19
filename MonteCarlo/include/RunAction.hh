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
// $Id: RunAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
//#include "G4Parameter.hh"
#include "globals.hh"
#include "SharedData.hh"
#include "TH2D.h"
#include "G4ThreeVector.hh"
#include <vector>


class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

class RunAction : public G4UserRunAction
{
  public:
    RunAction(SharedData *sd);
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  void ClearVectors() {
        fTrackID_v.clear();
	fModNb_v.clear();
        fRadNb_v.clear();
        fEdep_v.clear();
        fPid_v.clear();
        fNCherenkovs_v.clear();
        fX_v.clear();
        fY_v.clear();
        fZ_v.clear();
        fPx_v.clear();
        fPy_v.clear();
        fPz_v.clear();
        fEventNo_v.clear();  
        fEnergy_v.clear();
        fCharge_v.clear();
        fVelocity_v.clear();
        fBeta_v.clear();}
  inline SharedData* GetSharedData() {return fSharedData;}
  inline  void SetRadNo(G4int rNo) {fRadNb_v.push_back(rNo);}
  inline  void SetRodNo(G4int rNo) {fRodNb_v.push_back(rNo);}  
  inline  void SetEdep(G4double edep) {fEdep_v.push_back(edep);}
  inline  void SetModNb(G4int mNo) {fModNb_v.push_back(mNo);}
  inline  void SetTrackID(G4int trackid) {fTrackID_v.push_back(trackid);}
  void SetPosition(G4ThreeVector pos) {fX_v.push_back(pos.getX()); fY_v.push_back(pos.getY()); fZ_v.push_back(pos.getZ());}
  void SetMomentum(G4ThreeVector p) {fPx_v.push_back(p.getX()); fPy_v.push_back(p.getY()); fPz_v.push_back(p.getZ());}
  inline  void SetEnergy(G4double e) {fEnergy_v.push_back(e);}
  inline  void SetPid(G4int pid) {fPid_v.push_back(pid);}
  inline  void SetNCherenkovs(G4int cherenkovs) {fNCherenkovs_v.push_back(cherenkovs);}
  inline  void SetEventNo(G4int eventNo) {fEventNo_v.push_back(eventNo);}
  inline  void SetCharge(G4double charge) {fCharge_v.push_back(charge);}
 												   inline  void SetVelocity(G4double velocity) {fVelocity_v.push_back(velocity);}

												   inline  void SetBeta(G4double beta) {fBeta_v.push_back(beta);}
  
  private:
    SharedData *fSharedData;
  std::vector<int>         fTrackID_v;
  std::vector<int>         fModNb_v;
  std::vector<int>         fRadNb_v;
  std::vector<int>         fRodNb_v;  
  std::vector<double>      fEdep_v;
  std::vector<int>         fPid_v;
  std::vector<double>      fX_v;
  std::vector<double>      fY_v;
  std::vector<double>      fZ_v;
  std::vector<double>      fPx_v;
  std::vector<double>      fPy_v;
  std::vector<double>      fPz_v;
  std::vector<int>         fEventNo_v;
  std::vector<int>         fNCherenkovs_v;  
  std::vector<double>      fEnergy_v;
  std::vector<double>      fCharge_v;
  std::vector<double>    fVelocity_v;
  std::vector<double>      fBeta_v;
    
};

#endif

