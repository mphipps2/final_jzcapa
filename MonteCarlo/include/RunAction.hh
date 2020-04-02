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
#include "globals.hh"
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
    RunAction( );
    virtual ~RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

  void ClearVectors() {
        fTrackID_v.clear();
		    fModNb_v.clear();
        fRadNb_v.clear();
		fRodNb_v.clear();
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
        fBeta_v.clear();
        fGap_Cherenk_v.clear();

		    fTrackID_v2.clear();
		    fModNb_v2.clear();
        fRadNb_v2.clear();
		    fRodNb_v2.clear();
        fEdep_v2.clear();
        fPid_v2.clear();
        fNCherenkovs_v2.clear();
        fX_v2.clear();
        fY_v2.clear();
        fZ_v2.clear();
        fPx_v2.clear();
        fPy_v2.clear();
        fPz_v2.clear();
        fEventNo_v2.clear();
        fEnergy_v2.clear();
        fCharge_v2.clear();
        fVelocity_v2.clear();
        fBeta_v2.clear();

		fTrackID_v3.clear();
		fModNb_v3.clear();
        fRadNb_v3.clear();
		    fRodNb_v3.clear();
        fEdep_v3.clear();
        fPid_v3.clear();
        fNCherenkovs_v3.clear();
        fX_v3.clear();
        fY_v3.clear();
        fZ_v3.clear();
        fPx_v3.clear();
        fPy_v3.clear();
        fPz_v3.clear();
        fEventNo_v3.clear();
        fEnergy_v3.clear();
        fCharge_v3.clear();
        fVelocity_v3.clear();
        fBeta_v3.clear();
		}

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
  inline  void SetGapCherenkovs(std::vector<int> cherenkovs) {
    fGap_Cherenk_v.resize(24,0);
    for(int i=0;i<24;i++) fGap_Cherenk_v.at(i)=cherenkovs.at(i);}

  inline  void SetRadNo_rpd(G4int rNo) {fRadNb_v2.push_back(rNo);}
  inline  void SetRodNo_rpd(G4int rNo) {fRodNb_v2.push_back(rNo);}
  inline  void SetEdep_rpd(G4double edep) {fEdep_v2.push_back(edep);}
  inline  void SetModNb_rpd(G4int mNo) {fModNb_v2.push_back(mNo);}
  inline  void SetTrackID_rpd(G4int trackid) {fTrackID_v2.push_back(trackid);}
  void SetPosition_rpd(G4ThreeVector pos) {fX_v2.push_back(pos.getX()); fY_v2.push_back(pos.getY()); fZ_v2.push_back(pos.getZ());}
  void SetMomentum_rpd(G4ThreeVector p) {fPx_v2.push_back(p.getX()); fPy_v2.push_back(p.getY()); fPz_v2.push_back(p.getZ());}
  inline  void SetEnergy_rpd(G4double e) {fEnergy_v2.push_back(e);}
  inline  void SetPid_rpd(G4int pid) {fPid_v2.push_back(pid);}
  inline  void SetNCherenkovs_rpd(G4int cherenkovs) {fNCherenkovs_v2.push_back(cherenkovs);}
  inline  void SetEventNo_rpd(G4int eventNo) {fEventNo_v2.push_back(eventNo);}
  inline  void SetCharge_rpd(G4double charge) {fCharge_v2.push_back(charge);}
  inline  void SetVelocity_rpd(G4double velocity) {fVelocity_v2.push_back(velocity);}
  inline  void SetBeta_rpd(G4double beta) {fBeta_v2.push_back(beta);}

  inline  void SetRadNo_fiber(G4int rNo) {fRadNb_v3.push_back(rNo);}
  inline  void SetRodNo_fiber(G4int rNo) {fRodNb_v3.push_back(rNo);}
  inline  void SetEdep_fiber(G4double edep) {fEdep_v3.push_back(edep);}
  inline  void SetModNb_fiber(G4int mNo) {fModNb_v3.push_back(mNo);}
  inline  void SetTrackID_fiber(G4int trackid) {fTrackID_v3.push_back(trackid);}
  void SetPosition_fiber(G4ThreeVector pos) {fX_v3.push_back(pos.getX()); fY_v3.push_back(pos.getY()); fZ_v3.push_back(pos.getZ());}
  void SetMomentum_fiber(G4ThreeVector p) {fPx_v3.push_back(p.getX()); fPy_v3.push_back(p.getY()); fPz_v3.push_back(p.getZ());}
  inline  void SetEnergy_fiber(G4double e) {fEnergy_v3.push_back(e);}
  inline  void SetPid_fiber(G4int pid) {fPid_v3.push_back(pid);}
  inline  void SetNCherenkovs_fiber(G4int cherenkovs) {fNCherenkovs_v3.push_back(cherenkovs);}
  inline  void SetEventNo_fiber(G4int eventNo) {fEventNo_v3.push_back(eventNo);}
  inline  void SetCharge_fiber(G4double charge) {fCharge_v3.push_back(charge);}
  inline  void SetVelocity_fiber(G4double velocity) {fVelocity_v3.push_back(velocity);}
  inline  void SetBeta_fiber(G4double beta) {fBeta_v3.push_back(beta);}

  private:
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
  std::vector<double>      fVelocity_v;
  std::vector<double>      fBeta_v;
  std::vector<int>         fGap_Cherenk_v;


  std::vector<int>         fTrackID_v2;
  std::vector<int>         fModNb_v2;
  std::vector<int>         fRadNb_v2;
  std::vector<int>         fRodNb_v2;
  std::vector<double>      fEdep_v2;
  std::vector<int>         fPid_v2;
  std::vector<double>      fX_v2;
  std::vector<double>      fY_v2;
  std::vector<double>      fZ_v2;
  std::vector<double>      fPx_v2;
  std::vector<double>      fPy_v2;
  std::vector<double>      fPz_v2;
  std::vector<int>         fEventNo_v2;
  std::vector<int>         fNCherenkovs_v2;
  std::vector<double>      fEnergy_v2;
  std::vector<double>      fCharge_v2;
  std::vector<double>      fVelocity_v2;
  std::vector<double>      fBeta_v2;

  std::vector<int>         fTrackID_v3;
  std::vector<int>         fModNb_v3;
  std::vector<int>         fRadNb_v3;
  std::vector<int>         fRodNb_v3;
  std::vector<double>      fEdep_v3;
  std::vector<int>         fPid_v3;
  std::vector<double>      fX_v3;
  std::vector<double>      fY_v3;
  std::vector<double>      fZ_v3;
  std::vector<double>      fPx_v3;
  std::vector<double>      fPy_v3;
  std::vector<double>      fPz_v3;
  std::vector<int>         fEventNo_v3;
  std::vector<int>         fNCherenkovs_v3;
  std::vector<double>      fEnergy_v3;
  std::vector<double>      fCharge_v3;
  std::vector<double>      fVelocity_v3;
  std::vector<double>      fBeta_v3;



};

#endif
