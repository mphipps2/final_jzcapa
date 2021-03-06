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
// $Id: EventAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "AnalysisManager.hh"

#include "G4UserEventAction.hh"
#include "FiberHit.hh"
#include "globals.hh"

#include <vector>

/// Event action class
///

class EventAction : public G4UserEventAction
{
public:
  EventAction( );
  virtual ~EventAction();

  virtual void BeginOfEventAction          ( const G4Event* event );
  virtual void EndOfEventAction            ( const G4Event* event );
  virtual void ProcessHitCollection        ( FiberHitsCollection* HC );
  virtual void ProcessOpticalHitCollection ( FiberHitsCollection* HC );


private:
  // Indecies are [mod#][dataSet][dataValue]
  // Vectors are created in AnalysisManager.cc
  std::vector< std::vector< std::vector<double> > > *m_RPDdblVec, *m_ZDCdblVec;
  std::vector< std::vector< std::vector< int  > > > *m_RPDintVec, *m_ZDCintVec;

  AnalysisManager* m_analysisManager;
  G4int    fEventNo, fnZDCs, fnRPDs;
  unsigned int fEventSeed1, fEventSeed2;
  G4double gunPosX, gunPosY, gunPosZ;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
