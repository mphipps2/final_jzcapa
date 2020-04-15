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

#include "Analysis.hh"

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

    virtual void MakeZDCTree       ( G4int nTupleNo, G4int zdcNo );
    virtual void MakeZDCOpticalTree( G4int nTupleNo, G4int zdcNo );
    virtual void MakeRPDTree       ( G4int nTupleNo, G4int rpdNo );
    virtual void MakeRPDOpticalTree( G4int nTupleNo, G4int rpdNo );



  private:
    G4String m_fileName;
    G4bool   CLUSTER;
    G4AnalysisManager* m_analysisManager;
    std::vector< std::vector< std::vector<double> > > *m_ZDCdblVec, *m_RPDdblVec;
    std::vector< std::vector< std::vector< int  > > > *m_ZDCintVec, *m_RPDintVec;

};

#endif
