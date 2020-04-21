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
/// \file AnalysisManager.hh
/// \brief Selection of the analysis technology
/// \author Chad Lantz
/// \date 16 April 2020

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "g4root.hh"
//#include "g4cvs.hh"
//#include "g4xml.hh"

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AnalysisManager
{
  public:
    static AnalysisManager* getInstance(void);
   ~AnalysisManager();

    void Book();
    void Save();

    void FillNtuples( );

    void MakeZDCTree       ( G4int nTupleNo, G4int zdcNo );
    void MakeZDCOpticalTree( G4int nTupleNo, G4int zdcNo );
    void MakeRPDTree       ( G4int nTupleNo, G4int rpdNo );
    void MakeRPDOpticalTree( G4int nTupleNo, G4int rpdNo );

    inline G4bool GetClusterFlag(){ return CLUSTER; }
    inline G4bool GetOpticalFlag(){ return OPTICAL; }

    inline  std::vector< std::vector< std::vector<double> > >*  GetRPDdoubleVectors( ){return m_RPDdblVec;}
    inline  std::vector< std::vector< std::vector< int  > > >*  GetRPDintVectors   ( ){return m_RPDintVec;}
    inline  std::vector< std::vector< std::vector<double> > >*  GetZDCdoubleVectors( ){return m_ZDCdblVec;}
    inline  std::vector< std::vector< std::vector< int  > > >*  GetZDCintVectors   ( ){return m_ZDCintVec;}

  private:
    AnalysisManager();
    static AnalysisManager* analysisManager;
    G4bool m_FactoryOn;
    G4bool CLUSTER, OPTICAL;
    G4AnalysisManager* m_analysisManager;
    DetectorConstruction* m_detectorConstruction;
    std::vector< std::vector< std::vector<double> > > *m_ZDCdblVec, *m_RPDdblVec;
    std::vector< std::vector< std::vector< int  > > > *m_ZDCintVec, *m_RPDintVec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
