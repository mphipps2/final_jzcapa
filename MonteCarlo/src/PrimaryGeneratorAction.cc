
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
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"

#include "CRMCinterface.h"

#include <iostream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fBeamType("gps"),
  fGenModelStr(""),
  fVertXingAngle(0.),
  fHorizXingAngle(0.),
  fProjPlane(0.),
  fpsrCut(5.),
  fnPrimaries(0),
  fGenModelCode(-1),
  PROJECT(false),
  CRMC_INITIALIZED(false),
  fpos( new G4ThreeVector(0.,0.,0.) )
{
  runManager = G4RunManager::GetRunManager();
  fGeneratorMessenger = new PrimaryGeneratorMessenger(this);
  fParticleGun = new G4GeneralParticleSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGeneratorMessenger;
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // Initialize CRMC if necessary
  if(fGenModelStr != "" && !CRMC_INITIALIZED) InitializeCRMC();


  if(fBeamType == "gps")
    fParticleGun->GeneratePrimaryVertex(anEvent);
  else if(fBeamType == "lhc")
    GenerateLHCEvent(anEvent);
  else if(fBeamType == "sps")
    GenerateSPSEvent(anEvent);
  else if(fBeamType == "fnal")
    GenerateFNALEvent(anEvent);
  else{
    G4cerr << "\nInvalid beam type selection. Aborting event\n" << G4endl;
    anEvent->SetEventAborted();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateLHCEvent(G4Event* anEvent)
{

  if( CRMC_INITIALIZED ){

    GenerateCRMCEvent();

  }else{ // Generate some neutrons
    // These values are hard coded for accuracy and consistency
    // Though they need to be updated with the correct values
    G4double sigmaThetaXZ = 0.;//3.57e-6; // 3.57e-6 Corresponds to a 1mm beam diameter
    G4double sigmaThetaYZ = 0.;//3.57e-6;
    G4double sigmaE = 1.e-3;
    G4double energy = (2.5 + G4RandGauss::shoot(0.0,sigmaE) )*TeV;

    G4ParticleDefinition* particleDefinition=
        G4ParticleTable::GetParticleTable()->FindParticle("neutron");

    // Adjust the momentum for the crossing angle
    G4ThreeVector momentum(0.,0.,1.);
    momentum.rotateY(fHorizXingAngle + G4RandGauss::shoot(0.0,sigmaThetaXZ));
    momentum.rotateX(fVertXingAngle  + G4RandGauss::shoot(0.0,sigmaThetaYZ));

    //If nPrimaries is 0, generate a random number from the distribution (to be implemented)
    if( fnPrimaries == 0 ) fnPrimaries = 1;
    // int nNeutrons = some distribution dependent random number;
    for(int i = 0; i < fnPrimaries; i++){
      fPrimaryVec.push_back( new G4PrimaryParticle(particleDefinition) );
      fPrimaryVec.back()->SetMomentumDirection( momentum );
      fPrimaryVec.back()->SetKineticEnergy( energy );
    }
  }

  for(uint i = 0; i < fPrimaryVec.size(); i++){

    // Project the beam to the plane requested by the user if requested
    // Otherwise carry forward the beam position
    G4ThreeVector* position;
    if(PROJECT){
      G4ThreeVector momentum = fPrimaryVec[i]->GetMomentumDirection();
      G4double projDist = fabs( fProjPlane - fpos->z() );
      G4double projectedX = fpos->x() + momentum.x()*projDist;
      G4double projectedY = fpos->y() + momentum.y()*projDist;
      position = new G4ThreeVector( projectedX, projectedY, fProjPlane );
    }else{
      position = new G4ThreeVector( fpos->x(), fpos->y(), fpos->z() );
    }

    G4PrimaryVertex* vertex = new G4PrimaryVertex( *position, 0. );
    vertex->SetPrimary( fPrimaryVec[i] );
    anEvent->AddPrimaryVertex( vertex );
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateSPSEvent(G4Event* anEvent)
{
  (void)anEvent; //Silence the unused variable message until this function is implemented
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateFNALEvent(G4Event* anEvent)
{
  (void)anEvent; //Silence the unused variable message until this function is implemented
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::InitializeCRMC()
{

  //Determine avaialble models
  std::vector< bool > modelsAvail(12,false);

  // EPOS_LHC and EPOS_1.99 are always avaialble apparently
  modelsAvail[0] = true;
  modelsAvail[1] = true;
  G4String modelList = "EPOS_LHC, EPOS_1.99";

  #ifdef __QGSJET01__
  modelsAvail[2] = true;
  modelList += ", Qgsjet01";
  #endif
  #ifdef __GHEISHA__
  modelsAvail[3] = true;
  modelList += ", Gheisha";
  #endif
  #ifdef __PYTHIA__
  modelsAvail[4] = true;
  modelList += ", Pythia";
  #endif
  #ifdef __HIJING__
  modelsAvail[5] = true;
  modelList += ", Hijing";
  #endif
  #ifdef __SIBYLL__
  modelsAvail[6] = true;
  modelList += ", Sibyll";
  #endif
  #ifdef __QGSJETII04__
  modelsAvail[7] = true;
  modelList += ", QgsjetII04";
  #endif
  #ifdef __PHOJET__
  modelsAvail[8] = true;
  modelList += ", Phojet";
  #endif
  #ifdef __QGSJETII03__
  modelsAvail[11] = true;
  modelList += ", QgsjetII03";
  #endif
  #ifdef __DPMJET__
  modelsAvail[12] = true;
  modelList += ", Dpmjet";
  #endif


  fGenModelStr.toLower();

       if(fGenModelStr.contains( "epos_lhc"   ) && modelsAvail[0]  ) fGenModelCode = 0;
  else if(fGenModelStr.contains( "epos_1.99"  ) && modelsAvail[1]  ) fGenModelCode = 1;
  else if(fGenModelStr.contains( "qgsjet01"   ) && modelsAvail[2]  ) fGenModelCode = 2;
  else if(fGenModelStr.contains( "gheisha"    ) && modelsAvail[3]  ) fGenModelCode = 3;
  else if(fGenModelStr.contains( "pythia"     ) && modelsAvail[4]  ) fGenModelCode = 4;
  else if(fGenModelStr.contains( "hijing"     ) && modelsAvail[5]  ) fGenModelCode = 5;
  else if(fGenModelStr.contains( "sibyll"     ) && modelsAvail[6]  ) fGenModelCode = 6;
  else if(fGenModelStr.contains( "qgsjetii04" ) && modelsAvail[7]  ) fGenModelCode = 7;
  else if(fGenModelStr.contains( "phojet"     ) && modelsAvail[8]  ) fGenModelCode = 8;
  else if(fGenModelStr.contains( "qgsjetii03" ) && modelsAvail[11] ) fGenModelCode = 11;
  else if(fGenModelStr.contains( "dpmjet"     ) && modelsAvail[12] ) fGenModelCode = 12;
  else{
    G4cerr << "Invalid event generator model. Available options are (not case sensitive):\n" << modelList << G4endl;


    G4cerr << "Exiting or crashing, depending on if I figured out a good way to exit" << G4endl;


    return;
  }

  // It should work fine at this point, but double check the init to be safe
  if( fCRMCInterface.init( fGenModelCode ) != 1) return;

  // open FORTRAN IO at first call
  fCRMCInterface.crmc_set(runManager->GetNumberOfEventsToBeProcessed(),
                      G4Random::getTheSeed(),
                      2500,           // Projectile momentum
                      -2500,          // Target momentum
                      208,            // Projectile ID
                      208,            // Target ID
                      fGenModelCode,  // Model
                      false,          // Produce tables (for output?)
                      1,              // Output type. Not applicable for us
                      "crmc.param "); // Parameter file name

  //call here variable settings from c++ interface
  //init models with set variables
  // fCRMCInterface.crmc_init(fCfg.GetOutputFileName().c_str(),fCfg.GetOutputFileName().size());
  fCRMCInterface.crmc_init("I don't think this matters",5);

  m_analysisManager = G4AnalysisManager::Instance();

  fNtupleNum = m_analysisManager->GetNofNtuples();


  //Integer
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "nGenParticles"  );
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "nSpectators"  );
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "model"  );

  //Doubles
  m_analysisManager->CreateNtupleDColumn( fNtupleNum, "impactParameter"  );

  ////std::vector< double >
  fdblVec->resize(5);
  m_analysisManager->CreateNtupleDColumn( fNtupleNum, "px",   fdblVec->at(0) );
  m_analysisManager->CreateNtupleDColumn( fNtupleNum, "py",   fdblVec->at(1) );
  m_analysisManager->CreateNtupleDColumn( fNtupleNum, "pz",   fdblVec->at(2) );
  m_analysisManager->CreateNtupleDColumn( fNtupleNum,  "E",   fdblVec->at(3) );
  m_analysisManager->CreateNtupleDColumn( fNtupleNum,  "m",   fdblVec->at(4) );


  //vector< int > branches
  fintVec->resize(3);
  m_analysisManager->CreateNtupleIColumn( fNtupleNum,      "pdgid", fintVec->at(0) );
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "CRMCstatus", fintVec->at(1) );
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "keptStatus", fintVec->at(2) );

  m_analysisManager->FinishNtuple( );

  CRMC_INITIALIZED = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateCRMCEvent()
{
  //double RPinclination = G4RandGauss::shoot(0.0,sigmaThetaXZ));
  G4ThreeVector momentum(0.,0.,0.);
  // cleanup vectors
  gCRMC_data.Clean();

  fCRMCInterface.crmc_generate(1,runManager->GetCurrentEvent()->GetEventID()+1,
                           gCRMC_data.fNParticles,
                           gCRMC_data.fImpactParameter,
                           gCRMC_data.fPartId[0],
                           gCRMC_data.fPartPx[0],
                           gCRMC_data.fPartPy[0],
                           gCRMC_data.fPartPz[0],
                           gCRMC_data.fPartEnergy[0],
                           gCRMC_data.fPartMass[0],
                           gCRMC_data.fPartStatus[0]);

  //Add all final state particles to
  for(int part = 0; part < gCRMC_data.fNParticles; part++){
    // Check for fragments in final state particles
    if(gCRMC_data.fPartStatus[part] == 1 && !(gCRMC_data.fPartId[part] > 1000) ){
      momentum.set( gCRMC_data.fPartPx[part],
                    gCRMC_data.fPartPy[part],
                    gCRMC_data.fPartPz[part]);

      //Rotate momentum for crossing angle and reaction plane
      momentum.rotateY(fHorizXingAngle);
      momentum.rotateX(fVertXingAngle);
      //momentum.rotateZ(RPinclination);


      fdblVec->at(0).push_back(momentum.x()); //px
      fdblVec->at(1).push_back(momentum.y()); //py
      fdblVec->at(2).push_back(momentum.z()); //pz
      fdblVec->at(3).push_back(gCRMC_data.fPartEnergy[part]); // E
      fdblVec->at(4).push_back(gCRMC_data.fPartMass[part]); // m

      fintVec->at(1).push_back(gCRMC_data.fPartStatus[part]); // CRMCstatus

      if(momentum.pseudoRapidity() < fpsrCut){
        fPrimaryVec.push_back( new G4PrimaryParticle(gCRMC_data.fPartId[part] ) ),
        fPrimaryVec.back()->SetKineticEnergy( gCRMC_data.fPartEnergy[part] );

        fintVec->at(2).push_back(1); //kept status == 1
      }else{
        fintVec->at(2).push_back(0); //Not kept status == 0
      }// end if pseudoRapidity
    }// end if fragment
  }// end particle loop

  gCRMC_data.sigtot = double(hadr5_.sigtot);      // ........ h-p total cross section in mb
  gCRMC_data.sigine = double(hadr5_.sigine);      // ........ h-p inelastic cross section in mb (all inelastic processes=sigtot-sigela)
  gCRMC_data.sigela = double(hadr5_.sigela);      // ........ h-p elastic cross section in mb
  gCRMC_data.sigdd = double(hadr5_.sigdd);        // ........ h-p double diffractive cross section in mb (both side)
  gCRMC_data.sigsd = double(hadr5_.sigsd);        // ........ h-p single diffractive cross section in mb (both side)
  gCRMC_data.sloela = double(hadr5_.sloela);      // ........ h-p elastic slope
  gCRMC_data.sigtotaa = double(hadr5_.sigtotaa);  // ........ h-A or A-A total cross section in mb
  gCRMC_data.sigineaa = double(hadr5_.sigineaa);  // ........ h-A or A-A cross section in mb (inelastic cross section to be used as CURRENT EVENT XS for defined projectile and target, previous h-p xs are really for h-p even if the projectile/target were defined as a nuclei)
  gCRMC_data.sigelaaa = double(hadr5_.sigelaaa);  // ........ h-A or A-A elastic cross section in mb
  gCRMC_data.npjevt = cevt_.npjevt;               // ........ number of primary projectile participants
  gCRMC_data.ntgevt = cevt_.ntgevt;               // ........ number of primary target participants
  gCRMC_data.kolevt = cevt_.kolevt;               // ........ number of collisions
  gCRMC_data.kohevt = cevt_.kohevt;               // ........ number of inelastic hard collisions
  gCRMC_data.npnevt = cevt_.npnevt;               // ........ number of primary projectile neutron spectators
  gCRMC_data.ntnevt = cevt_.ntnevt;               // ........ number of primary target neutron spectators
  gCRMC_data.nppevt = cevt_.nppevt;               // ........ number of primary projectile proton spectators
  gCRMC_data.ntpevt = cevt_.ntpevt;               // ........ number of primary target proton spectators
  gCRMC_data.nglevt = cevt_.nglevt;               // ........ number of collisions acc to  Glauber
  gCRMC_data.ng1evt = c2evt_.ng1evt;              // ........ number of collisions acc to  Glauber
  gCRMC_data.ng2evt = c2evt_.ng2evt;              // ........ number of Glauber participants with at least two IAs
  gCRMC_data.bimevt = double(cevt_.bimevt);       // ........ absolute value of impact parameter
  gCRMC_data.phievt = double(cevt_.phievt);       // ........ angle of impact parameter
  gCRMC_data.fglevt = double(c2evt_.fglevt);      // ........
  gCRMC_data.typevt = int(c2evt_.typevt);         // ........ type of event (1=Non Diff, 2=Double Diff, 3=Central Diff, 4=AB->XB, -4=AB->AX)



  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "nGenParticles"  );
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "nSpectators"  );
  m_analysisManager->CreateNtupleIColumn( fNtupleNum, "model"  );

  //Doubles
  m_analysisManager->CreateNtupleDColumn( fNtupleNum, "impactParameter"  );



/* Unused vaiables

  hadr5_
      float sigcut;    // ........ h-p cut cross section in mb : in principle it is the non-diffractive xs but the definition depends on the model
      float sigdif;    // ........ h-p diffractive cross section in mb (SD+DD+DPE) (in principle sigdif+sigcut=sigine but it depends how DPE xs is counted (in EPOS 1.99 it is counted as elastic because nothing was produced but in EPOS LHC DPE are produced)
      float sigcutaa;  // ........ h-A or A-A ND xs or production xs mb

  cevt_
      int   nevt;   // ........ error code. 1=valid event, 0=invalid event
      int   koievt; // ........ number of inelastic collisions
      float pmxevt; // ........ reference momentum
      float egyevt; // ........ pp cm energy (hadron) or string energy (lepton)
      int   jpnevt; // ........ number of absolute projectile neutron spectators
      int   jppevt; // ........ number of absolute projectile proton spectators
      int   jtnevt; // ........ number of absolute target neutron spectators
      int   jtpevt; // ........ number of absolute target proton spectators
      float xbjevt; // ........ bjorken x for dis
      float qsqevt; // ........ q**2 for dis
      float zppevt; // ........ average Z-parton-proj
      float zptevt; // ........ average Z-parton-targ
      int   minfra; // ........
      int   maxfra; // ........

  c2evt_
      float rglevt; // ........
      float sglevt; // ........
      float eglevt; // ........
      int   ikoevt; // ........ number of elementary parton-parton scatterings
*/



}
