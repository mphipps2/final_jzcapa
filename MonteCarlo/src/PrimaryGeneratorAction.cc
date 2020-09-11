
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

#include "CRMCconfig.h"

#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fBeamType("gps"),
  fGenModelStr("dpmjet"),
  fVertXingAngle(0.),
  fHorizXingAngle(0.),
  fProjPlane(0.),
  fpsrCut(5.),
  fnPrimaries(0),
  fCurrentEvent(0),
  PROJECT(false),
  CRMC_INITIALIZED(false),
  fpos( new G4ThreeVector(0.,0.,0.) ),
  eventGenFile(0)
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

  if( eventGenFile && eventGenFile->IsOpen() ){
    eventGenFile->Close();
    delete eventGenFile;
    for(std::vector<int>* vec : fintVec) delete vec;
    for(std::vector<double>* vec : fdblVec) delete vec;

  }
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

    // TODO: Add charge cut here

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

  //################################ Determine avaialble models

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


  //################################ Determine requested model

  fGenModelStr.toLower();
  G4int genModelCode = -1;
  G4String libName = "";
  G4String outputName;

       if(fGenModelStr.contains( "epos_lhc"   ) && modelsAvail[0]  ){ genModelCode = 0;  libName = "libEpos.so";       outputName = "epos199";    }
  else if(fGenModelStr.contains( "epos_1.99"  ) && modelsAvail[1]  ){ genModelCode = 1;  libName = "libEpos.so";       outputName = "eposlhc";    }
  else if(fGenModelStr.contains( "qgsjet01"   ) && modelsAvail[2]  ){ genModelCode = 2;  libName = "libQgsjet01.so";   outputName = "qgejet01";   }
  else if(fGenModelStr.contains( "pythia"     ) && modelsAvail[4]  ){ genModelCode = 4;  libName = "libPythia.so";     outputName = "pythia";     }
  else if(fGenModelStr.contains( "gheisha"    ) && modelsAvail[3]  ){ genModelCode = 3;  libName = "libGheisha.so";    outputName = "gheisha";    }
  else if(fGenModelStr.contains( "hijing"     ) && modelsAvail[5]  ){ genModelCode = 5;  libName = "libHijing.so";     outputName = "hijing";     }
  else if(fGenModelStr.contains( "sibyll"     ) && modelsAvail[6]  ){ genModelCode = 6;  libName = "libSibyll.so";     outputName = "sibyll";     }
  else if(fGenModelStr.contains( "qgsjetii04" ) && modelsAvail[7]  ){ genModelCode = 7;  libName = "libQgsjetII04.so"; outputName = "qgsjetII04"; }
  else if(fGenModelStr.contains( "phojet"     ) && modelsAvail[8]  ){ genModelCode = 8;  libName = "libPhojet.so";     outputName = "phojet";     }
  else if(fGenModelStr.contains( "qgsjetii03" ) && modelsAvail[11] ){ genModelCode = 11; libName = "libQgsjetII03.so"; outputName = "qgsjetII03"; }
  else if(fGenModelStr.contains( "dpmjet"     ) && modelsAvail[12] ){ genModelCode = 12; libName = "libDpmjet.so";     outputName = "dpmjet";     }
  else{
    G4cerr << "Invalid event generator model. Available options are (not case sensitive):\n" << modelList << G4endl;


    G4cerr << "Exiting or crashing, depending on if I figured out a good way to exit" << G4endl;


    return;
  }

  //################################ Generate the events

  long seed = CLHEP::RandFlat::shootInt(100000000);
  char command[128];
  sprintf(command,"%s/bin/crmc -o root -p2500 -P-2500 -n%d -s %ld -i208 -I208 -m%d",
          std::getenv("CRMC_INSTALL"),
          runManager->GetNumberOfEventsToBeProcessed(),
          seed,
          genModelCode );

  system(command);

  //################################ Initialize event gen output tree

  m_analysisManager = AnalysisManager::getInstance();

  //Assign vectors with reference names and place them in a vector for export to AnalysisManager
  int maxNpart = 200000;
  fintVec.push_back( fCRMCpdgid = new std::vector<int>(maxNpart,0) );
  fintVec.push_back( fCRMCstatus = new std::vector<int>(maxNpart,0) );
  fintVec.push_back( fCRMCkeptStatus = new std::vector<int>(maxNpart,0) );

  fdblVec.push_back( fCRMCpx = new std::vector<double>(maxNpart,0) );
  fdblVec.push_back( fCRMCpy = new std::vector<double>(maxNpart,0) );
  fdblVec.push_back( fCRMCpz = new std::vector<double>(maxNpart,0) );
  fdblVec.push_back( fCRMCenergy = new std::vector<double>(maxNpart,0) );
  fdblVec.push_back( fCRMCm = new std::vector<double>(maxNpart,0) );

  // Send the vectors to the analysis manager to set them as output
  m_analysisManager->MakeEventGenTree( fintVec, fdblVec );

  //################################ Open event gen file as input

  char fileName[128];
  sprintf(fileName, "crmc_%s_%ld_Pb_Pb_2500.root", outputName.data(), seed );

  eventGenFile = new TFile( fileName, "read" );
  if(eventGenFile->IsZombie()){
    G4cout << "file didn't read" << G4endl;
  }

  eventGenParticleTree = (TTree*)eventGenFile->Get("Particle");

  eventGenParticleTree->SetBranchAddress("nPart",&fCRMCnPart);
  eventGenParticleTree->SetBranchAddress("pdgid", &fCRMCpdgid->at(0) );
  eventGenParticleTree->SetBranchAddress("status",&fCRMCstatus->at(0) );
  eventGenParticleTree->SetBranchAddress("ImpactParameter",&fCRMCimpactPar);
  eventGenParticleTree->SetBranchAddress("px",&fCRMCpx->at(0) );
  eventGenParticleTree->SetBranchAddress("py",&fCRMCpy->at(0) );
  eventGenParticleTree->SetBranchAddress("pz",&fCRMCpz->at(0) );
  eventGenParticleTree->SetBranchAddress("E", &fCRMCenergy->at(0) );
  eventGenParticleTree->SetBranchAddress("m", &fCRMCm->at(0) );

  CRMC_INITIALIZED = true;
}//end InitializeCRMC

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateCRMCEvent()
{
  //double RPinclination = G4RandGauss::shoot(0.0,sigmaThetaXZ));
  G4ThreeVector momentum(0.,0.,0.);
  eventGenParticleTree->GetEntry( fCurrentEvent++ );

  //Clear vectors
  fPrimaryVec.clear();
  for(auto vec : fdblVec) vec->resize(fCRMCnPart);
  for(auto vec : fintVec) vec->resize(fCRMCnPart);


  //Add all final state particles to the particle vector
  for(int part = 0; part < fCRMCnPart; part++){

    // Cut out fragments in final state particles
    if(fCRMCstatus->at(part) == 1 && !(fCRMCpdgid->at(part) > 1000) ) continue;

    momentum.set( fCRMCpx->at(part), //px
                  fCRMCpy->at(part), //py
                  fCRMCpz->at(part));//pz

    //Rotate momentum for crossing angle and reaction plane
    momentum.rotateY(fHorizXingAngle);
    momentum.rotateX(fVertXingAngle);
    //momentum.rotateZ(RPinclination);

    // Push into the output vectors
    fCRMCpx->at(part) = momentum.x();
    fCRMCpy->at(part) = momentum.y();
    fCRMCpz->at(part) = momentum.z();

    if(momentum.pseudoRapidity() < fpsrCut){
      fPrimaryVec.push_back( new G4PrimaryParticle(fCRMCpdgid->at(part) ) ),
      fPrimaryVec.back()->SetKineticEnergy( fCRMCenergy->at(part) );

      fCRMCkeptStatus->push_back(1); //kept status == 1
    }else{
      fCRMCkeptStatus->push_back(0); //Not kept status == 0
    }// end if pseudoRapidity
  }// end particle loop
}// end GenerateCRMCEvent
