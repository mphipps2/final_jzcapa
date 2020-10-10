
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

#ifdef CRMC
#include "CRMCconfig.h"
#endif

#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fBeamType("gps"),
  fGenInputFile(""),
  fGenModel(""),
  fPtDist(""),
  fMultDist(""),
  fVertXingAngle(0.),
  fHorizXingAngle(0.),
  fProjPlane(0.),
  fpsrCut(9.),
  fCollisionPt(0.15),
  fFragmentationPt(0.15),
  fnPrimaries(0),
  fCurrentEvent(0),
  fMinNspec(0),
  fMaxNspec(120),
  PROJECT(false),
  INPUT_INITIALIZED(false),
  GENERATE_CRMC_EVENTS(false),
  RANDOMIZE_RP(false),
  fpos( new G4ThreeVector(0.,0.,0.) ),
  eventGenFile(0)
{
  runManager = G4RunManager::GetRunManager();
  fGeneratorMessenger = new PrimaryGeneratorMessenger(this);
  fParticleGun = new G4GeneralParticleSource();
  m_analysisManager = AnalysisManager::getInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGeneratorMessenger;
  delete fParticleGun;

  if( INPUT_INITIALIZED ){
    eventGenFile->Close();
    delete eventGenFile;
    for(std::vector<int>* vec : fintVec) delete vec;
    for(std::vector<double>* vec : fdblVec) delete vec;

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if(GENERATE_CRMC_EVENTS && !INPUT_INITIALIZED) InitializeCRMC();

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

void PrimaryGeneratorAction::SetGeneratorModel( G4String model ){
  if(!model.contains("toy")) GENERATE_CRMC_EVENTS = true;
  fGenModel = model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateLHCEvent(G4Event* anEvent)
{

       if( INPUT_INITIALIZED ) ReadEvent();
  else if( fGenModel.contains("toyv1") ) GenerateToyV1();
  else{ // Generate some neutrons
    // These values are hard coded for accuracy and consistency
    // Though they need to be updated with the correct values
    G4double sigmaThetaXZ = 0.;//3.57e-6; // 3.57e-6 Corresponds to a 1mm beam diameter
    G4double sigmaThetaYZ = 0.;//3.57e-6;
    G4double sigmaE = 1.e-3;
    G4double energy = (2.5 + G4RandGauss::shoot(0.0,sigmaE) )*TeV;
    fPrimaryVec.clear();

    G4ParticleDefinition* particleDefinition=
        G4ParticleTable::GetParticleTable()->FindParticle("neutron");

    // Adjust the momentum for the crossing angle
    G4double psiRP = (RANDOMIZE_RP) ? CLHEP::RandFlat::shootInt( CLHEP::twopi ) : 0.0;
    G4ThreeVector momentum(0.,0.,1.);
    momentum.rotateY(fHorizXingAngle + G4RandGauss::shoot(0.0,sigmaThetaXZ));
    momentum.rotateX(fVertXingAngle  + G4RandGauss::shoot(0.0,sigmaThetaYZ));
    momentum.rotateZ( psiRP );

    // If nPrimaries is 0, generate a random number from the distribution (to be implemented)
    if( fnPrimaries == 0 ) fnPrimaries = 1;
    // int nNeutrons = some distribution dependent random number;
    for(int i = 0; i < fnPrimaries; i++){
      fPrimaryVec.push_back( new G4PrimaryParticle(particleDefinition) );
      fPrimaryVec.back()->SetMomentumDirection( momentum );
      fPrimaryVec.back()->SetKineticEnergy( energy );
    }
  }

  for(uint i = 0; i < fPrimaryVec.size(); i++){

    // If the particle is charged it will be swept away by the steering magnets.
    // Also, if particle is outside of the acceptance of the detector (set by user)
    // Change the kept status to 0 (not kept) and continue to the next particle.
    G4ThreeVector momentum = fPrimaryVec[i]->GetMomentumDirection();
    if( fPrimaryVec[i]->GetCharge() != 0.0 || momentum.pseudoRapidity() < fpsrCut ){
      fCRMCkeptStatus->at( fCRMCkeptIndex[i] ) = 0;
      continue;
    }

    // Project the beam to the plane requested by the user if requested
    // Otherwise carry forward the beam position
    G4ThreeVector* position;
    if(PROJECT){
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

  #ifndef CRMC
  G4cerr << "CRMC wasn't found. Check your configuration before trying again" << G4endl;
  return;
  #endif

  //################################ Determine avaialble models as defined in CRMCconfig.h

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

  G4int genModelCode = -1;
  G4String libName = "";
  G4String model;

       if(fGenModel.contains( "epos_lhc"   ) && modelsAvail[0]  ){ genModelCode = 0;  libName = "libEpos.so";       model = "epos199";    }
  else if(fGenModel.contains( "epos_1.99"  ) && modelsAvail[1]  ){ genModelCode = 1;  libName = "libEpos.so";       model = "eposlhc";    }
  else if(fGenModel.contains( "qgsjet01"   ) && modelsAvail[2]  ){ genModelCode = 2;  libName = "libQgsjet01.so";   model = "qgejet01";   }
  else if(fGenModel.contains( "pythia"     ) && modelsAvail[4]  ){ genModelCode = 4;  libName = "libPythia.so";     model = "pythia";     }
  else if(fGenModel.contains( "gheisha"    ) && modelsAvail[3]  ){ genModelCode = 3;  libName = "libGheisha.so";    model = "gheisha";    }
  else if(fGenModel.contains( "hijing"     ) && modelsAvail[5]  ){ genModelCode = 5;  libName = "libHijing.so";     model = "hijing";     }
  else if(fGenModel.contains( "sibyll"     ) && modelsAvail[6]  ){ genModelCode = 6;  libName = "libSibyll.so";     model = "sibyll";     }
  else if(fGenModel.contains( "qgsjetii04" ) && modelsAvail[7]  ){ genModelCode = 7;  libName = "libQgsjetII04.so"; model = "qgsjetII04"; }
  else if(fGenModel.contains( "phojet"     ) && modelsAvail[8]  ){ genModelCode = 8;  libName = "libPhojet.so";     model = "phojet";     }
  else if(fGenModel.contains( "qgsjetii03" ) && modelsAvail[11] ){ genModelCode = 11; libName = "libQgsjetII03.so"; model = "qgsjetII03"; }
  else if(fGenModel.contains( "dpmjet"     ) && modelsAvail[12] ){ genModelCode = 12; libName = "libDpmjet.so";     model = "dpmjet";     }
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

  char fileName[128];
  sprintf(fileName, "crmc_%s_%ld_Pb_Pb_2500.root", model.data(), seed );
  OpenInputFile( fileName );

}// end InitializeCRMC

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::OpenInputFile(G4String fileName)
{

  // If input has already been initialized inform the user and bail
  if( INPUT_INITIALIZED ){
    G4cerr << "WARNING: Multiple input files/methods selected." << G4endl;
    G4cerr << "Make sure /beam/eventGen and /beam/input aren't both selected in your run macro" << G4endl;
    return;
  }

  //################################ Initialize event gen output tree

  // Assign vectors with reference names and place them in a vector for export to AnalysisManager
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
  m_analysisManager->MakeEventGenTree( fintVec, fdblVec, 0 );


  //################################ Open input file
  eventGenFile = new TFile( fileName.data(), "read" );
  if(eventGenFile->IsZombie()){
    G4cout << "file didn't read" << G4endl;
    return;
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

  INPUT_INITIALIZED = true;

}// end OpenInputFile

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::ReadEvent()
{
  if( eventGenParticleTree->GetEntries() < runManager->GetNumberOfEventsToBeProcessed() ){
      G4cerr << "Not enough entries in input file to complete the run" << G4endl;

      G4cerr << "Exiting or crashing or resetting the number of events, depending on if I figured out a good way to exit" << G4endl;
  }
  G4ThreeVector momentum(0.,0.,0.);
  G4double psiRP = (RANDOMIZE_RP) ? CLHEP::RandFlat::shootInt( CLHEP::twopi ) : 0.0;
  eventGenParticleTree->GetEntry( fCurrentEvent++ );

  // Clear vectors
  fPrimaryVec.clear();
  fCRMCkeptIndex.clear();
  fCRMCkeptStatus->clear();
  fCRMCkeptStatus->resize(fCRMCnPart,0);
  for(auto vec : fdblVec) vec->resize(fCRMCnPart);
  for(auto vec : fintVec) vec->resize(fCRMCnPart);


  // Add all final state particles to the particle vector
  for(int part = 0; part < fCRMCnPart; part++){

    momentum.set( fCRMCpx->at(part), //px
                  fCRMCpy->at(part), //py
                  fCRMCpz->at(part));//pz


    // Rotate momentum for crossing angle and reaction plane
    momentum.rotateY(fHorizXingAngle);
    momentum.rotateX(fVertXingAngle);
    momentum.rotateZ(psiRP);

    // Replace the content of the output vectors with the rotated components
    fCRMCpx->at(part) = momentum.x();
    fCRMCpy->at(part) = momentum.y();
    fCRMCpz->at(part) = momentum.z();


    // Cut out fragments in final state particles
    if(fCRMCstatus->at(part) == 1 && fCRMCpdgid->at(part) > 999999999 ) continue;

    fPrimaryVec.push_back(
      new G4PrimaryParticle( fCRMCpdgid->at(part),
                             fCRMCpx->at(part)*GeV,
                             fCRMCpy->at(part)*GeV,
                             fCRMCpz->at(part)*GeV,
                             fCRMCenergy->at(part)*GeV) );

    fCRMCkeptIndex.push_back(part);
    fCRMCkeptStatus->at(part) = 1; //kept status == 1
  }// end particle loop
}// end ReadEvent

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateToyV1(){

  //################################
  // Set up output variables
  //################################
    std::vector< double > *Px, *Py, *Pz, *E;
    if( fCurrentEvent == 0 ){

      fdblVec.push_back( Px = new std::vector<double> );
      fdblVec.push_back( Py = new std::vector<double> );
      fdblVec.push_back( Pz = new std::vector<double> );
      fdblVec.push_back( E  = new std::vector<double> );

      m_analysisManager->MakeEventGenTree( fintVec, fdblVec, 1);
    }else{
      Px = fdblVec[0];
      Py = fdblVec[1];
      Pz = fdblVec[2];
      E  = fdblVec[3];

      for(uint i = 0; i < fdblVec.size(); i++){
        fdblVec[i]->clear();
      }
      fPrimaryVec.clear();
    }
    fCurrentEvent++;

  //################################
  // Set up pt generation method
  //################################

  // Selection of pT generation scheme
  G4int pTChoice = 0;
         if( fPtDist.contains("gaus") ) pTChoice = 0; //< pT generation accoding to a gaussian of given sigma
    else if( fPtDist.contains("hist") ) pTChoice = 1; //< generation according to pT from histo
    else if( fPtDist.contains("func") ) pTChoice = 2; //< generation according to pT from function

    //Pt generation pre-stages
    TH1D* pThist = NULL;
    TFile* pTfile = NULL;
    if(pTChoice == 1){
      // Try to open the file in the Utils directory then try locally
      pTfile = TFile::Open( (std::getenv("JZCaPA") + (std::string)"Utils/eposPTneut.root").c_str() );
      if( pTfile->IsZombie() ){
        delete pTfile;
        pTfile = TFile::Open("eposPTneut.root");
      }

      // Grab the histogram
      pThist = (TH1D*)pTfile->Get("npt");
    }

    TF1 *pTfun = NULL;
    if(pTChoice == 2){
      pTfun = new TF1("pTfun","gaus(0)+pol1(3)",0,1.5);
      pTfun->SetParameter(0,850);
      pTfun->SetParameter(1,0.10); //pT avg in GeV
      pTfun->SetParameter(2,0.04);
      pTfun->SetParameter(3,300.);
      pTfun->SetParameter(4,-600.);
    }


    //################################
    // Set up multiplicity generation method
    //################################

    //Selection of # of neutrons generation method
    int spectatorChoice = 0;
         if( fMultDist.contains("func") ) spectatorChoice = 0; //< random from poisson distribution
    else if( fMultDist.contains("hist") ) spectatorChoice = 1; //< from ATLAS ZDC multiplicity histogram

    //Spectators Number pre-stages
    TFile* nAnBfile = NULL;
    TH2D* nAnBhist = NULL;

    if(spectatorChoice == 1){
      // Try to open the file in the Utils directory then try locally
      nAnBfile = TFile::Open( (std::getenv("JZCaPA") + (std::string)"Utils/ZDC_AvsC.root").c_str() );
      if( nAnBfile->IsZombie() ){
        delete nAnBfile;
        nAnBfile = TFile::Open("ZDC_AvsC.root");
      }
      nAnBhist = (TH2D*)nAnBfile->Get("hZdcAZdcC_MB");
    }

    //Per nucleon energy in the HI collision, in GeV
    double Ptot = 2760.;

    //PARAMETER FOR GENERATION
    double pTNuclearComponents[2];
    double reactionPlaneAngle; //Defining the reaction plane angle along which the pT of v1 will be applied

    int particles  = -1;  //# of neutrons going to ZDC + (or A)
    int particles2 = -1;  //# of neutrons going to ZDC - (or C)

    double bufferPtGen;
    //Random number generator engine
    TRandom3* myRand = new TRandom3(0);
    //Containers
    std::vector < double > pTplus;         //Vector containing pT of neutrons going towards ZDC + (or A)
    std::vector < double > azimuth_plus;   //Vector containing azimuthal angle of neutrons going towards ZDC + (or A)
    std::vector < double > pTminus;        //Vector containing pT of neutrons going towards ZDC - (or C)
    std::vector < double > azimuth_minus;  //Vector containing azimuthal angle of neutrons going towards ZDC - (or C)

    //Buffer variables for kinematic computation
    double neutron_mass = 0.939565; // GeV
    G4ParticleDefinition* particleDefinition=
        G4ParticleTable::GetParticleTable()->FindParticle("neutron");

    //################################
    // Implementation
    //################################

    //===================================
    //Spectators Block - Here is defined how many spectators we have for the event
    if(spectatorChoice == 0){
      particles = -1.;
      while( particles < 0){
        particles = myRand->Poisson(fnPrimaries);
        particles2 = myRand->Poisson(fnPrimaries);
        if(particles < fMinNspec || particles > fMaxNspec) particles = -1;
      }
    }else if(spectatorChoice == 1){
      double buffer, buffer2;
      particles = -1;
      particles2 = -1;
      //TODO: make this also for particles2 - at the moment 1 arm implementation
      //How TODO: Make more particles and change the sign of their momentum components
      while( particles < 1){
        nAnBhist->GetRandom2(buffer, buffer2);
        particles = (int)buffer;
        particles2 = (int)buffer2;
        if(particles < fMinNspec || particles > fMaxNspec) particles = -1;
      }
    }
    //===================================
    //pT nuclear block - Here we extract a direction for the pT nuclear and we compute components
    reactionPlaneAngle = -TMath::Pi()+(myRand->Rndm()*TMath::TwoPi());
    pTNuclearComponents[0] = fCollisionPt*TMath::Cos(reactionPlaneAngle);
    pTNuclearComponents[1] = fCollisionPt*TMath::Sin(reactionPlaneAngle);

    for(int i = 0; i < particles; i++){
      //Extraction of pT for particle i
      if(pTChoice == 0) pTplus.push_back(TMath::Abs(myRand->Gaus(0,fFragmentationPt)));
      else if(pTChoice == 1) pTplus.push_back(pThist->GetRandom());
      else if(pTChoice == 2) {
        double fVal = -1.;
        while( fVal < 0){
           bufferPtGen = -1.;
           bufferPtGen = pTfun->GetRandom();
           fVal = pTfun->Eval(bufferPtGen);
         }
        pTplus.push_back(bufferPtGen);
      }
      //Now azimuth
      azimuth_plus.push_back(-TMath::Pi()+(myRand->Rndm()*TMath::TwoPi()));
      //Neutrons without nuclear pT kick
      Px->push_back( pTplus.back()*TMath::Cos(azimuth_plus.back()) );
      Py->push_back( pTplus.back()*TMath::Sin(azimuth_plus.back()) );
      Pz->push_back( sqrt(pow(Ptot,2)-pow(Px->back(),2)-pow(Py->back(),2)) );
      E ->push_back( sqrt(pow(Ptot,2)+pow(neutron_mass,2)) );
      //Now adding the contribution of pT nuclear
      Px->back() += pTNuclearComponents[0];
      Py->back() += pTNuclearComponents[1];
      Pz->back() = sqrt(pow(Ptot,2)-pow(Px->back(),2)-pow(Py->back(),2));

      // Push the neutron into the primary particle vector
      fPrimaryVec.push_back( new G4PrimaryParticle( particleDefinition,
                                                    Px->back()*GeV,
                                                    Py->back()*GeV,
                                                    Pz->back()*GeV,
                                                    E->back()*GeV ) );

  	}//End of loop on particles of ZDC + (A)

    m_analysisManager->FillEventGenTree( particles,           // Number of neutrons on side A
                                         fCollisionPt,        // Spectator pt from collision
                                         fFragmentationPt,    // Additional pt from fragmentation
                                         reactionPlaneAngle );// Reaction plane angle (Psi)

}
