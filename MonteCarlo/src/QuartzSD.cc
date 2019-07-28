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
// Author: Michael Phipps

#include "QuartzSD.hh"
#include "SharedData.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"

#include "TString.h"

#include <string>
#include <iostream>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QuartzSD::QuartzSD(G4String sdName, SharedData* sd, G4int modNum)
  :G4VSensitiveDetector(sdName), m_sd(sd), m_modNum(modNum) {
  collectionName.insert(sdName);
  HCID = -1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

QuartzSD::~QuartzSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QuartzSD::HistInitialize(){
  std::string name = GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QuartzSD::Initialize(G4HCofThisEvent* HCE){
  quartzCollection = new QuartzHitsCollection(SensitiveDetectorName,
					      m_modNum);

  std::string name = collectionName[0];

  //  static G4int HCID = -1;

  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID( name );}

  HCE->AddHitsCollection( HCID, quartzCollection );
  std::cout << " HCID " << HCID << " name " << name << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool QuartzSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

  TEnv* config = m_sd->GetConfig();

  G4int    modNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);


  char modName[256];
  sprintf(modName,"mod%dType",modNum+1);
  G4int    nModules            = config->GetValue( "nModules",4);
  G4int modType[5]; G4int modNStripsPerGap[5];
  //G4int modNRadiators[5];
  //  G4double coreDiameter[5]; G4double claddingThickness[5]; G4double modCladdingIndexRefraction[5];

  modType[0]            = config->GetValue( "modType1",5);
  modType[1]            = config->GetValue( "modType2",5);
  modType[2]            = config->GetValue( "modType3",3);
  modType[3]            = config->GetValue( "modType4",3);
  modType[4]            = config->GetValue( "modType5",3);
  for (int i = 0; i < nModules; ++i) {
    char variable[256];
    sprintf(variable,"mod%dCoreIndexRefraction",i+1);
    m_modCoreIndexRefraction[i] = config->GetValue( variable,1.46);
    sprintf(variable,"mod%dNStripsPerGap",i+1);
    modNStripsPerGap[i] = config->GetValue(variable,52);
    //    sprintf(variable,"mod%dNRadiators",i+1);
    //modNRadiators[i] = config->GetValue(variable,12);
    //    sprintf(variable,"mod%dCladdingIndexRefraction",i+1);
    //    modCladdingIndexRefraction[i] = config->GetValue( variable,1.43);
  }

  G4double eDep   = aStep->GetTotalEnergyDeposit();

  G4int    totalRodNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
  G4int    radNum;
  G4int    rodNum;
  if (modType[modNum] < 4) {
    radNum = totalRodNum / 52;
    rodNum = totalRodNum % 52 ;
  }
  else {
    radNum = totalRodNum / modNStripsPerGap[modNum];
    rodNum = totalRodNum % modNStripsPerGap[modNum];
  }

  G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();
  G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentum();
  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();
  G4double charge = aStep->GetPreStepPoint()->GetCharge();

  int capturedPhotons = CalculateCherenkovs(aStep,modNum);

  QuartzHit* newHit = new QuartzHit();

  newHit->setCharge        ( charge );
  newHit->setTrackID       (aStep->GetTrack()->GetTrackID() );
  newHit->setModNb         (modNum );
  newHit->setRadNb         (radNum );
  newHit->setRodNb         (rodNum );
  newHit->setEdep          (eDep );
  newHit->setPos           (aStep->GetPostStepPoint()->GetPosition() );
  newHit->setParticle      (particle);
  newHit->setEnergy        (energy);
  newHit->setMomentum      (momentum);
  newHit->setNCherenkovs   (capturedPhotons);
  //  if (capturedPhotons != 0 ) std::cout << " capturedPhotons " << capturedPhotons << std::endl;
  quartzCollection->insert (newHit );

  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int QuartzSD::CalculateCherenkovs(G4Step* aStep,int modNum) {


  const G4DynamicParticle* aParticle = aStep->GetTrack()->GetDynamicParticle();
  const G4double           charge    = aParticle->GetDefinition()->GetPDGCharge();

  if (charge==0) { return false; }

  const G4StepPoint*       pPreStepPoint  = aStep->GetPreStepPoint();
  const G4StepPoint*       pPostStepPoint = aStep->GetPostStepPoint();
  const G4double           beta           = (pPreStepPoint ->GetBeta() + pPostStepPoint->GetBeta())/2.;

  if (beta==0) { return false; }

  TEnv* config = m_sd->GetConfig();

  G4double minWavelength     = config->GetValue("cherenkovMinWavelength",250) * CLHEP::nanometer ;
  G4double maxWavelength     = config->GetValue("cherenkovMaxWavelength",600) * CLHEP::nanometer ;
  G4double inverseWLDiff     = (1./minWavelength) - (1./maxWavelength);


  const float step_length(aStep->GetStepLength());

  G4double MeanNumberOfPhotons = 2*TMath::Pi()*(1./137.)*step_length*inverseWLDiff*(charge)*(charge)*(1.0 - 1.0/(beta*m_modCoreIndexRefraction[modNum]*beta*m_modCoreIndexRefraction[modNum]));

  if (MeanNumberOfPhotons <= 0.0) { return false; }
  //  std::cout << " n photons " << MeanNumberOfPhotons << std::endl;


  const G4int NumPhotons = (G4int)G4Poisson(MeanNumberOfPhotons);

  if (NumPhotons <= 0) { return false; }

    const G4ThreeVector pos = pPreStepPoint->GetPosition();
  //  float yPos              = pos.y();
  const G4ThreeVector p0  = aStep->GetDeltaPosition().unit();

  const float BetaInverse = 1./beta;
  char modName[256];
  sprintf(modName,"mod%dType",modNum+1);
  int modType             = config->GetValue(modName,1) ;
  float coreIndexRefraction;
  float claddingIndexRefraction;
  std::string modCladding;
  bool cladding;
  if (modType < 4) {
    coreIndexRefraction = 1.45;  // fused quartz in current ZDC
    claddingIndexRefraction = 1; // air in current ZDC
  }
  else {
    char name[256];
    sprintf(name,"mod%dCoreIndexRefraction",modNum+1);
    coreIndexRefraction = config->GetValue(name,1.46) ;
    sprintf(name,"mod%dCladdingIndexRefraction",modNum+1);
    claddingIndexRefraction = config->GetValue(name,1.43) ;
    sprintf(name,"mod%dCladding",modNum+1);
    modCladding = config->GetValue(name,"true") ;
    cladding = modCladding == "true" ? true : false;
    if (!cladding) claddingIndexRefraction = 1.;
  }

  float criticalAngle = asin(claddingIndexRefraction/coreIndexRefraction);
  //  std::cout << " mod " << modNum  <<" core n " << coreIndexRefraction <<  " cladding n " << claddingIndexRefraction << " critical angle " << criticalAngle << std::endl;
  int photonCount = 0;

  for (G4int I = 0; I < NumPhotons; I++) {
    G4double rand;
    //    float sampledEnergy;
    float cosTheta, sin2Theta;

    // sample an energy for Photon
    rand = G4UniformRand();
    // sampledEnergy = Pmin + rand * dp;

    cosTheta  = BetaInverse / coreIndexRefraction;
    sin2Theta = (1.0 - cosTheta)*(1.0 + cosTheta);

    // Generate random position of photon on cone surface defined by Theta
    rand = G4UniformRand();
    const float phi = 2.*M_PI*rand;
    const float sinPhi = sin(phi);
    const float cosPhi = cos(phi);

    // calculate x,y, and z components of photon momentum
    // (in coord system with primary particle direction aligned with the z axis)
    // and Rotate momentum direction back to global reference system
    const float sinTheta = sqrt(sin2Theta);
    const float px = sinTheta*cosPhi;
    const float py = sinTheta*sinPhi;
    const float pz = cosTheta;
    G4ParticleMomentum photonMomentum(px, py, pz);
    photonMomentum.rotateUz(p0);
    bool Transmission=0;
    const float PT = sqrt(photonMomentum.getX()*photonMomentum.getX() + photonMomentum.getZ()*photonMomentum.getZ());
    const float PY = photonMomentum.getY();

    if      (PY<=0) Transmission = 0;         // If photon is travelling with -ve PY, Its not going to reach the PMT
    else if (PT==0) Transmission = 1; // if photon is along Y-axis it will be completely transmitted into the PMT
    else {
      const float Theta = M_PI/2.0-atan(PT/PY);
      if (Theta < criticalAngle) Transmission = 0;
      else Transmission=1.0;
    }
    if (Transmission == 1.0) ++photonCount;
  }
  //  std::cout << " totalPhotonsCaptured " << photonCount << " % " << (float) photonCount/NumPhotons << " critical angle " << criticalAngle <<  std::endl;
  return photonCount;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void QuartzSD::EndOfEvent(G4HCofThisEvent*)
{

  //  G4int NbHits = quartzCollection->entries();

  /*
  if(verboseLevel>0) {
      std::cout << " if verbose loop" << std::endl;
      std::cout << "\n-------->Hits Collection: in this event they are " << NbHits
		<< " hits in the calorimeter cells: " << std::endl;
      for (G4int i=0;i<NbHits;i++) {
	if (i %100 == 0) std::cout << " i " << i << std::endl;
	(*quartzCollection)[i]->Print();
      }
  }
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
