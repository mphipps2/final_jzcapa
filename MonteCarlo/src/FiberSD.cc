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

#include "FiberSD.hh"
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

FiberSD::FiberSD(G4String sdName, SharedData* sd, G4int modNum)
  :G4VSensitiveDetector(sdName), m_sd(sd), m_modNum(modNum) {
  collectionName.insert(sdName);
  HCID = -1; 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FiberSD::~FiberSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::HistInitialize(){
  std::string name = GetName();
  /*  
  // Add some histograms
  h2_radNum_eDep = new TH2D( Form("h2_radNum_eDep_%s", name.c_str() ),
				  ";rad number;eDep [keV]",
				  14,0,14,
				  50,0,25);
  m_sd->AddOutputHistogram( h2_radNum_eDep );
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::Initialize(G4HCofThisEvent* HCE){
	
  fiberCollection = new FiberHitsCollection(SensitiveDetectorName,
					      m_modNum); 

  std::string name = collectionName[0];					    
  
  //  static G4int HCID = -1;
  
  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID( name );}
  
  HCE->AddHitsCollection( HCID, fiberCollection );
  G4cout << " HCID " << HCID << " name " << name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FiberSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){
	

  TEnv* config = m_sd->GetConfig();

  
    G4int    modNum = m_modNum;
  
    //G4int modNStripsPerGap;
	//modNStripsPerGap = 64;
	
    char variable[256];
    sprintf(variable,"mod%dCoreIndexRefraction",6);
    m_modCoreIndexRefraction = config->GetValue( variable,1.46);
   
    

  G4double eDep   = aStep->GetTotalEnergyDeposit();
  
  G4int    totalRodNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
 
  G4String    radNum_s = "7" +  aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() + "0";
  
  G4int    rodNum;
  G4int    radNum;
   
    //radNum = totalRodNum / modNStripsPerGap;
    //rodNum = totalRodNum % modNStripsPerGap;
	

	rodNum = totalRodNum;
	
	radNum = (std::stoi (radNum_s)*100)+rodNum; 
	
	/*
	G4cout << G4endl << "TrackID: " << aStep->GetTrack()->GetTrackID() << G4endl;
	G4cout << "RodNum = " << rodNum << G4endl;
	G4cout << "physvol " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
	G4cout << "logvol  " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName() << G4endl;
	*/
  
  G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();
  G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentum();
  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();
  G4double charge = aStep->GetPreStepPoint()->GetCharge();
  

  int capturedPhotons = CalculateCherenkovs(aStep, modNum);

  FiberHit* newHit = new FiberHit();

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
  fiberCollection->insert (newHit );
  

  return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int FiberSD::CalculateCherenkovs(G4Step* aStep,__attribute__((unused)) int modNum) {
 
	
  const G4DynamicParticle* aParticle = aStep->GetTrack()->GetDynamicParticle();
  const G4double           charge    = aParticle->GetDefinition()->GetPDGCharge();
  //  LogStream<<MSG::INFO<<"ZDC stripCharge "<< charge << endreq;

  
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
  
  G4double MeanNumberOfPhotons = 2*TMath::Pi()*(1./137.)*step_length*inverseWLDiff*(charge)*(charge)*(1.0 - 1.0/(beta*m_modCoreIndexRefraction*beta*m_modCoreIndexRefraction));
 
  if (MeanNumberOfPhotons <= 0.0) { return false; }
  //  std::cout << " n photons " << MeanNumberOfPhotons << std::endl;


  const G4int NumPhotons = (G4int)G4Poisson(MeanNumberOfPhotons);

  if (NumPhotons <= 0) { return false; }  

    const G4ThreeVector pos = pPreStepPoint->GetPosition();
  //  float yPos              = pos.y();
  const G4ThreeVector p0  = aStep->GetDeltaPosition().unit();

  const float BetaInverse = 1./beta;
 
  float coreIndexRefraction;
  float claddingIndexRefraction;
  std::string modCladding;
  bool cladding;
  
    char name[256];
    sprintf(name,"mod%dCoreIndexRefraction",6);
    coreIndexRefraction = config->GetValue(name,1.46) ;
    sprintf(name,"mod%dCladdingIndexRefraction",6);
    claddingIndexRefraction = config->GetValue(name,1.43) ;
    sprintf(name,"mod%dCladding",6);
    modCladding = config->GetValue(name,"true") ;
    cladding = modCladding == "true" ? true : false;
    if (!cladding) claddingIndexRefraction = 1.; 
  

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

    if      (PY<=0) Transmission = 1;         // If photon is travelling with -ve PY, Its not going to reach the PMT
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

void FiberSD::EndOfEvent(G4HCofThisEvent*)
{

  //  G4int NbHits = fiberCollection->entries();

  /*
  if(verboseLevel>0) {
      std::cout << " if verbose loop" << std::endl;
      std::cout << "\n-------->Hits Collection: in this event they are " << NbHits 
		<< " hits in the calorimeter cells: " << std::endl;
      for (G4int i=0;i<NbHits;i++) {
	if (i %100 == 0) std::cout << " i " << i << std::endl;
	(*fiberCollection)[i]->Print();
      }
  }
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

