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
//
// $Id: CherenkovSD.cc,v 1.9 2006/06/29 17:48:27 gunter Exp $
// GEANT4 tag $Name: geant4-08-01-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CherenkovSD.hh"
#include "SharedData.hh"

#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "TString.h"

#include <string>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CherenkovSD::CherenkovSD(G4String sdName, SharedData* sd, int modNumber)
  :G4VSensitiveDetector(sdName), m_sd(sd){
  collectionName.insert(sdName);
  HCID = -1;
  fModNumber = modNumber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CherenkovSD::~CherenkovSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CherenkovSD::HistInitialize(){

  //  sd->AddOutputHistogram( h_radNum_nParticles );
  /*
  m_sd->AddIntToTree("ID",&v_fTrackID);
  m_sd->AddIntToTree("ModNb",&v_fModNb);
  m_sd->AddIntToTree("RadNb",&v_fRadNb);
  m_sd->AddDoubleToTree("EDep",&v_fEdep);
  m_sd->AddIntToTree("Pid",&v_fPid);
  m_sd->AddDoubleToTree("X",&v_fX);
  m_sd->AddDoubleToTree("Y",&v_fY);
  m_sd->AddDoubleToTree("Z",&v_fZ);
  m_sd->AddDoubleToTree("Px",&v_fPx);
  m_sd->AddDoubleToTree("Py",&v_fPy);
  m_sd->AddDoubleToTree("Pz",&v_fPz);
  //  m_sd->AddDoubleToTree< std::vector<G4int> >("EventNo",&v_fEventNo);
  m_sd->AddDoubleToTree("Energy",&v_fEnergy);
  m_sd->AddDoubleToTree("Charge",&v_fCharge);
  */

  /*  
  m_sd->AddOutputToTree("ID",&fTrackID);
  m_sd->AddOutputToTree("ModNb",&fModNb);
  m_sd->AddOutputToTree("RadNb",&fRadNb);
  m_sd->AddOutputToTree("EDep",&fEdep);
  m_sd->AddOutputToTree("Pid",&fPid);
  m_sd->AddOutputToTree("X",&fX);
  m_sd->AddOutputToTree("Y",&fY);
  m_sd->AddOutputToTree("Z",&fZ);
  m_sd->AddOutputToTree("Px",&fPx);
  m_sd->AddOutputToTree("Py",&fPy);
  m_sd->AddOutputToTree("Pz",&fPz);
  m_sd->AddOutputToTree("EventNo",&fEventNo);
  m_sd->AddOutputToTree("Energy",&fEnergy);
  m_sd->AddOutputToTree("Charge",&fCharge);
  */
  // G4String is strange...
  //  std::string name = GetName();
  /*  
  // Add some histograms
  h2_rodNum_eDep = new TH2D( Form("h2_rodNum_eDep_%s", name.c_str() ),
				  ";rod number;eDep [keV]",
				  14,0,14,
				  50,0,25);
  m_m_sd->AddOutputHistogram( h2_rodNum_eDep );
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Note one hit collection per module
void CherenkovSD::Initialize(G4HCofThisEvent* HCE){
  radiatorCollection = new CherenkovHitsCollection(SensitiveDetectorName,
					      fModNumber); 

  std::string name = collectionName[0];					    
  
  std::cout << "collection name " << fModNumber << " sd name " << SensitiveDetectorName << std::endl;
  // Get Config
  TEnv* config = m_sd->GetConfig();
  int nModules = config->GetValue("nModules",4);
  char variable[256];
  for (int i = 0; i < nModules; ++i) {   
    sprintf(variable,"mod%dNRadiators",i+1);
    fModNRadiators[i] = config->GetValue(variable,12);
  }
  
  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID( name );
      std::cout << " HCID " << HCID << std::endl;
    }
  HCE->AddHitsCollection( HCID, radiatorCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// ProcessHits called each step in scoring logical volume which the sensitive detector is set
G4bool CherenkovSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

  G4double eDep   = aStep->GetTotalEnergyDeposit();
  G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();
  G4double velocity = aStep->GetPreStepPoint()->GetVelocity();
  G4double beta = aStep->GetPreStepPoint()->GetBeta();
  G4ThreeVector momentum = aStep->GetPreStepPoint()->GetMomentum();
  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();
  //Note: copy number index refers to the generation it belongs to. 0 is child (radiator, quartz, tungsten), 1 is the chamber, 2 is the "module" 
  //  G4int    modNum =aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);
  //this refers to the physical modules. ie) if you set 2 modules in the config file fModNumber goes from 0-1
  G4int    modNum = fModNumber;
  G4int    radNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
  /*
    std::cout << " copy 0 " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0) << " copy 1 " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1) << " copy 2 " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2) << " copy 3 " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3) << " fmod " << fModNumber << std::endl; */
  G4double charge = aStep->GetPreStepPoint()->GetCharge();
  //  std::cout << "mod num " << modNum << " radNum " << radNum << std::endl;

  //  h2_radNum_eDep->Fill( radNum, eDep );
  
  CherenkovHit* newHit = new CherenkovHit();
  //  int trackID = aStep->GetTrack()->GetTrackID();
  newHit->setCharge    (  charge );
  newHit->setTrackID   ( aStep->GetTrack()->GetTrackID() );
  newHit->setModNb     ( modNum );
  int rNum = radNum + (fModNRadiators[modNum] * modNum);
  newHit->setRadNb     ( rNum );
  //  std::cout << "rNum " << rNum << std::endl;
  newHit->setEdep      ( eDep );
  newHit->setParticle  (particle);
  newHit->setEnergy    (energy);
  newHit->setVelocity    (velocity);
  newHit->setBeta    (beta);
  newHit->setMomentum  (momentum);
  newHit->setPos       ( aStep->GetPostStepPoint()->GetPosition() );
  newHit->setEventNo   (m_sd->GetEventNo());
  //    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
    /*
    G4cout << "  trackID: " << trackID
	 << "  mod: " << modNum << "  rad: " << rNum
	 << "  energy deposit: " << G4BestUnit(eDep,"Energy")
	 << "  position: " << G4BestUnit(position,"Length") << G4endl;
    */

  //  /*  if (position.getX() < 0.0001 && position.getX() > -0.0001)*/  	std::cout << "xpos " << position.getX() << std::endl; 
  //          G4cout << "---- Hit # " << i << G4endl;
  //  if (particle->GetPDGEncoding() == 22 ) {
  /*
          G4cout << " Cherenkov No " <<  radNum << G4endl;
          G4cout << " Mod No "   <<  modNum << G4endl;
          G4cout << " eDep "     <<  eDep << G4endl;
	  G4cout << " etotal "     <<  energy << G4endl;
	  G4cout << " pid "     <<  particle->GetPDGEncoding() << G4endl;
	  // }
	  */
  // Note: when crossing volume boundaries there the post step occurs at the edge of volume 1 (of step 1) and the subsequent prestep occurs at the edge of volume 2 (of step 2). In this way, the prestep of something in the SD function always occurs inside the SD. The poststep could then be in the SD or another volume
  //  std::cout << "pre " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() << " post " << aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() << std::endl;
  /*  if (aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() == "Emitter" && aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName() != "Emitter") {
   */

  radiatorCollection->insert( newHit );
    //}
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void CherenkovSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	
  
  //  G4String name = HCE->GetHC()->GetName();
  //if(HCID<0) return;
   std::cout << " HCID end of event " << HCID << std::endl;
  CherenkovHitsCollection* HC = 0;
	__attribute__((unused)) G4int nCollections =  HCE->GetNumberOfCollections();
  __attribute__((unused)) G4int hitsCollID = 0;

  
  if(HCE) {
    //    while (hitsCollID < nCollections) {
      HC = (CherenkovHitsCollection*)(HCE->GetHC(HCID));
      int prevTrackId = 0;
      int prevCherenkovNo = 0;
      //      G4int NbHits = radiatorCollection->entries();
      int n_hits = HC->entries();
            std::cout << " n_hits " << n_hits << std::endl;
      if(verboseLevel>0) { 
        G4cout << "\n-------->Hits Collection: in this event they are " << n_hits 
	   << " hits in the calorimeter cells: " << G4endl;
      }

      for (G4int i=0;i<n_hits;i++) {
        int         radiatorNo  = (*HC)[i]->getRadNb();
        double      eDep =  (*HC)[i]->getEdep();
        int         trackID =  (*HC)[i]->getTrackID();
        if (trackID != prevTrackId || radiatorNo != prevCherenkovNo || eDep != 0) {
	  // std::cout << " cond 1 " << (trackID != prevTrackId) << " cond 2 " << (radiatorNo != prevCherenkovNo) << " cond 3 " << (eDep != 0) << std::endl;
	 
          fRadNb = (*HC)[i]->getRadNb();
	  std::cout << " radNum post " << fRadNb << std::endl;
          fEdep = (*HC)[i]->getEdep();
          fModNb = (*HC)[i]->getModNb();
          fTrackID = (*HC)[i]->getTrackID();
          fX = (*HC)[i]->getPos().getX();
          fY = (*HC)[i]->getPos().getY();
          fZ = (*HC)[i]->getPos().getZ();
          fPx = (*HC)[i]->getMomentum().getX();
          fPy = (*HC)[i]->getMomentum().getY();
          fPz = (*HC)[i]->getMomentum().getZ();
          fEnergy = (*HC)[i]->getEnergy();
          fPid = (*HC)[i]->getParticle()->GetPDGEncoding();
          fCharge = (*HC)[i]->getCharge();	  
	  //	  std::cout << " print " << std::endl;
	  //      (*radiatorCollection)[i]->Print();
	  m_sd->GetTree()->Fill();
	}
        prevTrackId = trackID;
        prevCherenkovNo = radiatorNo;
      }
     
      //hitsCollID++;
      // }
  }
  //  std::cout << " vec size " << v_fRadNb.size() << std::endl;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

