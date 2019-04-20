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

#include "G4UnitsTable.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "CherenkovHit.hh"
#include "QuartzHit.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "SharedData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
  hitsCollID = -1;
  
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{

  fEventNo = evt->GetEventID();
  //  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  hitsCollID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  QuartzHitsCollection* HC = 0;
  G4int nCollections =  HCE->GetNumberOfCollections();
  int totalPhotons = 0;
  int IDholder = 0;
  if(HCE) {
    while (hitsCollID < nCollections) {
      HC = (QuartzHitsCollection*)(HCE->GetHC(hitsCollID));
      int n_hit = HC->entries();
	  IDholder  = hitsCollID;
      int prevTrackId = 0;
      int prevRadiatorNo = 0;
      std::cout  << " hitsCollId " << hitsCollID << " nHits " << n_hit << std::endl;
      for ( int i = 0 ; i < n_hit; i++){
     
        G4int         radiatorNo  = (*HC)[i]->getRadNb();
		G4int         rodNo  = (*HC)[i]->getRodNb();
        G4double      eDep =  (*HC)[i]->getEdep();
        G4int         modNb =  (*HC)[i]->getModNb();
        G4int         trackID =  (*HC)[i]->getTrackID();
        G4ThreeVector position = (*HC)[i]->getPos();
        G4ThreeVector momentum = (*HC)[i]->getMomentum();
        G4double      energy   = (*HC)[i]->getEnergy();
        G4int         pid      = (*HC)[i]->getParticle()->GetPDGEncoding();
		G4int         nCherenkovs = (*HC)[i]->getNCherenkovs();
        G4double      charge   = (*HC)[i]->getCharge();
        G4double      velocity   = (*HC)[i]->getVelocity();
        G4double      beta   = (*HC)[i]->getBeta();
      
	/*		
         G4cout << "  trackID: " << trackID
	 << "  mod: " << modNb << "  rad: " << radiatorNo
	 << "  energy deposit: " << G4BestUnit(eDep,"Energy")
	 << "  position: " << G4BestUnit(position,"Length") <<
	   " beta " << beta << " velocity " << velocity << G4endl;
	*/

        //Add energy from every step in scoring volume as well as every particle
	if (trackID != prevTrackId || radiatorNo != prevRadiatorNo || eDep != 0) {

	  /*
	    G4cout << "PASSED  CUT!!!!!!!!!!!!!!  trackID: " << trackID
	    << "  mod: " << modNb << "  rad: " << radiatorNo
	    << "  energy deposit: " << G4BestUnit(eDep,"Energy")
	    << "  position: " << G4BestUnit(position,"Length") << G4endl;
	  */
	  
	  if(IDholder==2){ //2 corresponds to the RPD hitsCollID
	  fRunAction->SetRadNo_rpd(radiatorNo);
	  fRunAction->SetRodNo_rpd(rodNo);
	  fRunAction->SetNCherenkovs_rpd(nCherenkovs);
	  fRunAction->SetEdep_rpd(eDep);
	  fRunAction->SetModNb_rpd(modNb);
	  fRunAction->SetTrackID_rpd(trackID);
	  fRunAction->SetPosition_rpd(position);	 
	  fRunAction->SetMomentum_rpd(momentum);
	  fRunAction->SetEnergy_rpd(energy);
	  fRunAction->SetPid_rpd(pid);
	  fRunAction->SetEventNo_rpd(fEventNo);
	  fRunAction->SetCharge_rpd(charge);
	  fRunAction->SetVelocity_rpd(velocity);
	  fRunAction->SetBeta_rpd(beta);}
	  else{
	  totalPhotons += nCherenkovs; //not being used currently
	  fRunAction->SetRadNo(radiatorNo);
	  fRunAction->SetRodNo(rodNo);
	  fRunAction->SetNCherenkovs(nCherenkovs);
	  fRunAction->SetEdep(eDep);
	  fRunAction->SetModNb(modNb);
	  fRunAction->SetTrackID(trackID);
	  fRunAction->SetPosition(position);	 
	  fRunAction->SetMomentum(momentum);
	  fRunAction->SetEnergy(energy);
	  fRunAction->SetPid(pid);
	  fRunAction->SetEventNo(fEventNo);
	  fRunAction->SetCharge(charge);
	  fRunAction->SetVelocity(velocity);
	  fRunAction->SetBeta(beta);}
	}
	prevTrackId = trackID;
	prevRadiatorNo = radiatorNo;
      }
      hitsCollID++;
    }
    
	
	fRunAction->GetSharedData()->GetRPDTree()->Fill();
	
	fRunAction->GetSharedData()->GetZDCTree()->Fill();
	
	
	
  }

  fRunAction->ClearVectors();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
