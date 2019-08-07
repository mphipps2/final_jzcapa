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
#include "G4PrimaryVertex.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
  hitsCollID = -1;
  fRunAction->GetSharedData()->AddOutputToRPDTree("gunPosX",&gunPosX);
  fRunAction->GetSharedData()->AddOutputToRPDTree("gunPosY",&gunPosY);
  fRunAction->GetSharedData()->AddOutputToRPDTree("gunPosZ",&gunPosZ);

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
  G4PrimaryVertex* pVert = evt->GetPrimaryVertex();
  gunPosX = pVert->GetX0();
  gunPosY = pVert->GetY0();
  gunPosZ = pVert->GetZ0();

  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  QuartzHitsCollection* HC = 0;
  G4int nCollections =  HCE->GetNumberOfCollections();


 //////////// NEED TO KNOW ORDER OF DETECTORS TO ASSIGN CORRECT HITCOLLECTION ////////////
	std::string detector[3];
	__attribute__((unused)) int ZDC1=-1;
	__attribute__((unused)) int ZDC2=-1;
	__attribute__((unused)) int RPD=-1;

	bool bzdc1flag=false;
	bool bzdc2flag=false;
	bool brpdflag=false;


  TEnv* config = fRunAction->GetSharedData()->GetConfig();
	int runNum = config->GetValue( "RunNumber", -1);
  bool OPTICAL = config->GetValue("OPTICAL_ON", false);
  bool CLUSTER = config->GetValue("CLUSTER_ON", false);

	fRunAction->GetSharedData()->LoadAlignmentFile(runNum);

	Alignment	*align_run 	= fRunAction->GetSharedData()->GetAlignment();

	detector[0]=align_run->upstream_Det;
	detector[1]=align_run->mid_Det;
	detector[2]=align_run->downstream_Det;

	for(int i=0; i<3; i++){
		if(detector[i]=="ZDC1") {
			bzdc1flag=true;}
		if(detector[i]=="ZDC2") {
			bzdc2flag=true;}
		if(detector[i]=="RPD") {
			brpdflag=true;}
	}

	if(bzdc1flag){
		ZDC1=0;
	}
	if(bzdc2flag){
		if(bzdc1flag){
		ZDC2=1;
		}
		else{
		ZDC2=0;
		}
	}
	if(brpdflag){
		if(bzdc1flag && bzdc2flag){
		RPD=2;
		}
		else if(bzdc1flag && !bzdc2flag){
		RPD=1;
		}
		else if(!bzdc1flag && bzdc2flag){
		RPD=1;
		}
		else RPD=0;
	}

  /////////////////////////////////////////////////////////////////////////

  int totalPhotons = 0;
  std::vector<int> Gap_Cherenk;
  Gap_Cherenk.resize(24,0);

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

        G4int         radiatorNo    = (*HC)[i]->getRadNb();
		    G4int         rodNo         = (*HC)[i]->getRodNb();
        G4double      eDep          = (*HC)[i]->getEdep();
        G4int         modNb         = (*HC)[i]->getModNb();
        G4int         trackID       = (*HC)[i]->getTrackID();
        G4ThreeVector position      = (*HC)[i]->getPos();
        G4ThreeVector momentum      = (*HC)[i]->getMomentum();
        G4double      energy        = (*HC)[i]->getEnergy();
        G4int         pid           = (*HC)[i]->getParticle()->GetPDGEncoding();
		    G4int         nCherenkovs   = (*HC)[i]->getNCherenkovs();
        G4double      charge        = (*HC)[i]->getCharge();
        G4double      velocity      = (*HC)[i]->getVelocity();
        G4double      beta          = (*HC)[i]->getBeta();

        //Add energy from every step in scoring volume as well as every particle
	if (trackID != prevTrackId || radiatorNo != prevRadiatorNo || eDep != 0) {

	  if(IDholder==RPD){ //corresponds to the RPD hitsCollID
        if(!CLUSTER){
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
  	  fRunAction->SetRodNo_rpd(rodNo);
  	  fRunAction->SetPosition_rpd(position);}
    }
	  else if(IDholder==RPD+1 && RPD!=-1 && !OPTICAL && !CLUSTER){ //RPD+1 corresponds to the Fiber hitsCollID
	  fRunAction->SetRadNo_fiber(radiatorNo);
	  fRunAction->SetRodNo_fiber(rodNo);
	  fRunAction->SetNCherenkovs_fiber(nCherenkovs);
	  fRunAction->SetEdep_fiber(eDep);
	  fRunAction->SetModNb_fiber(modNb);
	  fRunAction->SetTrackID_fiber(trackID);
	  fRunAction->SetPosition_fiber(position);
	  fRunAction->SetMomentum_fiber(momentum);
	  fRunAction->SetEnergy_fiber(energy);
	  fRunAction->SetPid_fiber(pid);
	  fRunAction->SetEventNo_fiber(fEventNo);
	  fRunAction->SetCharge_fiber(charge);
	  fRunAction->SetVelocity_fiber(velocity);
	  fRunAction->SetBeta_fiber(beta);  }
	  else if(IDholder==ZDC1 || IDholder==ZDC2 ){//ZDC hitsCollID
        if(!CLUSTER){
      totalPhotons += nCherenkovs;
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
        else{
          //std::cout << "radiatorNo = " << radiatorNo << ", modNb = " << modNb << ", radiatorNo+(modNb*12) = " << radiatorNo+(modNb*12) << std::endl;
          Gap_Cherenk.at(radiatorNo+(modNb*12)) += nCherenkovs;}
      }
}
	prevTrackId = trackID;
	prevRadiatorNo = radiatorNo;
      }
      hitsCollID++;
    }//end of hit loop



if((bzdc1flag || bzdc2flag) && CLUSTER) {
  fRunAction->SetGapCherenkovs(Gap_Cherenk);
  Gap_Cherenk.clear();
}

  fRunAction->GetSharedData()->GetRPDTree()->Fill();

	fRunAction->GetSharedData()->GetZDCTree()->Fill();

	fRunAction->GetSharedData()->GetFiberTree()->Fill();

  }

  fRunAction->ClearVectors();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
