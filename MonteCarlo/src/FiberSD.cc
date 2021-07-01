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
/// \ingroup mc
/// \file FiberSD.cc
/// \author Michael Phipps

#include "FiberSD.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"

#include <TRandom3.h>

#include <string>
#include <iostream>
#include <cmath>

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FiberSD::FiberSD(G4String sdName, G4int modNum, G4bool optical, G4bool fullOptical)
  :G4VSensitiveDetector(sdName),
  m_modNum(modNum),
  FULLOPTICAL(fullOptical),
  OPTICAL(optical),
  REDUCED_TREE(false),
  ML_REDUCED_TREE(false),
  ZDC(false),
  RPD(false)
 {
  collectionName.insert(sdName);
  HCID = -1;
  m_isTIR = 0;
  m_TIR_count = 0;
  m_prevTrackID = -1;
  m_lastFiveCore = 1;
  m_avgYStepLength = 0;
  m_avgStepLength = 0;
  m_avgTimeLength = 0;
  m_firstStepTime = 0;
  is_upward_light = 0;
  is_downward_light = 0;
  if( sdName.contains("R") ) RPD = true;
  if( sdName.contains("Z") ) ZDC = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FiberSD::~FiberSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::HistInitialize(){
  std::string name = GetName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::Initialize(G4HCofThisEvent* HCE){
  fiberCollection = new FiberHitsCollection(SensitiveDetectorName,
					      m_modNum);

  std::string name = collectionName[0];
  if(REDUCED_TREE){
    // Grab the vector from AnalysisManager if reduced tree mode is active
    AnalysisManager* analysisManager = AnalysisManager::getInstance();
    m_cherenkovVec = analysisManager->GetFiberVector(ZDC,RPD,m_modNum);
    m_timeVec      = analysisManager->GetTimeVector(ZDC,RPD,m_modNum);

    // Clear it and resize in case this isn't the first event
    m_cherenkovVec->clear();
    m_cherenkovVec->resize( m_nFibers, 0 );
    m_timeVec->clear();
    m_timeVec->resize( 128*m_nChannels, 0 );
  }

  else if(ML_REDUCED_TREE){
    // Grab the vector from AnalysisManager if reduced tree mode is active
    AnalysisManager* analysisManager = AnalysisManager::getInstance();
    m_cherenkovVec = analysisManager->GetFiberVector(ZDC,RPD,m_modNum);
    m_timeVec      = analysisManager->GetTimeVector(ZDC,RPD,m_modNum);
    // Clear it and resize in case this isn't the first event
    m_cherenkovVec->clear();
    m_cherenkovVec->resize( m_nChannels, 0 );
    if (RPD) {
      // Clear it and resize in case this isn't the first event
      m_incidenceAngleVec = analysisManager->GetIncidenceAngleVector(m_modNum);
      m_incidenceAngleVec->clear();
      // Clear it and resize in case this isn't the first event
      m_yOriginVec = analysisManager->GetYOriginVector(m_modNum);
      m_yOriginVec->clear();
      // Clear it and resize in case this isn't the first event
      m_genCherenkovVec = analysisManager->GetFiberGenVector(m_modNum);
      m_genCherenkovVec->clear();
      m_genCherenkovVec->resize( m_nChannels, 0 );
      m_timeVec->clear();
      m_timeVec->resize( 128*m_nChannels, 0 );
    }
  }
  else {
    if(m_cherenkovVec){
      AnalysisManager* analysisManager = AnalysisManager::getInstance();
      m_cherenkovVec = analysisManager->GetFiberVector(ZDC,RPD,m_modNum);

      m_cherenkovVec->clear();
      m_cherenkovVec->resize( m_nFibers, 0 );
    }
  }
  
  m_nCherenkovs = 0;
  m_nHits = 0;

  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID( name );}

  HCE->AddHitsCollection( HCID, fiberCollection );
  G4cout << " HCID " << HCID << " name " << name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FiberSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();

  //Get the number of Cherenkov photons created in this step
  int generatedPhotons = 0;
  const std::vector<const G4Track*>* secVec = aStep->GetSecondaryInCurrentStep();
  for(uint i = 0; i < secVec->size(); i++){
    if( secVec->at(i)->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
      generatedPhotons++;
    }//end if photon
  }//end secondary track loop
  m_nCherenkovs += generatedPhotons; 
  
  //Don't record hits that didn't produce cherenkov photons
  // if(capturedPhotons == 0) return true;

  G4int rodNum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);

  G4ThreeVector pos = aStep->GetTrack()->GetPosition();
  //  G4ParticleDefinition *particle = aStep->GetTrack()->GetDefinition();

  // introducing a generated Cherenkov vector to look at capture efficiency event-by-event
  if (ML_REDUCED_TREE && RPD) {
    G4int channelNum;
    channelNum = GetRPDChannelMapping(rodNum);
    m_genCherenkovVec->at(channelNum)+=generatedPhotons;
  }
  
  // If FULLOPTICAL is true, determine if the photon has reached the top of the topOfVolume
  // If OPTICAL is true, determine if the photon is captured (requires # of steps inside SD to exceed m_captureThreshold and for the final five steps to not include the Kapton buffer material)
  // and add the hit to the collection if it has
  G4bool KEEP_DOWNWARD_LIGHT = 1;
  G4bool KEEP_UPWARD_LIGHT = 1;
  int trackID = aStep->GetTrack()->GetTrackID();
  if (trackID != m_prevTrackID) {
    is_downward_light = 0;
    is_upward_light = 0;
    m_isTIR = 0;
  }

  if(OPTICAL){
    if (aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {

      G4ThreeVector p = aStep->GetTrack()->GetMomentumDirection();
      double pT = sqrt( pow( p.x(), 2.0 ) + pow( p.z(), 2.0 ));
      double theta = M_PI / 2.0 - atan( pT / p.y() );
      // std::cout << " track " << trackID << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " rodNum " << rodNum << " material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " m_TIR_count " << m_TIR_count  << std::endl;
      if (theta > M_PI/2) {
	is_downward_light = 1;
	// kill downward going photons
	if (!KEEP_DOWNWARD_LIGHT) { 
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	  return true;
	}
      }
      else {
	is_upward_light = 1;
	// kill upward going photons
	if (!KEEP_UPWARD_LIGHT) {
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	  return true;
	}
      }
      if (trackID == m_prevTrackID) {
	++m_TIR_count;
      }
      else {
	m_TIR_count = 1;
	m_lastFiveCore = 1;
	m_avgYStepLength = 0;
	m_avgStepLength = 0;
	m_avgStepTime = 0;
	m_firstStepTime = aStep->GetTrack()->GetGlobalTime();
      }
      if (m_TIR_count > m_captureThreshold - 5) {
	if (aStep->GetPostStepPoint()->GetMaterial()->GetName() == "Kapton_UI") {
	  m_lastFiveCore = 0;
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	  return true;
	}

	m_avgYStepLength += (aStep->GetPostStepPoint()->GetPosition().y() - aStep->GetPreStepPoint()->GetPosition().y());
	m_avgStepLength += aStep->GetStepLength();
	m_avgStepTime += aStep->GetDeltaTime();
	
      }
      // light considered captured if number of TIR reflections >= threshold number or light reaches top of readout
      //      if ((m_TIR_count >= m_captureThreshold  && m_lastFiveCore) || (pos.y() >= m_topOfVolume - 0.1*mm)) {
      if ((m_TIR_count >= m_captureThreshold  && m_lastFiveCore)) {

	// Get average vertical step length
	// here worked
	m_avgYStepLength /= 5;
	m_avgStepLength /= 5;
	m_avgStepTime /= 5;
	double nTotalSteps;
	if (m_avgYStepLength < 0) m_avgYStepLength *= -1;
	//	if (is_downward_light) std::cout << " M_AVGYSTEPLENGTH " << m_avgYStepLength << std::endl;
	// if downward going light, fand path length to go down and then back up
	if (is_downward_light) nTotalSteps = ((pos.y() - m_bottomOfVolume) / m_avgYStepLength) + ((m_topOfVolume - m_bottomOfVolume) / m_avgYStepLength);
	else nTotalSteps = (m_topOfVolume - pos.y()) / m_avgYStepLength;
	double totalPathLength = m_avgStepLength * nTotalSteps;
	double timeArrived = m_firstStepTime + (m_avgStepTime * nTotalSteps);
	double cherenkovEnergy = aStep->GetPreStepPoint()->GetTotalEnergy();
	// if above/below max/min value, use max/min value
	// if somewhere between, linear interpolation
	double absorpLength = aStep->GetPreStepPoint()->GetMaterial()->GetMaterialPropertiesTable()->GetProperty("ABSLENGTH")->Value(cherenkovEnergy);
	double survivalProb = exp((-1 * totalPathLength) / absorpLength);
	// initializing to 0 invokes random seed
	//	std::cout << " isDownwardLight " << is_downward_light << " totalPathLength " << totalPathLength << " timeArrived " << timeArrived << " nTotalSteps " << nTotalSteps << std::endl;
	TRandom3 *r3 = new TRandom3(0);
	G4bool survived = r3->Binomial(1,survivalProb);
	delete r3;
	//	std::cout << " avgYStep " << m_avgYStepLength << " avgStepLength " << m_avgStepLength << " nTotalSteps " << nTotalSteps << " totalPathLength " << totalPathLength << " energy " << cherenkovEnergy << " absorp Length " << absorpLength << " survivalProb " << survivalProb << " survived " << survived << std::endl;

	if (!survived) {
	  m_lastFiveCore = 0;
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	  return true;
	}

	if(REDUCED_TREE){
	  m_cherenkovVec->at(rodNum)++;
	  FillTimeVector( rodNum, timeArrived );
	  m_nHits++;
	  m_lastFiveCore = 0;
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
	  return true;
	}
	else if(ML_REDUCED_TREE){
	  G4int channelNum;
	  if (RPD) channelNum = GetRPDChannelMapping(rodNum);
	  else     channelNum = GetZDCChannelMapping(rodNum);
	  if (RPD) {
	    FillChannelTimeVector( channelNum, timeArrived );
	    m_yOriginVec->push_back(aStep->GetTrack()->GetVertexPosition().y());
	    
	    double incidenceAngle = GetIncidenceAngle(aStep);
	    // incidence angles between [0,180]. eg) if critical angle is 82 degrees, retained light between [82,98]
	    // normalize s.t. distribution is between [0,180]
	    if (incidenceAngle > 90) incidenceAngle = 180 - incidenceAngle;
	    // then flip the convention for downward going light so we can visualize the directional dependency
	    if (is_downward_light) incidenceAngle = 180 - incidenceAngle;
	    m_incidenceAngleVec->push_back(incidenceAngle);	 
	  }
	  m_cherenkovVec->at(channelNum)++;
	  m_nHits++;
	  m_lastFiveCore = 0;
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
	  return true;
	}
	// Full Tree Mode
	else { 
	  FiberHit* newHit = new FiberHit();
	  newHit->setPos      ( pos );
	  newHit->setOrigin   ( aStep->GetTrack()->GetVertexPosition() );
	  newHit->setMomentum ( aStep->GetPreStepPoint()->GetMomentum() );
	  newHit->setEnergy   ( aStep->GetPreStepPoint()->GetTotalEnergy() );
	  newHit->setTime     ( aStep->GetTrack()->GetGlobalTime() );
	  newHit->setRodNb    ( rodNum );

	  fiberCollection->insert ( newHit );
	  m_nHits++;

	  m_lastFiveCore = 0;
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
	  return true;
	}
      }
      // photon not captured (in some non-TIR reflection mode)
      else if (m_TIR_count >= m_captureThreshold && !m_lastFiveCore ) {
	aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
      }
      m_prevTrackID = trackID;;
    
    }
  }
  else if (FULLOPTICAL) {
    if( aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      
      //      int trackID = aStep->GetTrack()->GetTrackID();
      G4ThreeVector p = aStep->GetTrack()->GetMomentumDirection();
      double pT = sqrt( pow( p.x(), 2.0 ) + pow( p.z(), 2.0 ));
      double theta = M_PI / 2.0 - atan( pT / p.y() );
      //      std::cout << " track " << trackID << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " rodNum " << rodNum << " material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " m_TIR_count " << m_TIR_count  << std::endl;
      //      if (RPD) std::cout << " trackID "  << trackID << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " rodNum " << rodNum << " material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " m_TIR_count " << m_TIR_count  << std::endl;
      if (!is_downward_light && !is_upward_light) {
	if (theta > M_PI/2) {
	  is_downward_light = 1;
	  // kill downward going photons
	  if (!KEEP_DOWNWARD_LIGHT) {
	    aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	    return true;
	  }
	}
	else {
	  is_upward_light = 1;
	  // kill upward going photons
	  if (!KEEP_UPWARD_LIGHT) {
	    aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	    return true;
	  }
	}
      }
      /*
      double incidenceAngle = GetIncidenceAngle(aStep);
      // incidence angles between [0,180]. eg) if critical angle is 82 degrees, retained light between [82,98]
      // normalize s.t. distribution is between [0,180]
      if (incidenceAngle > 90) incidenceAngle = 180 - incidenceAngle;
      // then flip the convention for downward going light so we can visualize the directional dependency

      if (RPD) {
	//	if (aStep->GetPreStepPoint()->GetMaterial()->GetName() == "SilicaCore_UI"  &&   aStep->GetPostStepPoint()->GetMaterial()->GetName() == "SilicaClad_UI"  ) {
	if (aStep->GetPreStepPoint()->GetMaterial()->GetName() == "SilicaClad_UI"  &&   aStep->GetPostStepPoint()->GetMaterial()->GetName() == "SilicaCore_UI"  ) {
	  if (incidenceAngle >= 82.148) m_isTIR = 1;
	  if (incidenceAngle < 82.148) {
	    aStep->GetTrack()->SetTrackStatus( fStopAndKill );
	    return true;
	  }
	}
	if (m_isTIR) std::cout << " trackID "  << trackID << " isTIR " << m_isTIR << " incidence angle: " << incidenceAngle << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " prestep material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() <<  " poststep material " << aStep->GetPostStepPoint()->GetMaterial()->GetName() << " total energy " << aStep->GetTrack()->GetTotalEnergy() * eV << " m_TIR_count " << m_TIR_count  << std::endl;
      }
      */
	//if (is_downward_light) incidenceAngle = 180 - incidenceAngle;

      
      if (trackID == m_prevTrackID) {
	++m_TIR_count;
	/*
	if (m_TIR_count > 15 ) {

	  aStep->GetTrack()->SetTrackStatus( fStopAndKill );
	  return true;
	}
	*/
	//	if (RPD) std::cout << " trackID "  << trackID << " incidence angle: " << incidenceAngle << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " rodNum " << rodNum << " material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " m_TIR_count " << m_TIR_count  << std::endl;
      }
      else {
	m_TIR_count = 1;
      }

      
      if (pos.y() >= m_topOfVolume - 0.1*mm){




	
	/*
	if (RPD) { 
	  if (is_downward_light) std::cout << " DOWNWARD LIGHT MADE IT TO THE TOP! " << std::endl;
	  else std::cout << " upward light " << std::endl;
	}
	*/
	if(REDUCED_TREE){
	  m_cherenkovVec->at(rodNum)++;
	  FillTimeVector( rodNum, aStep->GetTrack()->GetGlobalTime() );
	  m_nHits++;
	  return true;
	}
	else if(ML_REDUCED_TREE){
	  G4int channelNum;
	  //	  std::cout << " filling vecs in FiberSD " << std::endl;
	  if (RPD) channelNum = GetRPDChannelMapping(rodNum);
	  else     channelNum = GetZDCChannelMapping(rodNum);
	  if (RPD) {
	    //	    if (aStep->GetPreStepPoint()->GetMaterial()->GetName() == "SilicaCore_UI"  &&   aStep->GetPostStepPoint()->GetMaterial()->GetName() == "SilicaClad_UI"  ) {
	    if (aStep->GetPreStepPoint()->GetMaterial()->GetName() == "SilicaClad_UI"  &&   aStep->GetPostStepPoint()->GetMaterial()->GetName() == "SilicaCore_UI"  ) {
	      double incidenceAngle = GetIncidenceAngle(aStep);
	      // incidence angles between [0,180]. eg) if critical angle is 82 degrees, retained light between [82,98]
	      // normalize s.t. distribution is between [0,180]
	      if (incidenceAngle > 90) incidenceAngle = 180 - incidenceAngle;
	      //	      if (is_downward_light) incidenceAngle = 180 - incidenceAngle;
	      m_incidenceAngleVec->push_back(incidenceAngle);	 
	      FillChannelTimeVector( channelNum, aStep->GetTrack()->GetGlobalTime());
	      m_yOriginVec->push_back(aStep->GetTrack()->GetVertexPosition().y() );
	      aStep->GetTrack()->SetTrackStatus( fStopAndKill );
	      return true;
	    }
	  }
	  m_cherenkovVec->at(channelNum)++;
	  m_nHits++;
	  m_lastFiveCore = 0;
	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
	  return true;
	}
	// Full Tree Mode
	else {
	  FiberHit* newHit = new FiberHit();
	  newHit->setPos      ( pos );
	  newHit->setOrigin   ( aStep->GetTrack()->GetVertexPosition() );
	  newHit->setMomentum ( aStep->GetPreStepPoint()->GetMomentum() );
	  newHit->setEnergy   ( aStep->GetPreStepPoint()->GetTotalEnergy() );
	  newHit->setTime     ( aStep->GetTrack()->GetGlobalTime() );
	  newHit->setRodNb    ( rodNum );

	  fiberCollection->insert ( newHit );
	  m_nHits++;

	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); //Kill the track so we only record it once
	  return true;
	}


	
      }
      m_prevTrackID = trackID;;
    }
  }
  // Optical OFF Mode
  else{ // Otherwise record all hits
    // don't record Cherenkovs in optical-off mode
    if (aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) return true;
    if(REDUCED_TREE){
      m_cherenkovVec->at(rodNum) += generatedPhotons;
      FillTimeVector( rodNum, aStep->GetTrack()->GetGlobalTime(), generatedPhotons );
      m_nHits++;
      return true;
    }
    else if(ML_REDUCED_TREE){
      G4int channelNum;
      if (RPD) channelNum = GetRPDChannelMapping(rodNum);
      else     channelNum = GetZDCChannelMapping(rodNum);
      if (RPD) FillTimeVector( rodNum, aStep->GetTrack()->GetGlobalTime(), generatedPhotons );
      m_cherenkovVec->at(channelNum) += generatedPhotons;
      m_nHits++;
      return true;
    }
    // Full Tree Mode
    else {
      FiberHit* newHit = new FiberHit();
      newHit->setCharge      ( aStep->GetPreStepPoint()->GetCharge() );
      newHit->setTrackID     ( aStep->GetTrack()->GetTrackID() );
      newHit->setModNb       ( m_modNum );
      newHit->setRodNb       ( rodNum );
      newHit->setEdep        ( aStep->GetTotalEnergyDeposit() );
      newHit->setOrigin      ( aStep->GetTrack()->GetVertexPosition() );
      newHit->setPos         ( pos );
      newHit->setParticle    ( particle );
      newHit->setEnergy      ( aStep->GetPreStepPoint()->GetTotalEnergy() );
      newHit->setMomentum    ( aStep->GetPreStepPoint()->GetMomentum() );
      newHit->setTime        ( aStep->GetTrack()->GetGlobalTime() );
      newHit->setNCherenkovs ( generatedPhotons );

      fiberCollection->insert ( newHit );
      m_nHits++;
      return true;
    }
  }
  return false; //Something failed if we got here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::EndOfEvent(G4HCofThisEvent*)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::SetReducedTree ( G4int _nFibers, G4int _nChannels ){
  REDUCED_TREE = true;
  m_nFibers = _nFibers;
  m_nChannels = _nChannels;
  m_nFiberPerChannel = _nFibers/_nChannels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::SetMLReducedTree ( G4int _nFibers, G4int _nChannels ){
  ML_REDUCED_TREE = true;
  m_nChannels = _nChannels;
  m_nFiberPerChannel = _nFibers/_nChannels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::FillTimeVector( G4int fiberNo, G4double time, G4int weight ){
  // This vector consists of nSegments number of "histograms" each 128 bins long, covering 64ns
  // concatenated in order of segment number (fiber number/m_nFiberPerChannel a.k.a. Channel)

  int bin = (time <= 64*ns) ? time/0.5 : 127;  // Determine the bin for this segments histo
  int offset = fiberNo/m_nFiberPerChannel;    // Do this because there was a rounding issue when just using an int cat
  bin += 128*offset;                           // Add the offset for the current segment (channel)
  m_timeVec->at(bin) += weight;                // Add the weight (nPhotons).
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FiberSD::FillChannelTimeVector( G4int channelNo, G4double time, G4int weight ){
  // This vector consists of nChannel number of "histograms" each 128 bins long, covering 64ns
  // concatenated in order of channel number (fiber number/m_nFiberPerChannel a.k.a. Channel)
  // Note: 128 is just an arbitrary buffer length (64 ns)

  int bin = (time <= 64*ns) ? time/0.5 : 127;  // Determine the bin for this segments histo
  //  std::cout << "time " << time << " time/0.5 " << time/0.5 << " bin " << bin << std::endl;
  int offset = channelNo;    // Do this because there was a rounding issue when just using an int cat
  bin += 128*offset;                           // Add the offset for the current segment (channel)
  m_timeVec->at(bin) += weight;                // Add the weight (nPhotons).
}

G4double FiberSD::GetIncidenceAngle(const G4Step *aStep)
{
  G4StepPoint *preStep = aStep->GetPreStepPoint();
  G4ThreeVector photonDirection = preStep->GetMomentum() / preStep->GetMomentum().mag();
  G4ThreeVector stepPos = preStep->GetPosition();

  const G4VTouchable *touchable = preStep->GetTouchable();

  const G4RotationMatrix *rotation = touchable->GetRotation();
  G4RotationMatrix rotation_inv = rotation->inverse();
  G4ThreeVector translation = touchable->GetTranslation();
  G4VSolid *sector = touchable->GetSolid();

  G4ThreeVector posLocal = *rotation * (stepPos - translation);
  G4ThreeVector normal =  - sector->SurfaceNormal(posLocal);

  G4ThreeVector photonDirectionLocal = *rotation * photonDirection;

  G4double incidenceAngle = acos( normal.dot(photonDirectionLocal) );
  G4double incidenceDegrees = (incidenceAngle / (2*M_PI)) * 360;
  
  return incidenceDegrees;
}
