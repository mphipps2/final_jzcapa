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
      m_energyVec = analysisManager->GetEnergyVector(m_modNum);
      m_energyVec->clear();
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


  //  std::cout << " processing track " << aStep->GetTrack()->GetTrackID() << std::endl;
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
  G4bool KEEP_DOWNWARD_LIGHT = true;
  G4bool KEEP_UPWARD_LIGHT = true;
  int trackID = aStep->GetTrack()->GetTrackID();
  if (trackID != m_prevTrackID) {
    is_downward_light = 0;
    is_upward_light = 0;
    m_isTIR = 0;
  }
  // 380
  //  std::cout << " TOPOFVOLUME " << m_topOfVolume << std::endl;
  if( aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      
    //      int trackID = aStep->GetTrack()->GetTrackID();
    G4ThreeVector p = aStep->GetTrack()->GetMomentumDirection();
    double pT = sqrt( pow( p.x(), 2.0 ) + pow( p.z(), 2.0 ));
    double theta = M_PI / 2.0 - atan( pT / p.y() );
    
    //    std::cout << " track " << trackID << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " rodNum " << rodNum << " prestep material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " poststep material " << aStep->GetPostStepPoint()->GetMaterial()->GetName() << " stepNum " << aStep->GetTrack()->GetCurrentStepNumber()  << " m_TIR_count " << m_TIR_count  << std::endl;
    //      if (RPD) std::cout << " trackID "  << trackID << " theta " << theta << " y " << aStep->GetTrack()->GetPosition().y() << " z " << aStep->GetTrack()->GetPosition().z() << " rodNum " << rodNum << " material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " m_TIR_count " << m_TIR_count  << std::endl;
    

    if (!is_downward_light && !is_upward_light) {
      if (theta > M_PI/2) {
	is_downward_light = 1;
	// kill downward going photons
	if (!KEEP_DOWNWARD_LIGHT) {
	  //	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	  //	  return true;
	}
      }
      else {
	is_upward_light = 1;
	// kill upward going photons
	if (!KEEP_UPWARD_LIGHT) {
	  //	  aStep->GetTrack()->SetTrackStatus( fStopAndKill ); 
	  //	  return true;
	  //	}
	}
      }
    } 

    double incidencePreStepAngle = GetIncidenceAngle(aStep);
    double incidencePreStepCorrectedAngle = incidencePreStepAngle;
    double incidencePostStepAngle = GetPostStepIncidenceAngle(aStep);
    double incidencePostStepCorrectedAngle = incidencePostStepAngle;
    double incidenceTrackAngle = GetTrackIncidenceAngle(aStep);

    // incidence angles between [0,180]. eg) if critical angle is 82 degrees, retained light between [82,98]
    // normalize s.t. distribution is between [0,180]

    //      if (incidenceAngle > 90) incidenceAngle = 180 - incidenceAngle;
    if (incidencePreStepAngle > 90) incidencePreStepCorrectedAngle = 180 - incidencePreStepAngle;
    if (incidencePostStepAngle > 90) incidencePostStepCorrectedAngle = 180 - incidencePostStepAngle;
    if (incidenceTrackAngle > 90) incidenceTrackAngle = 180 - incidenceTrackAngle;
    double incidencePreStepCorrectedAngleOriginal = incidencePreStepCorrectedAngle;
    double incidencePostStepCorrectedAngleOriginal = incidencePostStepCorrectedAngle;
    /*
    // then flip the convention for downward going light so we can visualize the directional dependency
    if (is_downward_light) {
      incidencePreStepCorrectedAngle = 180 - incidencePreStepCorrectedAngle;
      incidencePostStepCorrectedAngle = 180 - incidencePostStepCorrectedAngle;
    }
    */
    if (RPD) {
      
      
            std::cout << "FiberSD trackID "  << trackID << " topOfVolume " << m_topOfVolume << " isTIR " << m_isTIR << " incidence preStep angle: " << incidencePreStepCorrectedAngle << " incidence postStep angle: " << incidencePostStepCorrectedAngle << " incidence track angle: " << incidenceTrackAngle << " x " << aStep->GetTrack()->GetPosition().x() << " y " << aStep->GetTrack()->GetPosition().y() << " prestep material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() <<  " poststep material " << aStep->GetPostStepPoint()->GetMaterial()->GetName()  << " process " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()  << " stepStatus " << aStep->GetPostStepPoint()->GetStepStatus() << " stepNum "<< aStep->GetTrack()->GetCurrentStepNumber() << " originx " << aStep->GetTrack()->GetVertexPosition().x() << " originy " << aStep->GetTrack()->GetVertexPosition().y() << " originVolume: " << aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName() << std::endl;
      
      //      /*      if (is_downward_light)*/ std::cout << "FiberSD trackID "  << trackID << " incidence track angle: " << incidenceTrackAngle << " x " << aStep->GetTrack()->GetPosition().x() << " y " << aStep->GetTrack()->GetPosition().y() << " prestep material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() <<  " prestep index " << aStep->GetPreStepPoint()->GetMaterial()->GetIndex() << " poststep material " << aStep->GetPostStepPoint()->GetMaterial()->GetName()  << " poststep index " << aStep->GetPostStepPoint()->GetMaterial()->GetIndex() << " stepNum "<< aStep->GetTrack()->GetCurrentStepNumber()  << " stepLength " << aStep->GetTrack()->GetStepLength() << " Time " << aStep->GetTrack()->GetGlobalTime() << std::endl;
    }

    //if (is_downward_light) incidenceAngle = 180 - incidenceAngle;
    //      if (is_downward_light) std::cout << " DOWNWARD LIGHT! " << std::endl;
    //      else std::cout << " UPWARD LIGHT " << std::endl;      
    // m_topOfVolume = 379.2 ie >= 379.1      
    if (pos.y() >= m_topOfVolume - 0.2*mm){


	
    /*
      if (RPD) { 
      if (is_downward_light) std::cout << " DOWNWARD LIGHT MADE IT TO THE TOP! " << std::endl;
      else std::cout << " upward light " << std::endl;
      }
    */

	  //    if (aStep->GetPreStepPoint()->GetMaterial()->GetName() == "SilicaCore_UI"  &&   aStep->GetPostStepPoint()->GetMaterial()->GetName() == "SilicaClad_UI"  ) {
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
	  //	    if (aStep->GetPreStepPoint()->GetMaterial()->GetName() == "SilicaClad_UI"  &&   aStep->GetPostStepPoint()->GetMaterial()->GetName() == "SilicaCore_UI"  ) {


	  //	  double incidenceAngle = GetIncidenceAngle(aStep);
	  // incidence angles between [0,180]. eg) if critical angle is 82 degrees, retained light between [82,98]
	  // normalize s.t. distribution is between [0,180]
	  //	  if (incidenceAngle > 90) incidenceAngle = 180 - incidenceAngle;
	  //	      if (is_downward_light) incidenceAngle = 180 - incidenceAngle;
	  //	  m_incidenceAngleVec->push_back(incidenceAngle);	 
	  m_incidenceAngleVec->push_back(incidencePostStepCorrectedAngle);	 
	  FillChannelTimeVector( channelNum, aStep->GetTrack()->GetGlobalTime());
	  m_yOriginVec->push_back(aStep->GetTrack()->GetVertexPosition().y() );
	  m_energyVec->push_back(aStep->GetTrack()->GetTotalEnergy() * 1e+6 );
	  //	  if (is_downward_light) std::cout << " RECORDING DOWNWARD LIGHT! " << std::endl;
	  //	  else if (is_upward_light) std::cout << " RECORDING UPWARD LIGHT! " << std::endl;
	  //	  std::cout << " RECORDING TRACK!!!! " << " topOfVolume " << m_topOfVolume << " trackID "  << trackID << " isTIR " << m_isTIR <<  " x " << aStep->GetTrack()->GetPosition().x() << " y " << aStep->GetTrack()->GetPosition().y() << " prestep material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() <<  " poststep material " << aStep->GetPostStepPoint()->GetMaterial()->GetName()  << " m_TIR_count " << m_TIR_count  << " stepNum "<< aStep->GetTrack()->GetCurrentStepNumber() << " originx " << aStep->GetTrack()->GetVertexPosition().x() << " originy " << aStep->GetTrack()->GetVertexPosition().y() << " originVolume: " << aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName() << std::endl;
	  //	if (is_downward_light)   std::cout << " RECORDING DOWNWARD TRACK!!!! " << " prestep incidence angle: " << incidencePreStepCorrectedAngleOriginal << " poststep incidence angle " << incidencePostStepCorrectedAngleOriginal << " prestep material " << aStep->GetPreStepPoint()->GetMaterial()->GetName() <<  " poststep material " << aStep->GetPostStepPoint()->GetMaterial()->GetName()  << " m_TIR_count " << m_TIR_count  <<  std::endl;
	  //	    std::cout << " light reached readout; isTIR " << m_isTIR << std::endl;
	  //	    aStep->GetTrack()->SetTrackStatus( fStopAndKill );
	  //	    return true;
	  //	    }
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
	//	std::cout << " recorded Track totalEnergy " << aStep->GetTrack()->GetTotalEnergy() * 1e+6 << std::endl;
	newHit->setEnergy   ( aStep->GetTrack()->GetTotalEnergy() * 1e+6);
	newHit->setEnergy   ( aStep->GetTrack()->GetTotalEnergy());
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

G4double FiberSD::GetPostStepIncidenceAngle(const G4Step *aStep)
{
  G4StepPoint *postStep = aStep->GetPostStepPoint();
  G4ThreeVector photonDirection = postStep->GetMomentum() / postStep->GetMomentum().mag();
  G4ThreeVector stepPos = postStep->GetPosition();

  const G4VTouchable *touchable = postStep->GetTouchable();

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

G4double FiberSD::GetTrackIncidenceAngle(const G4Step *aStep)
{
  G4Track *track = aStep->GetTrack();
  //  G4ThreeVector photonDirection = postStep->GetMomentum() / postStep->GetMomentum().mag();
  G4ThreeVector photonDirection = track->GetMomentumDirection();
  G4ThreeVector stepPos = track->GetPosition();

  const G4VTouchable *touchable = track->GetTouchable();

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

