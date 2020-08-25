
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
  fVertXingAngle(0.),
  fHorizXingAngle(0.),
  fProjPlane(0.),
  fnPrimaries(0),
  PROJECT(false),
  fpos( new G4ThreeVector(0.,0.,0.) )
{
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
  // These values are hard coded for accuracy and consistency
  // Though they need to be updated with the correct values
  G4double sigmaThetaXZ = 0.;//3.57e-6; // 3.57e-6 Corresponds to a 1mm beam diameter
  G4double sigmaThetaYZ = 0.;//3.57e-6;
  G4double sigmaX = 0.*mm;
  G4double sigmaY = 0.*mm;
  G4double sigmaZ = 0.*mm;
  G4double sigmaE = 1.e-3;
  G4double energy = (2.5 + G4RandGauss::shoot(0.0,sigmaE) )*TeV;

  G4ParticleDefinition* particleDefinition=
      G4ParticleTable::GetParticleTable()->FindParticle("neutron");

  // Adjust the momentum for the crossing angle
  G4ThreeVector momentum(0.,0.,1.);
  momentum.rotateY(fHorizXingAngle + G4RandGauss::shoot(0.0,sigmaThetaXZ));
  momentum.rotateX(fVertXingAngle  + G4RandGauss::shoot(0.0,sigmaThetaYZ));

  // Project the beam to the plane requested by the user if requested
  // Otherwise carry forward the beam position
  G4ThreeVector* position;
  if(PROJECT){
    G4double projDist = fabs( fProjPlane - fpos->z() );
    G4double projectedX = fpos->x() + G4RandGauss::shoot(0.0,sigmaX) + momentum.x()*projDist;
    G4double projectedY = fpos->y() + G4RandGauss::shoot(0.0,sigmaY) + momentum.y()*projDist;
    position = new G4ThreeVector( projectedX, projectedY, fProjPlane + G4RandGauss::shoot(0.0,sigmaZ) );
  }else{
    position = new G4ThreeVector( fpos->x(), fpos->y(), fpos->z() );
  }

  //If nPrimaries is 0, generate a random number from the distribution (to be implemented)
  if( fnPrimaries == 0 ) fnPrimaries = 1;
  // int nNeutrons = some distribution dependent random number;
  for(int i = 0; i < fnPrimaries; i++){
    G4PrimaryParticle* particle = new G4PrimaryParticle(particleDefinition);
    particle->SetMomentumDirection( momentum );
    particle->SetKineticEnergy( energy );

    G4PrimaryVertex* vertex = new G4PrimaryVertex( *position, 0. );
    vertex->SetPrimary( particle );
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
