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
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  m_particleGun(0), 
  m_world(0)
{
  G4int n_particle = 1;
  m_particleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic
  /*
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
  */
  // G4IonTable* IonTable = G4IonTable::GetIonTable();
  // G4String IonName;
  // G4ParticleDefinition* particle
    //    = particleTable->FindParticle(particleName="proton");
  //   = ionTable->FindIon(82,207,0.);
  //      = particleTable->FindParticle(particleName="opticalphoton");
  // m_particleGun->SetParticleDefinition(particle);
  //  m_particleGun->SetParticleMomentumDirection(G4ThreeVector(.5,.5,.5));
      m_particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1.));

      //m_particleGun->SetParticleEnergy(3.* eV);
     m_particleGun->SetParticleEnergy(150.* GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete m_particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Module volume
  // from G4LogicalVolumeStore.
  
  //  G4double worldSizeX = 0;
  //  G4double worldSizeY = 0;
  G4double worldSizeZ = 0;

  if ( !m_world )
  {
    G4LogicalVolume* moduleLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    if ( moduleLV ) m_world = dynamic_cast<G4Box*>(moduleLV->GetSolid());
  }
  if ( m_world ) {
    //    worldSizeX = m_world->GetXHalfLength()*2.;
    //    worldSizeY = m_world->GetYHalfLength()*2.;
    worldSizeZ = m_world->GetZHalfLength()*2.;
  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Module volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }

  //  G4double size = 0.25; 
  //  G4double x0 = size * worldSizeX * (G4UniformRand()-0.5) * mm;
  //  G4double y0 = size * worldSizeY * (G4UniformRand()-0.5) * mm;
  G4double z0 = 0.5 * worldSizeZ * mm;
  //  x0 = 0.5;
  //  y0 = 0.5;
  std::cout << " particle gun position x 0. y 0. z " << z0 << std::endl;
  //  m_particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  m_particleGun->SetParticlePosition(G4ThreeVector(0.,0.,z0));

  m_particleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

