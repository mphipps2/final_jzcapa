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
// January 22, 2012
// Yakov Kulinich
// Stony Brook & BNL
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef CherenkovHit_h
#define CherenkovHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class CherenkovHit : public G4VHit
{
public:

  CherenkovHit();
  ~CherenkovHit();
  CherenkovHit(const CherenkovHit&);
  const CherenkovHit& operator=(const CherenkovHit&);
  G4int operator==(const CherenkovHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

public:
  void setTrackID    (G4int track)        { trackID = track; };
  void setModNb      (G4int mod)          { modNb = mod; };  
  void setRadNb      (G4int rad)          { radNb = rad; };  
  void setEdep       (G4double de)        { edep = de; };
  void setPos        (G4ThreeVector xyz)  { pos = xyz; };
  void setMomentum   (G4ThreeVector mom)  { momentum = mom;};
  void setParticle   (G4ParticleDefinition *part) {particle = part;};
  void setEnergy     (G4double e)         {energy = e;};
  void setCharge     (G4double c)         {charge = c;};
  void setEventNo    (G4int eventNo)      {fEventNo = eventNo;}
  void setVelocity    (G4double vel)      {velocity = vel;}
  void setBeta    (G4double b)      {beta = b;}
  
  G4int         getTrackID()    { return trackID; };
  G4int         getModNb()      { return modNb; };
  G4int         getRadNb()      { return radNb; };
  G4ParticleDefinition*  getParticle()      { return particle; };
  G4double      getEdep()       { return edep; };
  G4double      getEnergy()       { return energy; };      
  G4ThreeVector getPos()        { return pos; };
  G4ThreeVector getMomentum()        { return momentum; };
  G4double      getCharge()     {return charge;};
  G4int         getEventNo()    {return fEventNo;}
  G4double      getVelocity()   {return velocity;}
  G4double      getBeta()   {return beta;}
private:
  G4int         trackID;
  G4int         modNb;
  G4int         radNb;
  G4double      edep;
  G4double      velocity;
  G4double      beta;
  G4ParticleDefinition *particle;
  G4ThreeVector pos;
  G4ThreeVector momentum;
  G4double      energy;
  G4double      charge;
  G4int         fEventNo;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<CherenkovHit> CherenkovHitsCollection;

extern G4Allocator<CherenkovHit> CherenkovHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CherenkovHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) CherenkovHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void CherenkovHit::operator delete(void *aHit)
{
  CherenkovHitAllocator.FreeSingle((CherenkovHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
