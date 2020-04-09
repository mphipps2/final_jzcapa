#ifndef PhysicsList_HH_
#define PhysicsList_HH_

#include "PhysicsMessenger.hh"

#include "G4VUserPhysicsList.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"

#include "G4Decay.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ionIonisation.hh"

// optical
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LossTableManager.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpMieHG.hh"
#include "G4OpWLS.hh"

// particles
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// gflash
#include "G4FastSimulationManagerProcess.hh"
// dy
#include "G4PionMinusInelasticProcess.hh"

#include <assert.h>


// This macro will be missing in Geant4.10.3 (for multi-threaded mode, which is not used in TGEANT), so we define it here
#define theParticleIterator ((this->subInstanceManager.offset[this->g4vuplInstanceID])._theParticleIterator)


class PhysicsList : public G4VUserPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList(void);
    void ConstructParticle(void);
    void ConstructProcess(void);

    void SetCuts(void);

    void SetList( G4String arg ){ Hadronic_PL = arg; }

    G4ParticleDefinition* getParticleByID(G4int id);

  private:
    PhysicsMessenger *m_messenger;
    G4String Hadronic_PL;
    int verboseLevel;

    void constructEM(void);
    void constructDecay(void);
    void addHadronic(void);
    void addGflash(void);
    void constructOptical(void);

    G4Decay* decay;
    std::vector<G4VPhysicsConstructor*>  hadronPhysics;

    // optical
    G4Scintillation* scintillation;
    G4Cerenkov* cerenkov;
    G4OpAbsorption* absorption;
    G4OpRayleigh* rayleigh;
    G4OpBoundaryProcess* boundary;
    G4OpMieHG* opmiehg;
    G4OpWLS* wlshift;

    //G4FastSimulationManagerProcess* gflash;

    G4StepLimiter* stepLimiter;
    std::vector<G4UserSpecialCuts*> userSpecialCuts;
    std::vector<std::pair<G4VProcess*, G4ParticleDefinition*> > g4VProcess;

};



#endif /* PhysicsList_HH_ */
