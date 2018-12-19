#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "OpticalPhysics.hh"
#include "G4Version.hh"
#include "G4LossTableManager.hh"

#if G4VERSION_NUMBER > 999
static G4ParticleTable::G4PTblDicIterator* aParticleIterator;
#endif


OpticalPhysics::OpticalPhysics(G4bool toggle)
    : G4VPhysicsConstructor("Optical")
{
  fWLSProcess                = NULL;
  fScintProcess              = NULL;
  fCerenkovProcess           = NULL;
  fBoundaryProcess           = NULL;
  fAbsorptionProcess         = NULL;
  fRayleighScattering        = NULL;
  fMieHGScatteringProcess    = NULL;

  fAbsorptionOn              = toggle;
}

OpticalPhysics::~OpticalPhysics() { }

#include "G4OpticalPhoton.hh"

void OpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}

#include "G4ProcessManager.hh"

void OpticalPhysics::ConstructProcess()
{
    G4cout << "OpticalPhysics:: Add Optical Physics Processes"
           << G4endl;

    /*      
  fWLSProcess = new G4OpWLS();

  fScintProcess = new G4Scintillation();
  fScintProcess->SetScintillationYieldFactor(1.);
  fScintProcess->SetTrackSecondariesFirst(true);
    */
  fCerenkovProcess = new G4Cerenkov();
  fCerenkovProcess->SetMaxNumPhotonsPerStep(300);
  fCerenkovProcess->SetTrackSecondariesFirst(true);

  // Note: this G4OpAbsorption could be turned back on if attenuation length in quartz wants to be studied
   fAbsorptionProcess      = new G4OpAbsorption();
  //  fRayleighScattering     = new G4OpRayleigh();
  //  fMieHGScatteringProcess = new G4OpMieHG();
  fBoundaryProcess        = new G4OpBoundaryProcess();

  G4ProcessManager* pManager =
                G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  if (!pManager) {
     std::ostringstream o;
     o << "Optical Photon without a Process Manager";
     G4Exception("OpticalPhysics::ConstructProcess()","",
                  FatalException,o.str().c_str());
  }

  if (fAbsorptionOn) pManager->AddDiscreteProcess(fAbsorptionProcess);

  //pManager->AddDiscreteProcess(fRayleighScattering);
  //pManager->AddDiscreteProcess(fMieHGScatteringProcess);

  pManager->AddDiscreteProcess(fBoundaryProcess);

  //  fWLSProcess->UseTimeProfile("delta");
  //fWLSProcess->UseTimeProfile("exponential");

  //  pManager->AddDiscreteProcess(fWLSProcess);

  //  fScintProcess->SetScintillationYieldFactor(1.);
  //  fScintProcess->SetScintillationExcitationRatio(0.0);
  //  fScintProcess->SetTrackSecondariesFirst(true);

  // Use Birks Correction in the Scintillation process

  //  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  //  fScintProcess->AddSaturation(emSaturation);

  //aParticleIterator->reset();
#if G4VERSION_NUMBER > 999
      aParticleIterator->reset();
      while ((*aParticleIterator)()) {
        G4ParticleDefinition* particle = aParticleIterator->value();
# else
  theParticleIterator->reset();
       while ((*theParticleIterator)()) {
         G4ParticleDefinition* particle = theParticleIterator->value();
#endif

    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
       std::ostringstream o;
       o << "Particle " << particleName << "without a Process Manager";
       G4Exception("OpticalPhysics::ConstructProcess()","",
                    FatalException,o.str().c_str());
    }

    if(fCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(fCerenkovProcess);
      pManager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    }
    /*
    if(fScintProcess->IsApplicable(*particle)){
      pManager->AddProcess(fScintProcess);
      pManager->SetProcessOrderingToLast(fScintProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(fScintProcess,idxPostStep);
    }
    */
  }
}

void OpticalPhysics::SetNbOfPhotonsCerenkov(G4int maxNumber)
{
  fCerenkovProcess->SetMaxNumPhotonsPerStep(maxNumber);
}
