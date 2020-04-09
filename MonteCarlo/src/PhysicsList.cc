#include "PhysicsList.hh"
#include "PhysicsMessenger.hh"

PhysicsList::PhysicsList( )
{

  m_messenger = new PhysicsMessenger( this );

  //Set Physics list via messenger
  //Hadronic_PL = _physicsList;

  decay = NULL;
  stepLimiter = NULL;

  // optical
  scintillation = NULL;
  cerenkov = NULL;
  absorption = NULL;
  rayleigh = NULL;
  boundary = NULL;
  opmiehg  = NULL;
  wlshift  = NULL;

//  gflash = NULL;

  verboseLevel = 1;
}

PhysicsList::~PhysicsList(void)
{
  if (decay != NULL)
    delete decay;

  if (stepLimiter != NULL)
    delete stepLimiter;

  for (unsigned int i = 0; i < userSpecialCuts.size(); i++)
    delete userSpecialCuts.at(i);
  userSpecialCuts.clear();

  for (unsigned int i = 0; i < g4VProcess.size(); i++)
    delete g4VProcess.at(i).first;
  g4VProcess.clear();

  for (unsigned int i = 0; i < hadronPhysics.size(); i++)
    delete hadronPhysics.at(i);
  hadronPhysics.clear();

  // optical
  if (scintillation != NULL)
    delete scintillation;
  if (cerenkov != NULL)
    delete cerenkov;
  if (absorption != NULL)
    delete absorption;
  if (rayleigh != NULL)
    delete rayleigh;
  if (boundary != NULL)
    delete boundary;
  if (opmiehg != NULL)
    delete opmiehg;
  if (wlshift != NULL)
    delete wlshift;


}

void PhysicsList::ConstructParticle(void)
{
  G4LeptonConstructor leptonConstructor;
  leptonConstructor.ConstructParticle();
  G4BosonConstructor bosonConstructor;
  bosonConstructor.ConstructParticle();
  G4MesonConstructor mesonConstructor;
  mesonConstructor.ConstructParticle();
  G4BaryonConstructor baryonConstructor;
  baryonConstructor.ConstructParticle();
  G4IonConstructor ionConstructor;
  ionConstructor.ConstructParticle();
  G4ShortLivedConstructor shortLivedConstructor;
  shortLivedConstructor.ConstructParticle();
}

void PhysicsList::ConstructProcess(void)
{
  AddTransportation();
  constructEM();
  constructDecay();
  addHadronic();
  for (unsigned int i = 0; i < hadronPhysics.size(); i++)
      hadronPhysics[i]->ConstructProcess();


// make a bool for optical option
if (1){
  constructOptical();
  std::cout << "Optical Physics Turned ON" << std::endl;
}
else std::cout << "Optical Physics Tunred OFF" << std::endl;

  stepLimiter = new G4StepLimiter();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    particle->GetProcessManager()->AddDiscreteProcess(stepLimiter);
    userSpecialCuts.push_back(new G4UserSpecialCuts());
    particle->GetProcessManager()->AddProcess(userSpecialCuts.back(), -1, -1,  1);
  }
}

void PhysicsList::constructEM(void)
{
  // from example N03
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      g4VProcess.push_back(std::make_pair(new G4PhotoElectricEffect, particle));
      g4VProcess.push_back(std::make_pair(new G4ComptonScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4GammaConversion, particle));

    } else if (particleName == "e-") {
      //electron
      g4VProcess.push_back(std::make_pair(new G4eMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4eIonisation, particle));
      g4VProcess.push_back(std::make_pair(new G4eBremsstrahlung, particle));

    } else if (particleName == "e+") {
      //positron
      g4VProcess.push_back(std::make_pair(new G4eMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4eIonisation, particle));
      g4VProcess.push_back(std::make_pair(new G4eBremsstrahlung, particle));
      g4VProcess.push_back(std::make_pair(new G4eplusAnnihilation, particle));

    } else if (particleName == "mu+" || particleName == "mu-") {
      //muon
      g4VProcess.push_back(
          std::make_pair(new G4MuMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4MuIonisation, particle));
      g4VProcess.push_back(std::make_pair(new G4MuBremsstrahlung, particle));
      g4VProcess.push_back(std::make_pair(new G4MuPairProduction, particle));

    } else if (particleName == "proton" || particleName == "pi-"
        || particleName == "pi+") {
      //proton
      g4VProcess.push_back(std::make_pair(new G4hMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4hIonisation, particle));
      g4VProcess.push_back(std::make_pair(new G4hBremsstrahlung, particle));
      g4VProcess.push_back(std::make_pair(new G4hPairProduction, particle));

    } else if (particleName == "alpha" || particleName == "He3") {
      //alpha
      g4VProcess.push_back(std::make_pair(new G4hMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4ionIonisation, particle));

    } else if (particleName == "GenericIon") {
      //Ions
      g4VProcess.push_back(std::make_pair(new G4hMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4ionIonisation, particle));

    } else if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0)
        && (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      g4VProcess.push_back(std::make_pair(new G4hMultipleScattering, particle));
      g4VProcess.push_back(std::make_pair(new G4hIonisation, particle));
    }
  }

  for (unsigned int i = 0; i < g4VProcess.size(); i++)
    ph->RegisterProcess(g4VProcess.at(i).first, g4VProcess.at(i).second);
}

void PhysicsList::constructDecay(void)
{
  // from example N03
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // Add Decay Process
  decay = new G4Decay();
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (decay->IsApplicable(*particle)) {
      ph->RegisterProcess(decay, particle);
    }
  }
}

void PhysicsList::addHadronic(void)
{
  if(Hadronic_PL == "QGSP_BERT")
     hadronPhysics.push_back(new G4HadronPhysicsQGSP_BERT(verboseLevel));
  else
     hadronPhysics.push_back(new G4HadronPhysicsFTFP_BERT(verboseLevel));

  G4NeutronTrackingCut* input = new G4NeutronTrackingCut(verboseLevel);
  input->SetKineticEnergyLimit(10.0 * CLHEP::MeV);
  input->SetTimeLimit(0.1 * CLHEP::ms);
  hadronPhysics.push_back(input);
  hadronPhysics.push_back(new G4StoppingPhysics(verboseLevel));
  hadronPhysics.push_back(new G4EmExtraPhysics(verboseLevel));
  hadronPhysics.push_back(new G4HadronElasticPhysics(verboseLevel));
}

void PhysicsList::addGflash(void)
{
  std::cout << " GFlash not implemented (yet)" << std::endl;
}

void PhysicsList::constructOptical(void)
{
  scintillation = new G4Scintillation();
  cerenkov = new G4Cerenkov();
  absorption = new G4OpAbsorption();
  rayleigh = new G4OpRayleigh();
  boundary = new G4OpBoundaryProcess();
  opmiehg = new G4OpMieHG();
  wlshift = new G4OpWLS();

    ////////////////////////////
    ////////////////////////////
    /*
G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();

 opticalPhysics->SetWLSTimeProfile("delta");

 opticalPhysics->SetScintillationYieldFactor(1.0);
 opticalPhysics->SetScintillationExcitationRatio(0.0);

 opticalPhysics->SetMaxNumPhotonsPerStep(100);
 opticalPhysics->SetMaxBetaChangePerStep(10.0);

 opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
 opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);

  scintillation->AddSaturation(G4LossTableManager::Instance()->EmSaturation());
  */
    ////////////////////////////
    ////////////////////////////


  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (scintillation->IsApplicable(*particle))
      this->RegisterProcess(scintillation, particle);
    if (cerenkov->IsApplicable(*particle))
      this->RegisterProcess(cerenkov, particle);



    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
      this->RegisterProcess(absorption, particle);
      this->RegisterProcess(rayleigh, particle);
      this->RegisterProcess(boundary, particle);
      this->RegisterProcess(opmiehg, particle);
      this->RegisterProcess(wlshift, particle);
    }
  }
}

void PhysicsList::SetCuts(void)
{

  std::cout << "Set Cuts to be defined, otherwise default is loaded." << std::endl;

  SetCutsWithDefault();
}

G4ParticleDefinition* PhysicsList::getParticleByID(G4int id)
{
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (particle->GetPDGEncoding() == id) {
      return particle;
    }
  }

  std::cout <<  "PhysicsList::getParticleByID: Unknown particle ID. Return G4Geantino." << std::endl;

  return G4Geantino::GeantinoDefinition();
}
