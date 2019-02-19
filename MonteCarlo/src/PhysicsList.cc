#include "PhysicsList.hh"
#include "SharedData.hh"

//#include "ExtraPhysics.hh"
#include "OpticalPhysics.hh"

#include <TEnv.h>

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//#include "G4PhysListFactory.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "StepMax.hh"

#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4IonPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4NeutronCrossSectionXS.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronHElasticPhysics.hh"     
#include "G4RadioactiveDecayPhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4Version.hh"
#include "G4LossTableManager.hh"

// *ARIC COMMENTED OUT*
//#if G4VERSION_NUMBER > 999
//static G4ParticleTable::G4PTblDicIterator* aParticleIterator;
//#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(G4String physName, SharedData *sd) : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();
    m_sd = sd;
    
    defaultCutValue  = 1.*mm;
    fCutForGamma     = defaultCutValue;
    fCutForElectron  = defaultCutValue;
    fCutForPositron  = defaultCutValue;

//    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = NULL;
    if (physName == "QGSP_BERT") {
       phys = new QGSP_BERT;
    } else {
       phys = new FTFP_BERT;
    }
//    if (factory.IsReferencePhysList(physName)) {
//       phys = factory.GetReferencePhysList(physName);
//       if(!phys)G4Exception("PhysicsList::PhysicsList","InvalidSetup",
//                            FatalException,"PhysicsList does not exist");
    //  fMessenger = new PhysicsListMessenger(this);
//    }

    for (G4int i = 0; ; ++i) {
       G4VPhysicsConstructor* elem =
                  const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
       if (elem == NULL) break;
       G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
       RegisterPhysics(elem);
    }

    fAbsorptionOn = false;

    TEnv* config = m_sd->GetConfig();
    if (config->GetValue("simCherenkov",false) == 1) SetAbsorption(fAbsorptionOn);
    //This looks complex, but it is not:
    //Get from base-class the pointer of the phsyicsVector
    //to be used. Remember: G4VModularPhysicsList is now a split class.
    //Why G4VModularPhysicsList::RegisterPhysics method is not used instead?
    //If possible we can remove this...
    fPhysicsVector =
                GetSubInstanceManager().offset[GetInstanceID()].physicsVector;
    
    //   fPhysicsVector->push_back(new ExtraPhysics());
/*
    fPhysicsVector->push_back(new G4RadioactiveDecayPhysics());
    fPhysicsVector->push_back(new G4HadronElasticPhysicsXS());
    fPhysicsVector->push_back(new G4StoppingPhysics());
    fPhysicsVector->push_back(new G4IonPhysics());
*/
    //    fPhysicsVector->push_back(new G4GenericIon());

    fStepMaxProcess = new StepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  //    delete fMessenger;

    delete fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ClearPhysics()
{
    for (G4PhysConstVector::iterator p  = fPhysicsVector->begin();
                                     p != fPhysicsVector->end(); ++p) {
        delete (*p);
    }
    fPhysicsVector->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
 
    G4VModularPhysicsList::ConstructParticle();

    G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
    MuonPlusDecayTable -> Insert(new
                           G4MuonDecayChannelWithSpin("mu+",0.986));
    MuonPlusDecayTable -> Insert(new
                           G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
    G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);

    G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
    MuonMinusDecayTable -> Insert(new
                            G4MuonDecayChannelWithSpin("mu-",0.986));
    MuonMinusDecayTable -> Insert(new
                            G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
    G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    SetVerbose(0);

    G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();

    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

    G4ProcessManager* pManager;
    pManager = G4MuonPlus::MuonPlus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      pManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

    pManager = G4MuonMinus::MuonMinus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      pManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    // Ions
    //    If ( pName == “GenericIon” || pName == “alpha” || pName ==“He3”) {
    /*      pManager->AddProcess (new G4MultipleScattering, -1, 1, 1 );
      pManager->AddProcess (new G4ionIonisation, -1, 2, 2 );
    */
      // Hadrons
      //  } else if (particle->GetPDGCharge() != 0 && particle->GetPDGMass() > 130.*MeV) {
      //  pmanager->AddProcess (new G4MultipleScattering, -1, 1, 1 );
      // pmanager->AddProcess (new G4hIonisation, -1, 2, 2 );
      // }

    G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

    decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

    pManager = G4PionPlus::PionPlus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(poldecay, idxPostStep);
      pManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

    pManager = G4PionMinus::PionMinus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(poldecay, idxPostStep);
      pManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    AddStepMax();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RemoveFromPhysicsList(const G4String& name)
{
    G4bool success = false;
    for (G4PhysConstVector::iterator p  = fPhysicsVector->begin();
                                     p != fPhysicsVector->end(); ++p) {
        G4VPhysicsConstructor* e = (*p);
        if (e->GetPhysicsName() == name) {
           fPhysicsVector->erase(p);
           success = true;
           break;
        }
    }
    if (!success) {
       G4ExceptionDescription message;
       message << "PhysicsList::RemoveFromEMPhysicsList "<< name << "not found";
       G4Exception("example PhysicsList::RemoveFromPhysicsList()",
       "ExamPhysicsList01",FatalException,message);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetAbsorption(G4bool toggle)
{
       fAbsorptionOn = toggle;

       RemoveFromPhysicsList("Optical");
       TEnv* config = m_sd->GetConfig();
       if (config->GetValue("simCherenkov",false) == 1) {
	 fPhysicsVector->push_back(fOpticalPhysics = new OpticalPhysics(toggle));       
	 fOpticalPhysics->ConstructProcess();
       }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
    if (verboseLevel >0) {
        G4cout << "PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length")
               << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    //SetCutValue(fCutForGamma, "gamma");
    SetCutValue(0.05*eV, "gamma");
    SetCutValue(fCutForElectron, "e-");
    SetCutValue(fCutForPositron, "e+");

    if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
    fCutForGamma = cut;
    SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
    fCutForElectron = cut;
    SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
    fCutForPositron = cut;
    SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetStepMax(G4double step)
{
  fStepMaxProcess->SetStepMax(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StepMax* PhysicsList::GetStepMaxProcess()
{
  return fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
#if G4VERSION_NUMBER >= 1030
  auto theParticleIterator1 = GetParticleIterator();
   
  theParticleIterator1->reset();
       
  while ((*theParticleIterator1)()) {
	
    G4ParticleDefinition* particle = theParticleIterator1->value();
# else
	
    theParticleIterator->reset();
  
    while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
#endif
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (fStepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
	{
	  if (pmanager) pmanager ->AddDiscreteProcess(fStepMaxProcess);
	}
    }
  }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetNbOfPhotonsCerenkov(G4int maxNumber)
{

    TEnv* config = m_sd->GetConfig();
    if (config->GetValue("simCherenkov",false) == 1)   fOpticalPhysics->SetNbOfPhotonsCerenkov(maxNumber);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVerbose(G4int verbose)
{
  TEnv* config = m_sd->GetConfig();
  if (config->GetValue("simCherenkov",false) == 1) {
    fOpticalPhysics->GetCerenkovProcess()->SetVerboseLevel(verbose);
    //   fOpticalPhysics->GetScintillationProcess()->SetVerboseLevel(verbose);
    fOpticalPhysics->GetAbsorptionProcess()->SetVerboseLevel(verbose);
    //   fOpticalPhysics->GetRayleighScatteringProcess()->SetVerboseLevel(verbose);
    //   fOpticalPhysics->GetMieHGScatteringProcess()->SetVerboseLevel(verbose);
    fOpticalPhysics->GetBoundaryProcess()->SetVerboseLevel(verbose);
  }

}
