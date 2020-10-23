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
/// \file src/PhysicsMessenger.cc
//---------------------------------------------------------------
//
//  From G4UserPhysicsListMessenger.cc
// ------------------------------------------------------------
//	History
//        first version                   09 Jan. 1998 by H.Kurashige
//        add buildPhysicsTable command   13 Apr. 1999 by H.Kurashige
//        add setStoredInAscii command    12 Mar. 2001 by H.Kurashige
//        add dumpOrderingParam command    3 May. 2011 by H.Kurashige
//        add SelectList command             Jul. 2020 by C.Lantz
// ------------------------------------------------------------

#include "PhysicsMessenger.hh"
#include "PhysicsList.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"


PhysicsMessenger::PhysicsMessenger(PhysicsList * pList)
:G4UImessenger(),fPhysicsList(pList)
{
  G4UIparameter* param = 0;
  // /run/particle    directory
  fPhysicsDir = new G4UIdirectory("/run/particle/");
  fPhysicsDir->SetGuidance("PhysicsList options");

  fSelectListCmd = new G4UIcmdWithAString("/run/particle/SelectList", this);
  fSelectListCmd->SetGuidance("Choose a physics list to use");
  fSelectListCmd->SetParameterName("List Name",false);
  fSelectListCmd->SetDefaultValue("FTFP_BERT");

  // /run/particle/Verbose command
  fverboseCmd = new G4UIcmdWithAnInteger("/run/particle/verbose",this);
  fverboseCmd->SetGuidance("Set the Verbose level of G4VUserPhysicsList.");
  fverboseCmd->SetGuidance(" 0 : Silent (default)");
  fverboseCmd->SetGuidance(" 1 : Display warning messages");
  fverboseCmd->SetGuidance(" 2 : Display more");
  fverboseCmd->SetParameterName("level",true);
  fverboseCmd->SetDefaultValue(0);
  fverboseCmd->SetRange("level >=0 && level <=3");

  // /run/setCut command
  fsetCutCmd = new G4UIcmdWithADoubleAndUnit("/run/setCut",this);
  fsetCutCmd->SetGuidance("Set default cut value ");
  fsetCutCmd->SetParameterName("cut",false);
  fsetCutCmd->SetDefaultValue(1.0);
  fsetCutCmd->SetRange("cut >=0.0");
  fsetCutCmd->SetDefaultUnit("mm");
  fsetCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // /run/setCutForAGivenParticle command
  fsetCutForAGivenParticleCmd = new G4UIcommand("/run/setCutForAGivenParticle",this) ;
  fsetCutForAGivenParticleCmd->SetGuidance("Set a cut value to a specific particle ") ;
  fsetCutForAGivenParticleCmd->SetGuidance("Usage: /run/setCutForAGivenParticle  gamma  1. mm") ;
  param = new G4UIparameter("particleName",'s',false) ;
  param->SetParameterCandidates("e- e+ gamma proton");
  fsetCutForAGivenParticleCmd->SetParameter(param) ;
  param = new G4UIparameter("cut",'d',false) ;
  param->SetDefaultValue("1.") ;
  param->SetParameterRange("cut>=0.0") ;
  fsetCutForAGivenParticleCmd->SetParameter(param) ;
  param = new G4UIparameter("unit",'s',false) ;
  param->SetDefaultUnit("mm");
  fsetCutForAGivenParticleCmd->SetParameter(param) ;
  fsetCutForAGivenParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // /run/getCutForAGivenParticle command
  fgetCutForAGivenParticleCmd = new G4UIcmdWithAString("/run/getCutForAGivenParticle",this) ;
  fgetCutForAGivenParticleCmd->SetGuidance("Get a cut value to a specific particle ") ;
  fgetCutForAGivenParticleCmd->SetGuidance("Usage: /run/getCutForAGivenParticle  gamma ") ;
  fgetCutForAGivenParticleCmd->SetParameterName("particleName",false,false) ;
  fgetCutForAGivenParticleCmd->SetCandidates("e- e+ gamma proton");
  fgetCutForAGivenParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // /run/setCutForRegion command
  fsetCutRCmd = new G4UIcommand("/run/setCutForRegion",this);
  fsetCutRCmd->SetGuidance("Set cut value for a region");
  param = new G4UIparameter("Region",'s',false);
  fsetCutRCmd->SetParameter(param);
  param = new G4UIparameter("cut",'d',false);
  param->SetParameterRange("cut >=0.0");
  fsetCutRCmd->SetParameter(param);
  param = new G4UIparameter("Unit",'s',true);
  param->SetDefaultValue("mm");
  param->SetParameterCandidates(fsetCutRCmd->UnitsList(fsetCutRCmd->CategoryOf("mm")));
  fsetCutRCmd->SetParameter(param);
  fsetCutRCmd->AvailableForStates(G4State_Idle);

  // /run/particle/DumpList command
  fdumpListCmd = new G4UIcmdWithoutParameter("/run/particle/dumpList",this);
  fdumpListCmd->SetGuidance("Dump List of particles in G4VUserPhysicsList. ");

  // /run/particle/addProcManager command
  faddProcManCmd = new G4UIcmdWithAString("/run/particle/addProcManager", this);
  faddProcManCmd->SetToBeBroadcasted(false);
  faddProcManCmd->SetGuidance("add process manager to specified particle type");
  faddProcManCmd->SetParameterName("particleType", true);
  faddProcManCmd->SetDefaultValue("");
  faddProcManCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // /run/particle/buildPhysicsTable command
  fbuildPTCmd = new G4UIcmdWithAString("/run/particle/buildPhysicsTable", this);
  fbuildPTCmd->SetGuidance("build physics table of specified particle type");
  fbuildPTCmd->SetParameterName("particleType", true);
  fbuildPTCmd->SetDefaultValue("");
  fbuildPTCmd->AvailableForStates(G4State_Init,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  // /run/particle/storePhysicsTable command
  fstoreCmd = new G4UIcmdWithAString("/run/particle/storePhysicsTable",this);
  fstoreCmd->SetGuidance("Store Physics Table");
  fstoreCmd->SetGuidance("  Enter directory name");
  fstoreCmd->SetParameterName("dirName",true);
  fstoreCmd->SetDefaultValue("");
  fstoreCmd->AvailableForStates(G4State_Idle);

  //  /run/particle/retrievePhysicsTable command
  fretrieveCmd = new G4UIcmdWithAString("/run/particle/retrievePhysicsTable",this);
  fretrieveCmd->SetGuidance("Retrieve Physics Table");
  fretrieveCmd->SetGuidance("  Enter directory name or OFF to switch off");
  fretrieveCmd->SetParameterName("dirName",true);
  fretrieveCmd->SetDefaultValue("");
  fretrieveCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //  /run/particle/setStoredInAscii command
  fasciiCmd = new G4UIcmdWithAnInteger("/run/particle/setStoredInAscii",this);
  fasciiCmd->SetGuidance("Switch on/off ascii mode in store/retrieve Physics Table");
  fasciiCmd->SetGuidance("  Enter 0(binary) or 1(ascii)");
  fasciiCmd->SetParameterName("ascii",true);
  fasciiCmd->SetDefaultValue(0);
  fasciiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fasciiCmd->SetRange("ascii ==0 || ascii ==1");

  //Commnad    /run/particle/applyCuts command
  fapplyCutsCmd = new G4UIcommand("/run/particle/applyCuts",this);
  fapplyCutsCmd->SetGuidance("Set applyCuts flag for a particle.");
  fapplyCutsCmd->SetGuidance(" Some EM processes which do not have infrared divergence");
  fapplyCutsCmd->SetGuidance("may generate gamma, e- and/or e+ with kinetic energies");
  fapplyCutsCmd->SetGuidance("below the production threshold. By setting this flag,");
  fapplyCutsCmd->SetGuidance("such secondaries below threshold are eliminated and");
  fapplyCutsCmd->SetGuidance("kinetic energies of such secondaries are accumulated");
  fapplyCutsCmd->SetGuidance("to the energy deposition of their mother.");
  fapplyCutsCmd->SetGuidance(" Note that 'applyCuts' makes sense only for gamma,");
  fapplyCutsCmd->SetGuidance("e- and e+. If this command is issued for other particle,");
  fapplyCutsCmd->SetGuidance("a warning message is displayed and the command is");
  fapplyCutsCmd->SetGuidance("ignored.");
  fapplyCutsCmd->SetGuidance(" If particle name is 'all', this command affects on");
  fapplyCutsCmd->SetGuidance("gamma, e- and e+.");
  param = new G4UIparameter("Flag",'s',true);
  param->SetDefaultValue("true");
  fapplyCutsCmd->SetParameter(param);
  param = new G4UIparameter("Particle",'s',true);
  param->SetDefaultValue("all");
  fapplyCutsCmd->SetParameter(param);
  fapplyCutsCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);

  //  /run/particle/dumpCutValues command
  fdumpCutValuesCmd = new G4UIcmdWithAString("/run/particle/dumpCutValues",this);
  fdumpCutValuesCmd->SetGuidance("Dump a list of production threshold values in range and energy");
  fdumpCutValuesCmd->SetGuidance("for all registered material-cuts-couples.");
  fdumpCutValuesCmd->SetGuidance("Dumping a list takes place when you issue 'beamOn' and");
  fdumpCutValuesCmd->SetGuidance("actual conversion tables from range to energy are available.");
  fdumpCutValuesCmd->SetGuidance("If you want a list 'immediately', use '/run/dumpRegion' for threshold");
  fdumpCutValuesCmd->SetGuidance("list given in range only. Also, '/run/dumpCouples' gives you the");
  fdumpCutValuesCmd->SetGuidance("current list if you have already issued 'run/beamOn' at least once.");
  fdumpCutValuesCmd->SetParameterName("particle",true);
  fdumpCutValuesCmd->SetDefaultValue("all");
  fdumpCutValuesCmd->AvailableForStates(G4State_Idle);

  //  /run/particle/dumpCutValues command
  fdumpOrdParamCmd = new G4UIcmdWithAnInteger("/run/particle/dumpOrderingParam",this);
  fdumpOrdParamCmd->SetGuidance("Dump a list of ordering parameter ");
  fdumpOrdParamCmd->SetParameterName("subtype",true);
  fdumpOrdParamCmd->SetDefaultValue(-1);
  fdumpOrdParamCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);
}

/*
 *
 */
PhysicsMessenger::~PhysicsMessenger(){
  delete fSelectListCmd;
  delete fsetCutCmd;
  delete fsetCutRCmd;
  delete fsetCutForAGivenParticleCmd;
  delete fgetCutForAGivenParticleCmd;
  delete fverboseCmd;
  delete fdumpListCmd;
  delete faddProcManCmd;
  delete fbuildPTCmd;
  delete fstoreCmd;
  delete fretrieveCmd;
  delete fasciiCmd;
  delete fapplyCutsCmd;
  delete fdumpCutValuesCmd;
  delete fdumpOrdParamCmd;
}

/*
 *
 */
void PhysicsMessenger::SetNewValue(G4UIcommand* command,G4String newValue){
  G4ExceptionDescription ed;
  if (command == fSelectListCmd) {
    fPhysicsList->SetList(newValue);

  } else if( command==fsetCutCmd ){
    G4double newCut = fsetCutCmd->GetNewDoubleValue(newValue);
    fPhysicsList->SetDefaultCutValue(newCut);
    fPhysicsList->SetCuts();

  } else if( command==fsetCutForAGivenParticleCmd ){
    G4String particleName, unit ; G4double cut ;
    std::istringstream str (newValue) ;
    str >> particleName >> cut >> unit ;
    fPhysicsList->SetCutValue(cut*G4UIcommand::ValueOf(unit), particleName) ;

  } else if( command==fgetCutForAGivenParticleCmd ){
    G4cout << fPhysicsList->GetCutValue(newValue)/mm <<"[mm]" << G4endl ;

  } else if( command==fsetCutRCmd ){
    std::istringstream is(newValue);
    G4String regName;
    G4String uniName;
    G4double cVal = -1.0;
    is >> regName >> cVal >> uniName;
    if (is.fail()) {
      ed << "illegal arguments : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    fPhysicsList->SetCutsForRegion(cVal*(fsetCutRCmd->ValueOf(uniName)),regName);

  } else if( command==fverboseCmd ) {
    fPhysicsList->SetVerboseLevel(fverboseCmd->GetNewIntValue(newValue));

  } else if( command==fdumpListCmd ){
    fPhysicsList->DumpList();

  } else if( command==fdumpOrdParamCmd ){
    G4int stype = fdumpOrdParamCmd->GetNewIntValue(newValue);
    G4PhysicsListHelper::GetPhysicsListHelper()->DumpOrdingParameterTable(stype);

  }  else if( command == faddProcManCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == 0)
    {
      ed << " Particle is not found : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    else if (particle->GetProcessManager() != 0)
    {
      ed << " Particle is not initialized : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    fPhysicsList->AddProcessManager(particle);

  }  else if( command == fbuildPTCmd ){
    G4ParticleDefinition* particle = (G4ParticleTable::GetParticleTable())->FindParticle(newValue);
    if (particle == 0)
    {
      ed << " Particle is not found : " << newValue;
      command->CommandFailed(ed);
      return;
    }
    fPhysicsList->PreparePhysicsTable(particle);
    fPhysicsList->BuildPhysicsTable(particle);

  } else if ( command == fstoreCmd ){
    fPhysicsList->StorePhysicsTable(newValue);

  } else if( command == fretrieveCmd ) {
    if ((newValue == "OFF") || (newValue == "off") ){
      fPhysicsList->ResetPhysicsTableRetrieved();
    } else {
      fPhysicsList->SetPhysicsTableRetrieved(newValue);
    }

  } else if( command == fasciiCmd ) {
    if (fasciiCmd->GetNewIntValue(newValue) == 0) {
      fPhysicsList->ResetStoredInAscii();
    } else {
      fPhysicsList->SetStoredInAscii();
    }

  } else if( command == fapplyCutsCmd ) {
    G4Tokenizer next( newValue );

    // check 1st argument
    G4String temp = G4String(next());
    G4bool flag = (temp =="true" || temp=="TRUE");

    // check 2nd argument
    G4String name = G4String(next());

    fPhysicsList->SetApplyCuts(flag, name);

  } else if( command == fdumpCutValuesCmd ) {
    fPhysicsList->DumpCutValuesTable(1);

  }
}

G4String PhysicsMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  G4String candidates("none");
  G4ParticleTable::G4PTblDicIterator *piter = (G4ParticleTable::GetParticleTable())->GetIterator();

  if( command==fsetCutCmd ) {
    cv = fsetCutCmd->ConvertToString( fPhysicsList->GetDefaultCutValue(), "mm" );

  } else if( command==fverboseCmd ){
    cv = fverboseCmd->ConvertToString(fPhysicsList->GetVerboseLevel());

  }  else if( command== faddProcManCmd ){
    // set candidate list
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    faddProcManCmd->SetCandidates(candidates);
    cv = "";

  }  else if( command== fbuildPTCmd ){
    // set candidate list
    piter -> reset();
    while( (*piter)() ){
      G4ParticleDefinition *particle = piter->value();
      candidates += " " + particle->GetParticleName();
    }
    faddProcManCmd->SetCandidates(candidates);
    cv = "";

  } else if ( command == fstoreCmd ){
    cv = fPhysicsList->GetPhysicsTableDirectory();

  }else if( command == fretrieveCmd ) {
    if (fPhysicsList->IsPhysicsTableRetrieved()) {
      cv = fPhysicsList->GetPhysicsTableDirectory();
    } else {
      cv = "OFF";
    }

  } else if( command==fasciiCmd ){
    if (fPhysicsList->IsStoredInAscii()){
      cv = "1";
    } else {
      cv = "0";
    }

//  } else if( command == fapplyCutsCmd ) {
//   if (fPhysicsList->GetApplyCuts("gamma")){
//     cv =  "true";
//   } else {
//     cv =  "false";
//   }
  }

  return cv;
}
