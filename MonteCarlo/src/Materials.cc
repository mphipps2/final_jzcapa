#include "Materials.hh"
#include "G4SystemOfUnits.hh"

Materials* Materials::materials = NULL;

Materials* Materials::getInstance(void)
{
    if (materials == NULL) {
        materials = new Materials();
    }
    return materials;
}

/* Function that constructs all the materials - to be then retrieved via pointer */
Materials::Materials(void){

  H  = new G4Element("Hydrogen", "H", 1., 1.01 * g/mole);
  N  = new G4Element("Nitrogen", "N", 7., 14.01 * g/mole);
  O  = new G4Element("Oxygen"  , "O", 8., 16.00 * g/mole);
  Si = new G4Element("Silicon" , "Si", 14., 28.08 * g/mole);
  C  = new G4Element("Carbon" , "C", 6., 12.00 * g/mole);

  //Materials definition
  Xenon = new G4Material("XenonGas", 54., 131.29 * g/mole, 5.458 * mg/cm3,
                                     kStateGas, 93.15 * kelvin, 1 * atmosphere);

  pQuartz = new G4Material("Quartz", 2.200 * g/cm3, 2);
  pQuartz->AddElement(Si,1);
  pQuartz->AddElement(O,2);

  EM_Quartz = new G4Material("EMQuartz", 2.200 * g/cm3, 2);
  EM_Quartz->AddElement(Si,1);
  EM_Quartz->AddElement(O,2);

  nist_manager = G4NistManager::Instance();
  nist_manager->SetVerbose(1);

  Air = nist_manager->FindOrBuildMaterial("G4_AIR");
  Al  = nist_manager->FindOrBuildMaterial("G4_Al");   	//G4_Al
  Cu  = nist_manager->FindOrBuildMaterial("G4_Cu");
  Pb  = nist_manager->FindOrBuildMaterial("G4_Pb");
  W   = nist_manager->FindOrBuildMaterial("G4_W");

  Polyethylene = new G4Material("Polyethylene", 1200 * kg/m3, 2);
  Polyethylene->AddElement(H,2);
  Polyethylene->AddElement(C,4);

  PMMA =  new G4Material("Polymethylmethacrylate", 1.18 * g/cm3, 3);
  PMMA->AddElement(C,5);
  PMMA->AddElement(H,8);
  PMMA->AddElement(O,2);

  Grease =  new G4Material("Grease", 1.0 * g/cm3, 3);
  Grease->AddElement(C,1);
  Grease->AddElement(H,1);
  Grease->AddElement(O,1);

}

void Materials::DefineOpticalProperties(void){

  float grease_RI = 1.46;
  float clad_RI = 1.49;
  float core_RI = 1.6;
  float tile_RI = 1.4585;

  //quartz optical properties
  G4double quartz_RIND[nEntriesWLS], quartz_ABSL[nEntriesWLS],quartz_RFLT[nEntriesWLS],quartz_EFIC[nEntriesWLS], PhotonEnergy[nEntriesWLS];
     for(int i = 0; i < nEntriesWLS; i++){
         PhotonEnergy[i] = 2.00*eV + i*0.03*eV;
         quartz_RIND[i] =tile_RI; //Refractive Index - constants
         quartz_ABSL[i] = 300.00*cm + i*20*cm; //Attenuation length
         quartz_RFLT[i] = 0.5;
         quartz_EFIC[i] = 0.5;

     }

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,quartz_RIND,nEntriesWLS);//index of refraction
   MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,quartz_ABSL,nEntriesWLS);//absorption length
  // MPT_Array.back()->AddProperty("REFLECTIVITY",PhotonEnergy,quartz_RFLT,nEntriesWLS);//refelectivity
  // MPT_Array.back()->AddProperty("EFFICIENCY",PhotonEnergy,quartz_EFIC,nEntriesWLS);//efficiency
   pQuartz->SetMaterialPropertiesTable(MPT_Array.back());




   //Aluminium optical properties
   G4double AllPhotonEnergies[nEntriesWLS] = {2.00*eV, 2.04*eV, 2.07*eV, 2.11*eV, 2.15*eV, 2.18*eV, 2.22*eV, 2.26*eV, 2.29*eV,
                                              2.33*eV, 2.37*eV, 2.40*eV, 2.44*eV, 2.48*eV, 2.51*eV, 2.55*eV, 2.59*eV, 2.62*eV,
                                              2.66*eV, 2.70*eV, 2.74*eV, 2.77*eV, 2.81*eV, 2.85*eV, 2.88*eV, 2.92*eV, 2.96*eV,
                                              2.99*eV, 3.03*eV, 3.07*eV, 3.10*eV, 3.14*eV, 3.18*eV, 3.21*eV, 3.25*eV, 3.29*eV,
                                              3.32*eV, 3.36*eV, 3.40*eV, 3.43*eV, 3.47*eV, 5.0*eV,  25.0*eV,100.0*eV, 1000.00*eV,
                                              10000.0*eV, 25000.0*eV, 50000.0*eV, 250000.*eV, 1000000.*eV};
   G4double Al_refl[nEntriesWLS] = {0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89, 0.89,
                                    0.89, 0.89, 0.89, 0.89, 0.79, 0.69, 0.49, 0.29, 0.19, 0.09};	//"Hard cut"
   G4double Al_ABSL[nEntriesWLS] = {0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm, 0.01*mm, 0.01*mm,0.01*mm,
                              0.01*mm, 0.10*mm, 1.00*mm,0.01*m,  0.01*m,  0.10*m, 0.10*m,  1.00*m,  10.0*m, 10.0*m};

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("REFLECTIVITY",AllPhotonEnergies,Al_refl,nEntriesWLS);
   MPT_Array.back()->AddProperty("ABSLENGTH",AllPhotonEnergies,Al_ABSL,nEntriesWLS);
   Al->SetMaterialPropertiesTable(MPT_Array.back());

   //Air optical properties
   G4double RefractiveIndexAir[nEntriesWLS];
   for (int i = 0; i < nEntriesWLS; i++) RefractiveIndexAir[i] = 1.0;//

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexAir, nEntriesWLS);
   Air->SetMaterialPropertiesTable(MPT_Array.back());

   //Polyethilene optical properties

   G4double RefractiveIndexClad1[nEntriesWLS];
   G4double AbsClad[nEntriesWLS];
   for (int i = 0; i < nEntriesWLS; i++){
     RefractiveIndexClad1[i] = clad_RI;//
     AbsClad[i] = 20.0*m;
   }

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexClad1,nEntriesWLS);
   MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntriesWLS);
   Polyethylene->SetMaterialPropertiesTable(MPT_Array.back());

   //PMMA optical properties
   G4double RefractiveIndexWLSfiber[nEntriesWLS];
   for (int i = 0; i < nEntriesWLS; i++) RefractiveIndexWLSfiber[i] = core_RI;
   G4double AbsWLSfiber[nEntriesWLS] = { 5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
                                         5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
                                         5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
                                         1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
                                         1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};
   G4double EmissionWLSfiber[nEntriesWLS] = {0.05, 0.10, 0.30, 0.50, 0.75, 1.00, 1.50, 1.85, 2.30, 2.75,
                                             3.25, 3.80, 4.50, 5.20, 6.00, 7.00, 8.50, 9.50, 11.1, 12.4,
                                             12.9, 13.0, 12.8, 12.3, 11.1, 11.0, 12.0, 11.0, 17.0, 16.9,
                                             15.0, 9.00, 2.50, 1.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00,
                                             0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("WLSABSLENGTH", PhotonEnergy, AbsWLSfiber, nEntriesWLS);
   MPT_Array.back()->AddProperty("WLSCOMPONENT", PhotonEnergy, EmissionWLSfiber, nEntriesWLS);
   MPT_Array.back()->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexWLSfiber, nEntriesWLS);
   MPT_Array.back()->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
   PMMA->SetMaterialPropertiesTable(MPT_Array.back());

   //Grease (silicone) optical properties
   G4double RefractiveIndexGrease[nEntriesWLS];
   for (int i = 0; i < nEntriesWLS; i++) RefractiveIndexGrease[i] = grease_RI;

   MPT_Array.push_back(new G4MaterialPropertiesTable());
   MPT_Array.back()->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexGrease,nEntriesWLS);
   MPT_Array.back()->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntriesWLS);
   Grease->SetMaterialPropertiesTable(MPT_Array.back());








   // set the optical boundary properties

  AlSurface = new G4OpticalSurface("AlSurface",unified, polished, dielectric_metal, 1);

   G4MaterialPropertiesTable* AlSurfaceProperty = new G4MaterialPropertiesTable();
   //AlSurfaceProperty->AddProperty("RINDEX",AllPhotonEnergies,RefractiveIndexGrease,nEntriesWLS);
   AlSurfaceProperty->AddProperty("REFLECTIVITY",AllPhotonEnergies,Al_refl,nEntriesWLS);
   AlSurfaceProperty->AddProperty("EFFICIENCY",AllPhotonEnergies,quartz_EFIC,nEntriesWLS);
   AlSurface->SetMaterialPropertiesTable(AlSurfaceProperty);

  TileSurface = new G4OpticalSurface("TileSurface",unified, polished, dielectric_dielectric,1);

   G4MaterialPropertiesTable* TileSurfaceProperty = new G4MaterialPropertiesTable();
   TileSurfaceProperty->AddProperty("RINDEX",AllPhotonEnergies,quartz_RIND,nEntriesWLS);
   //TileSurfaceProperty->AddProperty("REFLECTIVITY",AllPhotonEnergies,reflectivity,nEntriesWLS);
   TileSurfaceProperty->AddProperty("EFFICIENCY",AllPhotonEnergies,quartz_EFIC,nEntriesWLS);
   TileSurface->SetMaterialPropertiesTable(TileSurfaceProperty);

   photonDetSurface = new G4OpticalSurface("PhotonDetSurface", glisur, ground, dielectric_metal,1);

   G4MaterialPropertiesTable* photonDetSurfaceProperty = new G4MaterialPropertiesTable();

   G4double p_mppc[] = {2.00*eV, 3.47*eV};
   const G4int nbins = sizeof(p_mppc)/sizeof(G4double);
   G4double refl_mppc[] = {1, 1};
   assert(sizeof(refl_mppc) == sizeof(p_mppc));
   G4double effi_mppc[] = {1, 1};
   assert(sizeof(effi_mppc) == sizeof(p_mppc));

   photonDetSurfaceProperty->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,nbins);
   photonDetSurfaceProperty->AddProperty("EFFICIENCY",p_mppc,effi_mppc,nbins);

   photonDetSurface->SetMaterialPropertiesTable(photonDetSurfaceProperty);



}

Materials::~Materials(void)
{
  delete H;
  delete N;
  delete O;
  delete Si;
  delete C;

  delete Xenon;
  delete pQuartz;
  delete EM_Quartz;
  delete Air;
  delete Al;
  delete Cu;
  delete Pb;
  delete W;
  delete Polyethylene;
  delete PMMA;
  delete Grease;

  for (unsigned int i = 0; i < MPT_Array.size(); i++)
      delete MPT_Array.at(i);

  MPT_Array.clear();
}
