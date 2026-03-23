/*===========================================================
  PhysicsList.cc  — Geant4 11.3-patch-02 compatible
===========================================================*/

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProductionCutsTable.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
    SetVerboseLevel(1);

    RegisterPhysics(new G4EmStandardPhysics(1));
    RegisterPhysics(new G4EmExtraPhysics(1));
    RegisterPhysics(new G4OpticalPhysics());
    RegisterPhysics(new G4DecayPhysics(1));
    RegisterPhysics(new G4RadioactiveDecayPhysics(1));
    RegisterPhysics(new G4HadronElasticPhysicsHP(1));
    RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP(1));
    RegisterPhysics(new G4StoppingPhysics(1));
    RegisterPhysics(new G4IonPhysics(1));
    RegisterPhysics(new G4NeutronTrackingCut(1));
}

PhysicsList::~PhysicsList() {}

void PhysicsList::ConstructParticle()
{
    G4BaryonConstructor  bC;  bC.ConstructParticle();
    G4BosonConstructor   boC; boC.ConstructParticle();
    G4LeptonConstructor  lC;  lC.ConstructParticle();
    G4MesonConstructor   mC;  mC.ConstructParticle();
    G4IonConstructor     iC;  iC.ConstructParticle();
    G4ShortLivedConstructor slC; slC.ConstructParticle();
    G4OpticalPhoton::OpticalPhotonDefinition();
    G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    // -------------------------------------------------------
    // Optical process configuration via G4OpticalParameters.
    //
    // Important: in Geant4 11.x there is NO per-step cap for
    // scintillation photons in the API — that is controlled
    // entirely by SCINTILLATIONYIELD in the material properties
    // table (already set to 200/MeV in DetectorConstruction,
    // which limits yield to ~200 × Edep photons per event).
    //
    // The Cerenkov cap prevents runaway photon production from
    // fast charged secondaries (electrons, protons) that can
    // otherwise emit thousands of Cerenkov photons per step.
    // -------------------------------------------------------
    auto* optParams = G4OpticalParameters::Instance();

    // Cap Cerenkov photons per tracking step
    optParams->SetCerenkovMaxPhotonsPerStep(100);

    // Limit beta change per step for Cerenkov (keeps stepping fast)
    optParams->SetCerenkovMaxBetaChange(10.0);

    // Track scintillation secondaries before continuing the primary —
    // ensures photons are created and recorded before the step moves on
    optParams->SetScintTrackSecondariesFirst(true);

    // Silence optical process verbose output
    optParams->SetVerboseLevel(0);
}

void PhysicsList::SetCuts()
{
    SetCutValue(1.0*mm, "gamma");
    SetCutValue(1.0*mm, "e-");
    SetCutValue(1.0*mm, "e+");
    SetCutValue(0.1*mm, "proton");
    G4ProductionCutsTable::GetProductionCutsTable()
        ->SetEnergyRange(100.*eV, 100.*GeV);
}
