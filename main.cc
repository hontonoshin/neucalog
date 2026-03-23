/*===========================================================
  main.cc  —  GAGG Calorimeter Neutron Simulation
  4×4×2 array of GAGG (Gd3Al2Ga3O12) crystals
===========================================================*/

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

int main(int argc, char** argv)
{
    // -------------------------------------------------------
    // Random engine
    // -------------------------------------------------------
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(nullptr));

    // -------------------------------------------------------
    // Run manager  
    // -------------------------------------------------------
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

    // -------------------------------------------------------
    // Mandatory user classes
    // -------------------------------------------------------
    auto* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList());
    runManager->SetUserInitialization(new ActionInitialization(detector));

    // -------------------------------------------------------
    // UI / visualisation
    // -------------------------------------------------------
    G4UIExecutive* ui = nullptr;
    if (argc == 1) {
        // Interactive session — no macro given
        ui = new G4UIExecutive(argc, argv, "tcsh");
    }

    // Vis manager 
    G4VisManager* visManager = new G4VisExecutive("quiet");
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (argc > 1) {
        // Batch mode
        G4String macroFile = argv[1];
        G4String command   = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
    } else {
      
        UImanager->ApplyCommand("/control/execute macros/init.mac");
        ui->SessionStart();
        delete ui;
    }

    // -------------------------------------------------------
    // Cleanup
    // -------------------------------------------------------
    delete visManager;
    delete runManager;

    return 0;
}
