#include "RunAction.hh"
#include "AnalysisManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"

RunAction::RunAction()
    : G4UserRunAction(),
      fNeutronFluence(0.),
      fPhotonYield(0.),
      fEdep(0.)
{
    // Use Register() — RegisterAccumulable() is deprecated in newer Geant4 11
    G4AccumulableManager* accMan = G4AccumulableManager::Instance();
    accMan->Register(fNeutronFluence);
    accMan->Register(fPhotonYield);
    accMan->Register(fEdep);
}

RunAction::~RunAction() {}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    G4AccumulableManager::Instance()->Reset();

    auto am = G4AnalysisManager::Instance();
    BookHistograms();
    am->SetDefaultFileType("root");
    am->OpenFile("GAGG_output.root");

    G4cout << "\n=== GAGG Calorimeter Simulation Starting ===" << G4endl;
    G4cout << "    4x4x2 GAGG crystal array" << G4endl;
    G4cout << "    Output: GAGG_output.root" << G4endl;
    G4cout << "============================================\n" << G4endl;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    G4int nEvents = run->GetNumberOfEvent();
    if (nEvents == 0) return;

    G4AccumulableManager::Instance()->Merge();

    auto am = G4AnalysisManager::Instance();
    am->Write();
    am->CloseFile();

    G4cout << "\n=== Run Summary ===" << G4endl;
    G4cout << "  Events processed:  " << nEvents << G4endl;
    G4cout << "  Total Edep:        "
           << fEdep.GetValue() / CLHEP::MeV << " MeV" << G4endl;
    G4cout << "  Total photons:     "
           << (G4int)fPhotonYield.GetValue() << G4endl;
    if (nEvents > 0) {
        G4cout << "  Avg Edep/event:    "
               << fEdep.GetValue() / nEvents / CLHEP::MeV << " MeV" << G4endl;
        G4cout << "  Avg photons/event: "
               << fPhotonYield.GetValue() / nEvents << G4endl;
    }
    G4cout << "===================\n" << G4endl;
}
