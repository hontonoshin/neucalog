#include "EventAction.hh"
#include "RunAction.hh"
#include "AnalysisManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

EventAction::EventAction(RunAction* runAction)
    : G4UserEventAction(),
      fRunAction(runAction),
      fEdep(0.),
      fPhotonCount(0),
      fNeutronSteps(0)
{}

EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    fEdep         = 0.;
    fPhotonCount  = 0;
    fNeutronSteps = 0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    fRunAction->AddEnergyDeposit(fEdep);
    fRunAction->AddPhotonYield((G4double)fPhotonCount);

    G4int evtID = event->GetEventID();
    if (evtID % 100 == 0) {
        // G4BestUnit needs the G4UnitsTable header — included above
        G4cout << ">>> Event " << evtID
               << "  Edep=" << fEdep / CLHEP::MeV << " MeV"
               << "  Photons=" << fPhotonCount
               << G4endl;
    }
}
