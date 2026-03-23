#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector)
    : G4VUserActionInitialization(),
      fDetector(detector)
{}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::BuildForMaster() const
{
    RunAction* runAction = new RunAction();
    SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction());

    RunAction*   runAction   = new RunAction();
    EventAction* eventAction = new EventAction(runAction);

    SetUserAction(runAction);
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(eventAction, fDetector));
}
