#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

class RunAction;

class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction* runAction);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event) override;
    virtual void EndOfEventAction(const G4Event* event) override;

    void AddEnergyDeposit(G4double edep)  { fEdep += edep; }
    void AddPhotonCount(G4int n)          { fPhotonCount += n; }
    void AddNeutronStep()                 { fNeutronSteps++; }

private:
    RunAction* fRunAction;
    G4double   fEdep;
    G4int      fPhotonCount;
    G4int      fNeutronSteps;
};

#endif
