#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*) override;
    virtual void EndOfRunAction(const G4Run*) override;

    void AddNeutronFluence(G4double val) { fNeutronFluence += val; }
    void AddPhotonYield(G4double val)    { fPhotonYield    += val; }
    void AddEnergyDeposit(G4double val)  { fEdep           += val; }

private:
    G4Accumulable<G4double> fNeutronFluence;
    G4Accumulable<G4double> fPhotonYield;
    G4Accumulable<G4double> fEdep;
};

#endif
