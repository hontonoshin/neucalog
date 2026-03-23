#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "AnalysisManager.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalPhoton.hh"
#include "G4Neutron.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

SteppingAction::SteppingAction(EventAction* eventAction,
                                const DetectorConstruction* detector)
    : G4UserSteppingAction(),
      fEventAction(eventAction),
      fDetector(detector)
{}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    G4ParticleDefinition* particle =
        step->GetTrack()->GetDefinition();

    // ---- Optical photons ----
    // Counting is handled exclusively by CrystalSD to avoid double-counting.
    // Here we only accumulate the per-event total from the SD hit collection,
    // which is done in CrystalSD::EndOfEvent -> EventAction is not used for
    // photon counting at all from SteppingAction.
    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
        // Count photons created inside crystals for EventAction display only.
        // Use post-step volume to correctly catch the birth step.
        G4VPhysicalVolume* postVol =
            step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        if (postVol &&
            postVol->GetLogicalVolume() == fDetector->GetCrystalLogical() &&
            step->GetTrack()->GetCurrentStepNumber() == 1) {
            fEventAction->AddPhotonCount(1);
        }
        return;
    }

    // ---- Determine if pre-step point is inside a crystal ----
    // Use post-step for entering tracks (step 1), pre-step for subsequent steps.
    G4LogicalVolume* preVol =
        step->GetPreStepPoint()->GetTouchableHandle()
            ->GetVolume()->GetLogicalVolume();

    // ---- Steps inside crystal ----
    if (preVol != fDetector->GetCrystalLogical()) return;

    // Energy deposit: accumulate in EventAction for per-event printout.
    // The authoritative histogram fill is done in CrystalSD::EndOfEvent.
    G4double edep = step->GetTotalEnergyDeposit();
    if (edep > 0.) {
        fEventAction->AddEnergyDeposit(edep);
    }

    // Neutron-specific
    if (particle == G4Neutron::NeutronDefinition()) {
        fEventAction->AddNeutronStep();

        // Record neutron energy and PRIMARY direction (at generator, not scatter).
        // We use track vertex info for the angle to get the true beam angle,
        // and the kinetic energy at first entry for the energy spectrum.
        if (step->GetTrack()->GetCurrentStepNumber() == 1) {
            G4double ekin = step->GetPreStepPoint()->GetKineticEnergy();
            GetAM()->FillH1(HistIDs::hNeutronEnergy, ekin / MeV);

            // Use vertex momentum direction = primary generator direction.
            // Reference axis is downward (−Z) for the cosmic-ray geometry:
            // a vertical neutron gives theta=0°, grazing track gives theta->90°.
            // cosTheta_zenith = -dir.z()  (dir.z is negative for downward tracks)
            G4ThreeVector dir = step->GetTrack()->GetVertexMomentumDirection();
            G4double cosTheta = std::min(1.0, std::max(-1.0, -dir.z()));
            G4double theta    = std::acos(cosTheta);
            GetAM()->FillH1(HistIDs::hNeutronAngleTheta,
                            theta / CLHEP::deg);
        }
    }

    // Secondary particle (non-optical) energies
    const G4TrackVector* secondaries = step->GetSecondary();
    if (secondaries) {
        for (auto* sec : *secondaries) {
            if (sec->GetDefinition() !=
                    G4OpticalPhoton::OpticalPhotonDefinition()) {
                G4double ekin = sec->GetKineticEnergy();
                if (ekin > 0.)
                    GetAM()->FillH1(HistIDs::hSecondaryParticles,
                                    ekin / MeV);
            }
        }
    }
}
