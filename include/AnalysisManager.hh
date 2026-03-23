#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "G4AnalysisManager.hh"
#include "globals.hh"

// Centralized histogram/ntuple IDs
namespace HistIDs {
    // 1D histograms
    constexpr G4int hNeutronEnergy        = 0;  // Incoming neutron energy
    constexpr G4int hNeutronAngleTheta    = 1;  // Neutron polar angle at entrance
    constexpr G4int hPhotonYieldPerEvent  = 2;  // Optical photons per event
    constexpr G4int hPhotonYieldPerCrystal= 3;  // Optical photons per crystal hit
    constexpr G4int hEdepPerEvent         = 4;  // Energy deposit per event
    constexpr G4int hEdepPerCrystal       = 5;  // Energy deposit per crystal
    constexpr G4int hPhotonAngleTheta     = 6;  // Optical photon polar angle
    constexpr G4int hPhotonAnglePhi       = 7;  // Optical photon azimuthal angle
    constexpr G4int hNeutronFluxZ         = 8;  // Neutron flux vs Z depth
    constexpr G4int hSecondaryParticles   = 9;  // Secondary particle types
    constexpr G4int hNeutronScatterAngle  = 10; // Neutron scattering angle
    constexpr G4int hPhotonEnergy         = 11; // Optical photon energy spectrum

    // 2D histograms
    constexpr G4int h2EdepXY              = 0;  // Edep map XY (front face)
    constexpr G4int h2PhotonYieldXY       = 1;  // Photon yield map XY
    constexpr G4int h2NeutronAngleEnergy  = 2;  // Neutron angle vs energy
    constexpr G4int h2PhotonAngle2D       = 3;  // Photon theta vs phi

    // Ntuples
    constexpr G4int ntEvent               = 0;  // Per-event summary
    constexpr G4int ntPhoton              = 1;  // Per-optical-photon data
    constexpr G4int ntNeutron             = 2;  // Neutron track data
}

// Convenience wrapper — just grab the singleton
inline G4AnalysisManager* GetAM() {
    return G4AnalysisManager::Instance();
}

void BookHistograms();   // Called in RunAction::BeginOfRunAction

#endif
