/*===========================================================
  CrystalSD.cc
  Sensitive detector for GAGG crystals.
  Records:
    - Energy deposits per crystal
    - Optical photon counts (from scintillation secondaries)
    - Particle type, position, momentum direction
===========================================================*/

#include "CrystalSD.hh"
#include "AnalysisManager.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhoton.hh"

CrystalSD::CrystalSD(const G4String& name,
                     const G4String& hitsCollectionName)
    : G4VSensitiveDetector(name),
      fHitsCollection(nullptr),
      fHCID(-1)
{
    collectionName.insert(hitsCollectionName);
}

CrystalSD::~CrystalSD() {}

void CrystalSD::Initialize(G4HCofThisEvent* hce)
{
    fHitsCollection = new CrystalHitsCollection(SensitiveDetectorName,
                                                collectionName[0]);
    if (fHCID < 0)
        fHCID = G4SDManager::GetSDMpointer()
                    ->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(fHCID, fHitsCollection);

    // Create one hit per crystal (32 total)
    for (G4int i = 0; i < 32; i++) {
        CrystalHit* hit = new CrystalHit();
        hit->SetCrystalID(i);
        fHitsCollection->insert(hit);
    }
}

G4bool CrystalSD::ProcessHits(G4Step* step,
                               G4TouchableHistory* /*history*/)
{
    G4ParticleDefinition* particle =
        step->GetTrack()->GetDefinition();

    // ---- Optical photons ----
    if (particle == G4OpticalPhoton::OpticalPhotonDefinition()) {
        if (step->GetTrack()->GetCurrentStepNumber() == 1) {
            // For a newly born photon, the pre-step touchable is the crystal
            // where it was created (scintillation birth point).
            G4int copyNum = step->GetPreStepPoint()
                                ->GetTouchableHandle()
                                ->GetCopyNumber();
            if (copyNum < 0 || copyNum >= 32) return false;
            CrystalHit* hit = (*fHitsCollection)[copyNum];
            if (!hit) return false;

            hit->AddPhoton();

            G4ThreeVector dir = step->GetPreStepPoint()->GetMomentumDirection();
            G4double energy   = step->GetPreStepPoint()->GetKineticEnergy();

            auto am = GetAM();
            G4double cosTheta = dir.z();
            G4double theta    = std::acos(std::min(1.0, std::max(-1.0, cosTheta)));
            G4double phi      = std::atan2(dir.y(), dir.x());

            am->FillH1(HistIDs::hPhotonAngleTheta, theta / CLHEP::deg);
            am->FillH1(HistIDs::hPhotonAnglePhi,   phi   / CLHEP::deg);
            am->FillH1(HistIDs::hPhotonEnergy,     energy / eV);
            am->FillH2(HistIDs::h2PhotonAngle2D,
                       theta / CLHEP::deg,
                       phi   / CLHEP::deg);

            am->FillNtupleDColumn(HistIDs::ntPhoton, 0, energy / eV);
            am->FillNtupleDColumn(HistIDs::ntPhoton, 1, theta  / CLHEP::deg);
            am->FillNtupleDColumn(HistIDs::ntPhoton, 2, phi    / CLHEP::deg);
            am->FillNtupleDColumn(HistIDs::ntPhoton, 3, (G4double)copyNum);
            am->AddNtupleRow(HistIDs::ntPhoton);
        }
        return true;
    }

    // ---- Other particles: record energy deposit ----
    // For a particle entering the crystal on its first step, the pre-step
    // point touchable is the MOTHER volume (air gap / calorimeter envelope),
    // not the crystal. We must use the post-step touchable for step 1,
    // and the pre-step touchable for subsequent steps inside the crystal.
    G4int copyNum = -1;
    G4int stepNum = step->GetTrack()->GetCurrentStepNumber();
    if (stepNum == 1) {
        // Entering step: post-step point is inside the crystal
        G4VPhysicalVolume* postPV =
            step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        if (!postPV) return false;
        copyNum = step->GetPostStepPoint()
                      ->GetTouchableHandle()
                      ->GetCopyNumber();
    } else {
        copyNum = step->GetPreStepPoint()
                      ->GetTouchableHandle()
                      ->GetCopyNumber();
    }

    if (copyNum < 0 || copyNum >= 32) return false;

    CrystalHit* hit = (*fHitsCollection)[copyNum];
    if (!hit) return false;

    G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return false;

    hit->AddEdep(edep);
    hit->SetPosition(step->GetPreStepPoint()->GetPosition());
    hit->SetParticleName(particle->GetParticleName());

    auto am = GetAM();
    G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();
    am->FillH2(HistIDs::h2EdepXY, pos.x()/cm, pos.y()/cm, edep/MeV);

    // Neutron-specific
    if (particle->GetParticleName() == "neutron") {
        G4ThreeVector dir = step->GetPreStepPoint()->GetMomentumDirection();
        G4double cosTheta = dir.z();
        G4double theta    = std::acos(std::min(1.0, std::max(-1.0, cosTheta)));
        am->FillH1(HistIDs::hNeutronScatterAngle, theta / CLHEP::deg);
        am->FillH1(HistIDs::hNeutronFluxZ,        pos.z()/cm);

        G4double ekin = step->GetPreStepPoint()->GetKineticEnergy();
        am->FillNtupleDColumn(HistIDs::ntNeutron, 0, ekin / MeV);
        am->FillNtupleDColumn(HistIDs::ntNeutron, 1, theta / CLHEP::deg);
        am->FillNtupleDColumn(HistIDs::ntNeutron, 2, pos.x()/cm);
        am->FillNtupleDColumn(HistIDs::ntNeutron, 3, pos.y()/cm);
        am->FillNtupleDColumn(HistIDs::ntNeutron, 4, pos.z()/cm);
        am->AddNtupleRow(HistIDs::ntNeutron);

        am->FillH2(HistIDs::h2NeutronAngleEnergy,
                   ekin / MeV, theta / CLHEP::deg);
    }

    return true;
}

void CrystalSD::EndOfEvent(G4HCofThisEvent* /*hce*/)
{
    auto am = GetAM();

    G4double totalEdep    = 0.;
    G4int    totalPhotons = 0;

    for (G4int i = 0; i < 32; i++) {
        CrystalHit* hit = (*fHitsCollection)[i];
        G4double    edep = hit->GetEdep();
        G4int       nph  = hit->GetPhotonCount();

        if (edep > 0. || nph > 0) {
            hit->Print();

            am->FillH1(HistIDs::hEdepPerCrystal,         edep / MeV);
            am->FillH1(HistIDs::hPhotonYieldPerCrystal,  (G4double)nph);

            G4ThreeVector pos = hit->GetPosition();
            am->FillH2(HistIDs::h2PhotonYieldXY,
                       pos.x()/cm, pos.y()/cm, (G4double)nph);
        }
        totalEdep    += edep;
        totalPhotons += nph;
    }

    am->FillH1(HistIDs::hEdepPerEvent,        totalEdep    / MeV);
    am->FillH1(HistIDs::hPhotonYieldPerEvent, (G4double)totalPhotons);

    // Per-event ntuple row
    am->FillNtupleDColumn(HistIDs::ntEvent, 0, totalEdep / MeV);
    am->FillNtupleDColumn(HistIDs::ntEvent, 1, (G4double)totalPhotons);
    am->AddNtupleRow(HistIDs::ntEvent);
}
