/*===========================================================
  AnalysisManager.cc
  Books all histograms and ntuples at BeginOfRun.
===========================================================*/

#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

void BookHistograms()
{
    auto am = G4AnalysisManager::Instance();
    am->SetVerboseLevel(1);
    am->SetNtupleMerging(true);

    // -------------------------------------------------------
    // 1D Histograms
    // -------------------------------------------------------

    // hNeutronEnergy [0]: Incoming neutron kinetic energy (1 keV – 20 MeV)
    // Starting at 1e-3 avoids log-scale spike from thermalized neutrons at 0.
    am->CreateH1("NeutronEnergy",
                 "Neutron kinetic energy at crystal entry;E_{kin} [MeV];Counts",
                 200, 1.e-3, 20.);

    // hNeutronAngleTheta [1]: Polar angle of primary neutron
    am->CreateH1("NeutronAngleTheta",
                 "Neutron polar angle (wrt Z);#theta [deg];Counts",
                 180, 0., 180.);

    // hPhotonYieldPerEvent [2]
    am->CreateH1("PhotonYieldPerEvent",
                 "Optical photon yield per event;N_{photons};Events",
                 500, 0., 50000.);

    // hPhotonYieldPerCrystal [3]
    am->CreateH1("PhotonYieldPerCrystal",
                 "Optical photon yield per crystal per event;N_{photons};Entries",
                 300, 0., 30000.);

    // hEdepPerEvent [4]
    am->CreateH1("EdepPerEvent",
                 "Total energy deposit per event;E_{dep} [MeV];Events",
                 200, 0., 20.);

    // hEdepPerCrystal [5]
    am->CreateH1("EdepPerCrystal",
                 "Energy deposit per crystal;E_{dep} [MeV];Entries",
                 200, 0., 10.);

    // hPhotonAngleTheta [6]
    am->CreateH1("PhotonAngleTheta",
                 "Optical photon polar angle;#theta [deg];Counts",
                 180, 0., 180.);

    // hPhotonAnglePhi [7]
    am->CreateH1("PhotonAnglePhi",
                 "Optical photon azimuthal angle;#phi [deg];Counts",
                 360, -180., 180.);

    // hNeutronFluxZ [8]: neutron steps vs Z position
    am->CreateH1("NeutronFluxZ",
                 "Neutron flux profile along Z;Z [cm];Steps",
                 120, -6., 6.);

    // hSecondaryParticles [9]: secondary particle IDs (dummy — text in ntuple)
    am->CreateH1("SecondaryEkin",
                 "Secondary particle kinetic energy;E_{kin} [MeV];Counts",
                 200, 0., 10.);

    // hNeutronScatterAngle [10]
    am->CreateH1("NeutronScatterAngle",
                 "Neutron step polar angle in detector;#theta [deg];Steps",
                 180, 0., 180.);

    // hPhotonEnergy [11]
    am->CreateH1("PhotonEnergy",
                 "Optical photon energy;E [eV];Counts",
                 200, 1., 5.);

    // -------------------------------------------------------
    // 2D Histograms
    // -------------------------------------------------------

    // h2EdepXY [0]: energy deposit map on XY face
    am->CreateH2("EdepMapXY",
                 "Energy deposit map (XY projection);X [cm];Y [cm]",
                 40, -5., 5.,
                 40, -5., 5.);

    // h2PhotonYieldXY [1]
    am->CreateH2("PhotonYieldMapXY",
                 "Optical photon yield map (XY);X [cm];Y [cm]",
                 40, -5., 5.,
                 40, -5., 5.);

    // h2NeutronAngleEnergy [2]
    am->CreateH2("NeutronAngleVsEnergy",
                 "Neutron angle vs kinetic energy;E_{kin} [MeV];#theta [deg]",
                 100, 0., 10.,
                 90,  0., 180.);

    // h2PhotonAngle2D [3]
    am->CreateH2("PhotonAngle2D",
                 "Optical photon angular distribution;#theta [deg];#phi [deg]",
                 90, 0., 180.,
                 72, -180., 180.);

    // -------------------------------------------------------
    // Ntuples
    // -------------------------------------------------------

    // ntEvent [0]: per-event summary
    am->CreateNtuple("EventSummary", "Per-event summary");
    am->CreateNtupleDColumn("TotalEdep_MeV");      // col 0
    am->CreateNtupleDColumn("TotalPhotons");       // col 1
    am->FinishNtuple();

    // ntPhoton [1]: per optical photon
    am->CreateNtuple("OpticalPhoton", "Per-optical-photon data");
    am->CreateNtupleDColumn("Energy_eV");          // col 0
    am->CreateNtupleDColumn("Theta_deg");          // col 1
    am->CreateNtupleDColumn("Phi_deg");            // col 2
    am->CreateNtupleDColumn("CrystalID");          // col 3
    am->FinishNtuple();

    // ntNeutron [2]: per neutron step in crystal
    am->CreateNtuple("NeutronStep", "Neutron step data in crystals");
    am->CreateNtupleDColumn("Ekin_MeV");           // col 0
    am->CreateNtupleDColumn("Theta_deg");          // col 1
    am->CreateNtupleDColumn("X_cm");               // col 2
    am->CreateNtupleDColumn("Y_cm");               // col 3
    am->CreateNtupleDColumn("Z_cm");               // col 4
    am->FinishNtuple();
}
