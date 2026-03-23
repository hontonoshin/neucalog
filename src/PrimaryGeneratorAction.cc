/*===========================================================
  PrimaryGeneratorAction.cc
  Cosmic-ray-like atmospheric neutron source

  Energy spectrum
  ---------------
  Parameterised Goldhagen / Gordon sea-level atmospheric
  neutron spectrum (New York, solar minimum).  The analytic
  form follows Eq. (5) of:
    Gordon et al., IEEE TNS 51(6), 3427 (2004)
  with four spectral components:

    Φ(E) = A₁·E·exp(−E/E_th)              [thermal, Maxwellian]
          + A₂ / E                         [epithermal, 1/E]
          + A₃·E·exp(−E/T_nuc)            [evaporation peak ~1 MeV]
          + A₄·E^(−γ)                      [high-energy power law]

  Component weight fractions (fraction of total flux):
    thermal      6 %   (25.3 meV peak)
    epithermal  21 %   (1/E from 0.5 eV – 100 keV)
    evaporation 57 %   (~1–5 MeV peak)
    high-energy 16 %   (power law up to 10 GeV)

  Sampled analytically (thermal+evaporation: Erlang k=2 trick;
  epithermal: log-uniform; high-energy: power-law inverse CDF).

  Angular distribution
  --------------------
  Secondary cosmic-ray neutrons at ground level follow a
  cos^n(θ) zenith distribution with n ≈ 2–4 (Heidbreder 1971).
  We use n = 2:

    p(cosθ) ∝ cos²θ,   cosθ ∈ [0,1]  (downward hemisphere only)

  Sampled analytically: cosθ = ξ^(1/3), ξ ∈ (0,1] uniform.
  Azimuth φ is uniform in [0, 2π).

  Spatial sampling
  ----------------
  Uniform random point on a disk at height fSourceZ (20 cm
  above the calorimeter top face).  Disk radius = half-diagonal
  of the calorimeter XY face + 2 cm margin, guaranteeing that
  oblique tracks still illuminate the full detector.
===========================================================*/

#include "PrimaryGeneratorAction.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

#include <cmath>

// ---------------------------------------------------------------------------
// Spectral constants (energies in MeV throughout)
// ---------------------------------------------------------------------------
namespace SpecPar {
    // Sub-range boundaries
    constexpr double E_thMax  =  0.5e-6;   // thermal top  = 0.5 eV
    constexpr double E_epMax  =  0.1;      // epithermal top = 100 keV
    constexpr double E_evMax  = 20.0;      // evaporation top = 20 MeV
    constexpr double E_hiMax  = 1.0e4;     // high-energy top = 10 GeV

    // Characteristic energies
    constexpr double E_th     =  25.3e-9;  // thermal kT = 25.3 meV
    constexpr double T_nuc    =   1.0;     // nuclear temperature = 1 MeV

    // High-energy spectral index
    constexpr double gamma_hi =   0.9;

    // Component weight fractions (must sum to 1)
    constexpr double w_th  = 0.06;
    constexpr double w_ep  = 0.21;
    constexpr double w_ev  = 0.57;
    // w_hi = 1 - w_th - w_ep - w_ev = 0.16 (used as the else branch)
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(new G4ParticleGun(1))
{
    G4ParticleTable*      table   = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* neutron = table->FindParticle("neutron");
    fParticleGun->SetParticleDefinition(neutron);

    // Calorimeter geometry (must match DetectorConstruction):
    //   NX=4, NY=4  crystals, pitch 2.505 cm → face ≈ 10.02 cm × 10.02 cm
    //   NZ=2 layers, crystal half-length = 3.0 cm, gap = 0.1 mm
    //   Total Z depth ≈ 2 × (2×3.0 cm) + gap ≈ 12.01 cm
    //   Top face at z ≈ +6.005 cm ≈ +6 cm
    //   Source plane: 20 cm above top face → z = +26 cm
    fSourceZ    =  26.0 * cm;

    // Disk: half-diagonal of 10 cm × 10 cm face = sqrt(50) ≈ 7.07 cm
    //       + 2 cm margin for oblique tracks
    fDiskRadius =   9.1 * cm;

    // Reasonable defaults (all overridden in GeneratePrimaries)
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., fSourceZ));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
    fParticleGun->SetParticleEnergy(1.0 * MeV);
}

// ---------------------------------------------------------------------------
// Destructor
// ---------------------------------------------------------------------------
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

// ---------------------------------------------------------------------------
// GeneratePrimaries — called once per event
// ---------------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    fParticleGun->SetParticleEnergy(SampleEnergy());
    fParticleGun->SetParticleMomentumDirection(SampleDirection());

    // Uniform sampling on source disk (rejection method)
    G4double x0, y0;
    do {
        x0 = (2.0*G4UniformRand() - 1.0) * fDiskRadius;
        y0 = (2.0*G4UniformRand() - 1.0) * fDiskRadius;
    } while (x0*x0 + y0*y0 > fDiskRadius*fDiskRadius);

    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, fSourceZ));
    fParticleGun->GeneratePrimaryVertex(event);
}

// ---------------------------------------------------------------------------
// SampleEnergy — four-component Gordon atmospheric neutron spectrum
// ---------------------------------------------------------------------------
G4double PrimaryGeneratorAction::SampleEnergy() const
{
    using namespace SpecPar;
    G4double r = G4UniformRand();

    // ------------------------------------------------------------------
    // 1. Thermal  f(E) ∝ E·exp(−E/E_th)  (Maxwellian, peak at E_th)
    //    Erlang k=2 trick: if U1,U2 ~ Exp(1/E_th), then E = U1+U2
    //    has pdf ∝ E·exp(−E/E_th).  Here Ui = −E_th·ln(ui).
    // ------------------------------------------------------------------
    if (r < w_th) {
        G4double u1 = std::max(G4UniformRand(), 1e-30);
        G4double u2 = std::max(G4UniformRand(), 1e-30);
        G4double E  = E_th * (-std::log(u1) - std::log(u2));
        return std::min(E, E_thMax);  // clamp to thermal range
    }

    // ------------------------------------------------------------------
    // 2. Epithermal  f(E) ∝ 1/E,  E ∈ [E_thMax, E_epMax]
    //    Inverse CDF: E = E_thMax · (E_epMax/E_thMax)^u
    // ------------------------------------------------------------------
    if (r < w_th + w_ep) {
        G4double u     = G4UniformRand();
        G4double ratio = E_epMax / E_thMax;
        return E_thMax * std::pow(ratio, u);
    }

    // ------------------------------------------------------------------
    // 3. Evaporation  f(E) ∝ E·exp(−E/T_nuc)  (nuclear evaporation)
    //    Same Erlang k=2 trick, rejection-clamped to [E_epMax, E_evMax].
    //    The bulk of this distribution sits around 1–2 MeV so acceptance
    //    is high (>90 %).
    // ------------------------------------------------------------------
    if (r < w_th + w_ep + w_ev) {
        G4double E;
        G4int    tries = 0;
        do {
            G4double u1 = std::max(G4UniformRand(), 1e-30);
            G4double u2 = std::max(G4UniformRand(), 1e-30);
            E = T_nuc * (-std::log(u1) - std::log(u2));
            ++tries;
        } while ((E < E_epMax || E > E_evMax) && tries < 500);
        // Clamp on unlikely loop exit
        return std::max(E_epMax, std::min(E, E_evMax));
    }

    // ------------------------------------------------------------------
    // 4. High-energy power law  f(E) ∝ E^(−γ),  E ∈ [E_evMax, E_hiMax]
    //    Inverse CDF: E = [Ea + u·(Eb - Ea)]^(1/α)
    //    where α = 1−γ and Ea = E_evMax^α, Eb = E_hiMax^α.
    // ------------------------------------------------------------------
    {
        G4double alpha = 1.0 - gamma_hi;   // = 0.1
        G4double Ea    = std::pow(E_evMax, alpha);
        G4double Eb    = std::pow(E_hiMax, alpha);
        G4double u     = G4UniformRand();
        return std::pow(Ea + u*(Eb - Ea), 1.0/alpha);
    }
}

// ---------------------------------------------------------------------------
// SampleDirection — cos²(θ) downward zenith distribution
// ---------------------------------------------------------------------------
G4ThreeVector PrimaryGeneratorAction::SampleDirection() const
{
    // p(cosθ) ∝ cos²θ  on  cosθ ∈ [0,1]
    // CDF: F(u) = u³  → inverse: cosθ = ξ^(1/3)
    G4double xi       = std::max(G4UniformRand(), 1e-30);
    G4double cosTheta = std::cbrt(xi);
    G4double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta*cosTheta));

    // Uniform azimuth
    G4double phi = CLHEP::twopi * G4UniformRand();

    // Neutrons travel DOWNWARD in our geometry (−Z direction)
    return G4ThreeVector( sinTheta * std::cos(phi),
                          sinTheta * std::sin(phi),
                         -cosTheta );
}
