#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

/*===========================================================
  PrimaryGeneratorAction.hh
  Cosmic-ray-like atmospheric neutron source.

  Physics model
  -------------
  Energy spectrum: parameterised Goldhagen / Gordon atmospheric
  neutron spectrum at sea level (New York, solar minimum).
  The spectrum has four components sampled via rejection:
    1. Thermal peak      ~25 meV   (Maxwellian)
    2. Epithermal slope  1/E       (eV – 100 keV)
    3. Evaporation peak  ~1 MeV   (Maxwell-Boltzmann in E)
    4. High-energy tail  power-law (1–10 GeV)
  Relative normalisation matches the tabulated flux ratios
  from Gordon et al. (2004) IEEE TNS 51(6) 3427.

  Angular distribution: downward-biased cosine distribution
  matching the cos²(θ) zenith dependence of cosmic-ray
  secondary neutrons at ground level.  θ is measured from
  the downward vertical (−Z axis in our geometry).
  Neutrons arrive from a hemisphere 30 cm above the detector.

  Spatial sampling: uniform random point on a disk of radius
  r_disk (= half-diagonal of calorimeter face + 2 cm margin)
  centred 20 cm above the calorimeter top face, so the full
  detector is illuminated uniformly from above.
===========================================================*/

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    G4ParticleGun* fParticleGun;

    // Source plane geometry (set in constructor from detector dimensions)
    G4double fSourceZ;    // Z of source plane (above detector)
    G4double fDiskRadius; // sampling disk radius on source plane

    // Sample neutron kinetic energy from the atmospheric spectrum
    G4double SampleEnergy() const;

    // Sample a downward direction with cos²(θ) zenith weighting
    // Returns a unit vector pointing downward into the detector
    G4ThreeVector SampleDirection() const;
};

#endif
