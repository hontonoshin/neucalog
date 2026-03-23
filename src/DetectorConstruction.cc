/*===========================================================
  DetectorConstruction.cc
  4x4x2 array of GAGG crystals (Gd3Al2Ga3O12)
  Crystal size: 2.5 x 2.5 x 3.0 cm
  Total face: 10x10 cm, depth: 6 cm
===========================================================*/

#include "DetectorConstruction.hh"
#include "CrystalSD.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fCrystalLogical(nullptr),
      fWorldLogical(nullptr)
{}

DetectorConstruction::~DetectorConstruction() {}

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nist = G4NistManager::Instance();

    // -------------------------------------------------------
    // GAGG: Gd3Al2Ga3O12  (Gadolinium Aluminium Gallium Garnet)
    // Density: ~6.63 g/cm3
    // -------------------------------------------------------
    G4Element* Gd = nist->FindOrBuildElement("Gd");
    G4Element* Al = nist->FindOrBuildElement("Al");
    G4Element* Ga = nist->FindOrBuildElement("Ga");
    G4Element* O  = nist->FindOrBuildElement("O");

    G4Material* GAGG = new G4Material("GAGG", 6.63*g/cm3, 4,
                                       kStateSolid);
    // Gd3Al2Ga3O12  -> 12 oxygen, 3 Gd, 2 Al, 3 Ga
    GAGG->AddElement(Gd,  3);
    GAGG->AddElement(Al,  2);
    GAGG->AddElement(Ga,  3);
    GAGG->AddElement(O,  12);

    // -------------------------------------------------------
    // Optical properties of GAGG
    // Emission peak ~520 nm (green), decay time ~150-200 ns
    // Photon yield ~40,000 ph/MeV (with Ce:GAGG ~46,000)
    // -------------------------------------------------------
    // Photon energies from 1.5 eV to 4.5 eV
    const G4int nEntries = 12;
    G4double photonEnergy[nEntries] = {
        1.50*eV, 1.80*eV, 2.00*eV, 2.20*eV, 2.40*eV, 2.48*eV,
        2.53*eV, 2.60*eV, 2.80*eV, 3.00*eV, 3.50*eV, 4.00*eV
    };

    // Refractive index of GAGG (~1.91)
    G4double refractiveIndex[nEntries] = {
        1.91, 1.91, 1.91, 1.91, 1.91, 1.92,
        1.92, 1.92, 1.93, 1.94, 1.95, 1.96
    };

    // Absorption length (cm -> mm)
    G4double absLength[nEntries] = {
        380.*mm, 380.*mm, 370.*mm, 360.*mm, 340.*mm, 300.*mm,
        280.*mm, 260.*mm, 200.*mm, 150.*mm, 80.*mm,  40.*mm
    };

    // Emission spectrum (normalised, peak at ~520nm = 2.38 eV)
    G4double emissionSpectrum[nEntries] = {
        0.000, 0.002, 0.010, 0.080, 0.400, 0.850,
        1.000, 0.900, 0.300, 0.060, 0.005, 0.000
    };

    G4MaterialPropertiesTable* mptGAGG = new G4MaterialPropertiesTable();
    mptGAGG->AddProperty("RINDEX",    photonEnergy, refractiveIndex, nEntries);
    mptGAGG->AddProperty("ABSLENGTH", photonEnergy, absLength,       nEntries);
    mptGAGG->AddProperty("SCINTILLATIONCOMPONENT1",
                          photonEnergy, emissionSpectrum, nEntries);

    // Scintillation yield: 40,000 photons/MeV
    // Reduced yield for RAM efficiency (scale factor 200x in analysis)
    // Real GAGG yield: ~40000 ph/MeV. Set low here to avoid OOM.
    mptGAGG->AddConstProperty("SCINTILLATIONYIELD",        200./MeV);
    mptGAGG->AddConstProperty("RESOLUTIONSCALE",           1.0);
    mptGAGG->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 150.*ns);
    mptGAGG->AddConstProperty("SCINTILLATIONYIELD1",        1.0);

    GAGG->SetMaterialPropertiesTable(mptGAGG);

    // -------------------------------------------------------
    // Air for world volume
    // -------------------------------------------------------
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

    // Refractive index of air for optical photons
    G4double airEnergy[2]     = {1.5*eV, 4.0*eV};
    G4double airRIndex[2]     = {1.0003, 1.0003};
    G4MaterialPropertiesTable* mptAir = new G4MaterialPropertiesTable();
    mptAir->AddProperty("RINDEX", airEnergy, airRIndex, 2);
    air->SetMaterialPropertiesTable(mptAir);
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

    G4NistManager* nist = G4NistManager::Instance();
    G4Material* air  = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* GAGG = G4Material::GetMaterial("GAGG");

    // Crystal dimensions
    const G4double crystalX = 2.5*cm / 2.0;   // half-lengths
    const G4double crystalY = 2.5*cm / 2.0;
    const G4double crystalZ = 3.0*cm / 2.0;

    const G4double gap       = 0.1*mm;   // small air gap between crystals
    const G4double pitchXY   = 2.0*crystalX + gap;
    const G4double pitchZ    = 2.0*crystalZ + gap;

    // Total calorimeter envelope
    const G4double caloX = NX * pitchXY / 2.0;
    const G4double caloY = NY * pitchXY / 2.0;
    const G4double caloZ = NZ * pitchZ  / 2.0;

    // World
    G4double worldX = 3.0 * caloX;
    G4double worldY = 3.0 * caloY;
    G4double worldZ = 3.0 * caloZ + 20.*cm;  // extra room for beam entrance

    G4Box* worldSolid = new G4Box("World", worldX, worldY, worldZ);
    fWorldLogical     = new G4LogicalVolume(worldSolid, air, "World");
    G4VPhysicalVolume* worldPhys =
        new G4PVPlacement(nullptr, G4ThreeVector(), fWorldLogical,
                          "World", nullptr, false, 0, true);

    // Calorimeter envelope (air mother)
    G4Box* caloSolid = new G4Box("Calorimeter", caloX, caloY, caloZ);
    G4LogicalVolume* caloLogical =
        new G4LogicalVolume(caloSolid, air, "Calorimeter");
    new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), caloLogical,
                      "Calorimeter", fWorldLogical, false, 0, true);

    // Individual crystal
    G4Box* crystalSolid  = new G4Box("Crystal", crystalX, crystalY, crystalZ);
    fCrystalLogical      = new G4LogicalVolume(crystalSolid, GAGG, "Crystal");

    // Reflective surface on crystal skin (TiO2-like reflector)
    G4OpticalSurface* reflSurface = new G4OpticalSurface("ReflectorSurface");
    reflSurface->SetType(dielectric_metal);
    reflSurface->SetFinish(polished);
    reflSurface->SetModel(glisur);

    G4double reflEnergy[2]    = {1.5*eV, 4.0*eV};
    G4double reflectivity[2]  = {0.98, 0.98};
    G4double efficiency[2]    = {0.0,  0.0};
    G4MaterialPropertiesTable* mptRefl = new G4MaterialPropertiesTable();
    mptRefl->AddProperty("REFLECTIVITY", reflEnergy, reflectivity, 2);
    mptRefl->AddProperty("EFFICIENCY",   reflEnergy, efficiency,   2);
    reflSurface->SetMaterialPropertiesTable(mptRefl);
    new G4LogicalSkinSurface("ReflectorSkin", fCrystalLogical, reflSurface);

    // Place 4x4x2 array
    G4int copyNum = 0;
    for (G4int iz = 0; iz < NZ; iz++) {
        for (G4int iy = 0; iy < NY; iy++) {
            for (G4int ix = 0; ix < NX; ix++) {
                G4double xPos = (ix - (NX-1)/2.0) * pitchXY;
                G4double yPos = (iy - (NY-1)/2.0) * pitchXY;
                G4double zPos = (iz - (NZ-1)/2.0) * pitchZ;

                new G4PVPlacement(nullptr,
                                  G4ThreeVector(xPos, yPos, zPos),
                                  fCrystalLogical,
                                  "Crystal",
                                  caloLogical,
                                  false,
                                  copyNum,
                                  true);
                copyNum++;
            }
        }
    }

    // -------------------------------------------------------
    // Visualization attributes  (2D/wireframe-friendly)
    // -------------------------------------------------------
    fWorldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VisAttributes* caloVis = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9, 0.1));
    caloVis->SetVisibility(true);
    caloLogical->SetVisAttributes(caloVis);

    G4VisAttributes* crystalVis = new G4VisAttributes(G4Colour(0.0, 0.8, 0.4, 0.4));
    crystalVis->SetVisibility(true);
    crystalVis->SetForceSolid(false);
    crystalVis->SetForceWireframe(true);  // wireframe — lighter on GPU
    fCrystalLogical->SetVisAttributes(crystalVis);

    return worldPhys;
}

void DetectorConstruction::ConstructSDandField()
{
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    sdManager->SetVerboseLevel(1);

    CrystalSD* crystalSD = new CrystalSD("CrystalSD", "CrystalHitsCollection");
    sdManager->AddNewDetector(crystalSD);
    SetSensitiveDetector("Crystal", crystalSD, true);
}
