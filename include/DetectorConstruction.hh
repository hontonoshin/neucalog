#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

// 4x4x2 matrix of GAGG crystals
// Each crystal: 2.5cm x 2.5cm x 3.0cm
// Total face: 10cm x 10cm, thickness: 6cm (2 layers)

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct() override;
    virtual void ConstructSDandField() override;

    G4LogicalVolume* GetCrystalLogical() const { return fCrystalLogical; }

    static const G4int NX = 4;   // columns
    static const G4int NY = 4;   // rows
    static const G4int NZ = 2;   // layers

private:
    void DefineMaterials();

    G4LogicalVolume* fCrystalLogical;
    G4LogicalVolume* fWorldLogical;
};

#endif
