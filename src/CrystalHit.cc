#include "CrystalHit.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

G4ThreadLocal G4Allocator<CrystalHit>* CrystalHitAllocator = nullptr;

CrystalHit::CrystalHit()
    : G4VHit(),
      fCrystalID(-1),
      fEdep(0.),
      fPhotonCount(0),
      fPos(G4ThreeVector()),
      fParticleName(""),
      fMomentumDir(G4ThreeVector())
{}

CrystalHit::~CrystalHit() {}

G4bool CrystalHit::operator==(const CrystalHit& right) const
{
    return (this == &right) ? true : false;
}

void CrystalHit::Print()
{
    G4cout << "  Crystal[" << fCrystalID << "]"
           << "  Edep=" << std::setw(7) << G4BestUnit(fEdep, "Energy")
           << "  Photons=" << fPhotonCount
           << "  Pos=" << fPos/mm << " mm"
           << G4endl;
}
