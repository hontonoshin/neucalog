#ifndef CrystalHit_h
#define CrystalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class CrystalHit : public G4VHit
{
public:
    CrystalHit();
    CrystalHit(const CrystalHit&) = default;
    virtual ~CrystalHit();

    CrystalHit& operator=(const CrystalHit&) = default;
    G4bool operator==(const CrystalHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    virtual void Draw()  override {}
    virtual void Print() override;

    // Setters
    void SetCrystalID(G4int id)             { fCrystalID = id; }
    void SetEdep(G4double edep)             { fEdep = edep; }
    void AddEdep(G4double edep)             { fEdep += edep; }
    void SetPhotonCount(G4int n)            { fPhotonCount = n; }
    void AddPhoton()                        { fPhotonCount++; }
    void SetPosition(G4ThreeVector pos)     { fPos = pos; }
    void SetParticleName(G4String name)     { fParticleName = name; }
    void SetMomentumDir(G4ThreeVector dir)  { fMomentumDir = dir; }

    // Getters
    G4int          GetCrystalID()     const { return fCrystalID; }
    G4double       GetEdep()          const { return fEdep; }
    G4int          GetPhotonCount()   const { return fPhotonCount; }
    G4ThreeVector  GetPosition()      const { return fPos; }
    G4String       GetParticleName()  const { return fParticleName; }
    G4ThreeVector  GetMomentumDir()   const { return fMomentumDir; }

private:
    G4int         fCrystalID;
    G4double      fEdep;
    G4int         fPhotonCount;
    G4ThreeVector fPos;
    G4String      fParticleName;
    G4ThreeVector fMomentumDir;
};

using CrystalHitsCollection = G4THitsCollection<CrystalHit>;

extern G4ThreadLocal G4Allocator<CrystalHit>* CrystalHitAllocator;

inline void* CrystalHit::operator new(size_t)
{
    if (!CrystalHitAllocator)
        CrystalHitAllocator = new G4Allocator<CrystalHit>;
    return (void*) CrystalHitAllocator->MallocSingle();
}

inline void CrystalHit::operator delete(void* hit)
{
    CrystalHitAllocator->FreeSingle((CrystalHit*) hit);
}

#endif
