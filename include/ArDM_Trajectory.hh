#ifndef ARDM_TRAJECTORY
#define ARDM_TRAJECTORY

#define DEFAULTVALUE -11111111

#include "G4VTrajectory.hh"
#include "G4TrajectoryPoint.hh"

#include "G4Track.hh"
#include "vector"

using namespace std;
typedef vector<G4TrajectoryPoint*> ArDM_TrajectoryPointContainer;

class ArDM_Trajectory : public G4VTrajectory {
public:
  ArDM_Trajectory();
  ArDM_Trajectory(const G4Track* track);
  ArDM_Trajectory(ArDM_Trajectory& trajectory);
  ~ArDM_Trajectory();

  G4int         GetTrackID()                   const { return fTrackID;                  };
  G4int         GetParentID()                  const { return fParentID;                 };
  G4String      GetParticleName()              const { return fParticleName;             };
  G4double      GetCharge()                    const { return fParticleCharge;           };
  G4int         GetPDGEncoding()               const { return fPDGEncoding;              };
  G4ThreeVector GetInitialMomentum()           const { return fInitialMomentum;          };
  G4ThreeVector GetInitialMomentumDirection()  const { return fInitialMomentumDirection; };
  G4int         GetPointEntries()              const { return fPositionRecord->size();   };
  G4String      GetVolumeName()                const { return fVolumeName;               };
  G4double      GetEkin()                      const { return fEkin;                     };
  G4double      GetSteplength()                const { return fSteplength;               };
  G4VTrajectoryPoint* GetPoint(G4int i)        const { return (*fPositionRecord)[i];     };
  ArDM_TrajectoryPointContainer* GetPositionRecord() { return fPositionRecord;      };

  void AppendStep(const G4Step* step);
  void MergeTrajectory(G4VTrajectory* secondTrajectory);

  void* operator new(size_t);
  void  operator delete(void*);
  G4bool operator==(const G4VTrajectory& right) const { return (this==&right); };

private:
  G4int         fTrackID;
  G4int         fParentID;
  G4String      fParticleName;
  G4double      fParticleCharge;
  G4int         fPDGEncoding;
  G4ThreeVector fInitialMomentum;
  G4ThreeVector fInitialMomentumDirection;
  G4String      fVolumeName;
  G4double      fEkin;
  G4double      fSteplength;
  ArDM_TrajectoryPointContainer* fPositionRecord;

};

#endif
