#ifndef _SENSITIVEPMT_
#define _SENSITIVEPMT_ 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "ArDM_PMTHit.hh"
#include "G4VSensitiveDetector.hh"
#include "ArDM_Digitizer.hh"

#include "vector"
using namespace std;


class ArDM_SensitivePMT : public G4VSensitiveDetector{

protected:
  G4int    fnpmt;
  G4String fDetName;
  G4String fHitsColName;
  G4String fTop_or_bottom;
  
  G4ThreeVector fNormal;
  vector<G4ThreeVector> fPMTvector;

  ArDM_PMTHitsCollection* fHitsCol;

  static const G4double fTimeSample;
  static const G4int    fNTimeSamples;

  G4double detectionProb(G4double phi);
  G4double detectionProb(G4ThreeVector hitpos,G4ThreeVector rPMT);
  G4double detectionProb(G4ThreeVector hitpos,G4int pmtID);

  G4int detected(G4ThreeVector hitpos,G4ThreeVector rPMT);
  G4int detected(G4ThreeVector hitpos,G4int pmtID);

  static ArDM_Digitizer* digitizer;

public:

  ArDM_SensitivePMT(G4int npmt,G4String detName,G4String hitsColName,G4String top_or_bottom);
  ~ArDM_SensitivePMT();

  void Initialize(G4HCofThisEvent* hc);
  G4bool ProcessHits(G4Step* astep, G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* hc);
  
};




#endif //_SENSITIVEPMT_
