#ifndef ARDM_SENSITIVEBOTTOMPMT
#define ARDM_SENSITIVEBOTTOMPMT

#include "ArDM_SensitivePMT.hh"

class ArDM_SensitiveBottomPMT : public ArDM_SensitivePMT {

public:
  ArDM_SensitiveBottomPMT(G4int,G4String,G4String,G4double,G4int);
  ~ArDM_SensitiveBottomPMT();
  static ArDM_SensitiveBottomPMT* getInstance();

  void Initialize(G4HCofThisEvent*);
  void EndOfEvent(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);

private:
  static ArDM_SensitiveBottomPMT* fInstance;

};

#endif
