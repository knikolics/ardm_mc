#ifndef ARDM_SENSITIVETOPPMT
#define ARDM_SENSITIVETOPPMT

#include "ArDM_SensitivePMT.hh"

class ArDM_SensitiveTopPMT : public ArDM_SensitivePMT {

public:
  ArDM_SensitiveTopPMT(G4int,G4String,G4String,G4double,G4int);
  ~ArDM_SensitiveTopPMT();
  static ArDM_SensitiveTopPMT* getInstance();

  void Initialize(G4HCofThisEvent*);
  void EndOfEvent(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);

private:
  static ArDM_SensitiveTopPMT* fInstance;
  ArDM_PMTHitsCollection* fHitsCol;

};

#endif
