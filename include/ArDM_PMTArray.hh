#ifndef ARDM_PMTARRAY
#define ARDM_PMTARRAY

#include "globals.hh"
#include "TArrayF.h"
#include "TArrayI.h"

class ArDM_PMTData {
public:
  ArDM_PMTData(G4int pmtID=0,G4int nTimeSamples=2048);
  ~ArDM_PMTData();

private:
  G4int fpmtID;
  G4int fNTimeSamples;
  TArrayF fData;
  TArrayI fNphotonPMT;
};

class ArDM_PMTArray {
public:
  ArDM_PMTArray(G4int npmt,G4String arrayName="PMTArray",G4int nTimeSample=2048);
  ~ArDM_PMTArray();
private:
  G4int fNpmt;
  G4int fNTimeSample;
  G4String fArrayName;
  ArDM_PMTData* fPMTArray;
};

#endif
