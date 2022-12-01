#ifndef ARDM_ANALYSIS
#define ARDM_ANALYSIS

#include "TH1D.h"
#include "globals.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "ArDM_AddOns.hh"
#include "ArDM_PMTArray.hh"

#include <vector>

using namespace std;
class TTree;
class TFile;

class ArDM_Analysis {

public:
  ~ArDM_Analysis();
  ArDM_Analysis();
  static ArDM_Analysis* getInstance();

  void BeginOfRun();
  void EndOfRun();
  void BeginOfEvent();
  void EndOfEvent();
  void bookTree();
  void SetOutputFile(G4String name);

  G4String filename;
  G4double fReflReflector;
  G4double fReflPMTSupportTop;
  G4double fReflPMTSupportBottom;
  G4double fReflSideReflectorUpper;
  G4double fReflSideReflectorLower;
  G4int fDetected;
  G4int fNphotonDetectedTop;
  G4int fNphotonDetectedBottom;
  G4int fTotNphotonPMTBottom[NBOTTOMPMT];
  G4int fTotNphotonPMTTop[NTOPPMT];
  G4double fEkinInitial;
  G4double fErecoil;
  G4int NElasticScatterings;
  G4double neutronPosX;
  G4double neutronPosY;
  G4double neutronPosZ;
  G4float ScatterPosX;
  G4float ScatterPosY;
  G4float ScatterPosZ;
  TH1D*  fSpecHist;

  G4float LArVol;
  G4float SStankVol;
  G4float PolyVol;

private:
  static ArDM_Analysis* fInstance;

  TFile* fOutputFile;
  TTree* fAnaTree;

  G4double fTimeSampleTop;
  G4double fTimeSampleBottom;

  G4int fNevents;
  G4int fNbottomPMT;
  G4int fNtopPMT;
  G4int fNtimeSampleBottom;
  G4int fNtimeSampleTop;

  vector<G4int> fNphotonPMTBottom;
  vector<G4int> fNphotonPMTTop;

  G4int fTotNscintPhoton;
  vector<G4int> fNscintPhoton;
  G4double fTotEdep;
};

#endif
