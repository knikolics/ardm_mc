#ifndef ARDM_RUNACTION
#define ARDM_RUNACTION

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4GeometryManager.hh"

class ArDM_Analysis;

class ArDM_RunAction : public G4UserRunAction {

public:
  ArDM_RunAction();
  ~ArDM_RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:
  ArDM_Analysis* fAna;
  G4GeometryManager* geoManager;
};

#endif
