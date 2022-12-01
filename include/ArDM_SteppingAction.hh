#ifndef ARDM_STEPPINGACTION
#define ARDM_STEPPINGACTION

#include "G4UserSteppingAction.hh"

class G4Step;
class ArDM_Analysis;

class ArDM_SteppingAction : public G4UserSteppingAction {
public:
  ArDM_SteppingAction();
  ~ArDM_SteppingAction();

  void UserSteppingAction(const G4Step*);
  ArDM_Analysis* ana;
};

#endif
