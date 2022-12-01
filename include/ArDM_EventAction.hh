#ifndef ARDM_EVENTACTION
#define ARDM_EVENTACTION

#include "G4UserEventAction.hh"

class G4Event;
class ArDM_Analysis;

class ArDM_EventAction : public G4UserEventAction {

public:
  ArDM_EventAction();
  ~ArDM_EventAction();

  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

private:
  ArDM_Analysis* ana;

};

#endif
