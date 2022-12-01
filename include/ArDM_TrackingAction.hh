#ifndef ARDM_TRACKINGACTION
#define ARDM_TRACKINGACTION

#include "G4UserTrackingAction.hh"

class G4Track;

class ArDM_TrackingAction : public G4UserTrackingAction {

public:
  ArDM_TrackingAction();
  ~ArDM_TrackingAction();

  void PreUserTrackingAction(const G4Track*);

};

#endif
