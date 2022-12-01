#ifndef ARDM_PHYSICSLIST
#define ARDM_PHYSICSLIST

#include "G4VUserPhysicsList.hh"

class G4ProcessManager;

class ArDM_PhysicsList : public G4VUserPhysicsList {

public:
  ArDM_PhysicsList();
  ~ArDM_PhysicsList();

  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  void AddScintillation();

private:
  void electronPhysics(G4ProcessManager*);
  void positronPhysics(G4ProcessManager*);
  void muonPhysics(G4ProcessManager*);
  void opPhotonPhysics(G4ProcessManager*);
  void gammaPhysics(G4ProcessManager*);
  void protonPhysics(G4ProcessManager*);
  void neutronPhysics(G4ProcessManager*);
  void ionPhysics(G4ProcessManager*);
  void hadronPhysics(G4ProcessManager*);
  void scintillationPhysics(G4ProcessManager*);
};

#endif
