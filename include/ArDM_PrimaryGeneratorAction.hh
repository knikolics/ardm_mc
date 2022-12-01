#ifndef PRIMARY_GENERATORACTION
#define PRIMARY_GENERATORACTION

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include <string>

//class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class TH1D;
class G4ParticleDefinition;

class ArDM_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:
//  ArDM_PrimaryGeneratorAction(char* infile, char* process);
//  ArDM_PrimaryGeneratorAction(char* particleType);
  ArDM_PrimaryGeneratorAction();
  ~ArDM_PrimaryGeneratorAction();
  //static ArDM_PrimaryGeneratorAction* getInstance(char* infile, char* process);
  static ArDM_PrimaryGeneratorAction* getInstance();
//  static ArDM_PrimaryGeneratorAction* getInstance(char* particleType, char* infile, char* process);
  G4ParticleDefinition* particle;
  
  void GeneratePrimaries(G4Event* event);
  char geo[50];
//  char spectrumOnOff[2];
  char infile[50];
  char process[10];
  G4double energy;

private:
  static ArDM_PrimaryGeneratorAction* fInstance;

  //G4ParticleGun* fParticleSource;
  G4GeneralParticleSource* fParticleSource;
  G4double fParticleEnergy;
  G4double neutronEnergy;

  TH1D *h;
  G4double NeutronEnergy();
  G4ThreeVector getRandomPolarization(G4ThreeVector particleMomentumDirection);
  G4ThreeVector position;
  G4ThreeVector SetGunPosition(char* volume);
};

#endif
