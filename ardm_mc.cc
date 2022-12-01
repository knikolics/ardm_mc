#include "G4RunManager.hh"
#include "ArDM_DetectorConstruction.hh"
#include "ArDM_PhysicsList.hh"
#include "ArDM_Analysis.hh"
#include "ArDM_PrimaryGeneratorAction.hh"
#include "ArDM_RunAction.hh"
#include "ArDM_EventAction.hh"
#include "ArDM_TrackingAction.hh"
#include "ArDM_SteppingAction.hh"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//#ifdef VIS
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
//#endif

#include "TStopwatch.h"

using namespace std;

int main(int argc,char* argv[]){

  int nevents = 0;
  //char particleType[100] = "";
  char inputFile[100] = "";

if(argv[1]==NULL){
  //G4cout << "Usage: ardm_mc N(particle) particle(e.g. 'e-')" << G4endl;
  G4cout << "Usage: ardm_mc particle(e.g. 'e-')" << G4endl;
  return 0;
  }
 else {
//  nevents = atoi(argv[1]);
  sprintf(inputFile,argv[1]);

  TStopwatch timer;
  timer.Start();

  G4String outputname = (string)inputFile + ".root";
  G4cout << "Creating output file: " << outputname << G4endl;
  ArDM_Analysis::getInstance()->SetOutputFile(outputname);

  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new ArDM_PhysicsList);
  runManager->SetUserInitialization(new ArDM_DetectorConstruction);
  //runManager->SetUserAction(new ArDM_PrimaryGeneratorAction(particleType));
  runManager->SetUserAction(new ArDM_PrimaryGeneratorAction);
  runManager->SetUserAction(new ArDM_RunAction);
  runManager->SetUserAction(new ArDM_EventAction);
  runManager->SetUserAction(new ArDM_TrackingAction);
  runManager->SetUserAction(new ArDM_SteppingAction);
  runManager->Initialize();

#ifdef VIS
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  delete visManager;
#endif
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
//  UImanager->ApplyCommand("/run/verbose 1");
//  UImanager->ApplyCommand("/event/verbose 1");
//  UImanager->ApplyCommand("/tracking/verbose 1");
//  UImanager->ApplyCommand("/control/execute vis.mac");
  UImanager->ApplyCommand("/control/execute myGPS.mac");
  runManager->BeamOn(nevents);
  //delete UImanager;
// delete visManager;
//#else
//  runManager->BeamOn(nevents);
//#endif

  delete runManager;
  timer.Stop();
  G4cout << "CPU time " << timer.CpuTime() << "\t real time " << timer.RealTime() << G4endl;
  return 0;
 }
}

