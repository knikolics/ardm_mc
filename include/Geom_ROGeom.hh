#ifndef GEOM_ROGEOM
#define GEOM_ROGEOM

#include "ArDM_AddOns.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Colour.hh"

extern "C++" G4LogicalVolume* constructTube(G4Material* mat,G4double innerR,G4double outerR,G4double halfheight,
    G4double startPhi,G4double deltaPhi,const char* name="", int attribute=0,G4Colour color=G4Colour::Grey());
extern "C++" G4VPhysicalVolume* constructWorld(G4Material* fWorldMat,G4double world_half_size=WORLD_HALF_SIZE);
extern "C++" G4VPhysicalVolume* constructTank(G4VPhysicalVolume* fMotherLog,G4Material* fTankMat, G4ThreeVector fTankPos, int attribute=0);
extern "C++" G4VPhysicalVolume* constructDetector(G4VPhysicalVolume* fMotherLog,G4Material* fDetMat, G4ThreeVector fDetPos, int attribute=0);
extern "C++" G4LogicalVolume*   constructPMT(G4Material* fPMTMat, const char* topBottomPMT="bottom", int attribute=0,
   G4double pmt_inner_radius = PMT_INNER_RADIUS, G4double pmt_outer_radius = PMT_OUTER_RADIUS, const char* name="PMT");
extern "C++" vector<G4VPhysicalVolume*> placePMT(G4LogicalVolume* fPMTLog,G4VPhysicalVolume* fMotherLog, vector<G4ThreeVector> rPMT,const char* name);
extern "C++" G4VPhysicalVolume* constructLArCol(G4VPhysicalVolume* fMotherLog,G4Material* fMat, G4ThreeVector fPos=G4ThreeVector(0.,0.,0.), int attribute=0);
extern "C++" G4VPhysicalVolume* constructGArCol(G4VPhysicalVolume* fMotherLog,G4Material* fMat, G4ThreeVector fPos=G4ThreeVector(0.,0.,0.), int attribute=0);
extern "C++" G4LogicalVolume* constructHemisphere(G4Material* fMat, const char* topBottom="bottom",
  int attribute=0,G4double inner_radius = PMT_INNER_RADIUS, G4double outer_radius = PMT_OUTER_RADIUS,const char* name="bottom");
extern "C++" G4double findPhi(G4ThreeVector emissionPos,G4ThreeVector colcenter,G4double radius);
#endif
