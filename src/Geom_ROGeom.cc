#include "Geom_ROGeom.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "algorithm"
#include "G4UnionSolid.hh"

G4LogicalVolume* constructTube(G4Material* mat,G4double innerR,G4double outerR,G4double halfheight,
			       G4double startPhi,G4double deltaPhi,const char* name,
			       int attribute,G4Colour color) {
  
  G4Tubs*          tubSolid = new G4Tubs(name,innerR,outerR,halfheight,startPhi,deltaPhi);
  G4LogicalVolume* tubLog   = new G4LogicalVolume(tubSolid,mat,name);

  G4VisAttributes* att = new G4VisAttributes(attribute);
  att->SetColour(color);
  att->SetForceAuxEdgeVisible(true);
  tubLog->SetVisAttributes(att);
  
  return tubLog;
}




G4VPhysicalVolume* constructWorld(G4Material* fWorldMat,G4double world_half_size) {
  
  //defining world
  G4Box* worldSolid = new G4Box("world",world_half_size,world_half_size,world_half_size);
  G4LogicalVolume* fWorldLog  = new G4LogicalVolume(worldSolid,fWorldMat,"world");

  G4VisAttributes* worldAtt   = new G4VisAttributes(false); //visibility = false
  fWorldLog->SetVisAttributes(worldAtt);
  G4VPhysicalVolume* fWorldPhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWorldLog,"world",0,false,0);
  return fWorldPhys;
}




G4VPhysicalVolume* constructTank(G4VPhysicalVolume* fMotherPhys,G4Material* fTankMat,G4ThreeVector fTankPos,int attribute) {

  if(!fMotherPhys) {
    G4cout<<"in Geom_ROGeom.cc, constructTank(): fMotherPhys = 0. exit."<<G4endl;
    return NULL;
  }

  //add tank to world
  G4double    tank_innerR      = TANK_INNER_RADIUS;
  G4double    tank_outerR      = TANK_OUTER_RADIUS;
  G4double    tank_half_height = TANK_HALF_HEIGHT;

  G4Tubs*          tankSolid = new G4Tubs("tank",tank_innerR,tank_outerR,tank_half_height,0.*deg,360.*deg);
  G4LogicalVolume* fTankLog  = new G4LogicalVolume(tankSolid,fTankMat,"tank");


  G4VisAttributes* tankAtt = new G4VisAttributes(attribute);
  // tankAtt->SetColour(.7,.4,.1); //brown
  tankAtt->SetColour(1.,1.,1.); //white
  tankAtt->SetForceAuxEdgeVisible(true);
  fTankLog->SetVisAttributes(tankAtt);


  //put tank into world
  G4VPhysicalVolume* fTankPhys = new G4PVPlacement(0,fTankPos,"tank",fTankLog,fMotherPhys,false,0);

  return fTankPhys;
}



G4VPhysicalVolume* constructDetector(G4VPhysicalVolume* fMotherPhys,G4Material* fDetMat,G4ThreeVector fDetPos,int attribute) {

  if(!fMotherPhys) {
    G4cout<<"in Geom_ROGeom.cc, constructDetector(): fMotherPhys = 0. exit."<<G4endl;
    return NULL;
  }

  //add detector (active volume) into tank
  G4double det_innerR = 0;
  G4double det_outerR = TANK_INNER_RADIUS;
  G4double det_half_height = TANK_HALF_HEIGHT;
  
  G4Tubs* detSolid          = new G4Tubs("detector",det_innerR,det_outerR,det_half_height,0.*deg,360.*deg);
  G4LogicalVolume* fDetLog  = new G4LogicalVolume(detSolid,fDetMat,"detector");

  G4VisAttributes* detAtt = new G4VisAttributes(attribute);
  detAtt->SetColour(1.0,1.0,0.0); //yellow
  detAtt->SetForceAuxEdgeVisible(true);
  fDetLog->SetVisAttributes(detAtt);

  G4VPhysicalVolume* fDetPhys = new G4PVPlacement(0,fDetPos,"detector",fDetLog,fMotherPhys,false,0);
   
  return fDetPhys;
}


G4LogicalVolume* constructPMT(G4Material* fPMTMat, const char* bottom_or_top, int attribute,
			      G4double pmt_inner_radius, G4double pmt_outer_radius,const char* name) {

  if(strcmp(bottom_or_top,"top") && strcmp(bottom_or_top,"bottom")) {
    G4cout<<"in Geom_ROGeom.cc, constructPMT(): argument bottom_or_tom is neither \"bottom\" nor \"top\". exit."
	  <<G4endl;
    return NULL;    
  }

  G4int topBottom;
  if(!strcmp(bottom_or_top,"bottom"))    topBottom = 0;
  else if (!strcmp(bottom_or_top,"top")) topBottom = 1;
  G4double startTheta = topBottom?(180.*deg - PMT_ACTIVE_RANGE):0;
  G4double dTheta     = PMT_ACTIVE_RANGE;

  G4Sphere* pmtSolid = new G4Sphere("pmt",pmt_inner_radius,pmt_outer_radius,0.*deg,360.*deg,startTheta,dTheta);
  G4LogicalVolume* fPMTLog = new G4LogicalVolume(pmtSolid,fPMTMat,name);

  G4VisAttributes* pmtAtt = new G4VisAttributes(attribute);
  pmtAtt->SetColour(0.0,1.0,0.0); //green
  pmtAtt->SetForceAuxEdgeVisible(true);
  fPMTLog->SetVisAttributes(pmtAtt);
  return fPMTLog;
}



vector<G4VPhysicalVolume*> placePMT(G4LogicalVolume* fPMTLog,G4VPhysicalVolume* fMotherPhys,
				    vector<G4ThreeVector> rPMT,const char* name) {

  vector<G4VPhysicalVolume*> physVolVec;
  for(unsigned int i=0;i<rPMT.size();i++) 
    physVolVec.push_back((G4VPhysicalVolume*)(new G4PVPlacement(0,rPMT[i],name,fPMTLog,fMotherPhys,false,0)));
  return physVolVec;
}



G4VPhysicalVolume* constructLArCol(G4VPhysicalVolume* fMotherPhys,G4Material* fMat,G4ThreeVector fPos, int attribute) {

  
  if(!fMotherPhys) {
    G4cout<<"in Geom_ROGeom.cc, constructLArCol(): fMotherPhys = 0. exit."<<G4endl;
    return NULL;
  }

  //add detector (active volume) into tank
  G4double det_innerR      = 0;
  G4double det_outerR      = TANK_INNER_RADIUS;
  G4double det_half_height = LAR_COLUMN_HALF_HEIGHT;
 
  G4Tubs*          detSolid   = new G4Tubs("LArCol",det_innerR,det_outerR,det_half_height,0.*deg,360.*deg);
  G4LogicalVolume* fLArColLog = new G4LogicalVolume(detSolid,fMat,"LArCol");
  
  G4VisAttributes* detAtt = new G4VisAttributes(attribute);
  detAtt->SetColour(1.0,1.0,0.0); //yellow
  detAtt->SetForceAuxEdgeVisible(true);
  fLArColLog->SetVisAttributes(detAtt);

  G4VPhysicalVolume* fLArColPhys = new G4PVPlacement(0,fPos,"LArCol",fLArColLog,fMotherPhys,false,0);

  return fLArColPhys;
}




G4VPhysicalVolume* constructGArCol(G4VPhysicalVolume* fMotherPhys,G4Material* fMat,G4ThreeVector fPos, int attribute) {
  if(!fMotherPhys) {
    G4cout<<"in Geom_ROGeom.cc, constructLArCol(): fMotherPhys = 0. exit."<<G4endl;
    return NULL;
  }

  G4double det_innerR = 0;
  G4double det_outerR = TANK_INNER_RADIUS;
  G4double det_half_height = GAR_COLUMN_HALF_HEIGHT;

  G4Tubs* detSolid            = new G4Tubs("GArCol",det_innerR,det_outerR,det_half_height,0.*deg,360.*deg);
  G4LogicalVolume* fGArColLog = new G4LogicalVolume(detSolid,fMat,"GArCol");

  G4VisAttributes* detAtt = new G4VisAttributes(attribute);
  detAtt->SetColour(1,1,1);
  detAtt->SetForceAuxEdgeVisible(true);
  fGArColLog->SetVisAttributes(detAtt);

  G4VPhysicalVolume* fGArColPhys = new G4PVPlacement(0,fPos,"GArCol",fGArColLog,fMotherPhys,false,0);

  return fGArColPhys;
}

G4LogicalVolume* constructHemisphere(G4Material* fHemisphereMat, const char* bottom_or_top, 
				     int attribute,G4double inner_radius,G4double outer_radius,const char* name) {

  if(strcmp(bottom_or_top,"top") && strcmp(bottom_or_top,"bottom")) {
    G4cout<<"in Geom_ROGeom.cc, constructHemisphere(): argument bottom_or_top is neither \"bottom\" nor \"top\". exit."
	  <<G4endl;
    return NULL;
  }

  G4int topBottom;
  if(!strcmp(bottom_or_top,"bottom"))    topBottom = 1;
  else if (!strcmp(bottom_or_top,"top")) topBottom = 0;
  G4double startTheta = topBottom?90*deg:0*deg;
  G4double dTheta     = 90*deg;

  G4Sphere* hemisphereSolid = new G4Sphere("hemisphere",inner_radius,outer_radius,
					   0.*deg,360.*deg,startTheta,dTheta);

  G4LogicalVolume* fHemisphereLog = new G4LogicalVolume(hemisphereSolid,fHemisphereMat,name);

  G4VisAttributes* hemisphereAtt = new G4VisAttributes(attribute);
  hemisphereAtt->SetColour(1,0,0);
  hemisphereAtt->SetForceAuxEdgeVisible(true);
  fHemisphereLog->SetVisAttributes(hemisphereAtt);
  return fHemisphereLog;
}

G4LogicalVolume* constructReflectorHolderRing(G4Material* fMat,int attribute,
					      G4double innerR,G4double outerR,G4double halfheight,
					      const char* name,G4Colour color) {

  G4Tubs*          holderRingSolid = new G4Tubs("name",innerR,outerR,halfheight,0.*deg,360.*deg);
  G4LogicalVolume* holderRingLog   = new G4LogicalVolume(holderRingSolid,fMat,name);

  G4VisAttributes* holderRingAtt = new G4VisAttributes(attribute);
  holderRingAtt->SetColour(color);
  holderRingAtt->SetForceAuxEdgeVisible(true);
  holderRingLog->SetVisAttributes(holderRingAtt);

  return holderRingLog;
}

G4double findPhi(G4ThreeVector emissionPos,G4ThreeVector colcenter,G4double radius) {
  G4ThreeVector connectLine =colcenter - emissionPos;
  G4ThreeVector connectLine2D = connectLine.perpPart();

  const int n = 4;
  G4ThreeVector pos[n];

  pos[0] = radius*connectLine2D.unit();
  pos[1] = -pos[0];

  G4ThreeVector connectLine_orthogonal = connectLine2D.orthogonal().perpPart().unit();
  pos[2] = radius*connectLine_orthogonal;
  pos[3] = -pos[2];

  G4double z = colcenter.z();
  G4double* phi = new G4double[n];

  for(int i=0;i<n;i++) {
    pos[i].setZ(z);
    phi[i] = (pos[i] - emissionPos).angle(connectLine);
  }

  G4double maxphi = *max_element(phi,phi+n);
  return maxphi;
}
