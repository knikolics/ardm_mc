#include "ArDM_DetectorConstruction.hh"
//#include "ArDM_test.cc"
#include "ArDM_SensitivePMT.hh"
#include "G4SDManager.hh"

//#include "ArDM_Analysis.hh"

#include "preparation.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4BooleanSolid.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4Element.hh"
#include "G4RotationMatrix.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Colour.hh"

//#include "Geom_ROGeom.hh"


#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4Torus.hh"


#include "fstream"
using namespace std;


ArDM_DetectorConstruction::ArDM_DetectorConstruction() 
  : G4VUserDetectorConstruction() {

  fNist = G4NistManager::Instance();
  fVacuum = new G4Material("vacuum",1.,1.01*g/mole,universe_mean_density,kStateGas,0.1*kelvin,1.e-19*pascal);  
  fWorldMat = fNist->FindOrBuildMaterial("G4_AIR"); fWorldMat->SetName("Air");
  fTankMat  = getTankMat(); //stainless steel
  fLidMat   = fTankMat;
  fDetMat   = fNist->FindOrBuildMaterial("G4_lAr");
  fLAr      = fNist->FindOrBuildMaterial("G4_lAr");
  fGAr      = fNist->FindOrBuildMaterial("G4_Ar");
  fPolyethylene  = getPolyethylene();
  fTeflon   = getTeflon();
  setMatPropTab_LAr(); 
  setMatPropTab_GAr(); 
  setMatPropTab_WLS();
  setMatPropTab_PMTMat();

  fWLSSupportMat = fTeflon;
  fBottomSideReflectorMat = fTeflon;
  fTopSideReflectorMat = fTeflon;
  fPMTSupportMat = fTeflon;
  fCathodeGridMat= fTankMat;

  fPMTCathodeMat = fPMTMat; //pmt cathode material is irrelevant ! since its physical properties will not be used in this MC !

  fPMTElectrodeMat = fTankMat; //stainless steel
  
  fPMTBaseMat = getFR4(); //fPolyethylene; //arbitrary for the time being
  //fPMTBaseMat = fNist->FindOrBuildMaterial("G10");
  //cout<<"fPMTBaseMat "<<fPMTBaseMat<<endl; getchar();

  fHVrMat = getHVrMat();
  fPerlite = getPerlite();
  fBorotron = getBorotron();

  
}




ArDM_DetectorConstruction::~ArDM_DetectorConstruction(){
  //deleting pointers 
  
  if(!fWorldPhys)        delete fWorldPhys;
  if(!fTankPhys)         delete fTankPhys;
  if(!fDetPhys)          delete fDetPhys;
  if(!fBottomPMT_ROPhys) delete fBottomPMT_ROPhys;
  if(!fTopPMT_ROPhys)    delete fTopPMT_ROPhys;

  if(fBottomPMTArrayPhys.size())  for(int i=0;i<fBottomPMTArrayPhys.size();i++) delete fBottomPMTArrayPhys[i];
  if(fTopPMTArrayPhys.size())     for(int i=0;i<fTopPMTArrayPhys.size();   i++) delete fTopPMTArrayPhys[i];
}



G4VPhysicalVolume* ArDM_DetectorConstruction::Construct(){

  //construct world
  G4Box* worldSolid = new G4Box("world",WORLD_HALF_SIZE,WORLD_HALF_SIZE,WORLD_HALF_SIZE);
  G4LogicalVolume* fWorldLog  = new G4LogicalVolume(worldSolid,fWorldMat,"world");

  G4VisAttributes* worldAtt   = new G4VisAttributes(false); //visibility = false
  fWorldLog->SetVisAttributes(worldAtt);
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWorldLog,"world",0,false,0);

  //add detector components
  
  addTank(); 
  addLArColumn();
  addGArColumn();
  addWLS();
  addWLSSupport();
  addPMT("bottom");    
  addPMT("top");

#if TURN_ON_SIDE_REFLECTOR
  addBottomSideReflector();
  addTopSideReflector();
#endif //TURN_ON_SIDE_REFLECTOR

  addPMTSupport("bottom");
  addPMTSupport("top");

#if TURN_ON_CATHODE_GRID
  addCathodeGrid();
#endif //TURN_ON_CATHODE_GRID



#if TURN_ON_PROTECTION_GRID
  addProtectionGrid();
#endif //TURN_ON_PROTECTION_GRID


#if NEUTRON_SHIELD
  addNeutronShield();

#endif //NEUTRON_SHIELD


  //addSource();
#if GAS_TEST_2013_MARCH
#if GAS_TEST_2013_MARCH_SOURCE_HOLDER
  addSourceHolder();

#if GAS_TEST_2013_MARCH_SOURCE_COATING
  addSourceCoating();
#endif //GAS_TEST_2013_MARCH_SOURCE_COATING

#endif //GAS_TEST_2013_MARCH_SOURCE_HOLDER

#endif //GAS_TEST_2013_MARCH



#if TEST_BRANCH24
  //neutron background
  addHVr();
  addPerlite();
#endif //TEST_BRANCH24


  build_surfaces();

//   //addEfield();


  return fWorldPhys;
}



void ArDM_DetectorConstruction::addTopFlange(){

  
  if(!fWorldPhys){
    G4cout<<"in ArDM_DetectorConstruction::addTopFlange(), fWorldPhys = 0. exit."<<G4endl;
    return;
  }



  //top flange = flange with some holes + stainless steel pillars
  //thickness of the flange = 4 cm
  //simplify the geometry by a simple stainless steel cylinder,
  //instead of using 4cm as the cylinder height, increase it a bit to compensate for the omitted pillars
  //--> increase by how much ??

  G4double innerR = TOP_FLANGE_INNER_RADIUS;
  G4double outerR = TOP_FLANGE_OUTER_RADIUS;
  G4double halfz  = TOP_FLANGE_HALF_HEIGHT_EFFECTIVE;
  
  G4Tubs* solid = new G4Tubs("topFlange",innerR,outerR,halfz,0*deg,360*deg);
  G4LogicalVolume* log = new G4LogicalVolume(solid,fTankMat,"topFlange");
  
  G4VisAttributes* lidAtt = new G4VisAttributes(true);
  lidAtt->SetColour(0.0,1.0,1.0); //cyan
  lidAtt->SetForceAuxEdgeVisible(true);
  log->SetVisAttributes(lidAtt);

  G4ThreeVector pos(TOP_FLANGE_POS_X,TOP_FLANGE_POS_Y,TOP_FLANGE_POS_Z);

  fTopFlangePhys = new G4PVPlacement(0,pos,"topFlange",log,fWorldPhys,false,0);
  
  return;
}








void ArDM_DetectorConstruction::addHVr(){

  //there are 2 x 27 high voltage resistors
  //in stead of simulating all these HVr,
  //we simulate them 'collectively' by grouping them into 2 groups,
  //we approximate each group with a "bar" running along the linear sector (along z-axis)
  //of the field shaper / main reflector .


  G4double innerR = HV_RESISTOR_BAR_INNER_RADIUS;
  G4double outerR = HV_RESISTOR_BAR_OUTER_RADIUS;
  G4double halfz  = HV_RESISTOR_BAR_HALF_HEIGHT;


  G4Tubs* solid = new G4Tubs("HVrBar",innerR,outerR,halfz,0*deg,360*deg);
  
  G4LogicalVolume* log = new G4LogicalVolume(solid,fHVrMat,"HVr");
  
  G4VPhysicalVolume* fMotherPhys;
  G4double z1,z2;
  if(fLArColPhys){ fMotherPhys = fLArColPhys; z1 = HV_RESISTOR_BAR_1_POS_Z_LARCOL; z2 = HV_RESISTOR_BAR_2_POS_Z_LARCOL; }
  else{            fMotherPhys = fWorldPhys;  z1 = HV_RESISTOR_BAR_1_POS_Z;        z2 = HV_RESISTOR_BAR_2_POS_Z;        }
  

  G4ThreeVector pos1(HV_RESISTOR_BAR_1_POS_X,HV_RESISTOR_BAR_1_POS_Y,z1);
  G4ThreeVector pos2(HV_RESISTOR_BAR_2_POS_X,HV_RESISTOR_BAR_2_POS_Y,z2);

  fHVrBar1Phys = new G4PVPlacement(0,pos1,"HVr",log,fMotherPhys,false,0);
  fHVrBar2Phys = new G4PVPlacement(0,pos2,"HVr",log,fMotherPhys,false,0);

  return;
}









void ArDM_DetectorConstruction::addTank(){

  if(!fWorldPhys){
    G4cout<<"in ArDM_DetectorConstruction::addTank(), fWorldPhys = 0. exit."<<G4endl;
    return;
  }
  


  //// try to simulate the stainless steel dewar more precisely
  //// the dewar is a "union" of :
  //// 1. a cylinder and
  //// 2. 3 spherical parts at the bottom
  ////    i.  one big arc at the (very) bottom of the dewar with radius R = 1010 mm --> here, called R1010_ARC
  ////    ii. 2 connecting arcs (each on one side) with radius r = 101 mm, which connect the cylinder with the R1010-arc --> here, called R101_ARC
  ////
  ////
  //// moreover, the wall of the dewar has the form of a sandwich of many alternating layers
  //// from inside to outside it's
  //// 1. 6 mm SS
  //// 2. 5 mm LAr
  //// 3. 3 mm SS
  //// 4. 10 mm LAr
  //// 5. 6 mm SS
  //// 6. 25 mm vacuum
  //// 7. 5 mm SS
  ////
  ////--> in order to simplify things :
  //// 1. combine all the stainless steel (SS) layers to one single one. --> total thickness 20 mm
  //// 2. add the LAr layers to the LAr volume inside the tank, which increases the volume of the LAr inside the tank by a negligible amount,
  ////    but still we can account for this by slightly increase the inner radius of the dewar from 500 mm to (500 + 5 + 10) mm
  //// 3. ignore the vacuum layer.
  ////
  //// this doesn't make any big difference for neutrons emitted from the dewar wall. since the LAr layers are too thin, 
  //// and the vaccum layer doesn't have any effect on neutrons.
  



  G4Tubs* cylinderSolid = new G4Tubs("tankCylinder",TANK_CYLINDER_INNER_RADIUS,TANK_CYLINDER_OUTER_RADIUS,TANK_CYLINDER_HALF_HEIGHT,0*deg,360*deg);

  G4Sphere* r1010_arc_solid = new G4Sphere("tankR1010Arc",TANK_BTM_PART_R1010_ARC_INNER_RADIUS,TANK_BTM_PART_R1010_ARC_OUTER_RADIUS,
					   0*deg,360*deg,180*deg - TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2,TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2);
					   //0*deg,360*deg,180*deg - TANK_BTM_PART_R1010_ARC_OPENING_ANGLE,TANK_BTM_PART_R1010_ARC_OPENING_ANGLE);


  // the connecting part (r101-arc) is complicated !!
  // it's a part of a torus.
  //
  // how to build it logically ?
  // 1. you have an ("solid") arc of a ring with certain opening angle theta, and innerRadius r1, outerRadius r2
  // 2. now you rotate this arc around an axis, which is perpendicular to the normal vector of the arc, you'll get the connecting part.
  // (yea, it's also a ring. or rather a part of a torus (which you'd get if you rotate the whole ring, and not just an arc of it, around a certain axis))

  //
  // in GEANT4, you can generate a whole torus, or a part of a torus with some certain swept angle.
  // but the cross section of the torus is always a ring, and not an arc of a ring. 
  // GEANT4 doesn't have any standard class for that shape.
  //
  //--> now, in order to get the connecting part, we have to cut the full torus with some boxes.


  G4double sweptR = (TANK_CYLINDER_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS);
  if(0){
    cout<<"TANK_CYLINDER_INNER_RADIUS "<<TANK_CYLINDER_INNER_RADIUS
	<<"\nTANK_BTM_PART_R101_ARC_INNER_RADIUS "<<TANK_BTM_PART_R101_ARC_INNER_RADIUS
	<<"\nsweptR "<<sweptR
	<<endl;
    getchar();
  }

  G4Torus* torus = new G4Torus("torus",TANK_BTM_PART_R101_ARC_INNER_RADIUS,TANK_BTM_PART_R101_ARC_OUTER_RADIUS,sweptR,0*deg,360*deg);
  //G4Torus* torus = new G4Torus("torus",0,TANK_BTM_PART_R101_ARC_OUTER_RADIUS,sweptR,0*deg,360*deg);
  
  //now define boxes to cut the torus
  G4Box* upperBox = new G4Box("upperBox",TANK_CYLINDER_OUTER_RADIUS,TANK_CYLINDER_OUTER_RADIUS,TANK_BTM_PART_R101_ARC_OUTER_RADIUS/2);

  // G4Cons* cone = new G4Cons("cone",TANK_CYLINDER_INNER_RADIUS-TANK_BTM_PART_R101_ARC_INNER_RADIUS - TANK_BTM_PART_R101_ARC_OUTER_RADIUS,
  // 			    TANK_CYLINDER_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS,
  // 			    TANK_CYLINDER_INNER_RADIUS-TANK_BTM_PART_R101_ARC_INNER_RADIUS - TANK_BTM_PART_R101_ARC_OUTER_RADIUS,
  // 			    TANK_CYLINDER_INNER_RADIUS-TANK_BTM_PART_R101_ARC_INNER_RADIUS + TANK_BTM_PART_R101_ARC_OUTER_RADIUS*sin(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2),
  // 			    TANK_BTM_PART_R101_ARC_OUTER_RADIUS*cos(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2)/2,
  // 			    0,360.*deg);




  G4Cons* cone = new G4Cons("cone",TANK_CYLINDER_INNER_RADIUS-TANK_BTM_PART_R101_ARC_INNER_RADIUS - TANK_BTM_PART_R101_ARC_OUTER_RADIUS,
			    TANK_CYLINDER_INNER_RADIUS-TANK_BTM_PART_R101_ARC_INNER_RADIUS + TANK_BTM_PART_R101_ARC_OUTER_RADIUS*sin(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2),
			    TANK_CYLINDER_INNER_RADIUS-TANK_BTM_PART_R101_ARC_INNER_RADIUS - TANK_BTM_PART_R101_ARC_OUTER_RADIUS,
			    TANK_CYLINDER_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS,
			    TANK_BTM_PART_R101_ARC_OUTER_RADIUS*cos(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2)/2,
			    0,360.*deg);



  G4ThreeVector upperBoxPos(0,0,upperBox->GetZHalfLength());
  G4ThreeVector conePos(0,0,-cone->GetZHalfLength() );
   
  G4SubtractionSolid* r101_arc_solid = new G4SubtractionSolid("tankR101Arc",torus,upperBox,0,upperBoxPos);
  r101_arc_solid = new G4SubtractionSolid("tankR101Arc",r101_arc_solid,cone,0,conePos);
  

  //now that we have all the parts:
  // 1. tank cylinder
  // 2. r101 connecting part
  // 3. r1010 bottom part
  //--> put these things together to get the shape for the dewar

  G4ThreeVector r101_pos (0,0,TANK_BTM_PART_R101_Z_RELATIVE_TO_TANK_CYLINDER);
  G4ThreeVector r1010_pos(0,0,TANK_BTM_PART_R1010_Z_RELATIVE_TO_TANK_CYLINDER);

  G4UnionSolid* tankSolid = new G4UnionSolid("tank",cylinderSolid,r101_arc_solid,0,r101_pos);
  tankSolid = new G4UnionSolid("tank",tankSolid,r1010_arc_solid,0,r1010_pos);
  //G4VSolid* tankSolid = r1010_arc_solid;
  
  G4LogicalVolume* fTankLog = new G4LogicalVolume(tankSolid,fTankMat,"tank");


  G4VisAttributes* tankAtt = new G4VisAttributes(1);
  // tankAtt->SetColour(.7,.4,.1); //brown
  tankAtt->SetColour(1.,1.,1.); //white
  tankAtt->SetForceAuxEdgeVisible(true);
  fTankLog->SetVisAttributes(tankAtt);


  //put tank into world
  G4ThreeVector tankPos(0.,0.,TANK_CYLINDER_POS_Z);
  fTankPhys = new G4PVPlacement(0,tankPos,"tank",fTankLog,fWorldPhys,false,0);

// #if TEST_BRANCH24
//   ArDM_Analysis::getInstance()->SStankVol = fTankPhys->GetLogicalVolume()->GetSolid()->GetCubicVolume();
// #endif //TEST_BRANCH24

  addTopFlange();

  return;
}






void ArDM_DetectorConstruction::addLArColumn(){

#if GAS_TEST_2013_MARCH
  fLAr = fGAr;
#endif //GAS_TEST_2013_MARCH
				       

  //in this function :
  //LArCol = LAr cylinder + LAr btm part
  //LAr cylinder <--> tank_cylinder (or cylinderSolid) in constructApproxTank
  //LAr btm part <--> LAr volume contained in the 'curved' btm part of the tank, corresponding to r101 + r1010 arc together in constructApproxTank
  //
  //the LAr btm part is a union of :
  //1. a part of a torus (<--> r101 arc)
  //2. a part of a cone  (<--> r1010 arc + a part of r101 arc)

  G4double LArCylinder_innerR = APPROX_LARCOL_CYLINDER_INNER_RADIUS;
  G4double LArCylinder_outerR = APPROX_LARCOL_CYLINDER_OUTER_RADIUS;
  G4double LArCylinder_halfz  = APPROX_LARCOL_CYLINDER_HALFHEIGHT;
  
  G4Tubs* LArCylinder_solid = new G4Tubs("LArCylinder",LArCylinder_innerR,LArCylinder_outerR,LArCylinder_halfz,0*deg,360*deg);
  

  //a part of a torus (r101 arc)
  G4double LAr_r101_arc_innerR = 0;
  G4double LAr_r101_arc_outerR = TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  G4double sweptR              = TANK_CYLINDER_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  G4Torus* torus = new G4Torus("LArTorus",LAr_r101_arc_innerR,LAr_r101_arc_outerR,sweptR,0*deg,360*deg);
  
  //now define boxes / cones to cut the torus
  G4Box* upperBox = new G4Box("upperBox",TANK_CYLINDER_INNER_RADIUS,TANK_CYLINDER_INNER_RADIUS,TANK_BTM_PART_R101_ARC_INNER_RADIUS/2);

  //the cone used to cut the torus doesn't have to be this big,
  //but just to simplify things, we can also take a big cone.
  G4double cone_upper_innerR = 0;
  G4double cone_upper_outerR = TANK_CYLINDER_INNER_RADIUS - TANK_BTM_PART_R101_ARC_INNER_RADIUS;
  G4double cone_lower_innerR = 0;
  G4double cone_lower_outerR = TANK_BTM_PART_R1010_ARC_OUTER_RADIUS*tan(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2);
  G4double cone_halfz        = (DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT + TANK_CYLINDER_THICKNESS)/2;

  G4Cons* cone    = new G4Cons("cone",cone_lower_innerR,cone_lower_outerR,cone_upper_innerR,cone_upper_outerR,cone_halfz,0*deg,360*deg);
  
  G4ThreeVector upperBoxPos(0,0,upperBox->GetZHalfLength());
  G4ThreeVector conePos(0,0,-cone->GetZHalfLength() );
  
  G4SubtractionSolid* LAr_r101_solid = new G4SubtractionSolid("LArR101",torus,upperBox,0,upperBoxPos);
  LAr_r101_solid = new G4SubtractionSolid("LArR101",LAr_r101_solid,cone,0,conePos);


  //r101 + r1010 joint part

  //actually innerR = 0 mm, but we will later subtract a cone from this sphere, it's better to approximate 0 mm by 0.001 mm
  G4double sphere_innerR = 0.001* mm; 
  G4double sphere_outerR = TANK_BTM_PART_R1010_ARC_INNER_RADIUS;
  G4Sphere* r1010_r101 = new G4Sphere("sphere",sphere_innerR,sphere_outerR,0,360*deg,
				      180*deg - TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2,TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2);

  
  G4double r1010_r101_cone_lower_innerR = 0;
  G4double r1010_r101_cone_lower_outerR =(TANK_BTM_PART_R1010_ARC_INNER_RADIUS - 
					  TANK_BTM_PART_R101_ARC_INNER_RADIUS)*sin(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2) ;

  G4double r1010_r101_cone_upper_innerR = 0;
  G4double r1010_r101_cone_upper_outerR = sphere_innerR;
  
  G4double r1010_r101_cone_halfz = (TANK_BTM_PART_R1010_ARC_INNER_RADIUS - sphere_innerR - 
				    TANK_BTM_PART_R101_ARC_INNER_RADIUS)*cos(TANK_BTM_PART_R1010_ARC_OPENING_ANGLE/2);
  

  G4Cons* r1010_r101_cone = new G4Cons("cone",r1010_r101_cone_lower_innerR,r1010_r101_cone_lower_outerR,
				       r1010_r101_cone_upper_innerR,r1010_r101_cone_upper_outerR,r1010_r101_cone_halfz,0*deg,360*deg);


  G4ThreeVector r1010_r101_cone_pos(0,0,-(DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER - r1010_r101_cone_halfz - sphere_innerR));
  G4SubtractionSolid* r1010_r101_solid = new G4SubtractionSolid("r1010_r101",r1010_r101,r1010_r101_cone,0,r1010_r101_cone_pos);
  

  //put r101 and r1010_r101 joint part together
  G4ThreeVector r101_rel_r1010_r101_pos(0,0,DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER);
  G4UnionSolid* btm_part_solid = new G4UnionSolid("btmPartSolid",LAr_r101_solid,r1010_r101_solid,0,r101_rel_r1010_r101_pos);

  //now put all pieces together to get the whole LAr volume
  //position vector relative to the center of LAr_cylinder_solid
  // G4ThreeVector LAr_r101_solid_pos(0,0,-APPROX_LARCOL_CYLINDER_HALFHEIGHT);
  // G4ThreeVector r1010_r101_solid_pos(0,0,DISTANCE_R1010_ARC_CENTER_TO_LOWER_EDGE_OF_TANK_CYLINDER - APPROX_LARCOL_CYLINDER_HALFHEIGHT);
  // G4UnionSolid* LArCol_solid = new G4UnionSolid("LArCol_solid",LArCylinder_solid,LAr_r101_solid,0,LAr_r101_solid_pos);
  // LArCol_solid = new G4UnionSolid("LArCol_solid",LArCol_solid,r1010_r101_solid,0,r1010_r101_solid_pos);

  G4ThreeVector btmPart_pos(0,0,-APPROX_LARCOL_CYLINDER_HALFHEIGHT);
  G4UnionSolid* LArCol_solid = new G4UnionSolid("LArCol_solid",LArCylinder_solid,btm_part_solid,0,btmPart_pos);

  G4LogicalVolume* fLArColLog = new G4LogicalVolume(LArCol_solid,fLAr,"LArCol");
  
  G4VisAttributes* detAtt = new G4VisAttributes(1);
  detAtt->SetColour(1.0,1.0,0.0); //yellow
  detAtt->SetForceAuxEdgeVisible(true);
  fLArColLog->SetVisAttributes(detAtt);

  G4double z = APPROX_LARCOL_CYLINDER_POS_Z;
  G4ThreeVector pos(0.,0.,z);
  fLArColPhys = new G4PVPlacement(0,pos,"LArCol",fLArColLog,fWorldPhys,false,0);


// #if TEST_BRANCH24		
//   ArDM_Analysis::getInstance()->LArVol = fLArColPhys->GetLogicalVolume()->GetSolid()->GetCubicVolume();
// #endif //TEST_BRANCH24	

  return;
}




void ArDM_DetectorConstruction::addGArColumn(){

  if(!fWorldPhys){
    G4cout<<"in ArDM_DetectorConstruction::addGArColumn(), fWorldPhys = 0. exit."<<G4endl;
    return;
  }

  //add detector (active volume) into tank
  G4double det_innerR = 0;
  G4double det_outerR = TANK_CYLINDER_INNER_RADIUS;
  G4double det_half_height = GAR_COLUMN_HALF_HEIGHT;
  
  G4Tubs* detSolid            = new G4Tubs("GArCol",det_innerR,det_outerR,det_half_height,0.*deg,360.*deg);
  G4LogicalVolume* fGArColLog = new G4LogicalVolume(detSolid,fGAr,"GArCol");
  
  G4VisAttributes* detAtt = new G4VisAttributes(1);
  detAtt->SetColour(1.0,1.0,1.0); //white
  detAtt->SetForceAuxEdgeVisible(true);
  fGArColLog->SetVisAttributes(detAtt);

  G4double z = GAR_COLUMN_Z;
  G4ThreeVector pos(0.,0.,z);

  fGArColPhys = new G4PVPlacement(0,pos,"GArCol",fGArColLog,fWorldPhys,false,0);

  return;
}



G4VSensitiveDetector* ArDM_DetectorConstruction::getSD(G4String bottom_or_top){


  if(bottom_or_top == "top" || bottom_or_top == "Top" || bottom_or_top == "TOP"                 ) bottom_or_top = "top";
  else if(bottom_or_top == "bottom" || bottom_or_top == "Bottom" || bottom_or_top == "BOTTOM" || 
	  bottom_or_top == "btm"    || bottom_or_top == "BTM"                                   ) bottom_or_top = "btm";
  else return NULL;


  ostringstream detNameStr;
  detNameStr << "PMTarray/"<<bottom_or_top;
  //detNameStr << bottom_or_top <<"Array";

  G4String detName = detNameStr.str();

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* sensDet = SDman->FindSensitiveDetector(detName);

  if(sensDet) return sensDet;

  ostringstream hitsColNameStr;
  hitsColNameStr << bottom_or_top << "PMTHits";
  G4String hitsColName = hitsColNameStr.str();

  int npmt = (bottom_or_top == "top")? NTOPPMT : NBOTTOMPMT;

  return new ArDM_SensitivePMT(npmt,detName,hitsColName,bottom_or_top); 
}



vector<G4VPhysicalVolume*> ArDM_DetectorConstruction::getPMTArrayPhys(G4String bottom_or_top){

  if(bottom_or_top == "bottom"){
    return fBottomPMTArrayPhys;
  }
  else if(bottom_or_top == "top"){
    return fTopPMTArrayPhys;
  }
  else{
    vector<G4VPhysicalVolume*> vec;
    return vec; 
  }
}





vector<G4VPhysicalVolume*> ArDM_DetectorConstruction::placePMT(G4LogicalVolume* fPMTLog,G4VPhysicalVolume* fMotherPhys,
							       vector<G4ThreeVector> rPMT,const char* bottom_or_top){


  G4RotationMatrix* rot = NULL;
  
  if(!strcmp(bottom_or_top,"top")){
    rot = new G4RotationMatrix;
    rot->rotateX(180.*deg);
    rot->rotateY(0.*deg);
    rot->rotateZ(0.*deg);
  }


  vector<G4VPhysicalVolume*> physVolVec;
  ostringstream pmtname;
  for(int i=0;i<rPMT.size();i++){
    pmtname.str("");
    if(!strcmp(bottom_or_top,"top")) pmtname<<fPMTLog->GetName()<<i;
    else if(!strcmp(bottom_or_top,"bottom")) pmtname<<fPMTLog->GetName()<<i+12;
    else pmtname<<fPMTLog->GetName();
    physVolVec.push_back((G4VPhysicalVolume*)(new G4PVPlacement(rot,rPMT[i],pmtname.str().c_str(),fPMTLog,fMotherPhys,false,0)));
  }

  return physVolVec;
}




vector<G4VPhysicalVolume*> ArDM_DetectorConstruction::placePMT(vector<G4LogicalVolume*> fPMTLog,G4VPhysicalVolume* fMotherPhys,
				    vector<G4ThreeVector> rPMT,const char* bottom_or_top){


  G4RotationMatrix* rot = NULL;
  
  if(!strcmp(bottom_or_top,"top")){
    rot = new G4RotationMatrix;
    rot->rotateX(180.*deg);
    rot->rotateY(0.*deg);
    rot->rotateZ(0.*deg);
  }

  vector<G4VPhysicalVolume*> physVolVec;
  ostringstream pmtname;
  for(int i=0;i<rPMT.size();i++){
    pmtname.str("");
    
    if(!strcmp(bottom_or_top,"top")) pmtname<<fPMTLog[i]->GetName()<<i;
    else if(!strcmp(bottom_or_top,"bottom")) pmtname<<fPMTLog[i]->GetName()<<i+12;
    else pmtname<<fPMTLog[i]->GetName();
    

    //test
    //pmtname<<fPMTLog[i]->GetName();
    //end test

    physVolVec.push_back((G4VPhysicalVolume*)(new G4PVPlacement(rot,rPMT[i],pmtname.str().c_str(),fPMTLog[i],fMotherPhys,false,0)));
  }
  return physVolVec;
}







void ArDM_DetectorConstruction::addPMT(G4String bottom_or_top){

  ostringstream name;
  name<<"pmt";

  //define the shape of the PMT
  //double deltaTheta = PMT_DELTA_THETA;
  //G4LogicalVolume* fPMTLog = constructPMT(fPMTMat,bottom_or_top,1,PMT_INNER_RADIUS,PMT_OUTER_RADIUS,name.str().c_str(),deltaTheta);  


  //middle cylinder
  G4double middle_cylinder_innerR = APPROX_PMT_MIDDLE_CYLINDER_INNER_RADIUS;
  G4double middle_cylinder_outerR = APPROX_PMT_MIDDLE_CYLINDER_OUTER_RADIUS;
  G4double middle_cylinder_halfz  = APPROX_PMT_MIDDLE_CYLINDER_HALF_HEIGHT;
  G4Tubs*  middle_cylinder_solid  = new G4Tubs("pmtMiddleCylinder",middle_cylinder_innerR,middle_cylinder_outerR,middle_cylinder_halfz,0*deg,360*deg);
  
  //spherical part : top
  G4double spherical_part_innerR_top  = APPROX_PMT_SPHERICAL_PART_INNER_RADIUS;
  G4double spherical_part_outerR_top  = APPROX_PMT_SPHERICAL_PART_OUTER_RADIUS;
  G4double spherical_part_start_theta_top = 0;
  G4double spherical_part_delta_theta_top = APPROX_PMT_SPHERICAL_PART_OPENING_ANGLE/2;

  G4Sphere* sph_part_solid_top = new G4Sphere("sphPartTop",spherical_part_innerR_top,spherical_part_outerR_top,
					      0*deg,360*deg,spherical_part_start_theta_top,spherical_part_delta_theta_top);



  //spherical part : btm
  G4double spherical_part_innerR_btm  = APPROX_PMT_SPHERICAL_PART_INNER_RADIUS;
  G4double spherical_part_outerR_btm  = APPROX_PMT_SPHERICAL_PART_OUTER_RADIUS;
  G4double spherical_part_start_theta_btm = 180*deg - APPROX_PMT_SPHERICAL_PART_OPENING_ANGLE/2;
  G4double spherical_part_delta_theta_btm = (APPROX_PMT_SPHERICAL_PART_OPENING_ANGLE-APPROX_PMT_BTM_SPHERE_HOLE_OPENING_ANGLE)/2;

  G4Sphere* sph_part_solid_btm = new G4Sphere("sphPartBtm",spherical_part_innerR_btm,spherical_part_outerR_btm,
					      0*deg,360*deg,spherical_part_start_theta_btm,spherical_part_delta_theta_btm);



  //btm tube
  //connecting the pmt with pmt base
  G4double btmTube_innerR = APPROX_PMT_BTM_TUBE_INNER_RADIUS;
  G4double btmTube_outerR = APPROX_PMT_BTM_TUBE_OUTER_RADIUS;
  G4double btmTube_halfz  = APPROX_PMT_BTM_TUBE_HALF_HEIGHT;
  G4Tubs*  btmTube_solid  = new G4Tubs("pmtBtmTube",btmTube_innerR,btmTube_outerR,btmTube_halfz,0*deg,360*deg);


  //relative positions of parts of pmt to pmt center 
  //pmt center = center of middle cylinder
  G4ThreeVector topSphPos(0,0, -APPROX_PMT_DISTANCE_PMT_CENTER_TO_TOP_SPHERE_CENTER); 
  G4ThreeVector btmSphPos(0,0,  APPROX_PMT_DISTANCE_PMT_CENTER_TO_BTM_SPHERE_CENTER);
  G4ThreeVector btmTubePos(0,0,-APPROX_PMT_DISTANCE_PMT_CENTER_TO_BTM_TUBE_CENTER  );


  G4UnionSolid* pmtSolid = new G4UnionSolid("pmt",middle_cylinder_solid,sph_part_solid_top,0,topSphPos);
  pmtSolid = new G4UnionSolid("pmt",pmtSolid,sph_part_solid_btm,0,btmSphPos);
  pmtSolid = new G4UnionSolid("pmt",pmtSolid,btmTube_solid,0,btmTubePos);

  G4LogicalVolume* fPMTLog = new G4LogicalVolume(pmtSolid,fPMTMat,name.str().c_str());

  G4VisAttributes* pmtAtt = new G4VisAttributes(1);
  pmtAtt->SetColour(0.0,1.0,0.0); //green
  pmtAtt->SetForceAuxEdgeVisible(true);
  fPMTLog->SetVisAttributes(pmtAtt);

  //place PMT in the detector
  G4VPhysicalVolume* fMotherPhys;
  G4double z;

  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = APPROX_TOP_PMT_Z_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = APPROX_TOP_PMT_Z;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = APPROX_BTM_PMT_Z_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = APPROX_BTM_PMT_Z;       }
  }

  vector<G4ThreeVector> rPMT = setPMTVector(bottom_or_top);
  for(int i=0;i<rPMT.size();i++) rPMT[i].setZ(z);

  vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(fPMTLog,fMotherPhys,rPMT,bottom_or_top);

  if(bottom_or_top == "bottom")   fBottomPMTArrayPhys = pmtPhysArray;
  else if(bottom_or_top == "top") fTopPMTArrayPhys    = pmtPhysArray;

  addPMTCoat(bottom_or_top);
  addPMTCathode(bottom_or_top);  //<-- this is the detector !
  
#if TEST_BRANCH24
  //neutron background simulation
  addPMTElectrode(bottom_or_top);
  addPMTBase(bottom_or_top);
#endif //TEST_BRANCH24
  return;
}









void ArDM_DetectorConstruction::addPMTElectrode(G4String bottom_or_top){

  ostringstream name;
  name<<"pmtElectrode";



  G4double innerR = APPROX_PMT_ELECTRODE_INNER_RADIUS;
  G4double outerR = APPROX_PMT_ELECTRODE_OUTER_RADIUS;
  G4double halfz  = APPROX_PMT_ELECTRODE_HALF_HEIGHT;

  G4Tubs* solid = new G4Tubs("pmtElectrode",innerR,outerR,halfz,0*deg,360*deg);
  G4LogicalVolume* log = new G4LogicalVolume(solid,fPMTElectrodeMat,"pmtElectrode");


  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = APPROX_PMT_ELECTRODE_Z_TOP_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = APPROX_PMT_ELECTRODE_Z_TOP;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = APPROX_PMT_ELECTRODE_Z_BTM_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = APPROX_PMT_ELECTRODE_Z_BTM;       }
  }
  
  vector<G4ThreeVector> rPMTBase = setPMTVector(bottom_or_top);
  for(int i=0; i < rPMTBase.size();i++) rPMTBase[i].setZ(z);

  vector<G4VPhysicalVolume*> pmtElectrodePhysArray = placePMT(log,fMotherPhys,rPMTBase,bottom_or_top);

  if(bottom_or_top == "bottom")   fBottomPMTElectrodeArrayPhys = pmtElectrodePhysArray;
  else if(bottom_or_top == "top") fTopPMTElectrodeArrayPhys    = pmtElectrodePhysArray;
  

  return;
}










void ArDM_DetectorConstruction::addPMTCathode(G4String bottom_or_top){

  ostringstream name;
  name<<"pmtCathode";

  //define the shape of the PMT
  
  //G4Sphere* pmtSolid = new G4Sphere("pmt",pmt_inner_radius,pmt_outer_radius,0.*deg,360.*deg,startTheta,dTheta);
  G4Sphere* pmtSolid = new G4Sphere(name.str().c_str(),PMT_CATHODE_INNER_RADIUS,PMT_CATHODE_OUTER_RADIUS,0.*deg,360.*deg,0*deg,PMT_CATHODE_ACTIVE_RANGE);
  G4LogicalVolume* fPMTCathodeLog = new G4LogicalVolume(pmtSolid,fPMTCathodeMat,name.str().c_str());

  G4VisAttributes* pmtAtt = new G4VisAttributes(1);
  pmtAtt->SetColour(0.0,1.0,0.0); //green
  pmtAtt->SetForceAuxEdgeVisible(true);
  fPMTCathodeLog->SetVisAttributes(pmtAtt);  

  //define the behaviour of the sensitive region of the PMT
  G4VSensitiveDetector* sensPMT = getSD(bottom_or_top);
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->AddNewDetector(sensPMT);

  //mark PMT as sensitive
  //if not the whole PMT's surface is sensitive --> use readout geometry to define the sensitive region !
  fPMTCathodeLog->SetSensitiveDetector(sensPMT);


  //place PMT in the detector
  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = TOP_PMT_CATHODE_Z_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = TOP_PMT_CATHODE_Z;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = BTM_PMT_CATHODE_Z_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = BTM_PMT_CATHODE_Z;       }
  }

  vector<G4ThreeVector> rPMT = setPMTVector(bottom_or_top);
  for(int i=0;i<rPMT.size();i++) rPMT[i].setZ(z);

  vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(fPMTCathodeLog,fMotherPhys,rPMT,bottom_or_top);
  //cout<<"fPMTCathodeLog->GetName() "<<fPMTCathodeLog->GetName()<<endl;getchar();

  if(bottom_or_top == "bottom")   fBottomPMTCathodeArrayPhys = pmtPhysArray;
  else if(bottom_or_top == "top") fTopPMTCathodeArrayPhys    = pmtPhysArray;

  return;
}









void ArDM_DetectorConstruction::addPMTCoat_universal_pmtCoat_thickness(G4String bottom_or_top){

  ostringstream name;
  name<<"pmtCoat";

  double PMT_COATING_INNER_RADIUS;
  double PMT_COATING_OUTER_RADIUS;
  double PMT_COATING_OPENING_ANGLE;

  if(bottom_or_top == "top"){
    PMT_COATING_INNER_RADIUS = TOP_PMT_COATING_INNER_RADIUS;
    PMT_COATING_OUTER_RADIUS = TOP_PMT_COATING_OUTER_RADIUS;
    PMT_COATING_OPENING_ANGLE= TOP_PMT_COATING_OPENING_ANGLE;
  }else{
    PMT_COATING_INNER_RADIUS = BOTTOM_PMT_COATING_INNER_RADIUS;
    PMT_COATING_OUTER_RADIUS = BOTTOM_PMT_COATING_OUTER_RADIUS;
    PMT_COATING_OPENING_ANGLE= BOTTOM_PMT_COATING_OPENING_ANGLE;
  }

  G4double deltaTheta = PMT_COATING_OPENING_ANGLE/2;

  G4Sphere* pmtSolid = new G4Sphere(name.str().c_str(),PMT_COATING_INNER_RADIUS,PMT_COATING_OUTER_RADIUS,0.*deg,360.*deg,0*deg,PMT_COATING_OPENING_ANGLE);
  G4LogicalVolume* fPMTLog = new G4LogicalVolume(pmtSolid,fPMTCathodeMat,name.str().c_str());

  G4VisAttributes* pmtAtt = new G4VisAttributes(1);
  pmtAtt->SetColour(0.0,1.0,0.0); //green
  pmtAtt->SetForceAuxEdgeVisible(true);
  fPMTLog->SetVisAttributes(pmtAtt);  

  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = TOP_PMT_COATING_Z_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = TOP_PMT_COATING_Z;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = BOTTOM_PMT_COATING_Z_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = BOTTOM_PMT_COATING_Z;       }
  }

  vector<G4ThreeVector> rPMT = setPMTVector(bottom_or_top);
  for(int i=0;i<rPMT.size();i++) rPMT[i].setZ(z);
 
  vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(fPMTLog,fMotherPhys,rPMT,bottom_or_top);

  if(bottom_or_top == "bottom")   fBottomPMTCoatArrayPhys = pmtPhysArray;
  else if(bottom_or_top == "top") fTopPMTCoatArrayPhys    = pmtPhysArray;

  return;
}












void ArDM_DetectorConstruction::addPMTCoat_individual_pmtCoat_thickness(G4String bottom_or_top){

  ostringstream name;
  name<<"pmtCoat";

  double PMT_COATING_INNER_RADIUS;
  double PMT_COATING_OUTER_RADIUS;
  double PMT_COATING_OPENING_ANGLE;

  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = TOP_PMT_COATING_Z_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = TOP_PMT_COATING_Z;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = BOTTOM_PMT_COATING_Z_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = BOTTOM_PMT_COATING_Z;       }
  }

  vector<G4ThreeVector> rPMT = setPMTVector(bottom_or_top);
  for(int i=0;i<rPMT.size();i++) rPMT[i].setZ(z);
 
  vector<G4LogicalVolume*> fPMTLog;

  for(int i=0;i<rPMT.size();i++){

    if(bottom_or_top == "top"){
      PMT_COATING_INNER_RADIUS  = TOP_PMT_COATING_INNER_RADIUS;
      double pmtCoatThickness   = getWLSThickness(TOP_CONV_EFF[i]);
      PMT_COATING_OUTER_RADIUS  = PMT_COATING_INNER_RADIUS + pmtCoatThickness;
      PMT_COATING_OPENING_ANGLE = TOP_PMT_COATING_OPENING_ANGLE;
    }else{
      PMT_COATING_INNER_RADIUS  = BOTTOM_PMT_COATING_INNER_RADIUS;
      double pmtCoatThickness   = getWLSThickness(BOTTOM_CONV_EFF[i]);
      PMT_COATING_OUTER_RADIUS  = PMT_COATING_INNER_RADIUS + pmtCoatThickness;
      PMT_COATING_OPENING_ANGLE = BOTTOM_PMT_COATING_OPENING_ANGLE;
    }

    G4double deltaTheta = PMT_COATING_OPENING_ANGLE/2;

    G4Sphere* pmtSolid = new G4Sphere(name.str().c_str(),PMT_COATING_INNER_RADIUS,PMT_COATING_OUTER_RADIUS,0.*deg,360.*deg,0*deg,PMT_COATING_OPENING_ANGLE);
    G4LogicalVolume* log = new G4LogicalVolume(pmtSolid,fPMTCathodeMat,name.str().c_str());
    fPMTLog.push_back(log);
  }


  vector<G4VPhysicalVolume*> pmtPhysArray;
  pmtPhysArray = placePMT(fPMTLog,fMotherPhys,rPMT,bottom_or_top);

  if(bottom_or_top == "bottom")   fBottomPMTCoatArrayPhys = pmtPhysArray;
  else if(bottom_or_top == "top") fTopPMTCoatArrayPhys    = pmtPhysArray;

  return;
}








void ArDM_DetectorConstruction::addPMTCoat(G4String bottom_or_top){

  addPMTCoat_universal_pmtCoat_thickness(bottom_or_top);
  //addPMTCoat_individual_pmtCoat_thickness(bottom_or_top);

  return;
}




void ArDM_DetectorConstruction::addPMTBase(G4String bottom_or_top){


  //simplistic :
  //pmt base ~ cylinder

  G4double innerR = 0;
  G4double outerR = APPROX_PMT_BASE_RADIUS;
  G4double halfz  = APPROX_PMT_BASE_HALF_THICKNESS;

  G4Tubs* solid = new G4Tubs("pmtBase",innerR,outerR,halfz,0*deg,360*deg);
  G4LogicalVolume* log = new G4LogicalVolume(solid,fPMTBaseMat,"pmtBase");
  

  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = APPROX_PMT_BASE_Z_TOP_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = APPROX_PMT_BASE_Z_TOP;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = APPROX_PMT_BASE_Z_BTM_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = APPROX_PMT_BASE_Z_BTM;       }
  }
  

  vector<G4ThreeVector> rPMTBase = setPMTVector(bottom_or_top);
  for(int i=0; i < rPMTBase.size();i++) rPMTBase[i].setZ(z);

  vector<G4VPhysicalVolume*> pmtBasePhysArray;
  pmtBasePhysArray = placePMT(log,fMotherPhys,rPMTBase,bottom_or_top);

  if(bottom_or_top == "bottom")   fBottomPMTBaseArrayPhys = pmtBasePhysArray;
  else if(bottom_or_top == "top") fTopPMTBaseArrayPhys    = pmtBasePhysArray;
  
  return;
}







G4double ArDM_DetectorConstruction::GArEpsilon(G4double wavelength){

  G4double inverseWavelength = 1./wavelength/wavelength;
  G4double epsilon = 1.2055e-2*(.2075/(91.012 - inverseWavelength) + 
				.0415/(87.892 - inverseWavelength) + 
				4.333/(214.02 - inverseWavelength));
  
  epsilon *= 8./12.; //--> why ?
  //<-- compare equation (1) and (2) in Bideau-Mehu's paper
  //--> (n2 - 1)/(n2 + 2) = (n-1)*8./12.
  epsilon  = (2*epsilon+1)/(1-epsilon);
  return epsilon;
}



G4double ArDM_DetectorConstruction::GArRefIndex(G4double wavelength){

  return sqrt(GArEpsilon(wavelength));
}




G4double ArDM_DetectorConstruction::LArEpsilon(G4double wavelength){

  G4double inverseWavelength = 1./wavelength/wavelength;
  G4double epsilon = 1.2055e-2*(.2075/(91.012 - inverseWavelength) + 
				.0415/(87.892 - inverseWavelength) + 
				4.333/(214.02 - inverseWavelength));
  
  epsilon *= 8./12.; //--> why ?
  //<-- compare equation (1) and (2) in Bideau-Mehu's paper
  //--> (n2 - 1)/(n2 + 2) = (n-1)*8./12.

  epsilon *= LAR_RHO/GAR_RHO; //density-correction
  epsilon  = (2*epsilon+1)/(1-epsilon);
  return epsilon;

}


G4double ArDM_DetectorConstruction::LArRefIndex(G4double wavelength){

  G4double nliquid = sqrt(LArEpsilon(wavelength));
  if(VERBOSE) G4cout<<"in ArDM_DetectorConstruction::LArRefIndex(..), wavelength "<<wavelength<<"\t nliquid "<<nliquid<<G4endl;
  return nliquid;
}



G4double ArDM_DetectorConstruction::RayleighAttenuationLength_LAr(G4double wavelength){
  //for underlying physics --> see ArDM_DetectorConstruction::setMatPropTab_detMat()
  //wavelength in micrometer !!
  G4double epsilon    = LArEpsilon(wavelength);
  G4double derivative = (epsilon-1)*(epsilon+2)/3/LAR_RHO;
  derivative *= derivative;
  wavelength *= micrometer/m/2/PI; //wavelength was in micrometer, convert it back to meter
  G4double h = KE_LAR/wavelength/wavelength/wavelength/wavelength*derivative;
  return 1./h*m; //returned value in meter ! (since KE_LAR,derivative,wavelength are all in SI units)
}





G4double ArDM_DetectorConstruction::RayleighAttenuationLength_GAr(G4double wavelength){
 //  //for underlying physics --> see ArDM_DetectorConstruction::setMatPropTab_detMat()
//   //wavelength in micrometer !!
//   G4double epsilon    = GArEpsilon(wavelength);
//   G4double derivative = (epsilon-1)*(epsilon+2)/3/GAR_RHO;
//   derivative *= derivative;
//   wavelength *= micrometer/m/2/PI; //wavelength was in micrometer, convert it back to meter
//   G4double h = KE_LAR/wavelength/wavelength/wavelength/wavelength*derivative;
//   return 1./h*m; //returned value in meter ! (since KE_LAR,derivative,wavelength are all in SI units)

  //for the time being, just set it to 100m for GAr
  return 100. *m;
}



G4double ArDM_DetectorConstruction::getPhotonWavelength(G4double photonEnergy){
  //the returned wavelength is measured in micrometer !!
  return PLANCK_H*LIGHTSPEED/photonEnergy/micrometer;
}



void ArDM_DetectorConstruction::setMatPropTab_LAr(){
  //detMat = liquid Argon

  //LAr's properties table: refractive index, rayleigh attenuation length, absorption length.
  //i. absorption length ~ detector dimension 
  //   (for the time being: photon will be killed when touching the detector's walls)
  //
  //
  //
  //for the other 2 quantities: read the following 2 papers:
  //
  //1. "rayleigh scattering in rare-gas liquids" (G.M.Seidel et al., 2002)
  //   link: http://arxiv.org/abs/hep-ex/0111054
  //
  //2. "measurement of refractive indices of Ne,Ar,Kr,Xe ..." (A.Bideau-Mehu et al., 1980)
  //   linke: http://adsabs.harvard.edu/abs/1981JQSRT..25..395B
  //
  //
  //
  //ii. refractive index n = sqrt(epsilon_r)
  //    (epsilon_r : permittivity)
  //
  //iii. rayleigh attenuation length
  //
  //
  //***************************************************
  //
  //from the 2nd paper:
  //for rare-gas: n(lambda) - 1 ~ Ne^2/(8pi^2*epsilon0*mc^2)*Sum_i(fi/(lambda_i^-2 - lambda^-2))
  //
  //lambda    : wavelength being in consideration
  //n(lambda) : wavelength-dependent refractive index of the medium
  //N         : number density of the medium
  //e,m       : charge/mass of electron
  //epsilon0  : permittivity of vacuum
  //lambda_i  : wavelength
  //fi        : absorption oscillator strength for transitions of wavelength lambda_i
  //
  // ==> epsilon_r = ... (n = sqrt(epsilon_r))
  //
  //numerical result for gaseous(!!) argon:
  //this formula holds at standard condition (273 K, 1 bar)
  //
  //n(lamb) - 1 = 1.2055e-2*(.2075/(91.012 - lamb^-2) + .0415/(87.892 - lamb^-2) + 4.3330/(214.02 - lamb^-2))
  //
  //
  //***ATTENTION***: this is for rare gas ! --> for liquid : correction needed ! (e.g. density correction)
  //
  //for liquid: (epsilon-1)/(epsilon+2) = A*rho 
  //with epsilon being dependent on lambda.
  //rho: mass density
  //A  : look at table 2 in Seidel's paper ! (for Ar at boiling point: A = .1045 cm3/g)
  //
  //also in the 2nd paper: in visible wavelength range, values of A for gases and liquids are so close to each other
  //--> intermolecular interactions have **very little** effect on the atomic polarizability at liquid-density.
  //--> assumption: dielectric constant epsilon at scintillation wavelength of noble liquids can be obtained
  //from those determined in gases !
  //
  //
  //***************************************************
  //
  //from the 1st paper:
  //rayleigh attenuation length:
  //
  //l = h^-1 = [1/(6pi*lambda^4)*k*T*rho^2*kappa_T*((d(epsilon)/d(rho))_T)^2]^-1
  //
  //k         : boltzmann const.
  //T         : temperature
  //kappa_T   : isothermal compressibility
  //rho       : liquid density
  //(d(epsilon)/d(rho))_T  : partial derivative of permittivity epsilon by density rho at constant temperature T
  //
  //
  //the partial derivaty can be approximated by: 
  //
  //(d(epsilon)/d(rho))_T = (epsilon - 1)(epsilon + 2) /3/rho
  //
  //***************************************************
  
  if(!fLAr){
    G4cout<<"in ArDM_DetectorConstruction::setMatPropTab_LAr(), fLAr = 0. exit."<<G4endl;
    return;
  }


  G4double lambdaMin = 115. *nm; //from SVN; why this number ? 
  G4double lambdaMax = 598. *nm; //from SVN; why this number ?
  const G4int nentries = 100;      //why this number ?

  
  G4double photonE[nentries],refIndex[nentries],rayleighLength[nentries],absorptionLength[nentries];
  
  G4double photonMaxE = PLANCK_H*LIGHTSPEED/lambdaMin;
  G4double photonMinE = PLANCK_H*LIGHTSPEED/lambdaMax;
  G4double dE         = (photonMaxE-photonMinE)/nentries;


  //G4double absLength  = 2. *m;//2*sqrt(DETECTOR_HALF_HEIGHT*DETECTOR_HALF_HEIGHT + DETECTOR_RADIUS*DETECTOR_RADIUS);
  G4double absLength  = 30. *m; //<-- from SVN. <-- why this number ?
  for(int i=0;i<nentries;i++){
    photonE[i]          = photonMinE + i*dE; 
    refIndex[i]         = LArRefIndex(getPhotonWavelength(photonE[i])); 
    rayleighLength[i]   = RayleighAttenuationLength_LAr(getPhotonWavelength(photonE[i]));
    absorptionLength[i] = absLength;
    if(VERBOSE) G4cout<<"in ArDM_DetectorConstruction::setMatPropTab_LAr(..), wavelength (nm): "
		      <<getPhotonWavelength(photonE[i])*micrometer/nm
		      <<"\t E (eV): "<<photonE[i]/eV
		      <<"\t refIndex "<<refIndex[i]<<"\t rayleigh length (mm) : "<<rayleighLength[i]/mm<<G4endl;

  }

  G4MaterialPropertiesTable* LAr_propTab = new G4MaterialPropertiesTable;
  LAr_propTab->AddProperty("RINDEX",   photonE,refIndex,        nentries);//->SetSpline(true);
  LAr_propTab->AddProperty("RAYLEIGH", photonE,rayleighLength,  nentries);//->SetSpline(true);
  LAr_propTab->AddProperty("ABSLENGTH",photonE,absorptionLength,nentries);//->SetSpline(true);

  setArScintProperty(LAr_propTab,"LAr");

  fLAr->SetName("LAr");
  fLAr->SetMaterialPropertiesTable(LAr_propTab);
  return;  
}




void ArDM_DetectorConstruction::setMatPropTab_GAr(){
  
  if(!fGAr){
    G4cout<<"in ArDM_DetectorConstruction::setMatPropTab_GAr(), fGAr = 0. exit."<<G4endl;
    return;
  }


  G4double lambdaMin = 115. *nm; //from SVN; why this number ? 
  G4double lambdaMax = 598. *nm; //from SVN; why this number ?
  const G4int nentries = 100;      //why this number ?
  
  G4double photonE[nentries],refIndex[nentries],rayleighLength[nentries],absorptionLength[nentries];
  
  G4double photonMaxE = PLANCK_H*LIGHTSPEED/lambdaMin;
  G4double photonMinE = PLANCK_H*LIGHTSPEED/lambdaMax;
  G4double dE         = (photonMaxE-photonMinE)/nentries;

  G4double absLength  = 100. *m;
  for(int i=0;i<nentries;i++){
    photonE[i]          = photonMinE + i*dE;
    refIndex[i]         = GArRefIndex(getPhotonWavelength(photonE[i]));
    rayleighLength[i]   = RayleighAttenuationLength_GAr(getPhotonWavelength(photonE[i]));
    absorptionLength[i] = absLength;
    if(VERBOSE) G4cout<<"in ArDM_DetectorConstruction::setMatPropTab_detMat(..), wavelength (nm): "
		      <<getPhotonWavelength(photonE[i])*micrometer/nm
		      <<"\t E (eV): "<<photonE[i]/eV
		      <<"\t refIndex "<<refIndex[i]<<"\t rayleigh length : "<<rayleighLength[i]/mm<<G4endl;
  }

  G4MaterialPropertiesTable* GAr_propTab = new G4MaterialPropertiesTable;
  GAr_propTab->AddProperty("RINDEX",   photonE,refIndex,        nentries);//->SetSpline(true);
  GAr_propTab->AddProperty("RAYLEIGH", photonE,rayleighLength,  nentries);//->SetSpline(true);
  GAr_propTab->AddProperty("ABSLENGTH",photonE,absorptionLength,nentries);//->SetSpline(true);

  setArScintProperty(GAr_propTab,"GAr");



  fGAr->SetName("GAr");
  fGAr->SetMaterialPropertiesTable(GAr_propTab);
  return;
}






void ArDM_DetectorConstruction::setArScintProperty(G4MaterialPropertiesTable* propTab,G4String medium){



  //absorption and emission spectrum is the same for GAr and LAr (??)
  //test for scint. process
  G4double meanE  = 9.68*eV;
  G4double sigmaE = .5*eV; //for the time being : arbitrary width !!
  G4double minScintE = meanE - 5*sigmaE;
  G4double maxScintE = meanE + 5*sigmaE;
  const int nScintE = 500;
  G4double delScintE = (maxScintE - minScintE)/nScintE;
  G4double scintE[nScintE];
  G4double scintFast[nScintE];
  G4double scintSlow[nScintE];
  for(int i=0;i<nScintE;i++){
    scintE[i] = minScintE + i*delScintE;
    scintFast[i] = exp(-(scintE[i]-meanE)*(scintE[i]-meanE)/sigmaE/sigmaE/2);
    scintSlow[i] = exp(-(scintE[i]-meanE)*(scintE[i]-meanE)/sigmaE/sigmaE/2);
  }
  
  propTab->AddProperty("FASTCOMPONENT",scintE,scintFast,nScintE);//->SetSpline(true);
  propTab->AddProperty("SLOWCOMPONENT",scintE,scintSlow,nScintE);//->SetSpline(true);



  if(!strcmp(medium,"LAr")){

    propTab->AddConstProperty("SCINTILLATIONYIELD",         LAR_SCINTILLATIONYIELD);
    propTab->AddConstProperty("RESOLUTIONSCALE",            LAR_RESOLUTIONSCALE);
    propTab->AddConstProperty("FASTTIMECONSTANT",           LAR_FASTTIMECONSTANT);
    propTab->AddConstProperty("SLOWTIMECONSTANT",           LAR_SLOWTIMECONSTANT);
    propTab->AddConstProperty("FASTSCINTILLATIONRISETIME",  LAR_FASTSCINTILLATIONRISETIME);
    propTab->AddConstProperty("SLOWSCINTILLATIONRISETIME",  LAR_SLOWSCINTILLATIONRISETIME);
    propTab->AddConstProperty("YIELDRATIO",                 LAR_YIELDRATIO);
    
    propTab->AddConstProperty("ELECTRONEXCITATIONRATIO",    LAR_EXCITATIONRATIO_ELECTRON);
    propTab->AddConstProperty("NEUTRONEXCITATIONRATIO",     LAR_EXCITATIONRATIO_NEUTRON);
    //different component ratio for different particles
    //--> use EXCITATIONRATIO  instead !
    
    //end test scint. proc.
    


  }else if(!strcmp(medium,"GAr")){
    
    propTab->AddConstProperty("SCINTILLATIONYIELD",         GAR_SCINTILLATIONYIELD);
    propTab->AddConstProperty("RESOLUTIONSCALE",            GAR_RESOLUTIONSCALE);
    propTab->AddConstProperty("FASTTIMECONSTANT",           GAR_FASTTIMECONSTANT);
    propTab->AddConstProperty("SLOWTIMECONSTANT",           GAR_SLOWTIMECONSTANT);
    propTab->AddConstProperty("FASTSCINTILLATIONRISETIME",  GAR_FASTSCINTILLATIONRISETIME);
    propTab->AddConstProperty("SLOWSCINTILLATIONRISETIME",  GAR_SLOWSCINTILLATIONRISETIME);
    propTab->AddConstProperty("YIELDRATIO",                 GAR_YIELDRATIO);
    
    propTab->AddConstProperty("ELECTRONEXCITATIONRATIO",    GAR_EXCITATIONRATIO_ELECTRON);
    propTab->AddConstProperty("NEUTRONEXCITATIONRATIO",     GAR_EXCITATIONRATIO_NEUTRON);
  
  }
  
  return;
}








G4Material* ArDM_DetectorConstruction::getTankMat(){
  //code from SVN
  
  G4Element* Si = new G4Element("Silicon" , "Si", 14., 28.09   *g/mole);   
  G4Element* Ni = new G4Element("Nickel"  , "Ni", 28., 58.6934 *g/mole);
  G4Element* Cr = new G4Element("Chromium", "Cr", 24., 51.9961 *g/mole);  
  G4Element* Fe = new G4Element("Iron"    , "Fe", 26., 55.845  *g/mole);  
  
  G4Material* stainless_steel = new G4Material("StainlessSteel",8.00 *g/cm3,4); //2.arg = density, 3.arg = nelements
  stainless_steel->AddElement(Si,1. *perCent);
  stainless_steel->AddElement(Ni,10.*perCent);
  stainless_steel->AddElement(Cr,19.*perCent);
  stainless_steel->AddElement(Fe,70.*perCent);
  
  return stainless_steel;
}








void ArDM_DetectorConstruction::addWLS(){
  
  //ring sector
  G4double wlsRing_innerR      = WLS_RINGSEC_INNER_RADIUS;
  G4double wlsRing_outerR      = WLS_RINGSEC_OUTER_RADIUS;
  G4double wlsRing_half_height = WLS_RINGSEC_HALF_HEIGHT;
  G4double wlsRing_start_phi   = WLS_RINGSEC_START_PHI;
  G4double wlsRing_delta_phi   = WLS_RINGSEC_DELTA_PHI;
  
  G4Tubs* wlsRingSolid = new G4Tubs("wlsRing",wlsRing_innerR,wlsRing_outerR,wlsRing_half_height,
				    wlsRing_start_phi,wlsRing_delta_phi);
  
  //linear sector
  G4double wls_linsec_halfx = WLS_LINSEC_HALF_X;
  G4double wls_linsec_halfy = WLS_LINSEC_HALF_Y;
  G4double wls_linsec_halfz = WLS_LINSEC_HALF_Z;
  
  G4Box* wlsLinSecSolid = new G4Box("wlsLinSec",wls_linsec_halfx,wls_linsec_halfy,wls_linsec_halfz);
  
  
  //position of the linear sector relative to the ring sector
  G4RotationMatrix* wls_linsec_rotMat = new G4RotationMatrix;
  wls_linsec_rotMat->rotateX(0);
  wls_linsec_rotMat->rotateY(0);
  wls_linsec_rotMat->rotateZ(90.*deg);
  
  //G4ThreeVector wlsLinSecPos(WLS_LINSEC_POS_X,WLS_LINSEC_POS_Y,WLS_LINSEC_POS_Z); //relative to tank !
  G4ThreeVector wlsLinSecPos(WLS_LINSEC_POS_X,WLS_LINSEC_POS_Y,0.); //relative to the ring sector
  
  //build unionSolid of ring and linear sector
  G4UnionSolid* wls = new G4UnionSolid("WLS",wlsRingSolid,wlsLinSecSolid,wls_linsec_rotMat,wlsLinSecPos);   
  G4LogicalVolume* fWLSLog = new G4LogicalVolume(wls,fWLSMat,"WLS",0,0,0);


  G4VisAttributes* wlsAtt = new G4VisAttributes(true);
  wlsAtt->SetColour(1.0,.0,1.0); //magenta
  wlsAtt->SetForceAuxEdgeVisible(true);
  fWLSLog->SetVisAttributes(wlsAtt);

  //place fWLSLog in fLArColLog
  
  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(fLArColPhys){ fMotherPhys = fLArColPhys; z = WLS_RINGSEC_POS_Z_LARCOL;}
  else{            fMotherPhys = fWorldPhys;  z = WLS_RINGSEC_POS_Z;       }
  G4ThreeVector wlsPos(WLS_RINGSEC_POS_X,WLS_RINGSEC_POS_Y,z);
  fWLSPhys = new G4PVPlacement(0,wlsPos,"WLS",fWLSLog,fMotherPhys,false,0);
  return;   
}



void ArDM_DetectorConstruction::addWLSSupport(){
  
  //ring sector
  G4double wlsRing_innerR      = WLS_SUPPORT_RINGSEC_INNER_RADIUS;
  G4double wlsRing_outerR      = WLS_SUPPORT_RINGSEC_OUTER_RADIUS;
  G4double wlsRing_half_height = WLS_SUPPORT_RINGSEC_HALF_HEIGHT;
  G4double wlsRing_start_phi   = WLS_SUPPORT_RINGSEC_START_PHI;
  G4double wlsRing_delta_phi   = WLS_SUPPORT_RINGSEC_DELTA_PHI;
  
  G4Tubs* wlsRingSolid = new G4Tubs("wlsSupportRing",wlsRing_innerR,wlsRing_outerR,wlsRing_half_height,
				    wlsRing_start_phi,wlsRing_delta_phi);
  
  //linear sector
  G4double wls_linsec_halfx = WLS_SUPPORT_LINSEC_HALF_X;
  G4double wls_linsec_halfy = WLS_SUPPORT_LINSEC_HALF_Y;
  G4double wls_linsec_halfz = WLS_SUPPORT_LINSEC_HALF_Z;
  
  G4Box* wlsLinSecSolid = new G4Box("wlsSupportLinSec",wls_linsec_halfx,wls_linsec_halfy,wls_linsec_halfz);
  
  G4RotationMatrix* wls_linsec_rotMat = new G4RotationMatrix;
  wls_linsec_rotMat->rotateX(0);
  wls_linsec_rotMat->rotateY(0);
  wls_linsec_rotMat->rotateZ(90.*deg);
    
  G4ThreeVector wlsLinSecPos(WLS_SUPPORT_LINSEC_POS_X,WLS_SUPPORT_LINSEC_POS_Y,0.); //relative to the ring sector
  
//   //build unionSolid of ring and linear sector
  G4UnionSolid* wls = new G4UnionSolid("WLSSupport",wlsRingSolid,wlsLinSecSolid,wls_linsec_rotMat,wlsLinSecPos);   
  G4LogicalVolume* fWLSLog = new G4LogicalVolume(wls,fWLSSupportMat,"WLSSupport",0,0,0);

//   //place fWLSLog in fLArColLog
  G4ThreeVector wlsPos(WLS_SUPPORT_RINGSEC_POS_X,WLS_SUPPORT_RINGSEC_POS_Y,WLS_SUPPORT_LINSEC_POS_Z_LARCOL);
  fWLSSupportPhys = new G4PVPlacement(0,wlsPos,"WLSSupport",fWLSLog,fLArColPhys,false,0);

  G4VisAttributes* wlsAtt = new G4VisAttributes(true);
  wlsAtt->SetColour(1.0,.0,1.0); //magenta
  wlsAtt->SetForceAuxEdgeVisible(true);
  fWLSLog->SetVisAttributes(wlsAtt);

  return;   
}







void ArDM_DetectorConstruction::addBottomSideReflector(){
  
  //new bottom side reflector
  //a cone spanned by the lower edge of the reflector and the bottom pmtSupport
  //during gas test in mar/apr 2013
  
  G4double innerR1     = BOTTOM_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE;
  G4double outerR1     = BOTTOM_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE;
  G4double innerR2     = BOTTOM_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE;
  G4double outerR2     = BOTTOM_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE;
  G4double half_height = BOTTOM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST;
  G4double start_phi   = BOTTOM_SIDE_REFLECTOR_START_PHI;
  G4double delta_phi   = BOTTOM_SIDE_REFLECTOR_DELTA_PHI;
  
  G4Cons* coneSolid = new G4Cons("bottomSideReflector",innerR1,outerR1,innerR2,outerR2,half_height,start_phi,delta_phi);
  
  
  G4double cynPart_innerR = BOTTOM_SIDE_REFLECTOR_CYNLIDRICAL_PART_INNER_RADIUS;
  G4double cynPart_outerR = BOTTOM_SIDE_REFLECTOR_CYNLIDRICAL_PART_OUTER_RADIUS;
  G4double cynPart_half_height = BOTTOM_SIDE_REFLECTOR_CYNLIDRICAL_PART_HALF_HEIGHT;
  
  G4Tubs* tubeSolid = new G4Tubs("bottomSideReflector",cynPart_innerR,cynPart_outerR,
				 cynPart_half_height,0*deg,360.*deg);
  
  
  G4ThreeVector pos_tube_rel_cone(0.,0.,
				  BOTTOM_SIDE_REFLECTOR_CYNLIDRICAL_PART_POS_Z-BOTTOM_SIDE_REFLECTOR_POS_Z_GASTEST);
  
  G4UnionSolid* solid = new G4UnionSolid("bottomSideReflector",coneSolid,tubeSolid,0,pos_tube_rel_cone);
  
  G4LogicalVolume* fBottomSideReflectorLog = new G4LogicalVolume(solid,fBottomSideReflectorMat,
								 "bottomSideReflector",0,0,0);
  
  //place fBottomSideReflectorLog in fLArColLog
  G4ThreeVector pos(BOTTOM_SIDE_REFLECTOR_POS_X,
		    BOTTOM_SIDE_REFLECTOR_POS_Y,
		    BOTTOM_SIDE_REFLECTOR_POS_Z_GASTEST_LARCOL);
  
  fBottomSideReflectorPhys = new G4PVPlacement(0,pos,"bottomSideReflector",fBottomSideReflectorLog,fLArColPhys,false,0);
  
  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  sideReflAtt->SetColour(G4Color::Green());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fBottomSideReflectorLog->SetVisAttributes(sideReflAtt);
  
  //     cout<<"BOTTOM_SIDE_REFLECTOR_POS_Z_GASTEST_LARCOL "<<BOTTOM_SIDE_REFLECTOR_POS_Z_GASTEST_LARCOL
  // 	<<"\t BOTTOM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST "<<BOTTOM_SIDE_REFLECTOR_HALF_HEIGHT_GASTEST
  // 	<<endl;getchar();

  
  
  
#if TURN_ON_SIDE_REFLECTOR_COATING
  addBottomSideReflectorCoat();
#endif //TURN_ON_SIDE_REFLECTOR_COATING
  
  return;   
}










void ArDM_DetectorConstruction::addTopSideReflector(){

  //top side refl. may span over GArCol-volume and LArCol-volume
  //hence, placing the top side refl. in only GArCol or LArCol will cause weird things
  //namely :
  //assuming we place the top side refl. as daughter volume of GArCol, then :
  //
  //1. GEANT4 doesn't complain when the top side reflector sticks outside of the GArCol-volume
  //   --> keyword "overlapping"
  //
  //2. but when we build the optical surface between GArCol / LArCol and the top side reflector, 
  //only the surface between GArCol and the top side reflector works as we want it to. 
  //the surface between LArCol and the top side reflector doesn't ! 
  //
  //3. moreover, the optical photon can somehow just simply go through the part of the top side reflector, which sticks into the LArCol volume 
  //   --> go figure !!
  //
  //--> solution : 
  //devide the topside reflector into 2 parts
  //i.  the part in GArCol, and
  //ii. the part in LArCol


  G4double innerR1_GAr = TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double outerR1_GAr = TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double innerR2_GAr = TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double outerR2_GAr = TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double halfheight_GAr = TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_GARCOL;
  
  G4double innerR1_LAr = TOP_SIDE_REFLECTOR_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double outerR1_LAr = TOP_SIDE_REFLECTOR_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double innerR2_LAr = TOP_SIDE_REFLECTOR_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double outerR2_LAr = TOP_SIDE_REFLECTOR_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double halfheight_LAr = TOP_SIDE_REFLECTOR_HALF_HEIGHT_IN_LARCOL;
  
  G4double start_phi   = TOP_SIDE_REFLECTOR_START_PHI;
  G4double delta_phi   = TOP_SIDE_REFLECTOR_DELTA_PHI;

  G4Cons* ringSolid_GAr = NULL;
  G4Cons* ringSolid_LAr = NULL;

  if(halfheight_GAr) ringSolid_GAr = new G4Cons("topSideReflector",innerR2_GAr,outerR2_GAr,innerR1_GAr,outerR1_GAr,
						halfheight_GAr,start_phi,delta_phi);

  if(halfheight_LAr) ringSolid_LAr = new G4Cons("topSideReflector",innerR2_LAr,outerR2_LAr,innerR1_LAr,outerR1_LAr,
						halfheight_LAr,start_phi,delta_phi);



  G4LogicalVolume* fTopSideReflectorLog_GAr = NULL;
  G4LogicalVolume* fTopSideReflectorLog_LAr = NULL;



  if(ringSolid_GAr) fTopSideReflectorLog_GAr = new G4LogicalVolume(ringSolid_GAr,fTopSideReflectorMat,
								   "topSideReflector",0,0,0);

  if(ringSolid_LAr) fTopSideReflectorLog_LAr = new G4LogicalVolume(ringSolid_LAr,fTopSideReflectorMat,
								   "topSideReflector",0,0,0);

  
  G4ThreeVector pos_GAr(TOP_SIDE_REFLECTOR_POS_X,TOP_SIDE_REFLECTOR_POS_Y,TOP_SIDE_REFLECTOR_IN_GARCOL_POS_Z_GARCOL);
  //G4ThreeVector pos_LAr(TOP_SIDE_REFLECTOR_POS_X,TOP_SIDE_REFLECTOR_POS_Y,TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z_LARCOL);
  G4ThreeVector pos_LAr(TOP_SIDE_REFLECTOR_POS_X,TOP_SIDE_REFLECTOR_POS_Y,TOP_SIDE_REFLECTOR_IN_LARCOL_POS_Z_APPROX_LARCOL);

  if(fTopSideReflectorLog_GAr && fGArColPhys) 
    fTopSideReflectorGArPhys = new G4PVPlacement(0,pos_GAr,"topSideReflector",fTopSideReflectorLog_GAr,fGArColPhys,false,0);

  if(fTopSideReflectorLog_LAr && fLArColPhys) 
    fTopSideReflectorLArPhys = new G4PVPlacement(0,pos_LAr,"topSideReflector",fTopSideReflectorLog_LAr,fLArColPhys,false,0);
    

  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  sideReflAtt->SetColour(G4Color::Green());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fTopSideReflectorLog_GAr->SetVisAttributes(sideReflAtt);
  fTopSideReflectorLog_LAr->SetVisAttributes(sideReflAtt);



#if TURN_ON_SIDE_REFLECTOR_COATING
  addTopSideReflectorCoat();
#endif //TURN_ON_SIDE_REFLECTOR_COATING


  return;   
}








void ArDM_DetectorConstruction::setMatPropTab_WLS(){
  //from SVN
  
  //material
  //for the time being: TPB
  
  G4Element* C = new G4Element("Carbon"  , "C", 6., 12.0107  *g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", 1.,  1.0079 *g/mole); 

  
  //informations about TPB can be found here
  //http://www.chemnet.com/cas/en/1450-63-1/1,1,4,4-Tetraphenyl-1,3-butadiene.html
  fWLSMat = new G4Material("TPB",1.079 *g/cm3,2);
  fWLSMat->AddElement(C,28);
  fWLSMat->AddElement(H,22);

  
  //for a certain thickness d of the WLS layer
  //the conversion efficiency p is determined by the thickness d 
  //and the absorption length (WLSABSLENGTH) l
  //
  //the conversion efficiency can be understood as 1 - transmittance
  //--> probability that the photon will go through the WLS layer without being absorbed is
  // 1 - p = exp(-d/l)
  //--> absorption length l = -d/ln(1-p)
  
  
  //test !!!!
  const G4int n=300;
  G4double Eemis[n],emission[n];
  
  G4double minTestE1 = 6.5 *eV;
  G4double maxTestE1 = 13.5 *eV;

  G4double minTestE2 = 0.1 *eV;
  G4double maxTestE2 = 6.  *eV;

#if (!ARDM) //if we use sosuke's test setup

#if _TESTSETUP_DETECTORCONSTRUCTION_

  conversionEfficiency  = TESTSETUP_WLS_CONVERSION_EFFICIENCY;
  thickness_of_WLSLayer = 2*TESTSETUP_WLS_HALF_THICKNESS;

#endif //_TESTSETUP_DETECTORCONSTRUCTION_

#endif //(!ARDM)

  G4double meanEmisE  = 2.85 *eV;
  G4double sigmaEmisE = .5   *eV;
  G4double minEmisE   = meanEmisE - 5*sigmaEmisE;
  G4double maxEmisE   = meanEmisE + 5*sigmaEmisE;
  //G4double delEmisE   = (maxEmisE - minEmisE)/n;
  G4double delEmisE   = (maxTestE2 - minTestE2)/n;
  
  for(G4int i=0;i<n;i++){
    //Eemis[i]    = minEmisE + i*delEmisE;
    Eemis[i]    = minTestE2 + i*delEmisE;
    if(Eemis[i] < minEmisE || Eemis[i] > maxEmisE) emission[i] = 0;
    else
      emission[i] = exp(-(Eemis[i] - meanEmisE)*(Eemis[i] - meanEmisE)/2/sigmaEmisE/sigmaEmisE);
  }


  //end test

  G4double ErefIndex[] = {1.85,10.69};
  G4double refIndex[]  = {1.635,1.635};
  G4double Eabs[]      = {.1*eV,6.0*eV,6.0001*eV,13.5*eV};
  G4double absLength[] = {10*m,10*m,WLS_MEAN_ABSORPTION_LENGTH,WLS_MEAN_ABSORPTION_LENGTH};

  
  G4MaterialPropertiesTable* wls_propTab = new G4MaterialPropertiesTable;
  wls_propTab->AddProperty("RINDEX",ErefIndex,refIndex,2);
  //wls_propTab->AddProperty("WLSABSLENGTH",Eabs,absLength,sizeof(Eabs)/sizeof(double))->SetSpline(true);
  //wls_propTab->AddProperty("WLSCOMPONENT",Eemis,emission,sizeof(Eemis)/sizeof(double))->SetSpline(true);
  //<-- SetPline(..) option sometimes causes problem ! e.g. negative absoprtion length !!!

  wls_propTab->AddProperty("WLSABSLENGTH",Eabs,absLength,sizeof(Eabs)/sizeof(double));
  wls_propTab->AddProperty("WLSCOMPONENT",Eemis,emission,sizeof(Eemis)/sizeof(double));
  wls_propTab->AddConstProperty("WLSTIMECONSTANT",2.*ns); //why 2. ns ??
  //wls_propTab->AddConstProperty("WLSMEANNUMBERPHOTONS",1);

  //end test !!


  fWLSMat->SetMaterialPropertiesTable(wls_propTab);
  fWLS_ringSec_Mat = fWLSMat;
  fWLS_linSec_Mat  = fWLS_ringSec_Mat;
  return;
}


void ArDM_DetectorConstruction::setMatPropTab_PMTMat(){

  //define PMT material
  //for the time being, let it be SiO2 (quartz)
  //with (artifitial !!! ) refractive index 2. (to avoid total reflection within WLS layer on top of PMT's surface !!)
  //WLS (TPB) has refractive index 1.635

  G4Element* Si = new G4Element("Silicon", "Si", 14., 28.09   *g/mole);
  G4Element* O  = new G4Element("Oxygen" , "O" ,  8., 15.9994 *g/mole);

  fPMTMat = new G4Material("PMTMat",2.200*g/cm3,2);
  fPMTMat->AddElement(Si, 1);
  fPMTMat->AddElement(O , 2);

  G4double E[]        = {1.85,10.69};
  G4double refIndex[] = {1.5,1.5}; //test !!

  G4MaterialPropertiesTable* pmtMat_propTab = new G4MaterialPropertiesTable;
  pmtMat_propTab->AddProperty("RINDEX",E,refIndex,sizeof(E)/sizeof(double));

  fPMTMat->SetName("Quartz");
  fPMTMat->SetMaterialPropertiesTable(pmtMat_propTab);
  return;  
}






/*

void ArDM_DetectorConstruction::unifyPhysVol(){

  //call this function when running simulation for neutron background.
  //why ? --> see function defination in .cc file for more explanation.
  //--> neutrons from contaminations in detector components are emitted from "dirty" detector components
  //
  //--> use /gps/pos/confine command in geant4 macro
  //
  //--> now, e.g., if we want to emit neutrons from PMT glass
  //--> there are 24 PMTs, if we set
  //
  // /gps/pos/confine pmtGlass1
  // /gps/pos/confine pmtGlass2
  //
  //--> this won't work, since the later command overrides the previous one.
  //
  //--> unify all phys. volumes for pmtGlass and set
  //   /gps/pos/confine unified_pmtGlass
  // 

//--> use multiple sources !!
//--> pay attention to the relative intensity between the sources !!
//



  //unify physical volumes for stainless steel tank
  // uni_tank = tank + toplid + btmlid
  
  G4ThreeVector toplidPos(0,0,TANK_HALF_HEIGHT+UPPER_LID_HALF_HEIGHT);    //position relative to tank center
  G4ThreeVector btmlidPos(0,0,-(TANK_HALF_HEIGHT+UPPER_LID_HALF_HEIGHT));
  G4VSolid* toplidSolid = fTopLidPhys->GetLogicalVolume()->GetSolid();
  G4VSolid* btmlidSolid = fBottomLidPhys->GetLogicalVolume()->GetSolid();
  G4VSolid* tankSolid   = fTankPhys->GetLogicalVolume()->GetSolid();

  G4UnionSolid* uniTankSolid = new G4UnionSolid("uniTank",tankSolid,toplidSolid,0,toplidPos);
  uniTankSolid = new G4UnionSolid("uniTank",uniTankSolid,btmlidSolid,0,btmlidPos);
  G4LogicalVolume* uniTankLog = new G4LogicalVolume(uniTankSolid,fTankMat,"uniTank",0,0,0);
  
  //not finished yet !!

  return;
}



*/




G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_tank_surf(){

  G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_tank",unified,ground,dielectric_metal); //why ground ?
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("LAr_tank",fLArColPhys,fTankPhys,opSurf);
  
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TANK,REFLECTIVITY_OF_TANK};
  //test
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};
  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test



  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;

}




G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_bottomLid_surf(){
  
  //why polished ?
  G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_bottomLid",glisur,polished,dielectric_metal);
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("LAr_bottomLid",fLArColPhys,fBottomLidPhys,opSurf);
  
//   G4double lambdaMin = 115. *nm; //from SVN; why this number ? 
//   G4double lambdaMax = 598. *nm; //from SVN; why this number ?
//   G4double photonMaxE = PLANCK_H*LIGHTSPEED/lambdaMin;
//   G4double photonMinE = PLANCK_H*LIGHTSPEED/lambdaMax;
//   G4double E[]={photonMinE,photonMaxE};
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_BOTTOM_LID,REFLECTIVITY_OF_BOTTOM_LID};


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;
}





G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_WLS_surf(G4String order){

  G4OpticalSurface* opSurf     = new G4OpticalSurface("LAr_WLS",unified,ground,dielectric_dielectric); 
  //G4OpticalSurface* opSurf     = new G4OpticalSurface("LAr_WLS",unified,polished,dielectric_dielectric); 
  G4double E[]={0.1*eV,13.5*eV};
  //G4double refIndex[]      = {1.5,1.5}; //test
  G4double reflectivity[]  = {1.,1.};   //test //no absorption
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  //propTab->AddProperty("RINDEX",E,refIndex,ne);
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  G4LogicalBorderSurface* borderSurf;
  
  if(!strcmp(order,"LAr_WLS"))
    borderSurf = new G4LogicalBorderSurface(order,fLArColPhys,fWLSPhys,opSurf);
  else if(!strcmp(order,"WLS_LAr"))
    borderSurf = new G4LogicalBorderSurface(order,fWLSPhys,fLArColPhys,opSurf);

  return borderSurf;
}


G4LogicalBorderSurface* ArDM_DetectorConstruction::build_GAr_tank_surf(){
  
  G4OpticalSurface*       opSurf     = new G4OpticalSurface("GAr_tank",unified,ground,dielectric_metal); //why ground ?
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("GAr_tank",fGArColPhys,fTankPhys,opSurf);
  
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TANK,REFLECTIVITY_OF_TANK};

  //test
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};
  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;
}




G4LogicalBorderSurface* ArDM_DetectorConstruction::build_GAr_topLid_surf(){
  
  G4OpticalSurface*       opSurf     = new G4OpticalSurface("GAr_topLid",glisur,polished,dielectric_metal); //why polished ?
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("GAr_topLid",fGArColPhys,fTopLidPhys,opSurf);
  
//   G4double lambdaMin = 115. *nm; //from SVN; why this number ? 
//   G4double lambdaMax = 598. *nm; //from SVN; why this number ?
//   G4double photonMaxE = PLANCK_H*LIGHTSPEED/lambdaMin;
//   G4double photonMinE = PLANCK_H*LIGHTSPEED/lambdaMax;
//   G4double E[]={photonMinE,photonMaxE};
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TOP_LID,REFLECTIVITY_OF_TOP_LID};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;
}



G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_GAr_surf(G4String order){
  //this method is not neccessary !
  //GEANT4 knows how to deal with gas-liquid surface !
  //if the material property tables for both are set !
  //why polished ??
  G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_GAr",unified,polished,dielectric_dielectric);
  G4LogicalBorderSurface* borderSurf;
  if(!strcmp(order,"LAr_GAr"))
    borderSurf = new G4LogicalBorderSurface(order,fLArColPhys,fGArColPhys,opSurf);
  else if(!strcmp(order,"GAr_LAr"))
    borderSurf = new G4LogicalBorderSurface(order,fGArColPhys,fLArColPhys,opSurf);


  G4double E[]={0.*eV,13.5*eV};
  //G4double reflectivity[] = {1.,1.}; //<-- this actually means "no absorption" ! transmission is still possible !!
  G4double reflectivity[] = {REFLECTIVITY_OF_LAR_SURFACE,REFLECTIVITY_OF_LAR_SURFACE};
  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  opSurf->SetMaterialPropertiesTable(propTab);
  return borderSurf;
}


vector<G4LogicalSkinSurface*> ArDM_DetectorConstruction::build_LAr_PMT_surf(){
  
  G4OpticalSurface* opSurf = new G4OpticalSurface("LAr_PMT",unified,groundbackpainted,dielectric_dielectric);
  vector<G4LogicalSkinSurface*> surfVec;
  for(int i=0;i<fBottomPMTArrayPhys.size();i++)
    surfVec.push_back(new G4LogicalSkinSurface("LAr_PMT",fBottomPMTArrayPhys[i]->GetLogicalVolume(),opSurf));
  
  G4double lambdaMin = 115. *nm; //from SVN; why this number ? 
  G4double lambdaMax = 598. *nm; //from SVN; why this number ?
  G4double photonMaxE = PLANCK_H*LIGHTSPEED/lambdaMin;
  G4double photonMinE = PLANCK_H*LIGHTSPEED/lambdaMax;
  G4double E[]={photonMinE,photonMaxE};
  G4double reflectivity[] = {0.,0.};
  //G4double reflectivity[] = {.0,.0};
  //G4double reflectivity[] = {1.0,1.0};//test !!!
  //G4double refIndex[] = {1.5,1.5};//test !!!
  G4double refIndex[] = {1.635,1.635};//test !!!

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("RINDEX",E,refIndex,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return surfVec;
}




vector<G4LogicalSkinSurface*> ArDM_DetectorConstruction::build_GAr_PMT_surf(){
  
  G4OpticalSurface* opSurf = new G4OpticalSurface("GAr_PMT",unified,groundbackpainted,dielectric_dielectric);
  vector<G4LogicalSkinSurface*> surfVec;
  for(int i=0;i<fBottomPMTArrayPhys.size();i++)
    surfVec.push_back(new G4LogicalSkinSurface("GAr_PMT",fTopPMTArrayPhys[i]->GetLogicalVolume(),opSurf));
  
  G4double lambdaMin = 115. *nm; //from SVN; why this number ? 
  G4double lambdaMax = 598. *nm; //from SVN; why this number ?
  G4double photonMaxE = PLANCK_H*LIGHTSPEED/lambdaMin;
  G4double photonMinE = PLANCK_H*LIGHTSPEED/lambdaMax;
  G4double E[]={photonMinE,photonMaxE};
  G4double reflectivity[] = {0.,0.};
  //G4double reflectivity[] = {1.0,1.0};//test !!!
  G4double refIndex[] = {1.5,1.5};//test !!!

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("RINDEX",E,refIndex,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return surfVec;
}





G4LogicalBorderSurface* ArDM_DetectorConstruction::build_WLS_WLSSupport_surf(){

  G4OpticalSurface*       opSurf     = new G4OpticalSurface("WLS_WLSSupport",unified,ground,dielectric_metal);

  ////dielectric_dielectric doesn't work <-- no reflection on reflector --> why !!?????
  ////this does work for LAr_WLS layer ! why doesn't it work for WLS_WLSSupport !!????? --> go figure !!
  ////--> maybe because the refractive index for teflon (<--> WLSSupport) is not available !
  ////--> for the time being, just use dielectric_metal type !
  //G4OpticalSurface*       opSurf     = new G4OpticalSurface("WLS_WLSSupport",unified,ground,dielectric_dielectric);

  opSurf->SetSigmaAlpha(MAIN_REFLECTOR_SIGMA_ALPHA*deg);


  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("WLS_WLSSupport",fWLSPhys,fWLSSupportPhys,opSurf);

  //test
  G4double Emin  = 0.1  *eV;
  G4double Emid1 = 6.0  *eV;
  G4double Emid2 = 6.0001  *eV;
  G4double Emax  = 13.5 *eV;

  //interval [0.1,6.0]*eV <--> +- 5*sigma (sigma=.5 eV) around 2.85 eV (<--> 420 nm)
  //reflectivity = 95% for [Emin,Emid], 0 for the rest.

  G4double E[]             = {Emin,Emid1,Emid2,Emax};
  G4double reflectivity[]  = {REFLECTIVITY_OF_REFLECTOR,REFLECTIVITY_OF_REFLECTOR,0,0};
  G4double specularspike[] = {MAIN_REFLECTOR_SPECULAR_SPIKE,MAIN_REFLECTOR_SPECULAR_SPIKE,MAIN_REFLECTOR_SPECULAR_SPIKE,MAIN_REFLECTOR_SPECULAR_SPIKE};
  G4double specularlobe[]  = {MAIN_REFLECTOR_SPECULAR_LOBE,MAIN_REFLECTOR_SPECULAR_LOBE,MAIN_REFLECTOR_SPECULAR_LOBE,MAIN_REFLECTOR_SPECULAR_LOBE};
  G4double backscatter[]   = {MAIN_REFLECTOR_BACK_SCATTER,MAIN_REFLECTOR_BACK_SCATTER,MAIN_REFLECTOR_BACK_SCATTER,MAIN_REFLECTOR_BACK_SCATTER};

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflReflector = REFLECTIVITY_OF_REFLECTOR;
// #endif //TEST_BRANCH0

  //the reflectivity given to the material properties table actually means (1 - absorption)
  //it will be used to determine whether the photon is absorbed or not
  //if not, further measures will be taken to determine whether the photon should undergo
  //
  //i.    total internal reflection
  //ii.   transmission
  //iii.  reflection <-- if reflection, determine the type of reflection by calling G4OpBoundaryProcess:ChooseReflection()
  //          a. specular spike
  //          b. specular lobe
  //          c. back scattering
  //          d. lambertian = 1 - a. - b. - c.
  //
  //see G4OpBoundaryProcess::PostStepDoIt(..) for more !


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);
  

  //end test


  return borderSurf;

}








G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_WLSSupport_surf(){

  G4OpticalSurface*       opSurf     = new G4OpticalSurface("WLS_WLSSupport",unified,ground,dielectric_metal);

  ////dielectric_dielectric doesn't work <-- no reflection on reflector --> why !!?????
  ////this does work for LAr_WLS layer ! why doesn't it work for WLS_WLSSupport !!????? --> go figure !!
  ////--> maybe because the refractive index for teflon (<--> WLSSupport) is not available !
  ////--> for the time being, just use dielectric_metal type !
  //G4OpticalSurface*       opSurf     = new G4OpticalSurface("WLS_WLSSupport",unified,ground,dielectric_dielectric);

  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("WLS_WLSSupport",fLArColPhys,fWLSSupportPhys,opSurf);

  //test
  G4double Emin  = 0.1  *eV;
  G4double Emid1 = 6.0  *eV;
  G4double Emid2 = 6.0001  *eV;
  G4double Emax  = 13.5 *eV;

  //interval [0.1,6.0]*eV <--> +- 5*sigma (sigma=.5 eV) around 2.85 eV (<--> 420 nm)
  //reflectivity = 95% for [Emin,Emid], 0 for the rest.

  G4double E[]             = {Emin,Emid1,Emid2,Emax};
  G4double reflectivity[]  = {0.,0.,0,0};


  if(!fWLSPhys){
    for(int i=0;i <= sizeof(reflectivity)/sizeof(G4double) / 2;i++){
      reflectivity[i] = REFLECTIVITY_OF_REFLECTOR;
    }
  }

  G4double specularspike[] = {0.,0.,0.,0.};
  G4double specularlobe[]  = {0.,0.,0.,0.};
  G4double backscatter[]   = {0.,0.,0.,0.};

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflReflector = REFLECTIVITY_OF_REFLECTOR;
// #endif //TEST_BRANCH0

  //the reflectivity given to the material properties table actually means (1 - absorption)
  //it will be used to determine whether the photon is absorbed or not
  //if not, further measures will be taken to determine whether the photon should undergo
  //
  //i.    total internal reflection
  //ii.   transmission
  //iii.  reflection <-- if reflection, determine the type of reflection by calling G4OpBoundaryProcess:ChooseReflection()
  //          a. specular spike
  //          b. specular lobe
  //          c. back scattering
  //          d. lambertian = 1 - a. - b. - c.
  //
  //see G4OpBoundaryProcess::PostStepDoIt(..) for more !


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);
  

  //end test


  return borderSurf;

}




vector<G4LogicalBorderSurface*> ArDM_DetectorConstruction::build_PMT_PMTCoat_surf(G4String topBottom,G4String order){

  G4OpticalSurface* opSurf = new G4OpticalSurface("PMT_PMTCoat",unified,ground,dielectric_dielectric);
  
//   G4double efficiency[] = {.18,.18};

  G4double E[]={.1*eV,6*eV,6.0001*eV,13.*eV}; 
  //[.1,6.]*eV ~ +- 5sigma (sigma = .5eV) around 2.85 eV (<--> 420 nm)
  //128 nm <--> 9.68 eV

  //G4double reflectivity[] = {1.,1.,0.,0.}; //<-- this actually means (1-absorption) !!
//   G4double reflectivity[] = {REFLECTIVITY_OF_PMT_FRONTSIDE,REFLECTIVITY_OF_PMT_FRONTSIDE,0,0};

//   if(!strcmp(order,"PMT_PMTCoat")) 
//     for(int i=0;i<sizeof(reflectivity)/sizeof(G4double);i++)
//       reflectivity[i] = 0.;

  G4double reflectivity[] = {1.,1.,0.,0.};

  //G4double efficiency[] = {.18,.18,0.,0.}; //quantum efficiency

  G4double specularspike[] = {0.,0.,0.,0.}; 
  G4double specularlobe[]  = {0.,0.,0.,0.};
  G4double backscatter[]   = {0.,0.,0.,0.};
  

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  //propTab->AddProperty("EFFICIENCY",E,efficiency,ne);//<-- check this again ! still not working !!!!
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);


  vector<G4LogicalBorderSurface*> pmt_pmtCoat_surf;
  vector<G4VPhysicalVolume*> medium1;
  vector<G4VPhysicalVolume*> medium2;

  if(!strcmp(topBottom,"bottom")){
    medium1 = fBottomPMTArrayPhys;
    medium2 = fBottomPMTCoatArrayPhys;
  }else if(!strcmp(topBottom,"top")){
    medium1 = fTopPMTArrayPhys;
    medium2 = fTopPMTCoatArrayPhys;
  }

  if(!strcmp(order,"PMTCoat_PMT")){
    vector<G4VPhysicalVolume*> medium = medium1;
    medium1 = medium2;
    medium2 = medium;
  }
    
  for(int i=0;i<medium1.size();i++)
    pmt_pmtCoat_surf.push_back(new G4LogicalBorderSurface(order,medium1[i],medium2[i],opSurf));

  return pmt_pmtCoat_surf;
}



vector<G4LogicalBorderSurface*> ArDM_DetectorConstruction::build_Ar_PMTCoat_surf(G4String topBottom,G4String order){


  G4OpticalSurface* opSurf = new G4OpticalSurface("Ar_PMTCoat",unified,ground,dielectric_dielectric);
  
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_PMTCOAT,REFLECTIVITY_OF_PMTCOAT};
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};
  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  //propTab->AddProperty("RINDEX",E,refIndex,ne); //<-- surfaceFinish = ground -->refIndex of the surface irrelevant !
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);



  G4VPhysicalVolume* mediumPhys;
  vector<G4VPhysicalVolume*> pmtCoatArrayPhys;
  ostringstream surfName;
  if(!strcmp(topBottom,"bottom")){
    mediumPhys       = fLArColPhys;
    pmtCoatArrayPhys = fBottomPMTCoatArrayPhys;
    surfName<<"LAr";
  }
  else if(!strcmp(topBottom,"top")){
    mediumPhys       = fGArColPhys;
    pmtCoatArrayPhys = fTopPMTCoatArrayPhys;
    surfName<<"GAr";
  }
  

  vector<G4LogicalBorderSurface*> ar_pmtCoat_surf;
  if(!strcmp(order,"Ar_PMTCoat")){

    surfName<<"_PMTCoat";
    for(int i=0;i<pmtCoatArrayPhys.size();i++)
      ar_pmtCoat_surf.push_back(new G4LogicalBorderSurface(surfName.str().c_str(),mediumPhys,
							   pmtCoatArrayPhys[i],opSurf));
  }else if(!strcmp(order,"PMTCoat_Ar")){
    surfName.str("PMTCoat_"+surfName.str());
    for(int i=0;i<pmtCoatArrayPhys.size();i++)
      ar_pmtCoat_surf.push_back(new G4LogicalBorderSurface(surfName.str().c_str(),pmtCoatArrayPhys[i],
							   mediumPhys,opSurf));
  }

  return ar_pmtCoat_surf;
}



vector<G4LogicalBorderSurface*> ArDM_DetectorConstruction::build_Ar_PMT_surf(G4String topBottom,G4String order){


  G4OpticalSurface* opSurf = new G4OpticalSurface("Ar_PMT",unified,ground,dielectric_dielectric);
  //G4double refIndex[]={1.5,1.5}; //test !!
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_PMT_BACKSIDE,REFLECTIVITY_OF_PMT_BACKSIDE};
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};
  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);



  G4VPhysicalVolume* mediumPhys;
  vector<G4VPhysicalVolume*> pmtArrayPhys;
  ostringstream surfName;
  if(!strcmp(topBottom,"bottom")){
    mediumPhys   = fLArColPhys;
    pmtArrayPhys = fBottomPMTArrayPhys;
    surfName<<"LAr";
  }
  else if(!strcmp(topBottom,"top")){
    mediumPhys   = fGArColPhys;
    pmtArrayPhys = fTopPMTArrayPhys;
    surfName<<"GAr";
  }
  

  vector<G4LogicalBorderSurface*> ar_pmt_surf;
  if(!strcmp(order,"Ar_PMT")){
    surfName<<"_PMT";
    for(int i=0;i< pmtArrayPhys.size();i++)
      ar_pmt_surf.push_back(new G4LogicalBorderSurface(surfName.str().c_str(),mediumPhys,
						       pmtArrayPhys[i],opSurf));

  }else if(!strcmp(order,"PMT_Ar")){
    surfName.str("PMT_"+surfName.str());
    for(int i=0;i<pmtArrayPhys.size();i++)
      ar_pmt_surf.push_back(new G4LogicalBorderSurface(surfName.str().c_str(),pmtArrayPhys[i],
						       mediumPhys,opSurf));
  }

  return ar_pmt_surf;
}







vector<G4LogicalBorderSurface*> ArDM_DetectorConstruction::build_PMTCathode_PMT_surf(G4String topBottom,G4String order){


  G4OpticalSurface* opSurf = new G4OpticalSurface("PMTCathode_PMT",unified,ground,dielectric_dielectric);
  //G4double refIndex[]={1.5,1.5}; //test !!
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_PMT_BACKSIDE,REFLECTIVITY_OF_PMT_BACKSIDE};
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};
  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);


  vector<G4LogicalBorderSurface*> pmtCathode_pmt_surf;
  vector<G4VPhysicalVolume*> medium1;
  vector<G4VPhysicalVolume*> medium2;

  if(!strcmp(topBottom,"bottom")){
    medium1 = fBottomPMTArrayPhys;
    medium2 = fBottomPMTCathodeArrayPhys;
  }else if(!strcmp(topBottom,"top")){
    medium1 = fTopPMTArrayPhys;
    medium2 = fTopPMTCathodeArrayPhys;
  }


  if(!strcmp(order,"PMTCathode_PMT")){
    vector<G4VPhysicalVolume*> medium = medium1;
    medium1 = medium2;
    medium2 = medium;
  }
    
   
  for(int i=0;i<medium1.size();i++)
    pmtCathode_pmt_surf.push_back(new G4LogicalBorderSurface(order,medium1[i],medium2[i],opSurf));

  return pmtCathode_pmt_surf;
}









G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_BottomSideReflector_surf(){

  G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_BottomSideReflector",
							    unified,ground,dielectric_metal); //why ground ?

//   G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_BottomSideReflector",
// 							    unified,ground,dielectric_dielectric); //why ground ?


  opSurf->SetSigmaAlpha(BOTTOM_SIDE_REFLECTOR_SIGMA_ALPHA*deg);

  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("LAr_BottomSideReflector",
								  fLArColPhys,fBottomSideReflectorPhys,opSurf);
  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_BOTTOM_SIDE_REFLECTOR,REFLECTIVITY_OF_BOTTOM_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {BOTTOM_SIDE_REFLECTOR_SPECULAR_SPIKE,BOTTOM_SIDE_REFLECTOR_SPECULAR_SPIKE,
			      BOTTOM_SIDE_REFLECTOR_SPECULAR_SPIKE,BOTTOM_SIDE_REFLECTOR_SPECULAR_SPIKE};

  G4double specularlobe[]  = {BOTTOM_SIDE_REFLECTOR_SPECULAR_LOBE,BOTTOM_SIDE_REFLECTOR_SPECULAR_LOBE,
			      BOTTOM_SIDE_REFLECTOR_SPECULAR_LOBE,BOTTOM_SIDE_REFLECTOR_SPECULAR_LOBE};

  G4double backscatter[]   = {BOTTOM_SIDE_REFLECTOR_BACK_SCATTER,BOTTOM_SIDE_REFLECTOR_BACK_SCATTER,
			      BOTTOM_SIDE_REFLECTOR_BACK_SCATTER,BOTTOM_SIDE_REFLECTOR_BACK_SCATTER};

  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflSideReflectorBottom = REFLECTIVITY_OF_BOTTOM_SIDE_REFLECTOR;
// #endif //TEST_BRANCH0

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;

}
















G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_TopSideReflector_surf(){

  //why ground ?
  //why dielectric_metal ??

  G4OpticalSurface* opSurf = new G4OpticalSurface("LAr_TopSideReflector",unified,ground,dielectric_metal); 


  ////dielectric_dielectric somehow doesn't work !?? --> go figure !!!
//   G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_BottomSideReflector",
// 							    unified,ground,dielectric_dielectric); //why ground ?


  opSurf->SetSigmaAlpha(TOP_SIDE_REFLECTOR_SIGMA_ALPHA*deg);

  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("LAr_TopSideReflector",
								  fLArColPhys,fTopSideReflectorLArPhys,opSurf);
  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,
			      TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE};

  G4double specularlobe[]  = {TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE,
			      TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE};

  G4double backscatter[]   = {TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER,
			      TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER};

  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflSideReflectorTop = REFLECTIVITY_OF_TOP_SIDE_REFLECTOR;
// #endif //TEST_BRANCH0

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;

}












G4LogicalBorderSurface* ArDM_DetectorConstruction::build_GAr_TopSideReflector_surf(){

  //why ground ?
  //why dielectric_metal ??

  G4OpticalSurface* opSurf = new G4OpticalSurface("GAr_TopSideReflector",unified,ground,dielectric_metal); 


  ////dielectric_dielectric somehow doesn't work !?? --> go figure !!!
//   G4OpticalSurface*       opSurf     = new G4OpticalSurface("LAr_BottomSideReflector",
// 							    unified,ground,dielectric_dielectric); //why ground ?

  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("GAr_TopSideReflector",
								  fGArColPhys,fTopSideReflectorGArPhys,opSurf);
  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,
			      TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE};
  
  G4double specularlobe[]  = {TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE,
			      TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE};

  G4double backscatter[]   = {TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER,
			      TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER};

  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflSideReflectorTop = REFLECTIVITY_OF_TOP_SIDE_REFLECTOR;
// #endif //TEST_BRANCH0

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;

}





////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////








void ArDM_DetectorConstruction::build_surfaces(){

  //100% absorption when photons fall onto the wall of the tank or onto the lids
  fLAr_tank_surf          = build_LAr_tank_surf();     
  fLAr_bottomLid_surf     = build_LAr_bottomLid_surf();
  fGAr_tank_surf          = build_GAr_tank_surf();
  fGAr_topLid_surf        = build_GAr_topLid_surf();


  //reflector, (1 - absorption) = 95%
  fWLS_WLSSupport_surf      = build_WLS_WLSSupport_surf();
  fLAr_WLSSupport_surf      = build_LAr_WLSSupport_surf();
  

  fLAr_BottomSideReflector_surf = build_LAr_BottomSideReflector_surf();

#if TURN_ON_SIDE_REFLECTOR_COATING
  fLAr_BottomSideReflectorCoat_surf = build_LAr_BottomSideReflectorCoat_surf("LAr_BottomSideReflectorCoat");
  fBottomSideReflectorCoat_LAr_surf = build_LAr_BottomSideReflectorCoat_surf("BottomSideReflectorCoat_LAr");

  fBottomSideReflector_BottomSideReflectorCoat_surf = build_BottomSideReflector_BottomSideReflectorCoat_surf();

  fGAr_TopSideReflectorCoat_surf = build_GAr_TopSideReflectorCoat_surf("GAr_TopSideReflectorCoat");
  fTopSideReflectorCoat_GAr_surf = build_GAr_TopSideReflectorCoat_surf("TopSideReflectorCoat_GAr");

  fTopSideReflectorGAr_TopSideReflectorCoat_surf = build_TopSideReflectorGAr_TopSideReflectorCoat_surf();

  fLAr_TopSideReflectorCoat_surf = build_LAr_TopSideReflectorCoat_surf("LAr_TopSideReflectorCoat");
  fTopSideReflectorCoat_LAr_surf = build_LAr_TopSideReflectorCoat_surf("TopSideReflectorCoat_LAr");

  fTopSideReflectorLAr_TopSideReflectorCoat_surf = build_TopSideReflectorLAr_TopSideReflectorCoat_surf();


#else 
  fLAr_TopSideReflector_surf    = build_LAr_TopSideReflector_surf();
  fGAr_TopSideReflector_surf    = build_GAr_TopSideReflector_surf();
 

#endif //TURN_ON_SIDE_REFLECTOR_COATING


  fLAr_WLS_surf           = build_LAr_WLS_surf("LAr_WLS");
  fWLS_LAr_surf           = build_LAr_WLS_surf("WLS_LAr");


  //PMTs' optical surface  

  //PMT -- PMTCoat
  //if a photon falls onto the PMT surface, it'll be absorbed (killed) regardless its energy
  //--> set reflectivity to 0. for all wavelengths
  fBottomPMT_PMTCoat_surf = build_PMT_PMTCoat_surf("bottom","PMT_PMTCoat");
  fTopPMT_PMTCoat_surf    = build_PMT_PMTCoat_surf("top","PMT_PMTCoat");
  fPMTCoat_BottomPMT_surf = build_PMT_PMTCoat_surf("bottom","PMTCoat_PMT");
  fPMTCoat_TopPMT_surf    = build_PMT_PMTCoat_surf("top","PMTCoat_PMT");

 
  // if(fBottomPMTCathodeArrayPhys.size()){ //if pmt cathode is implemented 
  //   fPMTCathode_bottomPMT_surf = build_PMTCathode_PMT_surf("bottom","PMTCathode_PMT");
  //   fBottomPMT_PMTCathode_surf = build_PMTCathode_PMT_surf("bottom","PMT_PMTCathode");

  // }else{
  //   //prevent total reflection within the glass-layers of the PMTs
  //   //PMT -- LAr/GAr
  //   //100% absorption
  //   fLAr_PMT_surf = build_Ar_PMT_surf("bottom","Ar_PMT");
  //   fPMT_LAr_surf = build_Ar_PMT_surf("bottom","PMT_Ar");
  // }


  
  // if(fTopPMTCathodeArrayPhys.size()){ //if pmt cathode is implemented 
  //   fPMTCathode_topPMT_surf    = build_PMTCathode_PMT_surf("top"   ,"PMTCathode_PMT");
  //   fTopPMT_PMTCathode_surf    = build_PMTCathode_PMT_surf("top"   ,"PMT_PMTCathode");

  // }else{
  //   //prevent total reflection within the glass-layers of the PMTs
  //   //PMT -- LAr/GAr
  //   //100% absorption
  //   fGAr_PMT_surf = build_Ar_PMT_surf("top","Ar_PMT");
  //   fPMT_GAr_surf = build_Ar_PMT_surf("top","PMT_Ar");
  // }


  //replace the 2 if-clauses above by the following block
  //on 2014, may 12th , 19:30

  if(1){
  fPMTCathode_bottomPMT_surf = build_PMTCathode_PMT_surf("bottom","PMTCathode_PMT");
  fBottomPMT_PMTCathode_surf = build_PMTCathode_PMT_surf("bottom","PMT_PMTCathode");

  fLAr_PMT_surf = build_Ar_PMT_surf("bottom","Ar_PMT");
  fPMT_LAr_surf = build_Ar_PMT_surf("bottom","PMT_Ar");

  fPMTCathode_topPMT_surf    = build_PMTCathode_PMT_surf("top"   ,"PMTCathode_PMT");
  fTopPMT_PMTCathode_surf    = build_PMTCathode_PMT_surf("top"   ,"PMT_PMTCathode");

  fGAr_PMT_surf = build_Ar_PMT_surf("top","Ar_PMT");
  fPMT_GAr_surf = build_Ar_PMT_surf("top","PMT_Ar");
  }

  //PMTCoat -- LAr/GAr
  fLAr_bottomPMTCoat_surf = build_Ar_PMTCoat_surf("bottom","Ar_PMTCoat");
  fGAr_topPMTCoat_surf    = build_Ar_PMTCoat_surf("top","Ar_PMTCoat");
  fBottomPMTCoat_LAr_surf = build_Ar_PMTCoat_surf("bottom","PMTCoat_Ar");
  fTopPMTCoat_GAr_surf    = build_Ar_PMTCoat_surf("top","PMTCoat_Ar");


  //interface liquid -- gas
  fLAr_GAr_surf = build_LAr_GAr_surf("LAr_GAr");
  fGAr_LAr_surf = build_LAr_GAr_surf("GAr_LAr");

  //LAr/GAr -- PMTSupport
  //<-- not really needed !
  //it's just the "hidden side" of the PMT support, which is highly likely not seen by optical photons
  fGAr_topPMTSupport_surf    = build_Ar_PMTSupport_surf("top");
  fLAr_bottomPMTSupport_surf = build_Ar_PMTSupport_surf("bottom");

  //PMTSupport -- PMTSupportCoat
  fPMTSupportCoat_topPMTSupport_surf    = build_PMTSupportCoat_PMTSupport_surf("top");
  fPMTSupportCoat_bottomPMTSupport_surf = build_PMTSupportCoat_PMTSupport_surf("bottom");

  //LAr/GAr -- PMTSupportCoat

  fTopPMTSupportCoat_GAr_surf    = build_Ar_PMTSupportCoat_surf("top"   ,"PMTSupportCoat_Ar");
  fBottomPMTSupportCoat_LAr_surf = build_Ar_PMTSupportCoat_surf("bottom","PMTSupportCoat_Ar");

  fGAr_topPMTSupportCoat_surf    = build_Ar_PMTSupportCoat_surf("top"   ,"Ar_PMTSupportCoat");
  fLAr_bottomPMTSupportCoat_surf = build_Ar_PMTSupportCoat_surf("bottom","Ar_PMTSupportCoat");



  if(1){
  //protection grid <--> LAr                                                                                                                          
  fLAr_protectionGrid_surf = build_LAr_protectionGrid_surf();

  //cathode grid <--> LAr                                                                                                                             
  fLAr_cathodeGrid_surf = build_LAr_cathodeGrid_surf();
  }

  return;
}













////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////





G4Material* ArDM_DetectorConstruction::getTetratex(){
  //not yet finished !
  //tetratex = CF2

  G4Element* C = new G4Element("Carbon"  , "C" , 6., 12.0107 *g/mole);
  G4Element* H = new G4Element("Hydrogen", "H" , 1.,  1.0079 *g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O" , 8., 15.9994 *g/mole);

  G4double density = POLYCARBONATE_DENSITY;
  G4Material* tetratex = new G4Material("Polycarbonate", density,3);
  tetratex->AddElement(C,26);
  tetratex->AddElement(H,14);
  tetratex->AddElement(O,3) ;


  return tetratex;
}






G4Material* ArDM_DetectorConstruction::getPolyethylene(){

//   G4Element* C = new G4Element("Carbon"  , "C" , 6., 12.0107 *g/mole);
//   G4Element* H = new G4Element("Hydrogen", "H" , 1.,  1.0079 *g/mole);

//   G4double density = POLYCARBONATE_DENSITY;
//   G4Material* polyethylene = new G4Material("Polyethylene", density,2);
//   polyethylene->AddElement(C,1);
//   polyethylene->AddElement(H,2);

  G4Material* polyethylene = fNist->FindOrBuildMaterial("G4_POLYETHYLENE");

  //cout<<"density [g/cm3] "<<polyethylene->GetDensity()/(g/cm3)<<endl; getchar();

  return polyethylene;
}


G4Material* ArDM_DetectorConstruction::getTeflon(){

  G4Material* teflon = fNist->FindOrBuildMaterial("G4_TEFLON");

  return teflon;
}







G4Material* ArDM_DetectorConstruction::getPerlite(){

  //according to wikipedia http://en.wikipedia.org/wiki/Perlite
  //perlite is to ~70% composed of Si02
  //so for the time being (2014, feb. 19th), 
  //just approximate perlite ~= Si02
  //but set the density to the actual density of perlite

  G4Element* Si = new G4Element("Silicon", "Si", 14., 28.085  *g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O" ,  8., 15.9994 *g/mole);

  G4double density = 1.1 * g/cm3;
  int nelements = 2;

  G4Material* perlite = new G4Material("perlite",density,nelements);
  perlite->AddElement(Si,1);
  perlite->AddElement(O,2);

  return perlite;
}






G4Material* ArDM_DetectorConstruction::getHVrMat(){

  //assuming that HVr is made of ceramic,
  //which is basically Al2O3

  G4Element* Al = new G4Element("Aluminium", "Al", 13., 26.9815 *g/mole);
  G4Element*  O = new G4Element("Oxygen"   ,  "O",  8., 15.9994 *g/mole);

  //G4double density = 101.96 *g/mole;
  G4double density = 3.95* g/cm3;
  G4Material* Al2O3 = new G4Material("Al2O3",density,2);
  Al2O3->AddElement(Al,2);
  Al2O3->AddElement(O ,3);
  
  return Al2O3;
}







G4Material* ArDM_DetectorConstruction::getBorotron(){

  G4Material* borotron = NULL;

  G4int elemental_percentage=1;

  if(elemental_percentage){
    
    G4Element* B = new G4Element("Boron"  , "B",  5., 10.811  *g/mole);
    G4Element* O = new G4Element("Oxygen" , "O",  8., 15.9994 *g/mole);
    G4Element* H = new G4Element("Hydrogen", "H", 1.,  1.0079 *g/mole); 
    G4Element* C = new G4Element("Carbon"  , "C", 6., 12.0107  *g/mole);


    G4double density = 1.005 *g/cm3;
    G4int nelements  = 4;
    
    //no oxygen in borotron ??

    borotron = new G4Material("borotron",density,nelements);
    borotron->AddElement(C, 71.84/100);
    borotron->AddElement(H, 12.06/100);
    borotron->AddElement(B,  5.00/100);
    borotron->AddElement(O, 11.10/100);

  }else{
    
    G4Material* polyethylene = fPolyethylene;
    
    G4Element* B = new G4Element("Boron"  , "B",  5., 10.811  *g/mole);
    G4Element* O = new G4Element("Oxygen" , "O",  8., 15.9994 *g/mole);

    G4double B2O3_density = 69.6182 *g/mole;
    G4Material* B2O3 = new G4Material("borontrioxide", B2O3_density,2);
    B2O3->AddElement(B, 2);
    B2O3->AddElement(O, 3);
    
    G4double borotron_density = 1.005 *g/cm3;
    G4int ncomponents=2;
    borotron = new G4Material("borotron",borotron_density,ncomponents);
    borotron->AddMaterial(polyethylene, 83.9 /100);
    borotron->AddMaterial(B2O3        , 16.1 /100);
  }

  return borotron;
}






G4Material* ArDM_DetectorConstruction::getFR4(){

  //material of pmt base = FR4
  //"typically" :
  //FR4 = 60% fiberglass + 40% epoxy

  //epoxy : C21H25ClO5
  //fiberglass : table 1 in http://www.agy.com/technical_info/graphics_PDFs/HighStrengthTechPaperEng.pdf
  //fiberglass "typically" consists of :
  //


  G4Element* O  = new G4Element("Oxygen"   ,  "O",  8., 15.9994 *g/mole);
  G4Element* H  = new G4Element("Hydrogen" ,  "H",  1.,  1.0079 *g/mole); 
  G4Element* C  = new G4Element("Carbon"   ,  "C",  6., 12.0107 *g/mole);
  G4Element* Cl = new G4Element("Chlorine" , "Cl", 17., 35.4530 *g/mole);
  G4Element* Si = new G4Element("Silicon"  , "Si", 14., 28.0900 *g/mole);   
  G4Element* Al = new G4Element("Aluminium", "Al", 13., 26.9815 *g/mole);
  G4Element* B  = new G4Element("Boron"    ,  "B",  5., 10.811  *g/mole);
  G4Element* Na = new G4Element("Natrium"  , "Na", 11., 22.9897 *g/mole);
  G4Element* K  = new G4Element("Kalium"   ,  "K", 19., 39.0983 *g/mole);
  G4Element* Ca = new G4Element("Calcium"  , "Ca", 20., 40.0780 *g/mole);
  G4Element* Mg = new G4Element("Magnesium", "Mg", 12., 24.3050 *g/mole);


 
  G4double epoxy_fraction = 40./100;
  G4double fb_fraction    = 60./100;

  //contribution of elements in epoxy :
  //epoxy = C21_H25_Cl_O5
  //--> whole molecule : 392.5u
  G4double epoxy_mass = 392.5; //u
  G4double epoxy_C_fr = 21 * 12.0 /epoxy_mass;
  G4double epoxy_H_fr = 25 *  1.0 /epoxy_mass;
  G4double epoxy_Cl_fr= 1  * 35.5 /epoxy_mass;
  G4double epoxy_O_fr = 5  * 16.0 /epoxy_mass;


  //contribution of elements in fiberglass :
  //"typical" fiberglass composition : in percent

  
  G4double fr_SiO2  = 68. /100;
  G4double fr_Al2O3 = 4. /100;
  G4double fr_B2O3  = 4. /100;
  G4double fr_CaO   = 8. /100;
  G4double fr_MgO   = 2. /100;
  G4double fr_Na2O  = 7. /100;
  G4double fr_K2O   = 7. /100;


  //molecular mass of the components in fiber glass : in [u]
  G4double m_SiO2  = 60;
  G4double m_Al2O3 = 102;
  G4double m_B2O3  = 70;
  G4double m_CaO   = 56;
  G4double m_MgO   = 40;
  G4double m_Na2O  = 62;
  G4double m_K2O   = 94;
  
  //fb = fiber glass
  G4double fb_Si_fr = 1. * 28 / m_SiO2  * fr_SiO2  ; 
  G4double fb_Al_fr = 2. * 27 / m_Al2O3 * fr_Al2O3 ;
  G4double fb_B_fr  = 2. * 11 / m_B2O3  * fr_B2O3  ;
  G4double fb_Ca_fr = 1. * 40 / m_CaO   * fr_CaO   ;
  G4double fb_Mg_fr = 1. * 24 / m_MgO   * fr_MgO   ;
  G4double fb_Na_fr = 2. * 23 / m_Na2O  * fr_Na2O  ;
  G4double fb_K_fr  = 2. * 39 / m_K2O   * fr_K2O   ;

  G4double fb_O_fr  = 16*(2. / m_SiO2  * fr_SiO2  +
			  3. / m_Al2O3 * fr_Al2O3 +
			  3. / m_B2O3  * fr_B2O3  +
			  1. / m_CaO   * fr_CaO   +
			  1. / m_MgO   * fr_MgO   +
			  1. / m_Na2O  * fr_Na2O  +
			  1. / m_K2O   * fr_K2O    );



  G4double fr_Si = fb_Si_fr * fb_fraction;
  G4double fr_Al = fb_Al_fr * fb_fraction;
  G4double fr_B  = fb_B_fr  * fb_fraction;
  G4double fr_Ca = fb_Ca_fr * fb_fraction;
  G4double fr_Mg = fb_Mg_fr * fb_fraction;
  G4double fr_Na = fb_Na_fr * fb_fraction;
  G4double fr_K  = fb_K_fr  * fb_fraction;

  G4double fr_C  = epoxy_C_fr * epoxy_fraction;
  G4double fr_H  = epoxy_H_fr * epoxy_fraction;
  G4double fr_Cl = epoxy_Cl_fr* epoxy_fraction;

  G4double fr_O  = epoxy_O_fr * epoxy_fraction + fb_O_fr * fb_fraction;

  G4double sum = fr_Si + fr_Al + fr_B + fr_Ca + fr_Mg + fr_Na + fr_K + fr_C + fr_H + fr_Cl + fr_O;
                  
  G4double density   = 1.850 *g/cm3;
  G4int    nelements = 11;
  G4Material* FR4 = new G4Material("FR4",density,nelements);

  FR4->AddElement(Si,fr_Si);
  FR4->AddElement(Al,fr_Al);
  FR4->AddElement(B ,fr_B );
  FR4->AddElement(Ca,fr_Ca);
  FR4->AddElement(Mg,fr_Mg);
  FR4->AddElement(Na,fr_Na);
  FR4->AddElement(K ,fr_K );
  FR4->AddElement(C ,fr_C );
  FR4->AddElement(H ,fr_H );
  FR4->AddElement(Cl,fr_Cl);
  FR4->AddElement(O ,fr_O );

  return FR4;
}


void ArDM_DetectorConstruction::addPMTSupport(G4String bottom_or_top){

  G4double pmtSupport_innerR = 0;
  G4double pmtSupport_outerR = PMT_SUPPORT_RADIUS;
  G4double pmtSupport_half_height = PMT_SUPPORT_HALF_HEIGHT;

  //G4Tubs* pmtSupportPlateSolid = new G4Tubs("pmtSupportPlate",pmtSupport_innerR,pmtSupport_outerR,0*deg,360.*deg);
  G4Tubs* pmtSupportPlateSolid = new G4Tubs("pmtSupportPlate",pmtSupport_innerR,pmtSupport_outerR,
					    pmtSupport_half_height,0*deg,360.*deg);

  //making holes for PMTs on the PMT support structure

  //size of the hole

  G4double pmtHole_innerR = 0;
  //G4double pmtHole_outerR = PMT_OUTER_RADIUS*sin(PMT_ACTIVE_RANGE);
  G4double pmtHole_outerR = APPROX_PMT_MIDDLE_CYLINDER_OUTER_RADIUS;
  G4double pmtHole_half_height = pmtSupport_half_height + 0.1*mm;
  //0.1*mm is just there so that we actually see a hole in the plate when calling G4SubtractionSolid

  G4Tubs* pmtHoleSolid = new G4Tubs("pmtHole",pmtHole_innerR,pmtHole_outerR,pmtHole_half_height,0*deg,360.*deg);


  vector<TVector2> rpmt2D = setPMTVector2D();
  
  G4BooleanSolid* pmtSupportSolid;
  
  pmtSupportSolid = new G4SubtractionSolid("pmtSupportSolid",pmtSupportPlateSolid,pmtHoleSolid,0,
					   G4ThreeVector(rpmt2D[0].X(),rpmt2D[0].Y(),0));

  for(int i=1;i<rpmt2D.size();i++)    
    pmtSupportSolid = new G4SubtractionSolid("pmtSupportSolid",pmtSupportSolid,pmtHoleSolid,0,
					     G4ThreeVector(rpmt2D[i].X(),rpmt2D[i].Y(),0));
  

  
  ostringstream pmtSupportName;
  pmtSupportName<<bottom_or_top<<"PMTSupport";
  G4LogicalVolume* fPMTSupportLog = new G4LogicalVolume(pmtSupportSolid,fPMTSupportMat,pmtSupportName.str().c_str());
  
  
  G4VisAttributes* pmtSupportAtt = new G4VisAttributes(true);
  //pmtSupportAtt->SetColour(0.8,.2,0.3); //dunno what colour this is ! who cares !
  pmtSupportAtt->SetColour(G4Colour::Yellow()); //dunno what colour this is ! who cares !
  pmtSupportAtt->SetForceAuxEdgeVisible(true);
  fPMTSupportLog->SetVisAttributes(pmtSupportAtt);
  
  
  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = TOP_PMT_SUPPORT_POS_Z_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = TOP_PMT_SUPPORT_POS_Z;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = BOTTOM_PMT_SUPPORT_POS_Z_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = BOTTOM_PMT_SUPPORT_POS_Z;       }
  }

  G4ThreeVector pmtSupportPos(BOTTOM_PMT_SUPPORT_POS_X,BOTTOM_PMT_SUPPORT_POS_Y,z);

  G4VPhysicalVolume* fPMTSupportPhys = new G4PVPlacement(0,pmtSupportPos,pmtSupportName.str().c_str(),
							 fPMTSupportLog,fMotherPhys,false,0);

  if(bottom_or_top == "bottom")   fBottomPMTSupportPhys = fPMTSupportPhys;
  else if(bottom_or_top == "top") fTopPMTSupportPhys    = fPMTSupportPhys;


  addPMTSupportCoat(bottom_or_top);

  return;
}






void ArDM_DetectorConstruction::addPMTSupportCoat(G4String bottom_or_top){

  G4double pmtSupportCoat_innerR = 0;
  G4double pmtSupportCoat_outerR = PMT_SUPPORT_COATING_RADIUS;
  G4double pmtSupportCoat_half_height;
  
  if(bottom_or_top == "top") pmtSupportCoat_half_height = TOP_PMT_SUPPORT_COATING_HALF_HEIGHT;
  else pmtSupportCoat_half_height = BOTTOM_PMT_SUPPORT_COATING_HALF_HEIGHT;

  G4Tubs* pmtSupportCoatSolid = new G4Tubs("pmtSupportCoat",pmtSupportCoat_innerR,pmtSupportCoat_outerR,
					   pmtSupportCoat_half_height,0*deg,360.*deg);

  //making holes for PMTs on the PMT support structure

  //size of the hole

  G4double pmtHole_innerR = 0;
  G4double pmtHole_outerR;
  
  if(bottom_or_top == "top") pmtHole_outerR = TOP_PMT_COATING_OUTER_RADIUS*sin(APPROX_PMT_ACTIVE_RANGE_OPENING_ANGLE/2);
  else pmtHole_outerR = BOTTOM_PMT_COATING_OUTER_RADIUS*sin(APPROX_PMT_ACTIVE_RANGE_OPENING_ANGLE/2);

  G4double pmtHole_half_height = pmtSupportCoat_half_height + 0.01*mm;
  //0.01*mm is just there so that we actually see a hole in the plate when calling G4SubtractionSolid


  G4Tubs* pmtHoleSolid = new G4Tubs("pmtHole",pmtHole_innerR,pmtHole_outerR,pmtHole_half_height,0*deg,360.*deg);


  vector<TVector2> rpmt2D = setPMTVector2D();
  
  G4BooleanSolid* pmtSupportSolid;
  
  pmtSupportSolid = new G4SubtractionSolid("pmtSupportCoatSolid",pmtSupportCoatSolid,pmtHoleSolid,0,
					   G4ThreeVector(rpmt2D[0].X(),rpmt2D[0].Y(),0));

  for(int i=1;i<rpmt2D.size();i++)    
    pmtSupportSolid = new G4SubtractionSolid("pmtSupportCoatSolid",pmtSupportSolid,pmtHoleSolid,0,
					     G4ThreeVector(rpmt2D[i].X(),rpmt2D[i].Y(),0));
  

  ostringstream pmtSupportCoatName;
  pmtSupportCoatName<<bottom_or_top<<"PMTSupportCoat";
  G4LogicalVolume* fPMTSupportCoatLog = new G4LogicalVolume(pmtSupportSolid,fWLSMat,pmtSupportCoatName.str().c_str());
  
  G4VisAttributes* pmtSupportCoatAtt = new G4VisAttributes(true);
  //pmtSupportCoatAtt->SetColour(0.8,.2,0.3); //dunno what colour this is ! who cares !
  //pmtSupportCoatAtt->SetColour(G4Colour::Yellow()); 
  pmtSupportCoatAtt->SetColour(G4Colour::Green()); 
  pmtSupportCoatAtt->SetForceAuxEdgeVisible(true);
  fPMTSupportCoatLog->SetVisAttributes(pmtSupportCoatAtt);
  
  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(bottom_or_top == "top"){
    if(fGArColPhys){ fMotherPhys = fGArColPhys; z = TOP_PMT_SUPPORT_COATING_POS_Z_GARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = TOP_PMT_SUPPORT_COATING_POS_Z;       }
  }else{
    if(fLArColPhys){ fMotherPhys = fLArColPhys; z = BOTTOM_PMT_SUPPORT_COATING_POS_Z_LARCOL;}
    else{            fMotherPhys = fWorldPhys;  z = BOTTOM_PMT_SUPPORT_COATING_POS_Z;       }
  }

  G4ThreeVector pmtSupportCoatPos(BOTTOM_PMT_SUPPORT_COATING_POS_X,
				  BOTTOM_PMT_SUPPORT_COATING_POS_Y,
				  z);

  G4VPhysicalVolume* fPMTSupportCoatPhys = new G4PVPlacement(0,pmtSupportCoatPos,pmtSupportCoatName.str().c_str(),
							     fPMTSupportCoatLog,fMotherPhys,false,0);


  if(bottom_or_top == "bottom")   fBottomPMTSupportCoatPhys = fPMTSupportCoatPhys;
  else if(bottom_or_top == "top") fTopPMTSupportCoatPhys    = fPMTSupportCoatPhys;


  return;
}





G4LogicalBorderSurface* ArDM_DetectorConstruction::build_Ar_PMTSupport_surf(G4String topBottom){


  G4VPhysicalVolume* mediumPhys;
  G4VPhysicalVolume* pmtSupportPhys;

  if(!strcmp(topBottom,"top")){
    mediumPhys     = fGArColPhys;
    pmtSupportPhys = fTopPMTSupportPhys;
  }else{
    mediumPhys     = fLArColPhys;
    pmtSupportPhys = fBottomPMTSupportPhys;    
  }


  //G4OpticalSurface*       opSurf     = new G4OpticalSurface("Ar_PMTSupport",unified,ground,dielectric_metal); 
  G4OpticalSurface*       opSurf     = new G4OpticalSurface("Ar_PMTSupport",unified,polished,dielectric_metal); 
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("Ar_PMTSupport",mediumPhys,pmtSupportPhys,opSurf);
  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_PMT_SUPPORT,REFLECTIVITY_OF_PMT_SUPPORT,0,0};

  //not needed for polished surface
  G4double specularspike[] = {0.,0.,0.,0.};
  G4double specularlobe[]  = {0.,0.,0.,0.};
  G4double backscatter[]   = {0.,0.,0.,0.};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;
}




G4LogicalBorderSurface* ArDM_DetectorConstruction::build_PMTSupportCoat_PMTSupport_surf(G4String topBottom){


  G4VPhysicalVolume* pmtSupportCoatPhys;
  G4VPhysicalVolume* pmtSupportPhys;

  if(!strcmp(topBottom,"top")){
    pmtSupportCoatPhys  = fTopPMTSupportCoatPhys;
    pmtSupportPhys      = fTopPMTSupportPhys;
  }else{
    pmtSupportCoatPhys  = fBottomPMTSupportCoatPhys;
    pmtSupportPhys      = fBottomPMTSupportPhys;    
  }

  ////dielectric_dielectric doesn't work <-- no reflection on reflector --> why !!?????
  ////this does work for LAr_WLS layer ! why doesn't it work for WLS_WLSSupport !!????? --> go figure !!
  //   G4OpticalSurface*       opSurf     = new G4OpticalSurface("PMTSupportCoat_PMTSupport",
  // 							    unified,ground,dielectric_dielectric); 

//   G4OpticalSurface*       opSurf     = new G4OpticalSurface("PMTSupportCoat_PMTSupport",
// 							    unified,polished,dielectric_metal); 

  G4OpticalSurface*       opSurf     = new G4OpticalSurface("PMTSupportCoat_PMTSupport",unified,ground,dielectric_metal); 
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("PMTSupportCoat_PMTSupport",
								  pmtSupportCoatPhys,pmtSupportPhys,opSurf);

  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};

  G4double pmtSupport_reflectivity=0;
  double specularspikeConst,specularlobeConst,backscatterConst,sigmaAlpha;

  if(!strcmp(topBottom,"top")){
    pmtSupport_reflectivity = REFLECTIVITY_OF_PMT_SUPPORT_TOP;
    specularspikeConst = TOP_PMT_SUPPORT_SPECULAR_SPIKE;
    specularlobeConst  = TOP_PMT_SUPPORT_SPECULAR_LOBE;
    backscatterConst   = TOP_PMT_SUPPORT_BACK_SCATTER;
    sigmaAlpha         = TOP_PMT_SUPPORT_REFLECTOR_SIGMA_ALPHA;
  }else if(!strcmp(topBottom,"bottom")){
    pmtSupport_reflectivity = REFLECTIVITY_OF_PMT_SUPPORT_BOTTOM;
    specularspikeConst = BOTTOM_PMT_SUPPORT_SPECULAR_SPIKE;
    specularlobeConst  = BOTTOM_PMT_SUPPORT_SPECULAR_LOBE;
    backscatterConst   = BOTTOM_PMT_SUPPORT_BACK_SCATTER;
    sigmaAlpha         = BOTTOM_PMT_SUPPORT_REFLECTOR_SIGMA_ALPHA;
  }

  opSurf->SetSigmaAlpha(sigmaAlpha*deg);

  G4double reflectivity[] = {pmtSupport_reflectivity,pmtSupport_reflectivity,0,0};

  //not needed for polished surface
  G4double specularspike[] = {specularspikeConst,specularspikeConst,specularspikeConst,specularspikeConst};
  G4double specularlobe[]  = {specularlobeConst,specularlobeConst,specularlobeConst,specularlobeConst};
  G4double backscatter[]   = {backscatterConst,backscatterConst,backscatterConst,backscatterConst};

// #if TEST_BRANCH0
//   if(!strcmp(topBottom,"top"))         ArDM_Analysis::getInstance()->fReflPMTSupportTop    = REFLECTIVITY_OF_PMT_SUPPORT;
//   else if(!strcmp(topBottom,"bottom")) ArDM_Analysis::getInstance()->fReflPMTSupportBottom = REFLECTIVITY_OF_PMT_SUPPORT;
// #endif //TEST_BRANCH0

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;
}






G4LogicalBorderSurface* ArDM_DetectorConstruction::build_Ar_PMTSupportCoat_surf(G4String topBottom,G4String order){


  G4VPhysicalVolume* medium1;
  G4VPhysicalVolume* medium2;

  ostringstream surfName;
  if(!strcmp(topBottom,"top")){
    medium1  = fGArColPhys;
    medium2  = fTopPMTSupportCoatPhys;
    surfName<<"GAr";
  }else{
    medium1  = fLArColPhys;
    medium2  = fBottomPMTSupportCoatPhys;    
    surfName<<"LAr";
  }

  

  if(!strcmp(order,"PMTSupportCoat_Ar")){
    G4VPhysicalVolume* medium = medium1;
    medium1 = medium2;
    medium2 = medium;
    surfName.str("PMTSupportCoat_"+surfName.str());

  }else{
    surfName<<"_PMTSupportCoat";
  }



//   G4OpticalSurface*       opSurf     = new G4OpticalSurface("PMTSupportCoat_PMTSupport",
// 							    unified,ground,dielectric_metal); 

  G4OpticalSurface*       opSurf     = new G4OpticalSurface("Ar_PMTSupportCoat",
							    unified,ground,dielectric_dielectric); 
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface(surfName.str().c_str(),medium1,medium2,opSurf);
  
  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_PMT_SUPPORT_COAT,REFLECTIVITY_OF_PMT_SUPPORT_COAT};

  //not needed for polished surface
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;
}




void ArDM_DetectorConstruction::addNeutronShield(){

  //neutrons can wander around before entering the shield
  //if we have a "one-piece" shield, it's difficult to get the information
  //whether neutrons enter the shield in the direction tank-towards-detector-axis
  //or the other way around.
  //
  //(neutron can manage to go into the shielded volume w/o passing through the shield,
  //and from there enter the shield from "inside-out", passing through the shield and
  //get out of the shielded volume)

  //
  //--> in order to obtain the information, from where neutrons enter the shield,
  //divide the shield into "innerShield" and "outerShield"

  //
  //since "innerShield" and "outerShield" are just labels, it doesn't matter how thick each shield is
  //as long as the total thickness is equal to the actual shield's thickness
  //for latter convenience, e.g. to shoot background-neutrons from the neutron shield itself using /gps/pos/confine ...,
  //make the outerShield really thin, e.g. 1 micron,
  //and the innerShield really thick, so that practically, innerShield = actually shield
  //

  //G4Material* mat = fPolycarbonate;
  //G4Material* mat = fPolyethylene;
  G4Material* mat = fBorotron;
  
  G4double innerR = NEUTRON_SHIELD_INNER_RADIUS;
  //G4double interR = innerR + NEUTRON_SHIELD_THICKNESS/2;
  G4double interR = innerR + NEUTRON_SHIELD_THICKNESS - 1e-6*m;
  G4double outerR = NEUTRON_SHIELD_OUTER_RADIUS;
  G4double halfheight = NEUTRON_SHIELD_HALF_HEIGHT;

  if(!NEUTRON_SHIELD_THICKNESS) return;

  G4Tubs*          tubSolid1 = new G4Tubs("neutronInnerShield",innerR,interR,halfheight,0.*deg,360.*deg);
  G4LogicalVolume* neutronShieldLog1   = new G4LogicalVolume(tubSolid1,mat,"neutronInnerShield");

  G4VisAttributes* att = new G4VisAttributes(1);
  att->SetColour(G4Colour::Grey());
  att->SetForceAuxEdgeVisible(true);
  neutronShieldLog1->SetVisAttributes(att);


  G4Tubs*          tubSolid2 = new G4Tubs("neutronInnerShield",interR,outerR,halfheight,0.*deg,360.*deg);
  G4LogicalVolume* neutronShieldLog2   = new G4LogicalVolume(tubSolid2,mat,"neutronInnerShield");

  neutronShieldLog2->SetVisAttributes(att);

  G4VPhysicalVolume* fMotherPhys;
  G4double z;
  if(fLArColPhys){ fMotherPhys = fLArColPhys; z = NEUTRON_SHIELD_POS_Z_LAR;   }
  else{            fMotherPhys = fWorldPhys;  z = NEUTRON_SHIELD_POS_Z;       }
  
  G4ThreeVector pos(NEUTRON_SHIELD_POS_X,NEUTRON_SHIELD_POS_Y,z);

  fNeutronShieldPhys = new G4PVPlacement(0,pos,"neutronInnerShield",neutronShieldLog1,fMotherPhys,false,0);
  new G4PVPlacement(0,pos,"neutronOuterShield",neutronShieldLog2,fMotherPhys,false,0);
  

// #if TEST_BRANCH24	
//   //ArDM_Analysis::getInstance()->PolyVol = fNeutronShieldPhys->GetLogicalVolume()->GetSolid()->GetCubicVolume();
//   ArDM_Analysis::getInstance()->PolyVol = neutronShieldLog1->GetSolid()->GetCubicVolume() +
//                                           neutronShieldLog1->GetSolid()->GetCubicVolume() ;


// #endif //TEST_BRANCH24	

  return;
}











void ArDM_DetectorConstruction::addEfield(){

  if(!fLArColPhys) return;


  //defining volume confined by the reflector


  //full cylinder
  G4double refVol_innerR      = 0*mm;
  G4double refVol_outerR      = WLS_RINGSEC_OUTER_RADIUS;
  G4double refVol_half_height = WLS_RINGSEC_HALF_HEIGHT;
  G4double refVol_start_phi   = 0.*deg;
  G4double refVol_delta_phi   = 360.*deg;
  
  G4Tubs* cylinder = new G4Tubs("cylinder",refVol_innerR,refVol_outerR,refVol_half_height,
				refVol_start_phi,refVol_delta_phi);
  

  //cuboid
  G4double x = refVol_outerR*(1-cos(40.*deg));
  G4double y = WLS_LINSEC_HALF_X;
  G4double z = refVol_half_height+1*mm;

  G4Box* box = new G4Box("box",x,y,z);
  
  G4ThreeVector boxPos(refVol_outerR,0,0);

  G4SubtractionSolid* refVolSolid = new G4SubtractionSolid("refVol",cylinder,box,0,boxPos);

  G4LogicalVolume* fRefVolLog = new G4LogicalVolume(refVolSolid,fLAr,"refVol",0,0,0);
  //G4LogicalVolume* fRefVolLog = new G4LogicalVolume(refVolSolid,fGAr,"refVol",0,0,0);
  //G4LogicalVolume* fRefVolLog = new G4LogicalVolume(refVolSolid,fVacuum,"refVol",0,0,0);
  
  G4VisAttributes* refVolAtt = new G4VisAttributes(false);
  refVolAtt->SetColour(G4Color::Green()); 
  refVolAtt->SetForceAuxEdgeVisible(true);
  fRefVolLog->SetVisAttributes(refVolAtt);
  
  G4ThreeVector refVolPos(WLS_RINGSEC_POS_X,WLS_RINGSEC_POS_Y,WLS_RINGSEC_POS_Z_LARCOL);
  
  //place refVol in the LAr column
  fReflectorVolume = new G4PVPlacement(0,refVolPos,"refVol",fRefVolLog,fLArColPhys,false,0);
  

  //assign an Efield to refVol

  G4double         EfieldStrength  = ELECTRIC_FIELD_STRENGTH;
  G4ElectricField* Efield          = new G4UniformElectricField(G4ThreeVector(0.,0.,-EfieldStrength));
  G4EqMagElectricField* equation   = new G4EqMagElectricField(Efield);
  G4MagIntegratorStepper* stepper  = new G4ClassicalRK4(equation,8);

  G4double minstep = 10.*um;
  G4MagInt_Driver* intgrDriver     = new G4MagInt_Driver(minstep,stepper,stepper->GetNumberOfVariables());
  G4ChordFinder*   chordFinder     = new G4ChordFinder(intgrDriver);
  

  G4FieldManager* fieldManager = new G4FieldManager();
  fieldManager->SetDetectorField(Efield);
  fieldManager->SetChordFinder(chordFinder);

  fRefVolLog->SetFieldManager(fieldManager,1);

  return;
}









void ArDM_DetectorConstruction::addPerlite(){

  if(!fWorldPhys){
    G4cout<<"in ArDM_DetectorConstruction::addPerlite(), fWorldPhys = 0. exit."<<G4endl;
    return ;  
  }


  G4double innerR = PERLITE_COLUMN_INNER_RADIUS;
  G4double outerR = PERLITE_COLUMN_OUTER_RADIUS;
  G4double halfz  = PERLITE_COLUMN_HALF_HEIGHT;
  G4Tubs* coatSolid = new G4Tubs("coatSolid",innerR,outerR,halfz,0*deg, 360*deg);
  
  
  //the perlite volume surrounds the tank / top flange
  //--> cut the coatSolid with appropriate objects
  //    to give the perlite volume the correct shape


  //dewar
  G4double dewar_innerR = 0;
  G4double dewar_outerR = TANK_CYLINDER_OUTER_RADIUS;
  G4double dewar_halfz  = TANK_CYLINDER_HALF_HEIGHT;
  G4Tubs* dewarSolid = new G4Tubs("dewarSolid",dewar_innerR,dewar_outerR,dewar_halfz,0*deg,360*deg);

  
  //top flange
  G4double tf_innerR = TOP_FLANGE_INNER_RADIUS;
  G4double tf_outerR = TOP_FLANGE_OUTER_RADIUS;
  G4double tf_halfz  = TOP_FLANGE_HALF_HEIGHT;
  G4Tubs* tfSolid = new G4Tubs("tfSolid",tf_innerR,tf_outerR,tf_halfz,0*deg,360*deg);

  
  //position relative to the center of coatSolid
  G4ThreeVector dewarPos(TANK_CYLINDER_POS_X - PERLITE_COLUMN_POS_X,
			 TANK_CYLINDER_POS_Y - PERLITE_COLUMN_POS_Y,
			 TANK_CYLINDER_POS_Z - PERLITE_COLUMN_POS_Z);


  G4ThreeVector tfPos(TOP_FLANGE_POS_X - PERLITE_COLUMN_POS_X,
		      TOP_FLANGE_POS_Y - PERLITE_COLUMN_POS_Y,
		      TOP_FLANGE_POS_Z - PERLITE_COLUMN_POS_Z);

  G4SubtractionSolid* perliteSolid = new G4SubtractionSolid("perlite",coatSolid,dewarSolid,0,dewarPos);
  perliteSolid = new G4SubtractionSolid("perlite",perliteSolid,tfSolid,0,tfPos);


  G4LogicalVolume* perliteLog = new G4LogicalVolume(perliteSolid,fPerlite,"perlite");

  G4ThreeVector perlitePos(PERLITE_COLUMN_POS_X,PERLITE_COLUMN_POS_Y,PERLITE_COLUMN_POS_Z);

  fPerlitePhys = new G4PVPlacement(0,perlitePos,"perlite",perliteLog,fWorldPhys,false,0);
  
  G4VisAttributes* att = new G4VisAttributes(false);
  att->SetColour(0.3,0.7,.1); //dunno what color this is. whatever !                                                                                             
  //at->SetColour(kRed); //dunno what color this is. whatever !                                                                                                 
  att->SetForceAuxEdgeVisible(true);
  perliteLog->SetVisAttributes(att);
  

  return;
}










void ArDM_DetectorConstruction::addSource(){

  
  if(!fLArColPhys){
    G4cout<<"in ArDM_DetectorConstruction::addSource(): fLArCol = 0. exit."<<G4endl;
    return;
  }

  G4double radius = SOURCE_RADIUS;
  G4double half_thickness = SOURCE_HALF_THICKNESS;

  G4Tubs* solid = new G4Tubs("",0,radius,half_thickness,0*deg,360*deg);

  G4LogicalVolume* log = new G4LogicalVolume(solid,fPolyethylene,"source");
  
  G4VisAttributes* att = new G4VisAttributes(false);
  log->SetVisAttributes(att);

  G4ThreeVector pos(SOURCE_POS_X,SOURCE_POS_Y,SOURCE_POS_Z_REL_TO_LARCOL);
  fSource = new G4PVPlacement(0,pos,"source",log,fLArColPhys,false,0);
  
  return;
}






void ArDM_DetectorConstruction::addSourceHolder(){
  
  if(!fLArColPhys){
    G4cout<<"in ArDM_DetectorConstruction::addSourceHolder(): fLArCol = 0. exit."<<G4endl;
    return;
  }

  G4double sourceHolder_halfx = SOURCE_HOLDER_HALF_X;
  G4double sourceHolder_halfy = SOURCE_HOLDER_HALF_Y;
  G4double sourceHolder_halfz = SOURCE_HOLDER_HALF_Z;

  G4double sourceHolder_hole_r           = SOURCE_HOLDER_HOLE_RADIUS;
  G4double sourceHolder_hole_half_height = SOURCE_HOLDER_HOLE_HALF_HEIGHT;

  
  G4Box*  sourceHolderbox  = new G4Box("",sourceHolder_halfx,sourceHolder_halfy,sourceHolder_halfz);
  G4Tubs* sourceHolderHole = new G4Tubs("",0,sourceHolder_hole_r,sourceHolder_hole_half_height,0,360);

  G4ThreeVector hole_box_relPos(SOURCE_HOLDER_HOLE_POS_RELATIVE_TO_SOURCE_HOLDER_CENTER_X,
				SOURCE_HOLDER_HOLE_POS_RELATIVE_TO_SOURCE_HOLDER_CENTER_Y,
				SOURCE_HOLDER_HOLE_POS_RELATIVE_TO_SOURCE_HOLDER_CENTER_Z);

  G4SubtractionSolid* sourceHolderSolid = new G4SubtractionSolid("sourceHolder",sourceHolderbox,sourceHolderHole,
								 0,hole_box_relPos);

  G4LogicalVolume* sourceHolderLog = new G4LogicalVolume(sourceHolderSolid,fTeflon,"sourceHolder");

  G4VisAttributes* sourceHolderAtt = new G4VisAttributes(true);
  sourceHolderAtt->SetColour(0.3,0.7,.1); //dunno what color this is. whatever !
  //sourceHolderAtt->SetColour(kRed); //dunno what color this is. whatever !
  sourceHolderAtt->SetForceAuxEdgeVisible(true);
  sourceHolderLog->SetVisAttributes(sourceHolderAtt);

  G4ThreeVector pos(SOURCE_HOLDER_POS_X,SOURCE_HOLDER_POS_Y,SOURCE_HOLDER_POS_Z_REL_TO_LARCOL);

  fSourceHolder = new G4PVPlacement(0,pos,"sourceHolder",sourceHolderLog,fLArColPhys,false,0);

  //rotate the source holder +- 10 degrees around the x axis
  
  G4RotationMatrix* rotMat = new G4RotationMatrix;
  rotMat->rotateX(HANDEDNESS*(0)*deg);
  rotMat->rotateY(HANDEDNESS*(0)*deg);
  rotMat->rotateZ(HANDEDNESS*0*deg);

  fSourceHolder->SetRotation(rotMat);


  return;
}





void ArDM_DetectorConstruction::addSourceCoating(){

  
  
  if(!fLArColPhys){
    G4cout<<"in ArDM_DetectorConstruction::addSourceHolder(): fLArCol = 0. exit."<<G4endl;
    return;
  }

  G4double sourceCoating_r = SOURCE_COATING_RADIUS;
  G4double sourceCoating_halfheight = SOURCE_COATING_HALF_HEIGHT;
  
  G4Tubs* sourceCoatingSolid = new G4Tubs("",0,sourceCoating_r,sourceCoating_halfheight,0,360);
  
  G4LogicalVolume* sourceCoatingLog = new G4LogicalVolume(sourceCoatingSolid,fPolyethylene,"sourceCoating");

  G4ThreeVector pos(SOURCE_COATING_POS_X,SOURCE_COATING_POS_X,SOURCE_COATING_POS_Z_REL_TO_LARCOL);

  fSourceCoating = new G4PVPlacement(0,pos,"sourceCoating",sourceCoatingLog,fLArColPhys,false,0);

  return;
}



G4LogicalVolume* ArDM_DetectorConstruction::constructGrid(G4Material* fMat,double plateInnerR,double plateOuterR,
							  double plateHalfThickness,double wireInnerR,double wireOuterR,
							  double wirePitch,const char* name,int attribute,G4Colour color){


  


  //plate
  G4Tubs* plateSolid = new G4Tubs("plate",plateInnerR,plateOuterR,plateHalfThickness,0.*deg,360.*deg);
  
  //wire
  vector<double> wireHalfLengthx,wireHalfLengthy;

  vector<G4ThreeVector> wirePosx = setGridWireVector("x",wireHalfLengthx,plateInnerR,wirePitch,wireOuterR);
  vector<G4ThreeVector> wirePosy = setGridWireVector("y",wireHalfLengthy,plateInnerR,wirePitch,wireOuterR);

  

  if(!wirePosx.size() && !wirePosy.size()) return NULL;

  
  G4RotationMatrix* cathodeWireRotMatx = new G4RotationMatrix;
  cathodeWireRotMatx->rotateX(90.*deg);
  cathodeWireRotMatx->rotateY(0);
  cathodeWireRotMatx->rotateZ(0);

  G4RotationMatrix* cathodeWireRotMaty = new G4RotationMatrix;
  cathodeWireRotMaty->rotateX(0);
  cathodeWireRotMaty->rotateY(90.*deg);
  cathodeWireRotMaty->rotateZ(0);


  G4Tubs* wireSolid0 = new G4Tubs("wire",wireInnerR,wireOuterR,wireHalfLengthx[0],0.*deg,360.*deg);
  
  
  G4UnionSolid* gridSolid = new G4UnionSolid("grid",plateSolid,wireSolid0,cathodeWireRotMatx,wirePosx[0]);
  
  for(int i=1;i<wirePosx.size();i++){
    G4Tubs* wireSolid = new G4Tubs("wire",wireInnerR,wireOuterR,wireHalfLengthx[i],0.*deg,360.*deg);
    gridSolid = new G4UnionSolid("grid",gridSolid,wireSolid,cathodeWireRotMatx,wirePosx[i]);
  }
  
  for(int i=0;i<wirePosy.size();i++){
    G4Tubs* wireSolid = new G4Tubs("wire",wireInnerR,wireOuterR,wireHalfLengthy[i],0.*deg,360.*deg);
    gridSolid = new G4UnionSolid("grid",gridSolid,wireSolid,cathodeWireRotMaty,wirePosy[i]);
  }
  


  G4LogicalVolume* gridLog = new G4LogicalVolume(gridSolid,fMat,name,0,0,0);
  
  G4VisAttributes* att = new G4VisAttributes(true);
  att->SetColour(color);
  att->SetForceAuxEdgeVisible(true);
  gridLog->SetVisAttributes(att);

  return gridLog;
}













void ArDM_DetectorConstruction::addCathodeGrid(){

  double plateInnerR = CATHODE_PLATE_INNER_RADIUS;
  double plateOuterR = CATHODE_PLATE_OUTER_RADIUS;
  double plateHalfThickness = CATHODE_PLATE_HALF_THICKNESS;

  double wireInnerR = CATHODE_WIRE_INNER_RADIUS;
  double wireOuterR = CATHODE_WIRE_OUTER_RADIUS;
  double wirePitch  = CATHODE_WIRE_PITCH;

  G4String name = "cathodeGrid";

  G4LogicalVolume* gridLog = constructGrid(fCathodeGridMat,plateInnerR,plateOuterR,plateHalfThickness,
					   wireInnerR,wireOuterR,wirePitch,name.c_str(),1,G4Colour::Green());

  G4ThreeVector pos(CATHODE_POS_X,CATHODE_POS_Y,CATHODE_POS_Z_LARCOL);
  fCathodeGridPhys = new G4PVPlacement(0,pos,name,gridLog,fLArColPhys,false,0);

  return;
}
















void ArDM_DetectorConstruction::addProtectionGrid(){

  double plateInnerR = BOTTOM_PROTECTION_GRID_PLATE_INNER_RADIUS;
  double plateOuterR = BOTTOM_PROTECTION_GRID_PLATE_OUTER_RADIUS;
  double plateHalfThickness = BOTTOM_PROTECTION_GRID_PLATE_HALF_THICKNESS;

  double wireInnerR = BOTTOM_PROTECTION_GRID_WIRE_INNER_RADIUS;
  double wireOuterR = BOTTOM_PROTECTION_GRID_WIRE_OUTER_RADIUS;
  double wirePitch  = BOTTOM_PROTECTION_GRID_WIRE_PITCH;

  G4String name = "protectionGrid";

  G4LogicalVolume* gridLog = constructGrid(fCathodeGridMat,plateInnerR,plateOuterR,plateHalfThickness,
					   wireInnerR,wireOuterR,wirePitch,name.c_str(),1,G4Colour::Red());

  G4ThreeVector pos(BOTTOM_PROTECTION_GRID_POS_X,BOTTOM_PROTECTION_GRID_POS_Y,BOTTOM_PROTECTION_GRID_POS_Z_LARCOL);
  fProtectionGridPhys = new G4PVPlacement(0,pos,name,gridLog,fLArColPhys,false,0);

  return;
}




































void ArDM_DetectorConstruction::addBottomSideReflectorCoat(){

  G4double innerR1     = BOTTOM_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE;
  G4double outerR1     = BOTTOM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE;
  G4double innerR2     = BOTTOM_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE;
  G4double outerR2     = BOTTOM_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE;
  G4double half_height = BOTTOM_SIDE_REFLECTOR_COAT_HALF_HEIGHT_GASTEST;
  G4double start_phi   = BOTTOM_SIDE_REFLECTOR_START_PHI;
  G4double delta_phi   = BOTTOM_SIDE_REFLECTOR_DELTA_PHI;
  
  G4Cons* coneSolid = new G4Cons("bottomSideReflector",innerR1,outerR1,innerR2,outerR2,half_height,start_phi,delta_phi);
  


  G4double cynPart_innerR = BOTTOM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_INNER_RADIUS;
  G4double cynPart_outerR = BOTTOM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_OUTER_RADIUS;
  G4double cynPart_half_height = BOTTOM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_HALF_HEIGHT;
  
  G4Tubs* tubeSolid = new G4Tubs("bottomSideReflector",cynPart_innerR,cynPart_outerR,
				 cynPart_half_height,0*deg,360.*deg);

  G4ThreeVector pos_tube_rel_cone(0.,0.,
				  BOTTOM_SIDE_REFLECTOR_COAT_CYNLIDRICAL_PART_POS_Z
				  -BOTTOM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST);
  

  G4UnionSolid* solid = new G4UnionSolid("bottomSideReflector",coneSolid,tubeSolid,0,pos_tube_rel_cone);


  G4LogicalVolume* fBottomSideReflectorLog = new G4LogicalVolume(solid,fWLSMat,
								 "bottomSideReflectorCoat",0,0,0);
  
  //place fBottomSideReflectorLog in fLArColLog
  G4ThreeVector pos(BOTTOM_SIDE_REFLECTOR_COAT_POS_X,
		    BOTTOM_SIDE_REFLECTOR_COAT_POS_Y,
		    BOTTOM_SIDE_REFLECTOR_COAT_POS_Z_GASTEST_LARCOL);
  
  fBottomSideReflectorCoatPhys = new G4PVPlacement(0,pos,"bottomSideReflectorCoat",
						   fBottomSideReflectorLog,fLArColPhys,false,0);
  
  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  sideReflAtt->SetColour(G4Color::Magenta());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fBottomSideReflectorLog->SetVisAttributes(sideReflAtt);
  
  return;
}




















void ArDM_DetectorConstruction::addTopSideReflectorCoat(){



  G4double innerR1_GAr = TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double outerR1_GAr = TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_GARCOL;
  G4double innerR2_GAr = TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double outerR2_GAr = TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_GARCOL;
  G4double halfheight_GAr = TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_GARCOL;
  
  G4double innerR1_LAr = TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double outerR1_LAr = TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_UPPER_EDGE_IN_LARCOL;
  G4double innerR2_LAr = TOP_SIDE_REFLECTOR_COAT_INNER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double outerR2_LAr = TOP_SIDE_REFLECTOR_COAT_OUTER_RADIUS_LOWER_EDGE_IN_LARCOL;
  G4double halfheight_LAr = TOP_SIDE_REFLECTOR_COAT_HALF_HEIGHT_IN_LARCOL;
  
  G4double start_phi   = TOP_SIDE_REFLECTOR_COAT_START_PHI;
  G4double delta_phi   = TOP_SIDE_REFLECTOR_COAT_DELTA_PHI;

  G4Cons* ringSolid_GAr = NULL;
  G4Cons* ringSolid_LAr = NULL;

  if(halfheight_GAr) ringSolid_GAr = new G4Cons("topSideReflectorCoat",innerR2_GAr,outerR2_GAr,innerR1_GAr,outerR1_GAr,
						halfheight_GAr,start_phi,delta_phi);

  if(halfheight_LAr) ringSolid_LAr = new G4Cons("topSideReflectorCoat",innerR2_LAr,outerR2_LAr,innerR1_LAr,outerR1_LAr,
						halfheight_LAr,start_phi,delta_phi);

 
  G4LogicalVolume* fTopSideReflectorLog_GAr = NULL;
  G4LogicalVolume* fTopSideReflectorLog_LAr = NULL;



  if(ringSolid_GAr) fTopSideReflectorLog_GAr = new G4LogicalVolume(ringSolid_GAr,fWLSMat,
								   "topSideReflectorCoat",0,0,0);

  if(ringSolid_LAr) fTopSideReflectorLog_LAr = new G4LogicalVolume(ringSolid_LAr,fWLSMat,
								   "topSideReflectorCoat",0,0,0);

  
  G4ThreeVector pos_GAr(TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_X,
			TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Y,
			TOP_SIDE_REFLECTOR_COAT_IN_GARCOL_POS_Z_GARCOL);

  G4ThreeVector pos_LAr(TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_X,
			TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Y,
			TOP_SIDE_REFLECTOR_COAT_IN_LARCOL_POS_Z_LARCOL);

  if(fTopSideReflectorLog_GAr && fGArColPhys) 
    fTopSideReflectorCoatGArPhys = new G4PVPlacement(0,pos_GAr,"topSideReflectorCoat",
						     fTopSideReflectorLog_GAr,fGArColPhys,false,0);

  if(fTopSideReflectorLog_LAr && fLArColPhys) 
    fTopSideReflectorCoatLArPhys = new G4PVPlacement(0,pos_LAr,"topSideReflectorCoat",
						     fTopSideReflectorLog_LAr,fLArColPhys,false,0);
    

  G4VisAttributes* sideReflAtt = new G4VisAttributes(true);
  sideReflAtt->SetColour(G4Color::Magenta());
  sideReflAtt->SetForceAuxEdgeVisible(true);
  fTopSideReflectorLog_GAr->SetVisAttributes(sideReflAtt);
  fTopSideReflectorLog_LAr->SetVisAttributes(sideReflAtt);

  return;   
  
}









G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_BottomSideReflectorCoat_surf(G4String order){

  
  G4OpticalSurface* opSurf     = new G4OpticalSurface("LAr_BottomSideReflectorCoat",
						      unified,ground,dielectric_dielectric); 

  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[]  = {1.,1.};   //test //no absorption
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  G4LogicalBorderSurface* borderSurf;
  
  if(!strcmp(order,"LAr_BottomSideReflectorCoat"))
    borderSurf = new G4LogicalBorderSurface(order,fLArColPhys,fBottomSideReflectorCoatPhys,opSurf);
  else if(!strcmp(order,"BottomSideReflectorCoat_LAr"))
    borderSurf = new G4LogicalBorderSurface(order,fBottomSideReflectorCoatPhys,fLArColPhys,opSurf);

  return borderSurf;
}











G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_TopSideReflectorCoat_surf(G4String order){

  
  G4OpticalSurface* opSurf     = new G4OpticalSurface("LAr_TopSideReflectorCoat",
						      unified,ground,dielectric_dielectric); 

  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[]  = {1.,1.};   //test //no absorption
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  G4LogicalBorderSurface* borderSurf;
  
  if(!strcmp(order,"LAr_TopSideReflectorCoat"))
    borderSurf = new G4LogicalBorderSurface(order,fLArColPhys,fTopSideReflectorCoatLArPhys,opSurf);
  else if(!strcmp(order,"TopSideReflectorCoat_LAr"))
    borderSurf = new G4LogicalBorderSurface(order,fTopSideReflectorCoatLArPhys,fLArColPhys,opSurf);

  return borderSurf;
}












G4LogicalBorderSurface* ArDM_DetectorConstruction::build_GAr_TopSideReflectorCoat_surf(G4String order){

  
  G4OpticalSurface* opSurf     = new G4OpticalSurface("GAr_TopSideReflectorCoat",
						      unified,ground,dielectric_dielectric); 

  G4double E[]={0.1*eV,13.5*eV};
  G4double reflectivity[]  = {1.,1.};   //test //no absorption
  G4double specularspike[] = {0.,0.};
  G4double specularlobe[]  = {0.,0.};
  G4double backscatter[]   = {0.,0.};

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  G4LogicalBorderSurface* borderSurf;
  
  if(!strcmp(order,"GAr_TopSideReflectorCoat"))
    borderSurf = new G4LogicalBorderSurface(order,fGArColPhys,fTopSideReflectorCoatGArPhys,opSurf);
  else if(!strcmp(order,"TopSideReflectorCoat_GAr"))
    borderSurf = new G4LogicalBorderSurface(order,fTopSideReflectorCoatGArPhys,fGArColPhys,opSurf);

  return borderSurf;
}
















G4LogicalBorderSurface* ArDM_DetectorConstruction::build_BottomSideReflector_BottomSideReflectorCoat_surf(){



  
  G4OpticalSurface* opSurf = new G4OpticalSurface("BottomSideReflector_BottomSideReflectorCoat",
						  unified,ground,dielectric_metal);



  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("BottomSideReflector_BottomSideReflectorCoat",
								  fBottomSideReflectorCoatPhys,
								  fBottomSideReflectorPhys,opSurf);

  //test
  G4double Emin  = 0.1  *eV;
  G4double Emid1 = 6.0  *eV;
  G4double Emid2 = 6.0001  *eV;
  G4double Emax  = 13.5 *eV;
  
  G4double E[]             = {Emin,Emid1,Emid2,Emax};
  G4double reflectivity[]  = {REFLECTIVITY_OF_BOTTOM_SIDE_REFLECTOR,REFLECTIVITY_OF_BOTTOM_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {0.,0.,0.,0.};
  G4double specularlobe[]  = {0.,0.,0.,0.};
  G4double backscatter[]   = {0.,0.,0.,0.};

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflSideReflectorBottom = REFLECTIVITY_OF_BOTTOM_SIDE_REFLECTOR;
// #endif //TEST_BRANCH0


  //the reflectivity given to the material properties table actually means (1 - absorption)
  //it will be used to determine whether the photon is absorbed or not
  //if not, further measures will be taken to determine whether the photon should undergo
  //
  //i.    total internal reflection
  //ii.   transmission
  //iii.  reflection <-- if reflection, determine the type of reflection by calling G4OpBoundaryProcess:ChooseReflection()
  //          a. specular spike
  //          b. specular lobe
  //          c. back scattering
  //          d. lambertian = 1 - a. - b. - c.
  //
  //see G4OpBoundaryProcess::PostStepDoIt(..) for more !


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);
  

  //end test


  return borderSurf;


}















G4LogicalBorderSurface* ArDM_DetectorConstruction::build_TopSideReflectorGAr_TopSideReflectorCoat_surf(){



  
  G4OpticalSurface* opSurf = new G4OpticalSurface("TopSideReflectorGAr_TopSideReflectorCoat",
						  unified,ground,dielectric_metal);


  //cout<<"fTopSideReflectorGArPhys "<<fTopSideReflectorGArPhys<<"\t fTopSideReflectorCoatGArPhys "<<fTopSideReflectorCoatGArPhys<<endl;getchar();
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("TopSideReflectorGAr_TopSideReflectorCoat",
								  fTopSideReflectorCoatGArPhys,
								  fTopSideReflectorGArPhys,opSurf);

  //test
  G4double Emin  = 0.1  *eV;
  G4double Emid1 = 6.0  *eV;
  G4double Emid2 = 6.0001  *eV;
  G4double Emax  = 13.5 *eV;
  
  G4double E[]             = {Emin,Emid1,Emid2,Emax};
  G4double reflectivity[]  = {REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {0.,0.,0.,0.};
  G4double specularlobe[]  = {0.,0.,0.,0.};
  G4double backscatter[]   = {0.,0.,0.,0.};

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflSideReflectorTop = REFLECTIVITY_OF_TOP_SIDE_REFLECTOR;
// #endif //TEST_BRANCH0

  //the reflectivity given to the material properties table actually means (1 - absorption)
  //it will be used to determine whether the photon is absorbed or not
  //if not, further measures will be taken to determine whether the photon should undergo
  //
  //i.    total internal reflection
  //ii.   transmission
  //iii.  reflection <-- if reflection, determine the type of reflection by calling G4OpBoundaryProcess:ChooseReflection()
  //          a. specular spike
  //          b. specular lobe
  //          c. back scattering
  //          d. lambertian = 1 - a. - b. - c.
  //
  //see G4OpBoundaryProcess::PostStepDoIt(..) for more !


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);
  

  //end test


  return borderSurf;


}
















G4LogicalBorderSurface* ArDM_DetectorConstruction::build_TopSideReflectorLAr_TopSideReflectorCoat_surf(){



  
  G4OpticalSurface* opSurf = new G4OpticalSurface("TopSideReflectorLAr_TopSideReflectorCoat",
						  unified,ground,dielectric_metal);


  //cout<<"fTopSideReflectorLArPhys "<<fTopSideReflectorLArPhys<<"\t fTopSideReflectorCoatLArPhys "<<fTopSideReflectorCoatLArPhys<<endl;getchar();
  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("TopSideReflectorLAr_TopSideReflectorCoat",
								  fTopSideReflectorCoatLArPhys,
								  fTopSideReflectorLArPhys,opSurf);

  //test
  G4double Emin  = 0.1  *eV;
  G4double Emid1 = 6.0  *eV;
  G4double Emid2 = 6.0001  *eV;
  G4double Emax  = 13.5 *eV;
  
  G4double E[]             = {Emin,Emid1,Emid2,Emax};
  G4double reflectivity[]  = {REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {0.,0.,0.,0.};
  G4double specularlobe[]  = {0.,0.,0.,0.};
  G4double backscatter[]   = {0.,0.,0.,0.};

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflSideReflectorTop = REFLECTIVITY_OF_TOP_SIDE_REFLECTOR;
// #endif //TEST_BRANCH0


  //the reflectivity given to the material properties table actually means (1 - absorption)
  //it will be used to determine whether the photon is absorbed or not
  //if not, further measures will be taken to determine whether the photon should undergo
  //
  //i.    total internal reflection
  //ii.   transmission
  //iii.  reflection <-- if reflection, determine the type of reflection by calling G4OpBoundaryProcess:ChooseReflection()
  //          a. specular spike
  //          b. specular lobe
  //          c. back scattering
  //          d. lambertian = 1 - a. - b. - c.
  //
  //see G4OpBoundaryProcess::PostStepDoIt(..) for more !


  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);
  

  //end test


  return borderSurf;


}

















G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_protectionGrid_surf(){

  //why ground ?
  //why dielectric_metal ??

  G4OpticalSurface* opSurf = new G4OpticalSurface("LAr_protectionGrid",unified,ground,dielectric_metal); 


  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("LAr_protectionGrid",
								  fLArColPhys,fProtectionGridPhys,opSurf);
  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,
			      TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE};
  
  G4double specularlobe[]  = {TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE,
			      TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE};

  G4double backscatter[]   = {TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER,
			      TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER};

  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflProtectionGrid = REFLECTIVITY_OF_PROTECTION_GRID;
// #endif //TEST_BRANCH0

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;

}









G4LogicalBorderSurface* ArDM_DetectorConstruction::build_LAr_cathodeGrid_surf(){

  //why ground ?
  //why dielectric_metal ??

  G4OpticalSurface* opSurf = new G4OpticalSurface("LAr_cathodeGrid",unified,ground,dielectric_metal); 


  G4LogicalBorderSurface* borderSurf = new G4LogicalBorderSurface("LAr_cathodeGrid",
								  fLArColPhys,fCathodeGridPhys,opSurf);
  
  G4double E[]={0.1*eV,6*eV,6.001*eV,13.5*eV};
  G4double reflectivity[] = {REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,REFLECTIVITY_OF_TOP_SIDE_REFLECTOR,0,0};
  G4double specularspike[] = {TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,
			      TOP_SIDE_REFLECTOR_SPECULAR_SPIKE,TOP_SIDE_REFLECTOR_SPECULAR_SPIKE};
  
  G4double specularlobe[]  = {TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE,
			      TOP_SIDE_REFLECTOR_SPECULAR_LOBE,TOP_SIDE_REFLECTOR_SPECULAR_LOBE};

  G4double backscatter[]   = {TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER,
			      TOP_SIDE_REFLECTOR_BACK_SCATTER,TOP_SIDE_REFLECTOR_BACK_SCATTER};

  //specularspike + specularlob + backscatter + lambertian = 1 !!
  //--> in this case: the surface is assumed to be lambertian !
  //end test

// #if TEST_BRANCH0
//   ArDM_Analysis::getInstance()->fReflCathodeGrid = REFLECTIVITY_OF_CATHODE_GRID;
// #endif //TEST_BRANCH0

  int ne = sizeof(E)/sizeof(G4double);
  G4MaterialPropertiesTable* propTab = new G4MaterialPropertiesTable;
  propTab->AddProperty("REFLECTIVITY",E,reflectivity,ne);
  propTab->AddProperty("SPECULARSPIKECONSTANT",E,specularspike,ne);
  propTab->AddProperty("SPECULARLOBECONSTANT",E,specularlobe,ne);
  propTab->AddProperty("BACKSCATTERCONSTANT",E,backscatter,ne);
  opSurf->SetMaterialPropertiesTable(propTab);

  return borderSurf;

}




















void ArDM_DetectorConstruction::verboseInfo(string detComp){


  //print informations about detector components,
  //especially : size, position coordinates 
  //--> useful for making envelope surrounding one particular detector component,
  //    and then shoot some particles from that detector component.



  if(detComp == "pmtxy"){

    vector<TVector2> rpmt = setPMTVector2D();
    double topz = TOP_PMT_Z;
    double btmz = BOTTOM_PMT_Z;
    
    cout<<"top array : topz = "<<topz<<"\t bottom array : btmz = "<<btmz<<endl<<endl;
    cout<<"pmtbase : topz = "<<APPROX_PMT_BASE_Z_TOP<<"\t btmz = "<<APPROX_PMT_BASE_Z_BTM<<endl<<endl;
    cout<<"xy position : "<<endl<<endl;
    
    for(int i=0;i<rpmt.size();i++){
      cout<<i<<"\t x : "<<rpmt[i].X()<<"\t\t y : "<<rpmt[i].Y()<<endl;
    }
    cout<<endl;
    
  }else if(detComp == ""){
    verboseInfo("tank");
    verboseInfo("topFlange");
    verboseInfo("neutronShield");
    verboseInfo("pmtGlass");
    verboseInfo("pmtBase");
    verboseInfo("pmtElectrode");
    verboseInfo("pmtxy");
    verboseInfo("HVr");
    verboseInfo("perlite");
  
  }else{

    G4double innerR=0;
    G4double outerR=0;
    G4double halfz =0;
    G4double topz  =DEFAULTVALUE;
    G4double btmz  =DEFAULTVALUE;
    G4double posz  =DEFAULTVALUE;
    G4double posx  = 0, posy = 0;
    
    G4String shape="";
    G4String volName="";
    
    
    if(detComp=="tank"){

      innerR = 0;
      outerR = TANK_CYLINDER_OUTER_RADIUS;
      halfz  = 2*TANK_CYLINDER_HALF_HEIGHT + DISTANCE_LOWER_EDGE_OF_TANK_CYLINDER_TO_TANK_LOWERMOST_POINT + TANK_CYLINDER_THICKNESS;
      halfz /= 2;
      
      posx   = 0.;
      posy   = 0.;
      posz   = TANK_CYLINDER_POS_Z + TANK_CYLINDER_HALF_HEIGHT - halfz;
      shape  = "cylinder";
      volName= fTankPhys->GetName();
      
      
    }else if(detComp=="topFlange"){
      innerR = TOP_FLANGE_INNER_RADIUS;
      outerR = TOP_FLANGE_OUTER_RADIUS;
      halfz  = TOP_FLANGE_HALF_HEIGHT_EFFECTIVE;
      posz   = TOP_FLANGE_POS_Z;
      shape  = "cylinder";
      volName= fTopFlangePhys->GetName();
      
    }else if(detComp=="neutronShield"){
      innerR = NEUTRON_SHIELD_INNER_RADIUS;
      outerR = innerR + 100*mm; //maximum outerR
      halfz  = NEUTRON_SHIELD_HALF_HEIGHT;
      posz   = NEUTRON_SHIELD_POS_Z;
      shape  = "cylinder";
      volName= "neutronInnerShield";
  
    }else if(detComp=="pmtGlass"){
      
      G4double pmtGlass_total_height = APPROX_PMT_SPHERICAL_PART_HEIGHT + 2*APPROX_PMT_MIDDLE_CYLINDER_HALF_HEIGHT;

      pmtGlass_total_height += APPROX_PMT_SPHERICAL_PART_OUTER_RADIUS*cos(APPROX_PMT_BTM_SPHERE_HOLE_OPENING_ANGLE/2) 
	                       - APPROX_PMT_SPHERICAL_PART_INNER_RADIUS*cos(APPROX_PMT_SPHERICAL_PART_OPENING_ANGLE/2);

      pmtGlass_total_height += 2*APPROX_PMT_BTM_TUBE_HALF_HEIGHT;
      
      innerR = 0;
      outerR = APPROX_PMT_MIDDLE_CYLINDER_OUTER_RADIUS;
      halfz  = pmtGlass_total_height/2;
      topz = APPROX_TOP_PMT_Z - APPROX_PMT_MIDDLE_CYLINDER_HALF_HEIGHT - APPROX_PMT_SPHERICAL_PART_HEIGHT + pmtGlass_total_height/2;
      btmz = APPROX_BTM_PMT_Z + APPROX_PMT_MIDDLE_CYLINDER_HALF_HEIGHT + APPROX_PMT_SPHERICAL_PART_HEIGHT - pmtGlass_total_height/2;
      shape="cylinder";
      volName="pmt{PMTid}";
      

    }else if(detComp=="pmtBase"){
      innerR = 0. *mm;
      outerR = APPROX_PMT_BASE_RADIUS;
      halfz  = APPROX_PMT_BASE_HALF_THICKNESS;
      topz   = APPROX_PMT_BASE_Z_TOP;
      btmz   = APPROX_PMT_BASE_Z_BTM;
      shape  = "cylinder";
      volName= "pmtBase{PMTid}";
      
    }else if(detComp=="pmtElectrode"){
      innerR = APPROX_PMT_ELECTRODE_INNER_RADIUS;
      outerR = APPROX_PMT_ELECTRODE_OUTER_RADIUS;
      halfz  = APPROX_PMT_ELECTRODE_HALF_HEIGHT;

      topz   = APPROX_PMT_ELECTRODE_Z_TOP;
      btmz   = APPROX_PMT_ELECTRODE_Z_BTM;

      shape  = "cylinder";
      volName= "pmtElectrode{PMTid}";
      
    }else if(detComp=="HVr"){
      
      innerR = HV_RESISTOR_BAR_INNER_RADIUS;
      outerR = HV_RESISTOR_BAR_OUTER_RADIUS;
      halfz  = HV_RESISTOR_BAR_HALF_HEIGHT;
      topz   = HV_RESISTOR_BAR_1_POS_Z;
      btmz   = HV_RESISTOR_BAR_2_POS_Z;
      
      posx   = HV_RESISTOR_BAR_1_POS_X;
      posy   = HV_RESISTOR_BAR_1_POS_Y;
      
      shape  = "cylinder";
      volName= "HVr";      
      
    }else if(detComp=="perlite"){

      innerR = PERLITE_COLUMN_INNER_RADIUS;
      outerR = PERLITE_COLUMN_OUTER_RADIUS;
      halfz  = PERLITE_COLUMN_HALF_HEIGHT;
      posx   = PERLITE_COLUMN_POS_X;
      posy   = PERLITE_COLUMN_POS_Y;
      posz   = PERLITE_COLUMN_POS_Z;
      shape  = "cylinder";
      volName=fPerlitePhys->GetName();

    }
    
    cout<<"envelope for "<< detComp <<endl
	<<"innerR "<<innerR <<endl
	<<"outerR "<<outerR <<endl
	<<"halfz "<<halfz <<endl
	<<"posx "<<posx <<"\t posy "<<posy <<"\t posz "<<posz <<"\t topz "<<topz<<"\t btmz "<<btmz<<endl
	<<"shape "<<shape<<endl
	<<"volName "<<volName <<endl
      
	<<endl<<endl;
  
  }

  
  return;
}
















