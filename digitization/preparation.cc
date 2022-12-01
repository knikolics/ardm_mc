#include "preparation.hh"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "G4Poisson.hh"
#include "istream"


G4ThreeVector randomPos2D(G4double radius,G4double z){
  G4double phi = 2*TMath::Pi()*drand48();
  G4double r   = radius*sqrt(drand48());
  G4ThreeVector vec(r*cos(phi),r*sin(phi),z);
  return vec;
}


G4ThreeVector sphericalDistribution(const char* direction){
  G4double cosTheta = 2*drand48() - 1;
  G4double sinTheta = sqrt(1-cosTheta*cosTheta);
  G4double phi = 2*PI*drand48();

  G4ThreeVector vector;

  if(!strcmp(direction,"negative"))
    vector.set(sinTheta*cos(phi),sinTheta*sin(phi),-abs(cosTheta));
  else if(!strcmp(direction,"positive"))
    vector.set(sinTheta*cos(phi),sinTheta*sin(phi),abs(cosTheta));
  else vector.set(sinTheta*cos(phi),sinTheta*sin(phi),cosTheta);

  return vector;
}


G4ThreeVector polarizationVec(G4ThreeVector momDir){
  //generate a random 3D vector perpendicular to momDir

  if(!momDir.mag2()) return sphericalDistribution();

  G4ThreeVector polarization = momDir.orthogonal();
  G4ThreeVector perp         = momDir.cross(polarization);

  G4double phi = twopi*G4UniformRand();
  polarization = std::cos(phi)*polarization + std::sin(phi)*perp;

  return polarization.unit();
}


vector<G4ThreeVector> setPMTVector_old_geometry(const char* top_or_bottom){
  
  vector<G4ThreeVector> rPMT;
  G4double z;

  if(!strcmp(top_or_bottom,"bottom"))   z = BOTTOM_PMT_Z;
  else if(!strcmp(top_or_bottom,"top")) z = TOP_PMT_Z;
  
  G4double x[] = {-103.,103.,-206.,0.,206.,-309.,-103.,103.,309.,-206.,0.,206.,-103.,103.};
  G4double y[] = {356.8,356.8,178.4,178.4,178.4,0.,0.,0.,0.,-178.4,-178.4,-178.4,-356.8,-356.8};

  int n = sizeof(x)/sizeof(G4double);
#ifdef SCALEFACTOR
  for(int i=0;i<n;i++) rPMT.push_back(G4ThreeVector(x[i]*SCALEFACTOR,y[i]*SCALEFACTOR,z));
#else
  for(int i=0;i<n;i++) rPMT.push_back(G4ThreeVector(x[i],y[i],z));
#endif
  return rPMT;
}



vector<TVector2> setPMTVector2D(){ 

#ifdef NPMT
#if (NPMT-12)
#undef NPMT
#define NPMT 12
#endif
#else
#define NPMT 12
#endif

G4double scalefactor = 1;

#ifdef SCALEFACTOR
  scalefactor = SCALEFACTOR;  
#endif //SCALEFACTOR

  double R1 = 122.4*mm  *scalefactor;
  double R2 = 244.8*mm  *scalefactor;
  double R3 = 323.84*mm *scalefactor;

  TVector2 r[NPMT];
  
  //first calculate the coordinates of the PMTs 
  //assuming that the x-axis is the symmetry axis of the linear sector of the reflector
  //the positive x-direction is from the center of the cross section of the reflector towards the linear sector

  //the 3 first PMTs are located on the innermost circle with radius R1
  //120 degrees shifted from each other


  //1st PMT
  r[0] = TVector2(-R1,0);

  //2nd and 3rd PMT
  r[1] = r[0].Rotate(120*deg);
  r[2] = r[0].Rotate(240*deg);


  //PMT 3,4,5 are located on the 2nd circle with radius R2
  //relative configuration is the same as for PMT 0,1,2
  //but rotated around the z-axis by 180 degrees
  r[3] = r[0].Rotate(180*deg)*R2/R1;
  r[4] = r[3].Rotate(120*deg);
  r[5] = r[3].Rotate(240*deg);

  
  //PMT 6,7,8,9,10,11 are on the outermost circle with radius R3
  //6, 7, 8  are 120   degrees shifted from each other
  //9,10,11  are 120   degrees shifted from each other
  //6 is -79.1 degrees away from 3
  //8 and 9 are 38.21 degrees away from each other
  r[6]  = r[3].Rotate(-79.1*deg)*R3/R2;
  r[7]  = r[6].Rotate(120*deg);
  r[8]  = r[6].Rotate(240*deg);

  r[9]  = r[8].Rotate(38.21*deg);
  r[10] = r[9].Rotate(120*deg);
  r[11] = r[9].Rotate(240*deg);

  vector<TVector2> rPMT;
  for(int i=0;i<NPMT;i++) rPMT.push_back(r[i]);

  //print xy coordinates of PMTs
  if(0){
    for(int i=0;i<NPMT;i++)
      cout<<"PMT "<<i<<"\t (x,y) [mm] = ("<<rPMT[i].X()<<" , "<<rPMT[i].Y()<<")"<<endl;
    
    getchar();
  }
  return rPMT;
}

vector<G4ThreeVector> setPMTVector(const char* top_or_bottom){
  //new geometry
  //12 PMTs

  vector<G4ThreeVector> rPMT;
  G4double z;

  if(!strcmp(top_or_bottom,"bottom"))   z = BOTTOM_PMT_Z; 
  else if(!strcmp(top_or_bottom,"top")) z = TOP_PMT_Z;


  vector<TVector2> rPMT2D = setPMTVector2D();

  for(int i=0;i<rPMT2D.size();i++){
    rPMT.push_back(G4ThreeVector(rPMT2D[i].X(),rPMT2D[i].Y(),z));
  }

  return rPMT;
}


G4int isInFiducialVolume(G4ThreeVector pos){


  //check for the ring section
  if(pos.perp() > WLS_RINGSEC_INNER_RADIUS) return 0;
  G4double z = pos.z();
  if(z > WLS_RINGSEC_POS_Z + WLS_RINGSEC_HALF_HEIGHT ||
     z < WLS_RINGSEC_POS_Z - WLS_RINGSEC_HALF_HEIGHT )  return 0;

  //check for the linear section
  //assuming the linear section is perpendicular to the x-axis
  G4double maxx = WLS_LINSEC_POS_X-WLS_LINSEC_HALF_Y;  

  if(pos.x() > maxx) return 0;

  return 1;
}


G4int canUseLightMap(G4ThreeVector pos){


  //check for the ring section
  if(pos.perp() > WLS_RINGSEC_INNER_RADIUS) return 0;
  G4double z = pos.z();
  if(z > 580 || z < -500 )  return 0;

  //check for the linear section
  //assuming the linear section is perpendicular to the x-axis
  G4double maxx = WLS_LINSEC_POS_X-WLS_LINSEC_HALF_Y;  
  if(pos.x() > maxx) return 0;

  return 1;
}


G4int isInShieldedVolume(G4ThreeVector pos){


  if(pos.perp() > NEUTRON_SHIELD_INNER_RADIUS) return 0;
  G4double z = pos.z();
  if(z > NEUTRON_SHIELD_POS_Z+NEUTRON_SHIELD_HALF_HEIGHT ||
     z < NEUTRON_SHIELD_POS_Z-NEUTRON_SHIELD_HALF_HEIGHT )  return 0;

  return 1;
}


vector<G4ThreeVector> setGridWireVector(const char* xyAxis,vector<double>& wireHalfLength,
					double plateInnerR, double wirePitch, double wireOuterR){


  int nwires = (int)((plateInnerR-wirePitch/2)/wirePitch)+1;
  double R = plateInnerR;

  vector<G4ThreeVector> pos;

  wireHalfLength.clear();

  double z;
  
  if(!strcmp(xyAxis,"x"))      z = 0;
  else if(!strcmp(xyAxis,"y")) z = 2*wireOuterR+.1*mm;
  for(int i=-nwires;i<nwires;i++){
    double coord = CATHODE_WIRE_PITCH*(i+1./2);
    double halflength = sqrt(R*R - coord*coord);

    if(!halflength) continue;

    if(!strcmp(xyAxis,"x"))      pos.push_back(G4ThreeVector(coord,0.,z));
    else if(!strcmp(xyAxis,"y")) pos.push_back(G4ThreeVector(0.,coord,z));

    wireHalfLength.push_back(halflength);
  }


  return pos;
}


double getWLSThickness(double convEff){

  return ((convEff) > 0 && (convEff) < 1)?(-TMath::Log(1-(convEff))*(WLS_MEAN_ABSORPTION_LENGTH)):(WLS_THICKNESS_100_PERCENT_CONVERSION_EFFICIENCY);

}



vector<vector<string> > readTextFile_string(string filename,char delimiter){


  ifstream inputFile(filename.c_str());
  
  vector<vector<string> > output;
  
  if(!inputFile) return output;

  string line;
  while(getline(inputFile,line)){
    vector<string> row;
    string value;
    istringstream iss(line,istringstream::in);
    if(iss.str().empty()) continue;
    while(getline(iss,value,delimiter)){
      row.push_back(value);
    }

    output.push_back(row);
  }


  return output;
}



vector<vector<double> > readTextFile_float(string filename,char delimiter){

  ifstream inputFile(filename.c_str());
  
  vector<vector<double> > output;
  
  if(!inputFile) return output;

  int count=0;
  
  string line;
  while(getline(inputFile,line)){
    vector<double> row;
    string value;
    istringstream iss(line,istringstream::in);
    if(iss.str().empty()) continue;
    while(getline(iss,value,delimiter)){
      row.push_back(atof(value.c_str()));
    }

    output.push_back(row);
  }
  
  return output;
}



