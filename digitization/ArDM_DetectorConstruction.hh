#ifndef _ARDM_DETECTORCONSTRUCTION_
#define _ARDM_DETECTORCONSTRUCTION_ 1


#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VSensitiveDetector.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Colour.hh"

using namespace std;

class ArDM_DetectorConstruction : public G4VUserDetectorConstruction{

protected:
  G4NistManager* fNist;
  G4Material*    fVacuum;
  G4Material*    fWorldMat;
  G4Material*    fTankMat;
  G4Material*    fLidMat;
  G4Material*    fDetMat;
  G4Material*    fLAr;
  G4Material*    fGAr;
  G4Material*    fPMTMat;
  G4Material*    fWLS_ringSec_Mat;
  G4Material*    fWLS_linSec_Mat;
  G4Material*    fWLSMat;
  G4Material*    fWLSSupportMat;
  G4Material*    fBottomSideReflectorMat;
  G4Material*    fTopSideReflectorMat;
  G4Material*    fPMTSupportMat;
  G4Material*    fTetratex;
  G4Material*    fPolyethylene;
  G4Material*    fTeflon;
  G4Material*    fCathodeGridMat;
  G4Material*    fPMTCathodeMat;
  G4Material*    fPMTElectrodeMat;
  G4Material*    fPMTBaseMat;
  G4Material*    fHVrMat;
  G4Material*    fPerlite;
  G4Material*    fBorotron;



  G4VPhysicalVolume* fWorldPhys;
  G4VPhysicalVolume* fTankPhys;
  G4VPhysicalVolume* fBottomLidPhys;
  G4VPhysicalVolume* fTopLidPhys;     //used in the approximation : tank = tank cylinder + btm lid + top lid

  //used in the approximation : tank = dewar + top flange, dewar = tank cylinder + curved btm part.
  //<-- this approximation is closer to reality than the above approximation (tank = tank cylinder + top lid + btm lid)
  //fTopFlangePhys should have the same meaning as fTopLidPhys, but for the sake of clarity, use fTopFlangePhys instead of fTopLidPhys.
  G4VPhysicalVolume* fTopFlangePhys;  
  G4VPhysicalVolume* fHVrBar1Phys;  
  G4VPhysicalVolume* fHVrBar2Phys;  
  G4VPhysicalVolume* fLArColPhys;
  G4VPhysicalVolume* fGArColPhys;
  G4VPhysicalVolume* fDetPhys;
  G4VPhysicalVolume* fWLS_ringSecPhys;
  G4VPhysicalVolume* fWLS_linSecPhys;
  G4VPhysicalVolume* fWLSPhys;
  G4VPhysicalVolume* fWLSSupportPhys;
  G4VPhysicalVolume* fBottomSideReflectorPhys; //cylinder cover above the protection grid, w/o wls-coating
  G4VPhysicalVolume* fTopSideReflectorGArPhys; 
  G4VPhysicalVolume* fTopSideReflectorLArPhys; 
  G4VPhysicalVolume* fBottomSideReflectorCoatPhys; //cylinder cover above the protection grid, w/o wls-coating
  G4VPhysicalVolume* fTopSideReflectorCoatGArPhys; 
  G4VPhysicalVolume* fTopSideReflectorCoatLArPhys; 
  G4VPhysicalVolume* fTopPMTSupportPhys;   //support tructure for PMT
  G4VPhysicalVolume* fBottomPMTSupportPhys;//support tructure for PMT
  G4VPhysicalVolume* fTopPMTSupportCoatPhys;   //coating of support tructure for PMT
  G4VPhysicalVolume* fBottomPMTSupportCoatPhys;//coating of support tructure for PMT
  G4VPhysicalVolume* fSourceHolder;
  G4VPhysicalVolume* fSource;
  G4VPhysicalVolume* fSourceCoating;
  G4VPhysicalVolume* fCathodeGridPhys;
  G4VPhysicalVolume* fProtectionGridPhys;
  G4VPhysicalVolume* fPerlitePhys;
  
  

  //use replica / G4VPVParameterisation instead of vector<...> !!
  //for the time being, just use vector<..> for simplicity
  vector<G4VPhysicalVolume*> fBottomPMTArrayPhys;
  vector<G4VPhysicalVolume*> fTopPMTArrayPhys;
  vector<G4VPhysicalVolume*> fBottomPMTCathodeArrayPhys;
  vector<G4VPhysicalVolume*> fTopPMTCathodeArrayPhys;
  vector<G4VPhysicalVolume*> fBottomPMTElectrodeArrayPhys;
  vector<G4VPhysicalVolume*> fTopPMTElectrodeArrayPhys;
  vector<G4VPhysicalVolume*> fBottomPMTCoatArrayPhys;
  vector<G4VPhysicalVolume*> fTopPMTCoatArrayPhys;
  vector<G4VPhysicalVolume*> fBottomPMTBaseArrayPhys;
  vector<G4VPhysicalVolume*> fTopPMTBaseArrayPhys;
  G4VPhysicalVolume* fBottomPMT_ROPhys;
  G4VPhysicalVolume* fTopPMT_ROPhys;  

  G4VPhysicalVolume* fNeutronShieldPhys;

  G4VPhysicalVolume* fReflectorVolume;     //volume confined by the reflector
  

  G4LogicalBorderSurface* fLAr_tank_surf;
  G4LogicalBorderSurface* fLAr_BottomSideReflector_surf; 
  G4LogicalBorderSurface* fLAr_TopSideReflector_surf; 
  G4LogicalBorderSurface* fGAr_TopSideReflector_surf; 
  G4LogicalBorderSurface* fLAr_BottomSideReflectorCoat_surf; 
  G4LogicalBorderSurface* fBottomSideReflectorCoat_LAr_surf; 
  G4LogicalBorderSurface* fLAr_TopSideReflectorCoat_surf; 
  G4LogicalBorderSurface* fGAr_TopSideReflectorCoat_surf; 
  G4LogicalBorderSurface* fTopSideReflectorCoat_LAr_surf; 
  G4LogicalBorderSurface* fTopSideReflectorCoat_GAr_surf; 
  G4LogicalBorderSurface* fBottomSideReflector_BottomSideReflectorCoat_surf; 
  G4LogicalBorderSurface* fTopSideReflectorLAr_TopSideReflectorCoat_surf; 
  G4LogicalBorderSurface* fTopSideReflectorGAr_TopSideReflectorCoat_surf; 
  G4LogicalBorderSurface* fLAr_bottomLid_surf;
  G4LogicalBorderSurface* fLAr_WLS_surf;
  G4LogicalBorderSurface* fWLS_LAr_surf;
  G4LogicalBorderSurface* fLAr_GAr_surf;
  G4LogicalBorderSurface* fGAr_LAr_surf;
  G4LogicalBorderSurface* fGAr_tank_surf;
  G4LogicalBorderSurface* fGAr_topLid_surf;
  G4LogicalBorderSurface* fWLS_WLSSupport_surf;
  G4LogicalBorderSurface* fLAr_WLSSupport_surf;
  G4LogicalBorderSurface* fGAr_topPMTSupport_surf;  
  G4LogicalBorderSurface* fLAr_bottomPMTSupport_surf;  
  G4LogicalBorderSurface* fPMTSupportCoat_topPMTSupport_surf;  
  G4LogicalBorderSurface* fPMTSupportCoat_bottomPMTSupport_surf; 
  G4LogicalBorderSurface* fTopPMTSupportCoat_GAr_surf;  
  G4LogicalBorderSurface* fBottomPMTSupportCoat_LAr_surf;  
  G4LogicalBorderSurface* fGAr_topPMTSupportCoat_surf;  
  G4LogicalBorderSurface* fLAr_bottomPMTSupportCoat_surf;  
  G4LogicalBorderSurface* fLAr_protectionGrid_surf;  
  G4LogicalBorderSurface* fLAr_cathodeGrid_surf;  

  vector<G4LogicalBorderSurface*> fBottomPMT_PMTCoat_surf;  
  vector<G4LogicalBorderSurface*> fPMTCoat_BottomPMT_surf;  
  vector<G4LogicalBorderSurface*> fTopPMT_PMTCoat_surf;
  vector<G4LogicalBorderSurface*> fPMTCoat_TopPMT_surf;
  vector<G4LogicalBorderSurface*> fLAr_bottomPMTCoat_surf;
  vector<G4LogicalBorderSurface*> fGAr_topPMTCoat_surf;
  vector<G4LogicalBorderSurface*> fBottomPMTCoat_LAr_surf;
  vector<G4LogicalBorderSurface*> fTopPMTCoat_GAr_surf;
  vector<G4LogicalBorderSurface*> fLAr_PMT_surf;
  vector<G4LogicalBorderSurface*> fGAr_PMT_surf;
  vector<G4LogicalBorderSurface*> fPMT_LAr_surf;
  vector<G4LogicalBorderSurface*> fPMT_GAr_surf;

  vector<G4LogicalBorderSurface*> fPMTCathode_bottomPMT_surf;
  vector<G4LogicalBorderSurface*> fPMTCathode_topPMT_surf;
  vector<G4LogicalBorderSurface*> fBottomPMT_PMTCathode_surf;
  vector<G4LogicalBorderSurface*> fTopPMT_PMTCathode_surf;

   

  void addTank();      //** new ** tank = cylinder + curved bottom part
  void addTopFlange();
  void addHVr();       //high voltage resistors
  void addLArColumn();
  void addGArColumn();  
  void addPMT(G4String bottom_or_top="bottom");
  void addPMTCoat(G4String bottom_or_top="bottom");
  void addPMTCoat_universal_pmtCoat_thickness(G4String bottom_or_top="bottom");
  void addPMTCoat_individual_pmtCoat_thickness(G4String bottom_or_top="bottom");
  void addPMTCathode(G4String bottom_or_top="bottom");
  void addPMTElectrode(G4String bottom_or_top="bottom");
  //void addROGeom(G4String bottom_or_top="bottom"); //just for visual demonstration !

  void addPMTBase(G4String bottom_or_top="bottom");

  void addWLS();
  void addWLSSupport();
  void addBottomSideReflector();
  void addTopSideReflector();
  void addBottomSideReflectorCoat();
  void addTopSideReflectorCoat();
  void addPMTSupport(G4String bottom_or_top="bottom");
  void addPMTSupportCoat(G4String bottom_or_top="bottom");
  void addNeutronShield();

  void addSource();
  void addSourceHolder();
  void addSourceCoating();

  void addCathodeGrid();
  void addProtectionGrid();

  void addEfield(); 

  void addPerlite();


  vector<G4VPhysicalVolume*> placePMT(G4LogicalVolume* fPMTLog,G4VPhysicalVolume* fMotherPhys,vector<G4ThreeVector> rPMT,const char* bottom_or_top);
  vector<G4VPhysicalVolume*> placePMT(vector<G4LogicalVolume*> fPMTLog,G4VPhysicalVolume* fMotherPhys,vector<G4ThreeVector> rPMT,const char* bottom_or_top);

  G4LogicalVolume* constructGrid(G4Material* fMat,double plateInnerR,double plateOuterR,
				 double plateHalfThickness,double wireInnerR,double wireOuterR,
				 double wirePitch,const char* name,int attribute,G4Colour color);
  


  G4LogicalBorderSurface* build_LAr_tank_surf();
  G4LogicalBorderSurface* build_LAr_BottomSideReflector_surf(); 
  G4LogicalBorderSurface* build_LAr_TopSideReflector_surf();   //top side refl. can be both in LAr
  G4LogicalBorderSurface* build_GAr_TopSideReflector_surf();   //and GAr

  G4LogicalBorderSurface* build_LAr_BottomSideReflectorCoat_surf(G4String order="LAr_BottomSideReflectorCoat"); 
  G4LogicalBorderSurface* build_LAr_TopSideReflectorCoat_surf(G4String order="LAr_TopSideReflectorCoat");   
  G4LogicalBorderSurface* build_GAr_TopSideReflectorCoat_surf(G4String order="GAr_TopSideReflectorCoat");   
  
  G4LogicalBorderSurface* build_BottomSideReflector_BottomSideReflectorCoat_surf(); 
  G4LogicalBorderSurface* build_TopSideReflectorGAr_TopSideReflectorCoat_surf();   
  G4LogicalBorderSurface* build_TopSideReflectorLAr_TopSideReflectorCoat_surf();   
  
  G4LogicalBorderSurface* build_LAr_bottomLid_surf();
  G4LogicalBorderSurface* build_LAr_WLS_surf(G4String order="LAr_WLS"); 
  G4LogicalBorderSurface* build_LAr_GAr_surf(G4String order="LAr_GAr"); 
  G4LogicalBorderSurface* build_GAr_tank_surf(); 
  G4LogicalBorderSurface* build_GAr_topLid_surf(); 
  G4LogicalBorderSurface* build_WLS_WLSSupport_surf();
  G4LogicalBorderSurface* build_LAr_WLSSupport_surf();
  vector<G4LogicalBorderSurface*> build_PMT_PMTCoat_surf(G4String topBottom="bottom",
							 G4String order="Ar_PMTCoat");
  vector<G4LogicalBorderSurface*> build_Ar_PMTCoat_surf(G4String topBottom="bottom",
							G4String order="PMT_PMTCoat");
  vector<G4LogicalBorderSurface*> build_Ar_PMT_surf(G4String topBottom="bottom",
						    G4String order="Ar_PMT");

  vector<G4LogicalBorderSurface*> build_PMTCathode_PMT_surf(G4String topBottom="bottom",
							    G4String order="PMTCathode_PMT");

  G4LogicalBorderSurface* build_Ar_PMTSupport_surf (G4String topBottom="bottom");  
  G4LogicalBorderSurface* build_PMTSupportCoat_PMTSupport_surf (G4String topBottom="bottom");  
  G4LogicalBorderSurface* build_Ar_PMTSupportCoat_surf (G4String topBottom="bottom",
							G4String order="Ar_PMTSupportCoat");
 
  G4LogicalBorderSurface* build_LAr_protectionGrid_surf ();
  G4LogicalBorderSurface* build_LAr_cathodeGrid_surf ();
 

  void build_surfaces();

  vector<G4LogicalSkinSurface*> build_LAr_PMT_surf();
  vector<G4LogicalSkinSurface*> build_GAr_PMT_surf();
  

  
  //void setMatPropTab_detMat(); //LAr
  G4Material* getTankMat(); //stainless steel
  G4Material* getTetratex(); 
  G4Material* getPolyethylene(); 
  G4Material* getTeflon();  
  G4Material* getHVrMat();
  G4Material* getPerlite(); 
  G4Material* getBorotron(); 
  G4Material* getFR4(); //material of PMT base

  void setMatPropTab_WLS(); //TPB
  void setMatPropTab_LAr(); //LAr
  void setMatPropTab_GAr(); //GAr
  //set the scintillation properties of argon (liquid and gas)
  void setArScintProperty(G4MaterialPropertiesTable* propTab,G4String medium="LAr");
  void setMatPropTab_PMTMat(); //PMTMat


  G4double GArRefIndex(G4double wavelength); //wavelength in micrometer !!
  G4double LArRefIndex(G4double wavelength); //wavelength in micrometer !!
  G4double GArEpsilon(G4double wavelength);  //wavelength in micrometer !!
  G4double LArEpsilon(G4double wavelength);  //wavelength in micrometer !!
  G4double RayleighAttenuationLength_LAr(G4double wavelength); //wavelength in micrometer !!
  G4double RayleighAttenuationLength_GAr(G4double wavelength); //wavelength in micrometer !!

  G4double getPhotonWavelength(G4double photonEnergy); //the returned wavelength is measured in micrometer !!

  G4VSensitiveDetector* getSD(G4String bottom_or_top);
  vector<G4VPhysicalVolume*> getPMTArrayPhys(G4String bottom_or_top);




  void verboseInfo(string detectorComponent="");





  //end test



public:
  ArDM_DetectorConstruction();
  ~ArDM_DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  
};







#endif // _ARDM_DETECTORCONSTRUCTION_
