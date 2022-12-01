#include "ArDM_SensitivePMT.hh"
#include "preparation.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"


const G4double ArDM_SensitivePMT :: fTimeSample   = TIMESAMPLE;
const G4int    ArDM_SensitivePMT :: fNTimeSamples = NTIMESAMPLE;
ArDM_Digitizer* ArDM_SensitivePMT:: digitizer     = NULL;

ArDM_SensitivePMT::ArDM_SensitivePMT(G4int npmt,G4String detName,G4String hitsColName,G4String top_or_bottom)
  :fnpmt(npmt),
   fDetName(detName),
   fHitsColName(hitsColName),
   fTop_or_bottom(top_or_bottom),
   G4VSensitiveDetector(detName){
  
  if(fTop_or_bottom == "top") fNormal.set(0.,0.,-1.);
  else fNormal.set(0.,0.,1.);
  
  
  collectionName.insert(hitsColName);
  fPMTvector=setPMTVector(fTop_or_bottom);

  double z = DEFAULTVALUE;
  if(fTop_or_bottom == "top") z = APPROX_TOP_PMT_CENTER_OF_TOP_SPHERE;
  else                        z = APPROX_BTM_PMT_CENTER_OF_TOP_SPHERE;


  //fPMTvector is not used to place the PMTarray in any physical volume,
  //it is rather used to calculate the angle, at which the photon arrives on the PMT surface,
  //so it should points to the center of the sensitive spherical part of the PMT,
  //and not to the center of the pmt as described in the approx_pmt_... part in preparation.hh,
  //or in preration.cc :: constructApproxPMT(...).
  //
  //the xy coordinate are the same for both center of sensitive part and center of pmt itself,
  //but z-coord. is different.
  //
  //--> so now, reset the z-coord.
  //
  for(int i=0;i<fPMTvector.size();i++) fPMTvector[i].setZ(z);

  digitizer = new ArDM_Digitizer("PMTDigitizer");
    
  
}  

  
  
ArDM_SensitivePMT::~ArDM_SensitivePMT(){
  fPMTvector.clear();
}

G4double ArDM_SensitivePMT::detectionProb(G4double phi){
  
  phi *= 180/PI; //convert from radian to degree.
  //cout<<"in ArDM_SensitivePMT::detectionProb(G4double phi) ... phi = "<<phi<<endl;

  if(phi > 46.5) return 0;
  return (2587 + 1134/(1 + (phi/27)*(phi/27)))/3722;  
}


G4double ArDM_SensitivePMT::detectionProb(G4ThreeVector hitpos, G4ThreeVector rPMT){
  return detectionProb(fNormal.angle(hitpos-rPMT));
}


G4double ArDM_SensitivePMT::detectionProb(G4ThreeVector hitpos,G4int pmtID){
  return detectionProb(hitpos,fPMTvector[ (pmtID < NTOPPMT) ? pmtID : (pmtID - NTOPPMT) ]);
}


G4int ArDM_SensitivePMT::detected(G4ThreeVector hitpos,G4ThreeVector rPMT){
  
  if(drand48() < detectionProb(hitpos,rPMT)) return 1;
  else return 0;
}


G4int ArDM_SensitivePMT::detected(G4ThreeVector hitpos,G4int pmtID){

#if VERBOSE
  cout<<"pmtID "<<pmtID<<"\t fPMTvector.size() "<<fPMTvector.size()<<endl;
#endif //VERBOSE

  return detected(hitpos,fPMTvector[(pmtID < NTOPPMT) ? pmtID : (pmtID-NTOPPMT) ]);
}





void ArDM_SensitivePMT::Initialize(G4HCofThisEvent* hc){
  fHitsCol = new ArDM_PMTHitsCollection(fDetName,fHitsColName);//"PMTarray","PMTHits");
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsColName);
  hc->AddHitsCollection(hcID,fHitsCol);
  verboseLevel = 0;

#if VERBOSE
  cout<<"hitsColName = "<<fHitsColName <<"\t id "<<hcID<<"\t nHitsColInThisEvt "<<G4SDManager::GetSDMpointer()->GetHCtable()->entries() <<endl;getchar();
#endif //VERBOSE

  return;
}




G4bool ArDM_SensitivePMT::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist){

  if(!aStep->GetTrack()->GetParticleDefinition()->GetParticleName().contains("opticalphoton")) return true;

  if(aStep->GetTrack()->GetDynamicParticle()->GetTotalEnergy() > PHOTON_ENERGY_THRESHOLD){
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);    
    return true;
  }

  G4int photonDetected = 0;
    
  ArDM_PMTHit* newhit = new ArDM_PMTHit;
  G4Track* track = aStep->GetTrack();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();

  //G4int globalTime = (int)(track->GetGlobalTime()/ns/fTimeSample);
  G4int globalTime = track->GetGlobalTime()/ns;
  G4String volName = track->GetVolume()->GetName();

  //10 = length of string "pmtCathode", take the last 2 character
  //e.g. : volName = "pmtCathode3"  --> pmtIDstr = volName.substr(10,2) = "3"  --> pmtID = 3
  //       volName = "pmtCathode"   --> pmtIDstr = volName.substr(10,2) = "17" --> pmtID = 17


  G4String pmtIDstr = volName(10,2); //<--> volName.substr(10,2)

  G4int pmtID = atoi(pmtIDstr.data());

  //check if the photon is detected
  //store only detected photons !
  photonDetected = detected(hitpos,pmtID);

  newhit->setHitGlobalTime(globalTime);
  newhit->setHitPMTid(pmtID);
  newhit->setIsDetected(photonDetected);

  fHitsCol->insert(newhit);


#if VERBOSE
  cout<<"in ArDM_SensitivePMT::ProcessHits(..) ... volName "<<volName<<"\t pmtID "<<pmtID<<"\t isDetected "<<photonDetected <<endl;
#endif //VERBOSE

  //photon falls onto PMT's surface --> kill it !
  aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

  return true;
}



void ArDM_SensitivePMT::EndOfEvent(G4HCofThisEvent* hc){

  digitizer->Digitize_simple(fHitsColName);

#if VERBOSE
  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
  ArDM_PMTHitsCollection* hitCol0 = (ArDM_PMTHitsCollection*) (digiMan->GetHitsCollection(0));
  ArDM_PMTHitsCollection* hitCol1 = (ArDM_PMTHitsCollection*) (digiMan->GetHitsCollection(1));

  cout<<hitCol0->GetName()<<"\t nentries "<<hitCol0->entries()
      <<"\t "<<hitCol1->GetName()<<"\t nentries "<<hitCol1->entries()
      <<endl;

  getchar();

#endif //VERBOSE
  

  return;
}




