
#include "preparation.hh"
#include "ArDM_Digitizer.hh"

#include "G4DigiManager.hh"
#include "ArDM_PMTHit.hh"
#include "ArDM_Analysis.hh"

#include "TRandom.h"

#include "assert.h"


#define ARDM_PMT_RISETIME (12*ns)
#define ARDM_TRANSIT_TIME (68*ns) 
#define ARDM_SIGMA_TRANSIT_TIME (2.8*ns) 
#define ARDM_MEAN_PEAK_WIDTH (28.72*ns)
#define ARDM_MEAN_FWHM (23.2 *ns)
#define ARDM_PMT_TAU (7.51*ns)


//taken from run1123
const double PMT_RESPONSE[] = {106.833, 104.107, 103.06,  102.068,
			       102.326, 103.905, 104.418, 101.413,
			       102.251, 102.041, 101.803, 101.762,
			       103.148, 103.382, 101.084, 100.6,
			       97.5557, 101.384, 68.9673, 103.344,
			       103.302, 101.969, 102.893, 102.24
                              };

const double PMT_RESPONSE_SIGMA[] = {41.9689, 36.0583, 33.6824, 34.4087,
				     34.154,  38.3709, 40.6972, 34.7465,
				     41.4877, 39.1481, 41.4032, 40.607,
				     37.0913, 37.4369, 33.3611, 33.6224,
				     33.5627, 37.5153, 13.1009, 34.1722,
				     34.2261, 35.532,  35.1569, 30.8242
                                     };





//from real data, run 7354
const double PMT_DARK_COUNT_RATE[] = {2.51932,0.393232,0.615367,0.379898,3.0342,0.543858,
				      0.284677,0.322333,10.4282,1.38244,0.19232,0.636074,
				      0.354449,0.252523,0.292142,0.214059,0.257819,0.293832,
				      0.20927,0.170534,0.34703,0.311346,0.168609,0.194339
                                     };




const double PMT_WHITE_NOISE_SIGMA[] = {0.873489,1.06833,0.828103,0.890891,0.824046,0.799382,
					0.922938,0.793321,0.790785,0.792672,0.847572,0.840875,
					0.921675,0.896091,0.868562,0.891789,0.785457,0.813048,
					0.809918,0.883621,0.797634,0.932994,0.9817,0.898223
                                       };



int ArDM_Digitizer::eNTimeSamples = NTIMESAMPLE;

ArDM_Digitizer::ArDM_Digitizer(G4String name)
  :G4VDigitizerModule(name){;}



ArDM_Digitizer::~ArDM_Digitizer(){;}






void ArDM_Digitizer::reset(float* array, int nentries){

  for(int i=0;i<nentries;i++) array[i] = 0;
  return;
}





void ArDM_Digitizer::Digitize(G4HCofThisEvent* hc){
  
  Digitize_simple(hc);
  //Digitize_pmtRes(hc);
  
  return; 
}





void ArDM_Digitizer::Digitize(){;}




void ArDM_Digitizer::Digitize_simple(ArDM_PMTHit* hit){
  /*
  if(!hit->getIsDetected()) return;


  int pmtid = hit->getHitPMTid();
  
  ArDM_Analysis* ana = ArDM_Analysis::getInstance();

#if TEST_BRANCH0
  
  if(pmtid < 12) ana->fNphotonDetectedTop++;
  else           ana->fNphotonDetectedBottom++;
  
#if TEST_BRANCH00
  
  int globalTime = (int)(hit->getHitGlobalTime()/TIMESAMPLE);
  
  ana->fTotNphotonPMT[pmtid]++;

  if(globalTime < NTIMESAMPLE) ana->fNphotonPMT[pmtid][globalTime]++;
    
    
#endif //TEST_BRANCH00
    
  

#endif //TEST_BRANCH0

  */

  return;
}









void ArDM_Digitizer::Digitize_simple(ArDM_PMTHitsCollection* hitsCol){

  int nhits = hitsCol->entries();
  for(int i=0; i < nhits; i ++ )   Digitize_simple((*hitsCol)[i]);
  
  return;
}







void ArDM_Digitizer::Digitize_simple(G4String hitsColName){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();

#if VERBOSE
  cout<<"hitsColName "<<hitsColName<<endl; //getchar();
#endif //VERBOSE

  G4int hitColID = digiMan->GetHitsCollectionID(hitsColName);

  ArDM_PMTHitsCollection* hitCol = (ArDM_PMTHitsCollection*) (digiMan->GetHitsCollection(hitColID));
  if(hitCol) Digitize_simple(hitCol);

  return; 
}





void ArDM_Digitizer::Digitize_simple(G4HCofThisEvent* hc){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
  
  for(int hci=0; hc->GetNumberOfCollections();hci++){
    ArDM_PMTHitsCollection* hitCol = (ArDM_PMTHitsCollection*) (digiMan->GetHitsCollection(hci));
    if(hitCol) Digitize_simple(hitCol);
  }

  return; 
}








void ArDM_Digitizer::Digitize_pmtRes_try1(ArDM_PMTHit* hit,float* rawDataArray){

  //simulting PMT (voltage) response to one single photon hit
  //right flank : exponentially decaying funtion,
  //left  flank : gaussian


  if(!hit->getIsDetected()) return;

  double peakStartTime = hit->getHitGlobalTime();
  double peakTime      = peakStartTime + ARDM_PMT_RISETIME;
  double peakStopTime  = peakStartTime + ARDM_MEAN_PEAK_WIDTH;

  int pmtid = hit->getHitPMTid();
  double pmtRes = gRandom->Gaus(PMT_RESPONSE[pmtid],PMT_RESPONSE_SIGMA[pmtid]);
  while(pmtRes <= 0  ) pmtRes = gRandom->Gaus(PMT_RESPONSE[pmtid],PMT_RESPONSE_SIGMA[pmtid]);


  double tau = ARDM_PMT_TAU;
  double Vo = pmtRes * TIMESAMPLE / (tau + sqrt(PI/2) * ARDM_SIGMA_TRANSIT_TIME ) ;

  int peakStartTimeInt = (int) (peakStartTime / TIMESAMPLE);
  int peakTimeInt      = (int) (peakTime      / TIMESAMPLE);
  int peakStopTimeInt  = (int) (peakStopTime  / TIMESAMPLE);

#if VERBOSE
  cout<<"peakStartTimeInt "<<peakStartTimeInt
      <<"\t peakTimeInt "<<peakTimeInt
      <<"\t peakStopTimeInt "<<peakStopTimeInt
      <<"\t Vo "<<Vo
      <<"\t pmtRes "<<pmtRes
      <<endl;
  
#endif //VERBOSE

  for(int i=peakStartTimeInt; i < peakTimeInt; i++) rawDataArray[i] += Vo * TMath::Gaus((i+1./2)*TIMESAMPLE, peakTime,ARDM_SIGMA_TRANSIT_TIME,1);
  for(int i=peakTimeInt; i <= peakStopTimeInt; i++) rawDataArray[i] += Vo * exp( - ((i+1./2)*TIMESAMPLE - peakTime)/tau);

  return ;

}






void ArDM_Digitizer::Digitize_pmtRes_try2(ArDM_PMTHit* hit,float* rawDataArray){

  //simulting PMT (voltage) response to one single photon hit
  //pmtResponse = log-normal distribution
  //V = Vo * exp( - log(t/tau) * log(t/tau) / 2 / sigma /sigma  )

  //or more precisely :
  //
  //V = Vo * exp( - log( (t-peakStartTime)/tau) * log( (t-peakStartTime)/tau) / 2 / sigma / sigma )
  //
  //can try : tau ~ FWHM ~ 4*5.8 ns , sigma ~ .25 (from real data fitting)
  //
  //Vo varies from pmt to pmt
  //tau, sigma can be assumed to be the same for all PMTs


  if(!hit->getIsDetected()) return;

  double peakStartTime = hit->getHitGlobalTime();
  double peakStopTime  = peakStartTime + ARDM_MEAN_PEAK_WIDTH;

  int pmtid = hit->getHitPMTid();
  double pmtRes = gRandom->Gaus(PMT_RESPONSE[pmtid],PMT_RESPONSE_SIGMA[pmtid]);
  while(pmtRes <= 0  ) pmtRes = gRandom->Gaus(PMT_RESPONSE[pmtid],PMT_RESPONSE_SIGMA[pmtid]);

  //double Vo = pmtRes * TIMESAMPLE / (tau + sqrt(PI/2) * ARDM_SIGMA_TRANSIT_TIME ) ;
  double tau = ARDM_MEAN_FWHM;
  double sigma = .25;

  double Vo = pmtRes * TIMESAMPLE / (tau * sigma * sqrt(2*PI) * exp(sigma*sigma/2));

  int peakStartTimeInt = (int) (peakStartTime / TIMESAMPLE);
  int peakStopTimeInt  = (int) (peakStopTime  / TIMESAMPLE);

#if VERBOSE
  cout<<"peakStartTimeInt "<<peakStartTimeInt
      <<"\t peakStopTimeInt "<<peakStopTimeInt
      <<"\t Vo "<<Vo
      <<"\t pmtRes "<<pmtRes
      <<endl;
  
#endif //VERBOSE

  for(int i=peakStartTimeInt; i <= peakStopTimeInt; i++) 
    rawDataArray[i] += Vo * TMath::Gaus(log( ( (i+1./2)*TIMESAMPLE - peakStartTime) / tau),0,sigma);


  return ;
}





void ArDM_Digitizer::Digitize_pmtRes(ArDM_PMTHit* hit,float* rawDataArray){

  //Digitize_pmtRes_try1(hit,rawDataArray);
  Digitize_pmtRes_try2(hit,rawDataArray);
  return;
}






void ArDM_Digitizer::Digitize_pmtRes(ArDM_PMTHit* hit){

  /*
#if TEST_BRANCH0
#if TEST_BRANCH00


  //the electron cloud from dynodes needs some transit time ts to reach the anode,
  //the transit time ts follows a gaussian distribution with mean mean_ts and width sigma_ts 
  //given by the technical sheet of the PMT.
  //
  //--> how to simulate the voltage signal of a single photon ? 
  //
  //i.   photon hits PMT surface at time t0
  //ii.  electron cloud hits the anode at time t1 = t0 + gauss(mean=transitTime,sigma=sigma_ts)
  //       --> peakTime = t1
  //iii. peakStartTime = peakTime - riseTime
  //iv.  stopTime = t2 = peakTime + mean(stopTime - peakTime) <-- mean(stopTime - peakTime) from real data
  //v.   get pmtRes = gauss(mean_pmtRes,sigma_pmtRes) <-- mean, sigma from real data for each PMT
  //vi.  now we have to throw random numbers to generate the amplitude of the voltage signal for each time sample,
  //     so that the peakIntegral we get from that is equal to the pmtRes in (v.)
  //
  //since the transit time of the electron cloud is small compared to the processing time of the amplifier after that,
  //in order to simplify things, instead of proceed like described above, we can do the following :
  //
  //i.   photon hits surface at time t0 =: startTime
  //ii.  t1 = t0 + rise time  =: peakTime
  //iii. t2 = t1 + mean(stopTime - peakTime) =: stopTime
  //iv.  pmtRes = gauss(mean_pmtRes,sigma_pmtRes)
  //v.   now we have to simulate the amplitude of the voltage signal :
  //      1. left  flank (startTime < t < peakTime) = (half-)gaussian with mean = peakTime, sigma = transit time spread
  //      2. right flank (peakTime  < t < stopTime) = exponential decay with tau = 1/RC <-- R,C from the circuit of the PMT, emperical value.
  //      
  //vi.  left blank  = Vo * exp( - (t-peakTime) * (t-peakTime) / 2/sigma_ts/sigma_s)
  //     right blank = Vo * exp( - (t-peakTime)/tau)
  //
  //vii. integral(left blank) = int(left blank)_peakStartTime ^peakTime ~= int(left blank)_minusInf ^peakTime = Vo * sqrt(2pi) * sigma_ts
  //     integral(left blank) = int(left blank)_peakTime ^peakStopTime ~= int(left blank)_peakTime ^inf = Vo * tau
  // 
  //     --> integral = pmtRes = Vo * (tau + sqrt(2pi) * sigma_ts)
  //     --> Vo = pmtRes / (tau + sqrt(2pi) * sigma_ts)
  //
  //


  //in class ArDM_Analysis :
  // float fRawData[NTOPPMT+NBOTTOMPMT][NTIMESAMPLEBOTTOM];


  int pmtid = hit->getHitPMTid();
  Digitize_pmtRes(hit,ArDM_Analysis::getInstance()->fRawData[pmtid]);

#endif //TEST_BRANCH0
#endif //TEST_BRANCH00

  */
  
  return;

}









void ArDM_Digitizer::Digitize_pmtRes(ArDM_PMTHitsCollection* hitcol){

  int nhits = hitcol->entries();
  for(int i=0;i<nhits;i++) Digitize_pmtRes((*hitcol)[i]);

  return;
}




void ArDM_Digitizer::Digitize_pmtRes(G4String hitsColName){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();

  G4int hitColID = digiMan->GetHitsCollectionID(hitsColName);
  ArDM_PMTHitsCollection* hitCol = (ArDM_PMTHitsCollection*) (digiMan->GetHitsCollection(hitColID));
  Digitize_pmtRes(hitCol);
  return;
}





void ArDM_Digitizer::Digitize_pmtRes(G4HCofThisEvent* hc){

  G4DigiManager* digiMan = G4DigiManager::GetDMpointer();
  
  for(int hci=0; hci < hc->GetNumberOfCollections(); hci++){
    ArDM_PMTHitsCollection* hitCol = (ArDM_PMTHitsCollection*) (digiMan->GetHitsCollection(hci));
    Digitize_pmtRes(hitCol);
  }

  return;
}









void ArDM_Digitizer::Digitize_darkCurrent(int pmtid,float* rawDataArray){

  assert(pmtid >=0 && pmtid < 24);

  int nelectrons = gRandom->Poisson(PMT_DARK_COUNT_RATE[pmtid]);

  //cout<<"nelectrons "<<nelectrons<<endl;

  if(nelectrons <= 0) return;

  double timespan = TIMESAMPLE * eNTimeSamples;
  for(int i=0; i < nelectrons;i++){
    double time = gRandom->Uniform(0,timespan);
    
    ArDM_PMTHit* hit;
    hit->setHitGlobalTime(time);
    hit->setHitPMTid(pmtid);
    hit->setIsDetected(1);

    Digitize_pmtRes(hit,rawDataArray);  
  }

  return;
}








void ArDM_Digitizer::Digitize_darkCurrent(int pmtid){

  //probability that an electron gets out of the PMT cathode by itself is given by dark count rate
  //the array PMT_DARK_COUNT_RATE above gives the mean number of "dark electrons"
  //in a time span, which corresponds to the timespan we would record for a normal event ( = 2048 samples * 4ns/sample)
  //--> from real data
  //
  //--> throw a poisson distribution around that mean value to get the number of dark electrons to be emitted.
  //
  //pulse shape for "dark electron" is the same for the one induced by an optical photon.
  //the probability that a dark electron appears is uniform over the whole time range of an optical-photon-event.


  assert(pmtid >=0 && pmtid < 24);

  int nelectrons = gRandom->Poisson(PMT_DARK_COUNT_RATE[pmtid]);

  if(nelectrons <= 0) return;

  double timespan = TIMESAMPLE * eNTimeSamples;
  for(int i=0; i < nelectrons;i++){
    double time = gRandom->Uniform(0,timespan);
    
    ArDM_PMTHit* hit;
    hit->setHitGlobalTime(time);
    hit->setHitPMTid(pmtid);
    hit->setIsDetected(1);

    Digitize_pmtRes(hit);  
  }

  return;
}









void ArDM_Digitizer::Digitize_whiteNoise(int pmtid,float* rawDataArray){

  assert(pmtid >=0 && pmtid < 24);

  for(int i=0;i < eNTimeSamples;i++){
    rawDataArray[i] += gRandom->Gaus(0,PMT_WHITE_NOISE_SIGMA[pmtid]);  
  }

  return;
}








void ArDM_Digitizer::Digitize_whiteNoise(int pmtid){

  /*
#if TEST_BRANCH0
#if TEST_BRANCH00

  assert(pmtid >=0 && pmtid < 24);

  Digitize_whiteNoise(pmtid,ArDM_Analysis::getInstance()->fRawData[pmtid]);

#endif //TEST_BRANCH0
#endif //TEST_BRANCH00
  */

  return;
}










void ArDM_Digitizer::Digitize_pmtRes(int* array,int pmtid,int nentries,float* rawDataArray){

  //array = int array[],
  //= number of photons detected in each time sample by PMT pmtid

  assert(pmtid >=0 && pmtid < 24);

  reset(rawDataArray,nentries);

  for(int i=0;i<nentries;i++){
    int nhits = array[i];
    //cout<<"pmt "<<pmtid<<"\t sample "<<i<<"\t nhist "<<nhits<<endl;
    for(int ii=0;ii<nhits;ii++){
      ArDM_PMTHit* hit;
      hit->setHitGlobalTime(i*TIMESAMPLE);
      hit->setHitPMTid(pmtid);
      hit->setIsDetected(1);
      Digitize_pmtRes(hit,rawDataArray);
    }  
  }


#if VERBOSE
  for(int i=0;i<nentries;i++){
    cout<<"pmtid "<<pmtid<<"\t sample "<< i << "\t rawData "<<rawDataArray[i]<<endl;
  }

#endif //VERBOSE


  return;
}









void ArDM_Digitizer::Digitize_pmtRes(int* array,int pmtid,int nentries){

  /*
#if TEST_BRANCH0
#if TEST_BRANCH00

  //array = int array[],
  //= number of photons detected in each time sample by PMT pmtid

  Digitize_pmtRes(array,pmtid,nentries,ArDM_Analysis::getInstance()->fRawData[pmtid]);

#endif //TEST_BRANCH0
#endif //TEST_BRANCH00

  */

  return;
}











void ArDM_Digitizer::Digitize_pmtRes_fromTree(string filepath){

  TFile* inputFile = new TFile(filepath.c_str());
  assert(inputFile->IsOpen());

  TTree* inputTree = (TTree*)inputFile->Get("tree");
  if(!inputTree) inputTree = (TTree*)inputFile->Get("Data");
  assert(inputTree);

  inputTree->SetBranchStatus("*RawData*",0);

  string outputFilename = filepath.substr(0,filepath.rfind(".root")) + "_digitized.root";

  TFile* outputFile = new TFile(outputFilename.c_str(),"recreate");
  TTree* outputTree = inputTree->CloneTree(0);


  float fRawData[NTOPPMT+NBOTTOMPMT][NTIMESAMPLEBOTTOM];

  ostringstream branchname,leafname;
  for(int i=0;i<NTOPPMT+NBOTTOMPMT;i++){
    branchname.str("");
    leafname.str("");

    //branchname = "eRawData" (instead of "fRawData") just to be consistent with the naming in realData.
    branchname<<"eRawData"<<i;
    leafname<<"eRawData["<< NTIMESAMPLEBOTTOM <<"]/F";
    outputTree->Branch(branchname.str().c_str(),&fRawData[i],leafname.str().c_str());
  }
  

  int fNphotonPMT[NTOPPMT+NBOTTOMPMT][NTIMESAMPLEBOTTOM];
  for(int i=0; i < NTOPPMT+NBOTTOMPMT; i++){
    branchname.str("");
    branchname << "fNphotonPMT" << i;
    inputTree->SetBranchAddress(branchname.str().c_str(),&fNphotonPMT[i]);
  }


  int nevents = inputTree->GetEntries();

  for(int evi=0;evi < nevents; evi++){
    if(!(evi%1000)) cout<<"event "<<evi<<endl;
    //if(evi > 0) break;
    inputTree->GetEntry(evi);


    for(int pmti=0;pmti < NTOPPMT+NBOTTOMPMT;pmti++){
      //reset(fRawData[pmti],NTIMESAMPLEBOTTOM);
      Digitize_pmtRes(fNphotonPMT[pmti],pmti,NTIMESAMPLEBOTTOM,fRawData[pmti]);

#if DARK_CURRENT_DIGITIZATION

      Digitize_darkCurrent(pmti,fRawData[pmti]);

#endif //DARK_CURRENT_DIGITIZATION


#if WHITE_NOSIE_DIGITIZATION

      Digitize_whiteNoise(pmti,fRawData[pmti]);

#endif //WHITE_NOISE_DIGITIZATION

    }


    outputTree->Fill();
  }
  
  outputFile->cd();
  outputTree->Write();
  outputFile->Close();

  delete inputFile;
  delete outputFile;

  return;
}
