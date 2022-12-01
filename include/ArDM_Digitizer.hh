#ifndef _ARDM_DIGITIZER_
#define _ARDM_DIGITIZER_

#include "G4VDigi.hh"
#include "G4VDigitizerModule.hh"
#include "ArDM_AddOns.hh"
#include "ArDM_PMTHit.hh"
#include "G4HCofThisEvent.hh"

using namespace std;

class ArDM_Digitizer : public G4VDigitizerModule {

  static int eNTimeSamples;

  void Digitize_simple(ArDM_PMTHit* hit); //simple digitization, just record global time of the hit
  void Digitize_simple(ArDM_PMTHitsCollection* hitcol); //simple digitization, just record global time of the hit

  void Digitize_pmtRes(ArDM_PMTHit* hit); //take pmtResponse into account. simulate analog signal for a hit.
  void Digitize_pmtRes_try1(ArDM_PMTHit* hit,float* array); //take pmtResponse into account. simulate analog signal for a hit.
  void Digitize_pmtRes_try2(ArDM_PMTHit* hit,float* array); //take pmtResponse into account. simulate analog signal for a hit.
  void Digitize_pmtRes(ArDM_PMTHit* hit,float* array); //take pmtResponse into account. simulate analog signal for a hit.

  void Digitize_pmtRes(ArDM_PMTHitsCollection* hitscol); //take pmtResponse into account. simulate analog signal for all hits.

  void Digitize_darkCurrent(int pmtid);
  void Digitize_darkCurrent(int pmtid,float* rawDataArray);

  void Digitize_whiteNoise(int pmtid);
  void Digitize_whiteNoise(int pmtid,float* rawDataArray);

  void Digitize_pmtRes(int* array, int pmtid, int nentries = NTIMESAMPLE); //arry = int array[] = number of photons detected in each timesample by 1 PMT

  void Digitize_pmtRes(int* array, int pmtid, int nentries, float* rawDataArray); //arry = int array[] = number of photons detected in each timesample by 1 PMT

public :

  ArDM_Digitizer(G4String name="PMTDigitizer");
  ~ArDM_Digitizer();

  void reset(float* array, int nentries);
  //void reset();

  void Digitize_simple(G4String hitsColName);
  void Digitize_pmtRes(G4String hitsColName);

  void Digitize_simple(G4HCofThisEvent* hc);
  void Digitize_pmtRes(G4HCofThisEvent* hc);
  
  void Digitize_pmtRes_fromTree(string filepath);
  void Digitize(G4HCofThisEvent* hc);

  void Digitize();

  static void setNTimeSamples(int ntimeSamples){ eNTimeSamples = ntimeSamples;};
};





#endif // _ARDM_DIGITIZER_
