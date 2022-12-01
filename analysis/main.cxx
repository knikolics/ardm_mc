#include <TROOT.h>
#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <Riostream.h>
#include "TSystem.h"
//#include "neutron.h"
#include "singleScatterAna.h"
#include <string>

using namespace std;
//-----------------------------
int main(int argc,char **argv) {
  if (argc != 2) {
    cout << "Give a ROOT file as input !" << endl;
    return 0;
  }
  //string name = (string)argv[1];
  string output = "histos_" + (string)argv[1];
  cout << "Create output file " << output << endl;
  TFile* myfile = new TFile(output.c_str(),"RECREATE");
  
  TChain* chain = new TChain("tree");
  //TString infile;
  string infile = (string)argv[1];
  cout << "Adding file " << infile << " to TChain" << endl;
  chain->Add(infile.c_str());

  cout << "Entries in analysis chain: " << chain->GetEntries() << endl;
  //neutron* ardm = new neutron(chain);
  singleScatterAna* ardm = new singleScatterAna(chain);
  ardm->Loop();

  myfile->Write();
  myfile->Close();
 
  return 0;
}

