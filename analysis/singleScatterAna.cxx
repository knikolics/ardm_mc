#define singleScatterAna_cxx
#include "singleScatterAna.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream> 
#include <vector>
#include <math.h>

using namespace std;

void singleScatterAna::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L singleScatterAna.C
//      Root > singleScatterAna t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

/* neutron bkg histos */
   const int TANK_INNER_RADIUS = 480; //mm
   const int TankRadiusSquare = 2340; //cm^2
   const float TANK_HALF_HEIGHT = 1044.7; //mm
   const int nbinsR = 120;
   const int nbinsE = 250;
   const int nbinsZ = 570;

   TH3F* singleScatter = new TH3F("singleScatter","",nbinsR,0,TANK_INNER_RADIUS,nbinsR,0,TANK_INNER_RADIUS,nbinsZ,-1055,655);
   TH2D* phaseSpace = new TH2D("phaseSpace","",nbinsE,0,250,nbinsR,0,nbinsR);
   TH1D* singleErecoil = new TH1D("singleErecoil","",nbinsE,0,250);
   TH1D* singleErecoil_Cut = new TH1D("singleErecoil_Cut","",nbinsE,0,250);
   TH1D* initialE = new TH1D("initialE","",100,0,10);
   TH1D* scatterRadius = new TH1D("scatterRadius","",nbinsR,0,TANK_INNER_RADIUS);
   TH2D* radialSqFid = new TH2D("radialSqFid","",nbinsE,0,250,nbinsR,0,TankRadiusSquare);
   TH2D* radialFid = new TH2D("radialFid","",nbinsE,0,250,nbinsR,0,TANK_INNER_RADIUS);
   TH2D* zFid = new TH2D("zFid","",nbinsE,0,250,nbinsZ,-1055,655);
   TH2D* fidVol = new TH2D("fidVol","",nbinsR,0,2500,nbinsZ,-1055,655);
   TH1I* NumScatter_moreThan0 = new TH1I("NumScatter_moreThan0","",25,0,25);
   TH1I* NumScatter_moreThan1 = new TH1I("NumScatter_moreThan1","",25,0,25);
/* end neutron bkg */

   vector<double> *myEdep;
//   myEdep->clear();
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      int totEdep = 0;

      if ( fErecoil < 10 ) continue;
      if ( NElasticScatterings>0 )    NumScatter_moreThan0->Fill(NElasticScatterings);
      if ( NElasticScatterings>1 )    NumScatter_moreThan1->Fill(NElasticScatterings);
      if ( NElasticScatterings==1 ) { 
          initialE->Fill(fEkinInitial);
          singleErecoil->Fill(fErecoil);
          double xPos = fIntPosx->at(0);
          double yPos = fIntPosy->at(0);
          double radius = sqrt(xPos*xPos+yPos*yPos);
          double radiusSq = xPos*xPos+yPos*yPos;
          radialFid->Fill(fErecoil,radius);
          radialSqFid->Fill(fErecoil,radiusSq);
          double zPos = fIntPosz->at(0);
          zFid->Fill(fErecoil,zPos);
      // if (Cut(ientry) < 0) continue;
//      myEdep->clear();
      myEdep = fEdep;
      double dummy;
      unsigned int nScatterings = myEdep->size();
      if ( nScatterings > 3) continue;
      for (unsigned int iScatter = 0; iScatter<nScatterings; iScatter++) {
          dummy = 0;
          dummy = myEdep->at(iScatter);
          if ( dummy < 10 ) { 
             totEdep = 0;
             break;
          }
          else totEdep += myEdep->at(iScatter);
      }
//      cout << "totEdep " << totEdep << endl;
      if ( totEdep<30 || totEdep>100 ) continue;
      singleErecoil_Cut->Fill(fErecoil);
    }
  }
}

