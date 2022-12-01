#define neutron_cxx
#include "neutron.h"
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TMath.h>
#include <iostream>

using namespace std;

void neutron::Loop()
{
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
   const int TANK_INNER_RADIUS = 500; //mm
   const float TANK_HALF_HEIGHT = 1044.7; //mm
   const int nbinsR = 500;
   const int nbinsE = 400;

   TH3F* singleScatter = new TH3F("singleScatter","",nbinsR,0,TANK_INNER_RADIUS,nbinsR,0,TANK_INNER_RADIUS,1045,-1045,1045);
   TH2D* phaseSpace = new TH2D("phaseSpace","",nbinsE,0,400,nbinsR,0,nbinsR);
   TH1D* singleErecoil = new TH1D("singleErecoil","",nbinsE,0,400);
   TH1D* initialE = new TH1D("initialE","",100,0,10);
   TH1D* radialFid = new TH1D("radialFid","",nbinsR,0,TANK_INNER_RADIUS);
   TH1D* zFid = new TH1D("zFid","",1045,-1045,1045);
   TH2D* fidVol = new TH2D("fidVol","",nbinsR,0,2500,1045,-1045,1045);
   TH2D* fidVol_30keV = new TH2D("fidVol_30keV","",500,0,2500,1045,-1045,1045);
   TH1I* NumScatter = new TH1I("NumScatter","",20,0,20);
/* end neutron bkg */

   const int nbinsLog = 200;

   TH2D* signalEnergy_noCuts = new TH2D("signalEnergy_noCuts","",nbinsE/2,0,200,nbinsLog,-1,4);
   TH2D* signalEnergy_Cuts = new TH2D("signalEnergy_Cuts","",nbinsE/2,0,200,nbinsLog,-1,4);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if ( NElasticScatterings>0 && fErecoil>2)    NumScatter->Fill(NElasticScatterings);     

//      if ( fErecoil < 30 ) continue;  // set energy threshold for acceptance
      /* 0 or more than 1 elastic scattering? */
      if ( NElasticScatterings==1 ) {
         double radius = sqrt(ScatterPosX*ScatterPosX + ScatterPosY*ScatterPosY);
         singleScatter->Fill(ScatterPosX,ScatterPosY,ScatterPosZ);
         singleErecoil->Fill(fErecoil);
         initialE->Fill(fEkinInitial);
      
         //double radius = sqrt(ScatterPosX*ScatterPosX + ScatterPosY*ScatterPosY);
         phaseSpace->Fill(fErecoil,radius);
         radialFid->Fill(radius);
         zFid->Fill(ScatterPosZ);
         double radSquare = radius*radius/100; // in cm^2
         fidVol->Fill(radius,ScatterPosZ);
        // if ( fErecoil > 30 )  fidVol_30keV->Fill(radius,ScatterPosZ);
      }
      if ( fNphotonDetectedTop==0 || fNphotonDetectedBottom==0 ) continue;
//      if ( fNphotonDetectedBottom==0 ) continue;
      double logSignal = log10 (fNphotonDetectedTop/fNphotonDetectedBottom);
      signalEnergy_noCuts->Fill(fErecoil,logSignal);

      //signalEnergy_Cuts->Fill(fErecoil,logSignal);
      // if (Cut(ientry) < 0) continue;
   }
   
}


