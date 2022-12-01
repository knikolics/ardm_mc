//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 11 09:11:38 2014 by ROOT version 5.34/01
// from TTree tree//cluster/home04/phys/khnguyen/scratch_xl/neutron_background_borotron_shield/170/ArDM_Analysis_tree_neutronBkg_SS_238U_ant_neutronShield_10cm.root
// found on file: neutronBkg/SS_238U_ant_neutronShield_10cm.root
//////////////////////////////////////////////////////////

#ifndef singleScatterAna_h
#define singleScatterAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class singleScatterAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           fNIntInLAr;
   Int_t           fNIntInFV;
   Int_t           fNElasIntInLAr;
   Int_t           fNElasIntInFV;
   Int_t           NElasticScatterings;
   Double_t        fEkinInitial;
   Double_t        fEkinEnterLAr;
   Double_t        fPosx;
   Double_t        fPosy;
   Double_t        fPosz;
   Double_t        fKillPosx;
   Double_t        fKillPosy;
   Double_t        fKillPosz;
   vector<double>  *fEdep;
   Double_t        fErecoil;
   vector<double>  *fIntPosx;
   vector<double>  *fIntPosy;
   vector<double>  *fIntPosz;

   // List of branches
   TBranch        *b_fNIntInLAr;   //!
   TBranch        *b_fNIntInFV;   //!
   TBranch        *b_fNElasIntInLAr;   //!
   TBranch        *b_fNElasIntInFV;   //!
   TBranch        *b_NElasticScatterings;   //!
   TBranch        *b_fEkinInitial;   //!
   TBranch        *b_fEkinEnterLAr;   //!
   TBranch        *b_fPosx;   //!
   TBranch        *b_fPosy;   //!
   TBranch        *b_fPosz;   //!
   TBranch        *b_fKillPosx;   //!
   TBranch        *b_fKillPosy;   //!
   TBranch        *b_fKillPosz;   //!
   TBranch        *b_fEdep;   //!
   TBranch        *b_fErecoil;   //!
   TBranch        *b_fIntPosx;   //!
   TBranch        *b_fIntPosy;   //!
   TBranch        *b_fIntPosz;   //!

   singleScatterAna(TTree *tree=0);
   virtual ~singleScatterAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef singleScatterAna_cxx
singleScatterAna::singleScatterAna(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("neutronBkg/SS_238U_ant_neutronShield_10cm.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("neutronBkg/SS_238U_ant_neutronShield_10cm.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

singleScatterAna::~singleScatterAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t singleScatterAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t singleScatterAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void singleScatterAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   fEdep = 0;
   fIntPosx = 0;
   fIntPosy = 0;
   fIntPosz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fNIntInLAr", &fNIntInLAr, &b_fNIntInLAr);
   fChain->SetBranchAddress("fNIntInFV", &fNIntInFV, &b_fNIntInFV);
   fChain->SetBranchAddress("fNElasIntInLAr", &fNElasIntInLAr, &b_fNElasIntInLAr);
   fChain->SetBranchAddress("fNElasIntInFV", &fNElasIntInFV, &b_fNElasIntInFV);
   fChain->SetBranchAddress("NElasticScatterings", &NElasticScatterings, &b_NElasticScatterings);
   fChain->SetBranchAddress("fEkinInitial", &fEkinInitial, &b_fEkinInitial);
   fChain->SetBranchAddress("fEkinEnterLAr", &fEkinEnterLAr, &b_fEkinEnterLAr);
   fChain->SetBranchAddress("fPosx", &fPosx, &b_fPosx);
   fChain->SetBranchAddress("fPosy", &fPosy, &b_fPosy);
   fChain->SetBranchAddress("fPosz", &fPosz, &b_fPosz);
   fChain->SetBranchAddress("fKillPosx", &fKillPosx, &b_fKillPosx);
   fChain->SetBranchAddress("fKillPosy", &fKillPosy, &b_fKillPosy);
   fChain->SetBranchAddress("fKillPosz", &fKillPosz, &b_fKillPosz);
   fChain->SetBranchAddress("fEdep", &fEdep, &b_fEdep);
   fChain->SetBranchAddress("fErecoil", &fErecoil, &b_fErecoil);
   fChain->SetBranchAddress("fIntPosx", &fIntPosx, &b_fIntPosx);
   fChain->SetBranchAddress("fIntPosy", &fIntPosy, &b_fIntPosy);
   fChain->SetBranchAddress("fIntPosz", &fIntPosz, &b_fIntPosz);
   Notify();
}

Bool_t singleScatterAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void singleScatterAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t singleScatterAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef singleScatterAna_cxx
