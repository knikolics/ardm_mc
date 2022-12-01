//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  2 14:44:30 2013 by ROOT version 5.34/01
// from TTree tree/SS_238U_alphaN.root
// found on file: SS_238U_alphaN.root
//////////////////////////////////////////////////////////

#ifndef neutron_h
#define neutron_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class neutron {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           fNtopPMT;
   Int_t           fNbottomPMT;
   Int_t           fNphotonDetectedTop;
   Int_t           fNphotonDetectedBottom;
   UInt_t          fTotNphotonPMTTop[12];   //[fNtopPMT]
   UInt_t          fTotNphotonPMTBottom[12];   //[fNbottomPMT]
   UInt_t          NElasticScatterings;
   Double_t        fEkinInitial;
   vector<double>  *fEdep;
   Double_t        fErecoil;
   Float_t         ScatterPosX;
   Float_t         ScatterPosY;
   Float_t         ScatterPosZ;

   // List of branches
   TBranch        *b_fNtopPMT;   //!
   TBranch        *b_fNbottomPMT;   //!
   TBranch        *b_fNphotonDetectedTop;   //!
   TBranch        *b_fNphotonDetectedBottom;   //!
   TBranch        *b_fTotNphotonPMTTop;   //!
   TBranch        *b_fTotNphotonPMTBottom;   //!
   TBranch        *b_NElasticScatterings;   //!
   TBranch        *b_fEkinInitial;   //!
   TBranch        *b_fErecoil;   //!
   TBranch        *b_ScatterPosX;   //!
   TBranch        *b_ScatterPosY;   //!
   TBranch        *b_ScatterPosZ;   //!

   neutron(TTree *tree=0);
   virtual ~neutron();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef neutron_cxx
neutron::neutron(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SS_238U_alphaN.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SS_238U_alphaN.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

neutron::~neutron()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t neutron::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t neutron::LoadTree(Long64_t entry)
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

void neutron::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fNtopPMT", &fNtopPMT, &b_fNtopPMT);
   fChain->SetBranchAddress("fNbottomPMT", &fNbottomPMT, &b_fNbottomPMT);
   fChain->SetBranchAddress("fNphotonDetectedTop", &fNphotonDetectedTop, &b_fNphotonDetectedTop);
   fChain->SetBranchAddress("fNphotonDetectedBottom", &fNphotonDetectedBottom, &b_fNphotonDetectedBottom);
   fChain->SetBranchAddress("fTotNphotonPMTTop", fTotNphotonPMTTop, &b_fTotNphotonPMTTop);
   fChain->SetBranchAddress("fTotNphotonPMTBottom", fTotNphotonPMTBottom, &b_fTotNphotonPMTBottom);
   fChain->SetBranchAddress("NElasticScatterings", &NElasticScatterings, &b_NElasticScatterings);
   fChain->SetBranchAddress("fEkinInitial", &fEkinInitial, &b_fEkinInitial);
   fChain->SetBranchAddress("fErecoil", &fErecoil, &b_fErecoil);
   fChain->SetBranchAddress("ScatterPosX", &ScatterPosX, &b_ScatterPosX);
   fChain->SetBranchAddress("ScatterPosY", &ScatterPosY, &b_ScatterPosY);
   fChain->SetBranchAddress("ScatterPosZ", &ScatterPosZ, &b_ScatterPosZ);
   Notify();
}

Bool_t neutron::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void neutron::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t neutron::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef neutron_cxx
