//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 17 06:54:34 2015 by ROOT version 5.34/18
// from TTree validationT/validationT
// found on file: inputEfficiency/data/defaultLUT/no_dr/267638/user.mtanaka.test.20150916151842.00267638.r6818_p2358_p2361_hist-output.root.43107335/user.mtanaka.6491787._000001.hist-output.root
//////////////////////////////////////////////////////////

#ifndef validEfficiency_h
#define validEfficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class validEfficiency {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         offline_pt;
   Float_t         offline_eta;
   Float_t         offline_phi;
   Float_t         offline_charge;
   Float_t         offline_ext_eta;
   Float_t         offline_ext_phi;
   Float_t         roi_eta;
   Float_t         roi_phi;

   // List of branches
   TBranch        *b_offline_pt;   //!
   TBranch        *b_offline_eta;   //!
   TBranch        *b_offline_phi;   //!
   TBranch        *b_offline_charge;   //!
   TBranch        *b_offline_ext_eta;   //!
   TBranch        *b_offline_ext_phi;   //!
   TBranch        *b_roi_eta;   //!
   TBranch        *b_roi_phi;   //!

   validEfficiency(TTree *tree=0);
   virtual ~validEfficiency();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef validEfficiency_cxx
validEfficiency::validEfficiency(TTree *tree) : fChain(0) 
{
  TChain *chain = new TChain("validationT");
  const char* rec = "datalist/inputEfficiencyDefaultLUTNoDr.list";

  ifstream finlist(rec);
  string file_rec;
  cout << "input sample" << endl;
  while(finlist>>file_rec) {
    cout << file_rec.c_str() << endl;
    chain->Add(file_rec.c_str());
  }
  tree = static_cast<TTree*>(chain);

  Init(tree);
}

validEfficiency::~validEfficiency()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t validEfficiency::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t validEfficiency::LoadTree(Long64_t entry)
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

void validEfficiency::Init(TTree *tree)
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

   fChain->SetBranchAddress("offline_pt", &offline_pt, &b_offline_pt);
   fChain->SetBranchAddress("offline_eta", &offline_eta, &b_offline_eta);
   fChain->SetBranchAddress("offline_phi", &offline_phi, &b_offline_phi);
   fChain->SetBranchAddress("offline_charge", &offline_charge, &b_offline_charge);
   fChain->SetBranchAddress("offline_ext_eta", &offline_ext_eta, &b_offline_ext_eta);
   fChain->SetBranchAddress("offline_ext_phi", &offline_ext_phi, &b_offline_ext_phi);
   fChain->SetBranchAddress("roi_eta", &roi_eta, &b_roi_eta);
   fChain->SetBranchAddress("roi_phi", &roi_phi, &b_roi_phi);
   Notify();
}

Bool_t validEfficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void validEfficiency::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t validEfficiency::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef validEfficiency_cxx
