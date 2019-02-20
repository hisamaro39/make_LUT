//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 27 10:06:00 2016 by ROOT version 5.34/18
// from TTree t_tap/TrigMuonTagAndProbe
// found on file: user.mtanaka.147407.PowhegPythia8_AZNLO_Zmumu.DKTAP.e3099_s2578_r7514_check_efficiency_zmumu_r7514_v01_EXT0.63107565/user.mtanaka.7536451.EXT0._000001.test.root
//////////////////////////////////////////////////////////

#ifndef CheckEff_h
#define CheckEff_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class CheckEff {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        probe_offline_pt;
   Double_t        probe_offline_eta;
   Double_t        probe_offline_phi;
   Int_t           probe_offline_charge;
   Double_t        probe_sa_pt;
   Double_t        probe_sa_eta;
   Double_t        probe_sa_phi;
   Int_t           probe_sa_charge;
   Double_t        probe_comb_pt;
   Double_t        probe_comb_eta;
   Double_t        probe_comb_phi;
   Int_t           probe_comb_charge;
   Double_t        probe_ef_pt;
   Double_t        probe_ef_eta;
   Double_t        probe_ef_phi;
   Int_t           probe_ef_charge;
   Double_t        probe_roi_eta;
   Double_t        probe_roi_phi;
   Bool_t          probe_pass_L1;
   Bool_t          probe_pass_SA;
   Bool_t          probe_pass_Comb;
   Bool_t          probe_pass_EF;
   Bool_t          probe_threshold;

   // List of branches
   TBranch        *b_probe_offline_pt;   //!
   TBranch        *b_probe_offline_eta;   //!
   TBranch        *b_probe_offline_phi;   //!
   TBranch        *b_probe_offline_charge;   //!
   TBranch        *b_probe_sa_pt;   //!
   TBranch        *b_probe_sa_eta;   //!
   TBranch        *b_probe_sa_phi;   //!
   TBranch        *b_probe_sa_charge;   //!
   TBranch        *b_probe_comb_pt;   //!
   TBranch        *b_probe_comb_eta;   //!
   TBranch        *b_probe_comb_phi;   //!
   TBranch        *b_probe_comb_charge;   //!
   TBranch        *b_probe_ef_pt;   //!
   TBranch        *b_probe_ef_eta;   //!
   TBranch        *b_probe_ef_phi;   //!
   TBranch        *b_probe_ef_charge;   //!
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_probe_pass_L1;   //!
   TBranch        *b_probe_pass_SA;   //!
   TBranch        *b_probe_pass_Comb;   //!
   TBranch        *b_probe_pass_EF;   //!
   TBranch        *b_probe_threshold;   //!

   CheckEff(TTree *tree=0);
   virtual ~CheckEff();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef CheckEff_cxx

CheckEff::CheckEff(TTree *tree) : fChain(0) 
{
  TChain *chain = new TChain("t_tap");
  //const char* rec = "datalist/inputCheckEfficiencyZmumur7447.list";
  //const char* rec = "datalist/inputCheckEfficiencyZmumur7463.list";
  //const char* rec = "datalist/inputCheckEfficiencyZmumuR7514.list";
  //const char* rec = "datalist/inputCheckEfficiencyZmumuR7534.list";
  //const char* rec = "datalist/inputEfficiencyNoTandPr7507.list";
  //const char* rec = "datalist/inputEfficiencyNoTandPr7540.list";
  
  //const char* rec = "datalist/inputEfficiencyZmumuUseMiddleEtaDefault.list";
  //const char* rec = "datalist/inputEfficiencyJpsimu4mu4UseMiddleEtaDefault.list";
    
  //const char* rec = "datalist/inputEfficiencyJpsimu2p5mu15.list";
  
  //const char* rec = "datalist/inputEfficiencyJpsiThondaDefault.list";
  const char* rec = "datalist/inputEfficiencyJpsiThondaModified.list";

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

CheckEff::~CheckEff()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CheckEff::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CheckEff::LoadTree(Long64_t entry)
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

void CheckEff::Init(TTree *tree)
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

   fChain->SetBranchAddress("probe_offline_pt", &probe_offline_pt, &b_probe_offline_pt);
   fChain->SetBranchAddress("probe_offline_eta", &probe_offline_eta, &b_probe_offline_eta);
   fChain->SetBranchAddress("probe_offline_phi", &probe_offline_phi, &b_probe_offline_phi);
   fChain->SetBranchAddress("probe_offline_charge", &probe_offline_charge, &b_probe_offline_charge);
   fChain->SetBranchAddress("probe_sa_pt", &probe_sa_pt, &b_probe_sa_pt);
   fChain->SetBranchAddress("probe_sa_eta", &probe_sa_eta, &b_probe_sa_eta);
   fChain->SetBranchAddress("probe_sa_phi", &probe_sa_phi, &b_probe_sa_phi);
   fChain->SetBranchAddress("probe_sa_charge", &probe_sa_charge, &b_probe_sa_charge);
   fChain->SetBranchAddress("probe_comb_pt", &probe_comb_pt, &b_probe_comb_pt);
   fChain->SetBranchAddress("probe_comb_eta", &probe_comb_eta, &b_probe_comb_eta);
   fChain->SetBranchAddress("probe_comb_phi", &probe_comb_phi, &b_probe_comb_phi);
   fChain->SetBranchAddress("probe_comb_charge", &probe_comb_charge, &b_probe_comb_charge);
   fChain->SetBranchAddress("probe_ef_pt", &probe_ef_pt, &b_probe_ef_pt);
   fChain->SetBranchAddress("probe_ef_eta", &probe_ef_eta, &b_probe_ef_eta);
   fChain->SetBranchAddress("probe_ef_phi", &probe_ef_phi, &b_probe_ef_phi);
   fChain->SetBranchAddress("probe_ef_charge", &probe_ef_charge, &b_probe_ef_charge);
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("probe_pass_L1", &probe_pass_L1, &b_probe_pass_L1);
   fChain->SetBranchAddress("probe_pass_SA", &probe_pass_SA, &b_probe_pass_SA);
   fChain->SetBranchAddress("probe_pass_Comb", &probe_pass_Comb, &b_probe_pass_Comb);
   fChain->SetBranchAddress("probe_pass_EF", &probe_pass_EF, &b_probe_pass_EF);
   fChain->SetBranchAddress("probe_threshold", &probe_threshold, &b_probe_threshold);
   Notify();
}

Bool_t CheckEff::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CheckEff::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CheckEff::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CheckEff_cxx
