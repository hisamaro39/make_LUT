//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 21 01:13:53 2016 by ROOT version 5.34/18
// from TTree t_tap/TrigMuonTagAndProbe
// found on file: user.mtanaka.300004.Pythia8BPhotospp_A14_CTEQ6L1_pp_Jpsimu4mu20.DKTAP.e4397_s2608_r6869_2muon_jpsi_mu4_mu20_again4_EXT0.62147120/user.mtanaka.7485178.EXT0._000001.test.root
//////////////////////////////////////////////////////////

#ifndef check2muon_h
#define check2muon_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class check2muon {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        probe_offline_pt;
   Double_t        probe_offline_eta;
   Double_t        probe_offline_phi;
   Int_t           probe_offline_charge;
   Double_t        tag_offline_pt;
   Double_t        tag_offline_eta;
   Double_t        tag_offline_phi;
   Int_t           tag_offline_charge;
   Double_t        probe_roi_eta;
   Double_t        probe_roi_phi;
   Double_t        tag_roi_eta;
   Double_t        tag_roi_phi;
   Double_t        probe_sa_pt;
   Double_t        probe_sa_eta;
   Double_t        probe_sa_phi;
   Int_t           probe_sa_charge;
   Double_t        tag_sa_pt;
   Double_t        tag_sa_eta;
   Double_t        tag_sa_phi;
   Int_t           tag_sa_charge;
   Double_t        probe_comb_pt;
   Double_t        probe_comb_eta;
   Double_t        probe_comb_phi;
   Int_t           probe_comb_charge;
   Double_t        tag_comb_pt;
   Double_t        tag_comb_eta;
   Double_t        tag_comb_phi;
   Int_t           tag_comb_charge;
   Double_t        probe_ef_pt;
   Double_t        probe_ef_eta;
   Double_t        probe_ef_phi;
   Int_t           probe_ef_charge;
   Double_t        tag_ef_pt;
   Double_t        tag_ef_eta;
   Double_t        tag_ef_phi;
   Int_t           tag_ef_charge;
   Bool_t          tag_passSA;
   Bool_t          tag_passComb;
   Bool_t          tag_passEF;
   Bool_t          probe_passL1;
   Bool_t          probe_passSA;
   Bool_t          probe_passComb;
   Bool_t          probe_passEF;

   // List of branches
   TBranch        *b_probe_offline_pt;   //!
   TBranch        *b_probe_offline_eta;   //!
   TBranch        *b_probe_offline_phi;   //!
   TBranch        *b_probe_offline_charge;   //!
   TBranch        *b_tag_offline_pt;   //!
   TBranch        *b_tag_offline_eta;   //!
   TBranch        *b_tag_offline_phi;   //!
   TBranch        *b_tag_offline_charge;   //!
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_tag_roi_eta;   //!
   TBranch        *b_tag_roi_phi;   //!
   TBranch        *b_probe_sa_pt;   //!
   TBranch        *b_probe_sa_eta;   //!
   TBranch        *b_probe_sa_phi;   //!
   TBranch        *b_probe_sa_charge;   //!
   TBranch        *b_tag_sa_pt;   //!
   TBranch        *b_tag_sa_eta;   //!
   TBranch        *b_tag_sa_phi;   //!
   TBranch        *b_tag_sa_charge;   //!
   TBranch        *b_probe_comb_pt;   //!
   TBranch        *b_probe_comb_eta;   //!
   TBranch        *b_probe_comb_phi;   //!
   TBranch        *b_probe_comb_charge;   //!
   TBranch        *b_tag_comb_pt;   //!
   TBranch        *b_tag_comb_eta;   //!
   TBranch        *b_tag_comb_phi;   //!
   TBranch        *b_tag_comb_charge;   //!
   TBranch        *b_probe_ef_pt;   //!
   TBranch        *b_probe_ef_eta;   //!
   TBranch        *b_probe_ef_phi;   //!
   TBranch        *b_probe_ef_charge;   //!
   TBranch        *b_tag_ef_pt;   //!
   TBranch        *b_tag_ef_eta;   //!
   TBranch        *b_tag_ef_phi;   //!
   TBranch        *b_tag_ef_charge;   //!
   TBranch        *b_tag_passSA;   //!
   TBranch        *b_tag_passComb;   //!
   TBranch        *b_tag_passEF;   //!
   TBranch        *b_probe_passL1;   //!
   TBranch        *b_probe_passSA;   //!
   TBranch        *b_probe_passComb;   //!
   TBranch        *b_probe_passEF;   //!

   check2muon(TTree *tree=0);
   virtual ~check2muon();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef check2muon_cxx
check2muon::check2muon(TTree *tree) : fChain(0) 
{
  TChain *chain = new TChain("t_tap");
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20.list";
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLTmu4.list";
  const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLTmu4msonly.list";

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

check2muon::~check2muon()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t check2muon::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t check2muon::LoadTree(Long64_t entry)
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

void check2muon::Init(TTree *tree)
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
   fChain->SetBranchAddress("tag_offline_pt", &tag_offline_pt, &b_tag_offline_pt);
   fChain->SetBranchAddress("tag_offline_eta", &tag_offline_eta, &b_tag_offline_eta);
   fChain->SetBranchAddress("tag_offline_phi", &tag_offline_phi, &b_tag_offline_phi);
   fChain->SetBranchAddress("tag_offline_charge", &tag_offline_charge, &b_tag_offline_charge);
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("tag_roi_eta", &tag_roi_eta, &b_tag_roi_eta);
   fChain->SetBranchAddress("tag_roi_phi", &tag_roi_phi, &b_tag_roi_phi);
   fChain->SetBranchAddress("probe_sa_pt", &probe_sa_pt, &b_probe_sa_pt);
   fChain->SetBranchAddress("probe_sa_eta", &probe_sa_eta, &b_probe_sa_eta);
   fChain->SetBranchAddress("probe_sa_phi", &probe_sa_phi, &b_probe_sa_phi);
   fChain->SetBranchAddress("probe_sa_charge", &probe_sa_charge, &b_probe_sa_charge);
   fChain->SetBranchAddress("tag_sa_pt", &tag_sa_pt, &b_tag_sa_pt);
   fChain->SetBranchAddress("tag_sa_eta", &tag_sa_eta, &b_tag_sa_eta);
   fChain->SetBranchAddress("tag_sa_phi", &tag_sa_phi, &b_tag_sa_phi);
   fChain->SetBranchAddress("tag_sa_charge", &tag_sa_charge, &b_tag_sa_charge);
   fChain->SetBranchAddress("probe_comb_pt", &probe_comb_pt, &b_probe_comb_pt);
   fChain->SetBranchAddress("probe_comb_eta", &probe_comb_eta, &b_probe_comb_eta);
   fChain->SetBranchAddress("probe_comb_phi", &probe_comb_phi, &b_probe_comb_phi);
   fChain->SetBranchAddress("probe_comb_charge", &probe_comb_charge, &b_probe_comb_charge);
   fChain->SetBranchAddress("tag_comb_pt", &tag_comb_pt, &b_tag_comb_pt);
   fChain->SetBranchAddress("tag_comb_eta", &tag_comb_eta, &b_tag_comb_eta);
   fChain->SetBranchAddress("tag_comb_phi", &tag_comb_phi, &b_tag_comb_phi);
   fChain->SetBranchAddress("tag_comb_charge", &tag_comb_charge, &b_tag_comb_charge);
   fChain->SetBranchAddress("probe_ef_pt", &probe_ef_pt, &b_probe_ef_pt);
   fChain->SetBranchAddress("probe_ef_eta", &probe_ef_eta, &b_probe_ef_eta);
   fChain->SetBranchAddress("probe_ef_phi", &probe_ef_phi, &b_probe_ef_phi);
   fChain->SetBranchAddress("probe_ef_charge", &probe_ef_charge, &b_probe_ef_charge);
   fChain->SetBranchAddress("tag_ef_pt", &tag_ef_pt, &b_tag_ef_pt);
   fChain->SetBranchAddress("tag_ef_eta", &tag_ef_eta, &b_tag_ef_eta);
   fChain->SetBranchAddress("tag_ef_phi", &tag_ef_phi, &b_tag_ef_phi);
   fChain->SetBranchAddress("tag_ef_charge", &tag_ef_charge, &b_tag_ef_charge);
   fChain->SetBranchAddress("tag_passSA", &tag_passSA, &b_tag_passSA);
   fChain->SetBranchAddress("tag_passComb", &tag_passComb, &b_tag_passComb);
   fChain->SetBranchAddress("tag_passEF", &tag_passEF, &b_tag_passEF);
   fChain->SetBranchAddress("probe_passL1", &probe_passL1, &b_probe_passL1);
   fChain->SetBranchAddress("probe_passSA", &probe_passSA, &b_probe_passSA);
   fChain->SetBranchAddress("probe_passComb", &probe_passComb, &b_probe_passComb);
   fChain->SetBranchAddress("probe_passEF", &probe_passEF, &b_probe_passEF);
   Notify();
}

Bool_t check2muon::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void check2muon::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t check2muon::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef check2muon_cxx
