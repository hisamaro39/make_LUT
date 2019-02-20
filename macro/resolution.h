//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  2 21:13:24 2015 by ROOT version 5.34/18
// from TTree validationT/validationT
// found on file: user.mtanaka.6405422._000001.hist-output.root
//////////////////////////////////////////////////////////

#ifndef resolution_h
#define resolution_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class resolution {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         tag_offline_pt;
   Float_t         tag_offline_eta;
   Float_t         tag_offline_phi;
   Float_t         tag_offline_charge;
   Float_t         probe_offline_pt;
   Float_t         probe_offline_eta;
   Float_t         probe_offline_phi;
   Float_t         probe_offline_charge;
   Float_t         probe_sa_pt;
   Float_t         probe_sa_pt_tgc;
   Float_t         probe_sa_pt_alpha;
   Float_t         probe_sa_pt_beta;
   Float_t         probe_sa_pt_ec_radius;
   Float_t         probe_sa_eta;
   Float_t         probe_sa_phi;
   Int_t           probe_sa_charge;
   Int_t           probe_sa_saddress;
   Float_t         probe_sa_alpha;
   Float_t         probe_sa_beta;
   Float_t         probe_sa_ec_radius;
   Float_t         probe_sa_br_radius;
   Int_t           probe_sa_eta_bin;
   Int_t           probe_sa_phi_bin;
   Float_t         probe_roi_eta;
   Float_t         probe_roi_phi;
   Bool_t          probe_passL1_1;
   Bool_t          probe_passL1_2;
   Bool_t          probe_passSA_1;
   Bool_t          probe_passSA_2;

   // List of branches
   TBranch        *b_tag_offline_pt;   //!
   TBranch        *b_tag_offline_eta;   //!
   TBranch        *b_tag_offline_phi;   //!
   TBranch        *b_tag_offline_charge;   //!
   TBranch        *b_probe_offline_pt;   //!
   TBranch        *b_probe_offline_eta;   //!
   TBranch        *b_probe_offline_phi;   //!
   TBranch        *b_probe_offline_charge;   //!
   TBranch        *b_probe_sa_pt;   //!
   TBranch        *b_probe_sa_pt_tgc;   //!
   TBranch        *b_probe_sa_pt_alpha;   //!
   TBranch        *b_probe_sa_pt_beta;   //!
   TBranch        *b_probe_sa_pt_ec_radius;   //!
   TBranch        *b_probe_sa_eta;   //!
   TBranch        *b_probe_sa_phi;   //!
   TBranch        *b_probe_sa_charge;   //!
   TBranch        *b_probe_sa_saddress;   //!
   TBranch        *b_probe_sa_alpha;   //!
   TBranch        *b_probe_sa_beta;   //!
   TBranch        *b_probe_sa_ec_radius;   //!
   TBranch        *b_probe_sa_br_radius;   //!
   TBranch        *b_probe_sa_eta_bin;   //!
   TBranch        *b_probe_sa_phi_bin;   //!
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_probe_passL1_1;   //!
   TBranch        *b_probe_passL1_2;   //!
   TBranch        *b_probe_passSA_1;   //!
   TBranch        *b_probe_passSA_2;   //!

   resolution(TTree *tree=0);
   virtual ~resolution();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef resolution_cxx
resolution::resolution(TTree *tree) : fChain(0) 
{
  TChain *chain = new TChain("validationT");
  //const char* rec = "datalist/inputEfficiencyOnline.list";
  const char* rec = "datalist/inputEfficiencyDefaultLUT.list";
  //const char* rec = "datalist/inputEfficiencyNewLUT.list";
  //const char* rec = "datalist/inputEfficiencyNewLUTFinal2.list";

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

resolution::~resolution()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t resolution::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t resolution::LoadTree(Long64_t entry)
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

void resolution::Init(TTree *tree)
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

   fChain->SetBranchAddress("tag_offline_pt", &tag_offline_pt, &b_tag_offline_pt);
   fChain->SetBranchAddress("tag_offline_eta", &tag_offline_eta, &b_tag_offline_eta);
   fChain->SetBranchAddress("tag_offline_phi", &tag_offline_phi, &b_tag_offline_phi);
   fChain->SetBranchAddress("tag_offline_charge", &tag_offline_charge, &b_tag_offline_charge);
   fChain->SetBranchAddress("probe_offline_pt", &probe_offline_pt, &b_probe_offline_pt);
   fChain->SetBranchAddress("probe_offline_eta", &probe_offline_eta, &b_probe_offline_eta);
   fChain->SetBranchAddress("probe_offline_phi", &probe_offline_phi, &b_probe_offline_phi);
   fChain->SetBranchAddress("probe_offline_charge", &probe_offline_charge, &b_probe_offline_charge);
   fChain->SetBranchAddress("probe_sa_pt", &probe_sa_pt, &b_probe_sa_pt);
   fChain->SetBranchAddress("probe_sa_pt_tgc", &probe_sa_pt_tgc, &b_probe_sa_pt_tgc);
   fChain->SetBranchAddress("probe_sa_pt_alpha", &probe_sa_pt_alpha, &b_probe_sa_pt_alpha);
   fChain->SetBranchAddress("probe_sa_pt_beta", &probe_sa_pt_beta, &b_probe_sa_pt_beta);
   fChain->SetBranchAddress("probe_sa_pt_ec_radius", &probe_sa_pt_ec_radius, &b_probe_sa_pt_ec_radius);
   fChain->SetBranchAddress("probe_sa_eta", &probe_sa_eta, &b_probe_sa_eta);
   fChain->SetBranchAddress("probe_sa_phi", &probe_sa_phi, &b_probe_sa_phi);
   fChain->SetBranchAddress("probe_sa_charge", &probe_sa_charge, &b_probe_sa_charge);
   fChain->SetBranchAddress("probe_sa_saddress", &probe_sa_saddress, &b_probe_sa_saddress);
   fChain->SetBranchAddress("probe_sa_alpha", &probe_sa_alpha, &b_probe_sa_alpha);
   fChain->SetBranchAddress("probe_sa_beta", &probe_sa_beta, &b_probe_sa_beta);
   fChain->SetBranchAddress("probe_sa_ec_radius", &probe_sa_ec_radius, &b_probe_sa_ec_radius);
   fChain->SetBranchAddress("probe_sa_br_radius", &probe_sa_br_radius, &b_probe_sa_br_radius);
   fChain->SetBranchAddress("probe_sa_eta_bin", &probe_sa_eta_bin, &b_probe_sa_eta_bin);
   fChain->SetBranchAddress("probe_sa_phi_bin", &probe_sa_phi_bin, &b_probe_sa_phi_bin);
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("probe_passL1_1", &probe_passL1_1, &b_probe_passL1_1);
   fChain->SetBranchAddress("probe_passL1_2", &probe_passL1_2, &b_probe_passL1_2);
   fChain->SetBranchAddress("probe_passSA_1", &probe_passSA_1, &b_probe_passSA_1);
   fChain->SetBranchAddress("probe_passSA_2", &probe_passSA_2, &b_probe_passSA_2);
   Notify();
}

Bool_t resolution::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void resolution::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t resolution::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef resolution_cxx
