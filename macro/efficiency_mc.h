//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  2 21:13:24 2015 by ROOT version 5.34/18
// from TTree validationT/validationT
// found on file: user.mtanaka.6405422._000001.hist-output.root
//////////////////////////////////////////////////////////

#ifndef efficiency_mc_h
#define efficiency_mc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class efficiency_mc {
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
   Float_t         probe_sa_delta_pt;
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
   Float_t         probe_roi_eta;
   Float_t         probe_roi_phi;
   Bool_t          probe_passL1mu4;
   Bool_t          probe_passSAmu4;
   Bool_t          probe_passSAmu4new;
   Bool_t          probe_passSAmu4newee;
   Bool_t          probe_passL1mu6;
   Bool_t          probe_passSAmu6;
   Bool_t          probe_passSAmu6new;
   Bool_t          probe_passSAmu6newee;

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
   TBranch        *b_probe_sa_delta_pt;   //!
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
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_probe_passL1mu4;   //!
   TBranch        *b_probe_passSAmu4;   //!
   TBranch        *b_probe_passSAmu4new;   //!
   TBranch        *b_probe_passSAmu4newee;   //!
   TBranch        *b_probe_passL1mu6;   //!
   TBranch        *b_probe_passSAmu6;   //!
   TBranch        *b_probe_passSAmu6new;   //!
   TBranch        *b_probe_passSAmu6newee;   //!

   efficiency_mc(TTree *tree=0);
   virtual ~efficiency_mc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef efficiency_mc_cxx
efficiency_mc::efficiency_mc(TTree *tree) : fChain(0) 
{
  TChain *chain = new TChain("validationT");
  //const char* rec = "datalist/inputEfficiencyMCZmumuDefault2.list";
  //const char* rec = "datalist/inputEfficiencyDataZmumu.list";
  //const char* rec = "datalist/inputEfficiencyDataZmumuMu6.list";
  //const char* rec = "datalist/inputEfficiencyDataZmumuDeltaPt.list";
  const char* rec = "datalist/inputEfficiencyZTandPVersion4.list";

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

efficiency_mc::~efficiency_mc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t efficiency_mc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t efficiency_mc::LoadTree(Long64_t entry)
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

void efficiency_mc::Init(TTree *tree)
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
   fChain->SetBranchAddress("probe_sa_delta_pt", &probe_sa_delta_pt, &b_probe_sa_delta_pt);
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
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("probe_passL1mu4", &probe_passL1mu4, &b_probe_passL1mu4);
   fChain->SetBranchAddress("probe_passSAmu4", &probe_passSAmu4, &b_probe_passSAmu4);
   fChain->SetBranchAddress("probe_passSAmu4new", &probe_passSAmu4new, &b_probe_passSAmu4new);
   fChain->SetBranchAddress("probe_passSAmu4newee", &probe_passSAmu4newee, &b_probe_passSAmu4newee);
   fChain->SetBranchAddress("probe_passL1mu6", &probe_passL1mu6, &b_probe_passL1mu6);
   fChain->SetBranchAddress("probe_passSAmu6", &probe_passSAmu6, &b_probe_passSAmu6);
   fChain->SetBranchAddress("probe_passSAmu6new", &probe_passSAmu6new, &b_probe_passSAmu6new);
   fChain->SetBranchAddress("probe_passSAmu6newee", &probe_passSAmu6newee, &b_probe_passSAmu6newee);
   Notify();
}

Bool_t efficiency_mc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void efficiency_mc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t efficiency_mc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef efficiency_mc_cxx
