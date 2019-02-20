//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 27 16:33:24 2015 by ROOT version 5.34/18
// from TTree t_tap/TrigMuonTagAndProbe
// found on file: inputEfficiency/data/defaultLUT/testJpsi/user.mtanaka.00279284.physics_Main.merge.DKTAP.f628_m1497_p2419_jpsi_mu4_test_EXT0.49528419/user.mtanaka.6806545.EXT0._000011.test.root
//////////////////////////////////////////////////////////

#ifndef effJpsi_h
#define effJpsi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class effJpsi {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        probe_offline_pt;
   Double_t        probe_offline_eta;
   Double_t        probe_offline_phi;
   Int_t           probe_offline_charge;
   //Double_t        tag_offline_pt;
   //Double_t        tag_offline_eta;
   //Double_t        tag_offline_phi;
   //Int_t           tag_offline_charge;
   Double_t        probe_sa_pt;
   Double_t        probe_sa_pt_alpha;
   Double_t        probe_sa_pt_beta;
   Double_t        probe_sa_pt_tgc;
   Double_t        probe_sa_pt_ee;
   Double_t        probe_sa_pt_ec_radius;
   Double_t        probe_sa_pt_new;
   Double_t        probe_sa_pt_alpha_new;
   Double_t        probe_sa_pt_beta_new;
   Double_t        probe_sa_pt_tgc_new;
   Double_t        probe_sa_alpha;
   Double_t        probe_sa_beta;
   Double_t        probe_sa_ec_radius;
   Double_t        probe_sa_br_radius;
   Double_t        probe_sa_eta;
   Double_t        probe_sa_phi;
   Double_t        probe_sa_deltaPt;
   Int_t           probe_sa_charge;
   Int_t           probe_sa_saddress;
   Int_t           probe_sa_ec_eta_bin;
   Int_t           probe_sa_ec_phi_bin;
   Double_t        probe_roi_eta;
   Double_t        probe_roi_phi;
   Bool_t          probe_pass_L1mu4;
   Bool_t          probe_pass_SAmu4;
   Bool_t          probe_pass_SAmu4_new;
   Bool_t          probe_pass_SAmu4_ee;
   Bool_t          probe_pass_L1mu6;
   Bool_t          probe_pass_SAmu6;
   Bool_t          probe_pass_SAmu6_new;
   Bool_t          probe_pass_SAmu6_ee;
   //vector<Float_t> *sp_r;
   //vector<Float_t> *sp_z;
   //vector<Float_t> *segment_r;
   //vector<Float_t> *segment_z;
   //vector<Float_t> *segment_chamber;

   // List of branches
   TBranch        *b_probe_offline_pt;   //!
   TBranch        *b_probe_offline_eta;   //!
   TBranch        *b_probe_offline_phi;   //!
   TBranch        *b_probe_offline_charge;   //!
   //TBranch        *b_tag_offline_pt;   //!
   //TBranch        *b_tag_offline_eta;   //!
   //TBranch        *b_tag_offline_phi;   //!
   //TBranch        *b_tag_offline_charge;   //!
   TBranch        *b_probe_sa_pt;   //!
   TBranch        *b_probe_sa_pt_alpha;   //!
   TBranch        *b_probe_sa_pt_beta;   //!
   TBranch        *b_probe_sa_pt_tgc;   //!
   TBranch        *b_probe_sa_pt_ee;   //!
   TBranch        *b_probe_sa_pt_ec_radius;   //!
   TBranch        *b_probe_sa_pt_new;   //!
   TBranch        *b_probe_sa_pt_alpha_new;   //!
   TBranch        *b_probe_sa_pt_beta_new;   //!
   TBranch        *b_probe_sa_pt_tgc_new;   //!
   TBranch        *b_probe_sa_alpha;   //!
   TBranch        *b_probe_sa_beta;   //!
   TBranch        *b_probe_sa_ec_radius;   //!
   TBranch        *b_probe_sa_br_radius;   //!
   TBranch        *b_probe_sa_eta;   //!
   TBranch        *b_probe_sa_phi;   //!
   TBranch        *b_probe_sa_deltaPt;   //!
   TBranch        *b_probe_sa_charge;   //!
   TBranch        *b_probe_sa_saddress;   //!
   TBranch        *b_probe_sa_ec_eta_bin;   //!
   TBranch        *b_probe_sa_ec_phi_bin;   //!
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_probe_pass_L1mu4;   //!
   TBranch        *b_probe_pass_SAmu4;   //!
   TBranch        *b_probe_pass_SAmu4_new;   //!
   TBranch        *b_probe_pass_SAmu4_ee;   //!
   TBranch        *b_probe_pass_L1mu6;   //!
   TBranch        *b_probe_pass_SAmu6;   //!
   TBranch        *b_probe_pass_SAmu6_new;   //!
   TBranch        *b_probe_pass_SAmu6_ee;   //!
   //TBranch        *b_sp_r;   //!
   //TBranch        *b_sp_z;   //!
   //TBranch        *b_segment_r;   //!
   //TBranch        *b_segment_z;   //!
   //TBranch        *b_segment_chamber;   //!

   effJpsi(TTree *tree=0);
   virtual ~effJpsi();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef effJpsi_cxx

effJpsi::effJpsi(TTree *tree) : fChain(0) 
{
  TChain *chain = new TChain("t_tap");
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9Selected2.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9All.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9Again.list";
  const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9Separate.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9EESym.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9SP.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9Seg.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9EEDistance.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9EEDistanceAll.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9EEAll.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9Final.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9Selected.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9NoWindow.list";
  //const char* rec = "datalist/inputEfficiencyJpsiTandPVersion9OfflineCharge.list";

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

effJpsi::~effJpsi()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t effJpsi::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t effJpsi::LoadTree(Long64_t entry)
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

void effJpsi::Init(TTree *tree)
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

   //sp_r=0;
   //sp_z=0;
   //segment_r=0;
   //segment_z=0;
   //segment_chamber=0;
   
   fChain->SetBranchAddress("probe_offline_pt", &probe_offline_pt, &b_probe_offline_pt);
   fChain->SetBranchAddress("probe_offline_eta", &probe_offline_eta, &b_probe_offline_eta);
   fChain->SetBranchAddress("probe_offline_phi", &probe_offline_phi, &b_probe_offline_phi);
   fChain->SetBranchAddress("probe_offline_charge", &probe_offline_charge, &b_probe_offline_charge);
   //fChain->SetBranchAddress("tag_offline_pt", &tag_offline_pt, &b_tag_offline_pt);
   //fChain->SetBranchAddress("tag_offline_eta", &tag_offline_eta, &b_tag_offline_eta);
   //fChain->SetBranchAddress("tag_offline_phi", &tag_offline_phi, &b_tag_offline_phi);
   //fChain->SetBranchAddress("tag_offline_charge", &tag_offline_charge, &b_tag_offline_charge);
   fChain->SetBranchAddress("probe_sa_pt", &probe_sa_pt, &b_probe_sa_pt);
   fChain->SetBranchAddress("probe_sa_pt_alpha", &probe_sa_pt_alpha, &b_probe_sa_pt_alpha);
   fChain->SetBranchAddress("probe_sa_pt_beta", &probe_sa_pt_beta, &b_probe_sa_pt_beta);
   fChain->SetBranchAddress("probe_sa_pt_tgc", &probe_sa_pt_tgc, &b_probe_sa_pt_tgc);
   fChain->SetBranchAddress("probe_sa_pt_ec_radius", &probe_sa_pt_ec_radius, &b_probe_sa_pt_ec_radius);
   fChain->SetBranchAddress("probe_sa_pt_ee", &probe_sa_pt_ee, &b_probe_sa_pt_ee);
   fChain->SetBranchAddress("probe_sa_pt_new", &probe_sa_pt_new, &b_probe_sa_pt_new);
   fChain->SetBranchAddress("probe_sa_pt_alpha_new", &probe_sa_pt_alpha_new, &b_probe_sa_pt_alpha_new);
   fChain->SetBranchAddress("probe_sa_pt_beta_new", &probe_sa_pt_beta_new, &b_probe_sa_pt_beta_new);
   fChain->SetBranchAddress("probe_sa_pt_tgc_new", &probe_sa_pt_tgc_new, &b_probe_sa_pt_tgc_new);
   fChain->SetBranchAddress("probe_sa_alpha", &probe_sa_alpha, &b_probe_sa_alpha);
   fChain->SetBranchAddress("probe_sa_beta", &probe_sa_beta, &b_probe_sa_beta);
   fChain->SetBranchAddress("probe_sa_ec_radius", &probe_sa_ec_radius, &b_probe_sa_ec_radius);
   fChain->SetBranchAddress("probe_sa_br_radius", &probe_sa_br_radius, &b_probe_sa_br_radius);
   fChain->SetBranchAddress("probe_sa_eta", &probe_sa_eta, &b_probe_sa_eta);
   fChain->SetBranchAddress("probe_sa_phi", &probe_sa_phi, &b_probe_sa_phi);
   fChain->SetBranchAddress("probe_sa_deltaPt", &probe_sa_deltaPt, &b_probe_sa_deltaPt);
   fChain->SetBranchAddress("probe_sa_charge", &probe_sa_charge, &b_probe_sa_charge);
   fChain->SetBranchAddress("probe_sa_saddress", &probe_sa_saddress, &b_probe_sa_saddress);
   fChain->SetBranchAddress("probe_sa_ec_eta_bin", &probe_sa_ec_eta_bin, &b_probe_sa_ec_eta_bin);
   fChain->SetBranchAddress("probe_sa_ec_phi_bin", &probe_sa_ec_phi_bin, &b_probe_sa_ec_phi_bin);
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("probe_pass_L1mu4", &probe_pass_L1mu4, &b_probe_pass_L1mu4);
   fChain->SetBranchAddress("probe_pass_SAmu4", &probe_pass_SAmu4, &b_probe_pass_SAmu4);
   fChain->SetBranchAddress("probe_pass_SAmu4_new", &probe_pass_SAmu4_new, &b_probe_pass_SAmu4_new);
   fChain->SetBranchAddress("probe_pass_SAmu4_ee", &probe_pass_SAmu4_ee, &b_probe_pass_SAmu4_ee);
   fChain->SetBranchAddress("probe_pass_L1mu6", &probe_pass_L1mu6, &b_probe_pass_L1mu6);
   fChain->SetBranchAddress("probe_pass_SAmu6", &probe_pass_SAmu6, &b_probe_pass_SAmu6);
   fChain->SetBranchAddress("probe_pass_SAmu6_new", &probe_pass_SAmu6_new, &b_probe_pass_SAmu6_new);
   fChain->SetBranchAddress("probe_pass_SAmu6_ee", &probe_pass_SAmu6_ee, &b_probe_pass_SAmu6_ee);
   //fChain->SetBranchAddress("sp_r", &sp_r, &b_sp_r);
   //fChain->SetBranchAddress("sp_z", &sp_z, &b_sp_z);
   //fChain->SetBranchAddress("segment_r", &segment_r, &b_segment_r);
   //fChain->SetBranchAddress("segment_z", &segment_z, &b_segment_z);
   //fChain->SetBranchAddress("segment_chamber", &segment_chamber, &b_segment_chamber);
   Notify();
}

Bool_t effJpsi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void effJpsi::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t effJpsi::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef effJpsi_cxx
