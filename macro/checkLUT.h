//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 27 22:46:21 2015 by ROOT version 6.02/05
// from TTree validationT/validationT
// found on file: hist-newEndcapLUT.root
//////////////////////////////////////////////////////////

#ifndef checkLUT_h
#define checkLUT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
using namespace std;

// Header file for the classes stored in the TTree if any.

class checkLUT {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         L1_eta;
   Float_t         L1_phi;
   Int_t           event;
   Float_t         L2_pt;
   Float_t         L2_pt_alpha;
   Float_t         L2_pt_beta;
   Float_t         L2_pt_ec_radius;
   Float_t         L2_pt_barrel;
   Float_t         L2_pt_alpha_rerun;
   Float_t         L2_pt_beta_rerun;
   Float_t         L2_pt_ec_radius_rerun;
   Float_t         L2_pt_barrel_rerun;
   Float_t         L2_eta;
   Float_t         L2_phi;
   Float_t         L2_etaMap;
   Float_t         L2_phiMap;
   Float_t         L2_charge;
   Float_t         L2_ec_alpha;
   Float_t         L2_ec_beta;
   Float_t         L2_ec_radius;
   Float_t         L2_br_radius;
   Int_t           L2_saddress;
   vector<Float_t> *L2_spr;
   vector<Float_t> *L2_spz;
   Float_t         tgcPt;
   Float_t         tgcMid1_r;
   Float_t         tgcMid1_z;
   Float_t         tgcMid2_r;
   Float_t         tgcMid2_z;
   Float_t         offline_pt;
   Float_t         offline_eta;
   Float_t         offline_phi;
   vector<Float_t> *offline_segment_r;
   vector<Float_t> *offline_segment_z;
   vector<Float_t> *offline_segment_eta;
   vector<Float_t> *offline_segment_phi;
   vector<Int_t> *offline_segment_chamber;

   // List of branches
   TBranch        *b_L1_eta;   //!
   TBranch        *b_L1_phi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_L2_pt;   //!
   TBranch        *b_L2_pt_alpha;   //!
   TBranch        *b_L2_pt_beta;   //!
   TBranch        *b_L2_pt_ec_radius;   //!
   TBranch        *b_L2_pt_barrel;   //!
   TBranch        *b_L2_pt_alpha_rerun;   //!
   TBranch        *b_L2_pt_beta_rerun;   //!
   TBranch        *b_L2_pt_ec_radius_rerun;   //!
   TBranch        *b_L2_pt_barrel_rerun;   //!
   TBranch        *b_L2_eta;   //!
   TBranch        *b_L2_phi;   //!
   TBranch        *b_L2_etaMap;   //!
   TBranch        *b_L2_phiMap;   //!
   TBranch        *b_L2_charge;   //!
   TBranch        *b_L2_ec_alpha;   //!
   TBranch        *b_L2_ec_beta;   //!
   TBranch        *b_L2_ec_radius;   //!
   TBranch        *b_L2_br_radius;   //!
   TBranch        *b_L2_saddress;   //!
   TBranch        *b_L2_spr;   //!
   TBranch        *b_L2_spz;   //!
   TBranch        *b_tgcPt;   //!
   TBranch        *b_tgcMid1_r;   //!
   TBranch        *b_tgcMid1_z;   //!
   TBranch        *b_tgcMid2_r;   //!
   TBranch        *b_tgcMid2_z;   //!
   TBranch        *b_offline_pt;   //!
   TBranch        *b_offline_eta;   //!
   TBranch        *b_offline_phi;   //!
   TBranch        *b_offline_segment_r;   //!
   TBranch        *b_offline_segment_z;   //!
   TBranch        *b_offline_segment_eta;   //!
   TBranch        *b_offline_segment_phi;   //!
   TBranch        *b_offline_segment_chamber;   //!

   checkLUT(TTree *tree=0);
   virtual ~checkLUT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef checkLUT_cxx
checkLUT::checkLUT(TTree *tree) : fChain(0) 
{

  TChain *chain = new TChain("validationT");
  const char* rec = "datalist/inputCheckLUT.list";
  //const char* rec = "datalist/test.list";
  //const char* rec = "datalist/inputRerunDefault.list";
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

checkLUT::~checkLUT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t checkLUT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t checkLUT::LoadTree(Long64_t entry)
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

void checkLUT::Init(TTree *tree)
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

   L2_spr=0;
   L2_spz=0;
   offline_segment_r=0;
   offline_segment_z=0;
   offline_segment_eta=0;
   offline_segment_phi=0;
   offline_segment_chamber=0;


   fChain->SetBranchAddress("L1_eta", &L1_eta, &b_L1_eta);
   fChain->SetBranchAddress("L1_phi", &L1_phi, &b_L1_phi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("L2_pt", &L2_pt, &b_L2_pt);
   fChain->SetBranchAddress("L2_pt_alpha", &L2_pt_alpha, &b_L2_pt_alpha);
   fChain->SetBranchAddress("L2_pt_beta", &L2_pt_beta, &b_L2_pt_beta);
   fChain->SetBranchAddress("L2_pt_ec_radius", &L2_pt_ec_radius, &b_L2_pt_ec_radius);
   fChain->SetBranchAddress("L2_pt_barrel", &L2_pt_barrel, &b_L2_pt_barrel);
   fChain->SetBranchAddress("L2_pt_alpha_rerun", &L2_pt_alpha_rerun, &b_L2_pt_alpha_rerun);
   fChain->SetBranchAddress("L2_pt_beta_rerun", &L2_pt_beta_rerun, &b_L2_pt_beta_rerun);
   fChain->SetBranchAddress("L2_pt_ec_radius_rerun", &L2_pt_ec_radius_rerun, &b_L2_pt_ec_radius_rerun);
   fChain->SetBranchAddress("L2_pt_barrel_rerun", &L2_pt_barrel_rerun, &b_L2_pt_barrel_rerun);
   fChain->SetBranchAddress("L2_eta", &L2_eta, &b_L2_eta);
   fChain->SetBranchAddress("L2_phi", &L2_phi, &b_L2_phi);
   fChain->SetBranchAddress("L2_etaMap", &L2_etaMap, &b_L2_etaMap);
   fChain->SetBranchAddress("L2_phiMap", &L2_phiMap, &b_L2_phiMap);
   fChain->SetBranchAddress("L2_charge", &L2_charge, &b_L2_charge);
   fChain->SetBranchAddress("L2_ec_alpha", &L2_ec_alpha, &b_L2_ec_alpha);
   fChain->SetBranchAddress("L2_ec_beta", &L2_ec_beta, &b_L2_ec_beta);
   fChain->SetBranchAddress("L2_ec_radius", &L2_ec_radius, &b_L2_ec_radius);
   fChain->SetBranchAddress("L2_br_radius", &L2_br_radius, &b_L2_br_radius);
   fChain->SetBranchAddress("L2_saddress", &L2_saddress, &b_L2_saddress);
   fChain->SetBranchAddress("L2_spr", &L2_spr, &b_L2_spr);
   fChain->SetBranchAddress("L2_spz", &L2_spz, &b_L2_spz);
   fChain->SetBranchAddress("tgcPt", &tgcPt, &b_tgcPt);
   fChain->SetBranchAddress("tgcMid1_r", &tgcMid1_r, &b_tgcMid1_r);
   fChain->SetBranchAddress("tgcMid1_z", &tgcMid1_z, &b_tgcMid1_z);
   fChain->SetBranchAddress("tgcMid2_r", &tgcMid2_r, &b_tgcMid2_r);
   fChain->SetBranchAddress("tgcMid2_z", &tgcMid2_z, &b_tgcMid2_z);
   fChain->SetBranchAddress("offline_pt", &offline_pt, &b_offline_pt);
   fChain->SetBranchAddress("offline_eta", &offline_eta, &b_offline_eta);
   fChain->SetBranchAddress("offline_phi", &offline_phi, &b_offline_phi);
   fChain->SetBranchAddress("offline_segment_r", &offline_segment_r, &b_offline_segment_r);
   fChain->SetBranchAddress("offline_segment_z", &offline_segment_z, &b_offline_segment_z);
   fChain->SetBranchAddress("offline_segment_eta", &offline_segment_eta, &b_offline_segment_eta);
   fChain->SetBranchAddress("offline_segment_phi", &offline_segment_phi, &b_offline_segment_phi);
   fChain->SetBranchAddress("offline_segment_chamber", &offline_segment_chamber, &b_offline_segment_chamber);
   Notify();
}

Bool_t checkLUT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void checkLUT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t checkLUT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef checkLUT_cxx
