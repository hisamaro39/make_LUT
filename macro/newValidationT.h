//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 25 20:30:28 2015 by ROOT version 5.34/18
// from TTree newValidationT/newValidationT
// found on file: hist-valid3.147407.PowhegPythia8_AZNLO_Zmumu.merge.DAOD_MUON0.e3099_s2579_r6164_p2324_tid05261900_00.root
//////////////////////////////////////////////////////////

#ifndef newValidationT_h
#define newValidationT_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include "macro/readLUT.h"
#include "macro/util.h"
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class newValidationT {
public :
  readLUT m_readLUT;
  util m_util;


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         L2_pt;
   //Float_t         L2_pt_alpha;
   Float_t         L2_eta;
   Float_t         L2_phi;
   Float_t         L1_eta;
   Float_t         L1_phi;
   Float_t         L2_etaMap;
   Float_t         L2_phiMap;
   Int_t         L2_charge;
   Float_t         L2_ec_alpha;
   Float_t         L2_ec_beta;
   Float_t         L2_ec_radius;
   Float_t         L2_br_radius;
   //Float_t         L2_road_aw0;
   //Float_t         L2_road_aw1;
   //Float_t         L2_road_aw2;
   //Float_t         L2_road_bw0;
   //Float_t         L2_road_bw1;
   //Float_t         L2_road_bw2;
   Int_t           L2_saddress;
   //Int_t           L2_SL;
   //Float_t         tgcPt;
   Float_t         tgcMid1_r;
   Float_t         tgcMid1_z;
   Float_t         tgcMid1_eta;
   Float_t         tgcMid1_phi;
   Float_t         tgcMid2_r;
   Float_t         tgcMid2_z;
   Float_t         offline_pt;
   Float_t         offline_eta;
   Float_t           offline_phi;
   Int_t         offline_charge;
   vector<Float_t> *sp_r;
   vector<Float_t> *sp_z;
   //vector<Float_t> *mdt_hit_r;
   //vector<Float_t> *mdt_hit_z;
   //vector<Float_t> *mdt_hit_residual;
   vector<Float_t> *segment_r;
   vector<Float_t> *segment_z;
   vector<Float_t> *segment_eta;
   vector<Float_t> *segment_phi;
   vector<Int_t> *segment_chamber;

   // List of branches
   TBranch        *b_L2_pt;   //!
   //TBranch        *b_L2_pt_alpha;   //!
   TBranch        *b_L2_eta;   //!
   TBranch        *b_L2_phi;   //!
   TBranch        *b_L1_eta;   //!
   TBranch        *b_L1_phi;   //!
   TBranch        *b_L2_etaMap;   //!
   TBranch        *b_L2_phiMap;   //!
   TBranch        *b_L2_charge;   //!
   TBranch        *b_L2_ec_alpha;   //!
   TBranch        *b_L2_ec_beta;   //!
   TBranch        *b_L2_ec_radius;   //!
   TBranch        *b_L2_br_radius;   //!
   //TBranch        *b_L2_road_aw0;   //!
   //TBranch        *b_L2_road_aw1;   //!
   //TBranch        *b_L2_road_aw2;   //!
   //TBranch        *b_L2_road_bw0;   //!
   //TBranch        *b_L2_road_bw1;   //!
   //TBranch        *b_L2_road_bw2;   //!
   TBranch        *b_L2_saddress;   //!
   //TBranch        *b_L2_SL;   //!
   //TBranch        *b_tgcPt;   //!
   TBranch        *b_tgcMid1_r;   //!
   TBranch        *b_tgcMid1_z;   //!
   TBranch        *b_tgcMid1_eta;   //!
   TBranch        *b_tgcMid1_phi;   //!
   TBranch        *b_tgcMid2_r;   //!
   TBranch        *b_tgcMid2_z;   //!
   TBranch        *b_offline_pt;   //!
   TBranch        *b_offline_eta;   //!
   TBranch        *b_offline_phi;   //!
   TBranch        *b_offline_charge;   //!
   TBranch        *b_sp_r;   //!
   TBranch        *b_sp_z;   //!
   //TBranch        *b_mdt_hit_r;   //!
   //TBranch        *b_mdt_hit_z;   //!
   //TBranch        *b_mdt_hit_residual;   //!
   TBranch        *b_segment_r;   //!
   TBranch        *b_segment_z;   //!
   TBranch        *b_segment_eta;   //!
   TBranch        *b_segment_phi;   //!
   TBranch        *b_segment_chamber;   //!

   newValidationT();
   virtual ~newValidationT();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(string list,string output);
   virtual void     makePtResidual(string list,string output);
   virtual void     fitResidual(string input, string output); 
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef newValidationT_cxx
newValidationT::newValidationT(){}

newValidationT::~newValidationT()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t newValidationT::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t newValidationT::LoadTree(Long64_t entry)
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

void newValidationT::Init(TTree *tree)
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

   sp_r=0;
   sp_z=0;
   //mdt_hit_r=0;
   //mdt_hit_z=0;
   //mdt_hit_residual=0;
   segment_r=0;
   segment_z=0;
   segment_eta=0;
   segment_phi=0;
   segment_chamber=0;

   fChain->SetBranchAddress("L2_pt", &L2_pt, &b_L2_pt);
   //fChain->SetBranchAddress("L2_pt_alpha", &L2_pt_alpha, &b_L2_pt_alpha);
   fChain->SetBranchAddress("L2_eta", &L2_eta, &b_L2_eta);
   fChain->SetBranchAddress("L2_phi", &L2_phi, &b_L2_phi);
   fChain->SetBranchAddress("L1_eta", &L1_eta, &b_L1_eta);
   fChain->SetBranchAddress("L1_phi", &L1_phi, &b_L1_phi);
   fChain->SetBranchAddress("L2_etaMap", &L2_etaMap, &b_L2_etaMap);
   fChain->SetBranchAddress("L2_phiMap", &L2_phiMap, &b_L2_phiMap);
   fChain->SetBranchAddress("L2_charge", &L2_charge, &b_L2_charge);
   fChain->SetBranchAddress("L2_ec_alpha", &L2_ec_alpha, &b_L2_ec_alpha);
   fChain->SetBranchAddress("L2_ec_beta", &L2_ec_beta, &b_L2_ec_beta);
   fChain->SetBranchAddress("L2_ec_radius", &L2_ec_radius, &b_L2_ec_radius);
   fChain->SetBranchAddress("L2_br_radius", &L2_br_radius, &b_L2_br_radius);
   //fChain->SetBranchAddress("L2_road_aw0", &L2_road_aw0, &b_L2_road_aw0);
   //fChain->SetBranchAddress("L2_road_aw1", &L2_road_aw1, &b_L2_road_aw1);
   //fChain->SetBranchAddress("L2_road_aw2", &L2_road_aw2, &b_L2_road_aw2);
   //fChain->SetBranchAddress("L2_road_bw0", &L2_road_bw0, &b_L2_road_bw0);
   //fChain->SetBranchAddress("L2_road_bw1", &L2_road_bw1, &b_L2_road_bw1);
   //fChain->SetBranchAddress("L2_road_bw2", &L2_road_bw2, &b_L2_road_bw2);
   fChain->SetBranchAddress("L2_saddress", &L2_saddress, &b_L2_saddress);
   //fChain->SetBranchAddress("L2_SL", &L2_SL, &b_L2_SL);
   //fChain->SetBranchAddress("tgcPt", &tgcPt, &b_tgcPt);
   fChain->SetBranchAddress("tgcMid1_r", &tgcMid1_r, &b_tgcMid1_r);
   fChain->SetBranchAddress("tgcMid1_z", &tgcMid1_z, &b_tgcMid1_z);
   fChain->SetBranchAddress("tgcMid1_eta", &tgcMid1_eta, &b_tgcMid1_eta);
   fChain->SetBranchAddress("tgcMid1_phi", &tgcMid1_phi, &b_tgcMid1_phi);
   fChain->SetBranchAddress("tgcMid2_r", &tgcMid2_r, &b_tgcMid2_r);
   fChain->SetBranchAddress("tgcMid2_z", &tgcMid2_z, &b_tgcMid2_z);
   fChain->SetBranchAddress("offline_pt", &offline_pt, &b_offline_pt);
   fChain->SetBranchAddress("offline_eta", &offline_eta, &b_offline_eta);
   fChain->SetBranchAddress("offline_phi", &offline_phi, &b_offline_phi);
   fChain->SetBranchAddress("offline_charge", &offline_charge, &b_offline_charge);
   fChain->SetBranchAddress("sp_r", &sp_r, &b_sp_r);
   fChain->SetBranchAddress("sp_z", &sp_z, &b_sp_z);
   //fChain->SetBranchAddress("mdt_hit_r", &mdt_hit_r, &b_mdt_hit_r);
   //fChain->SetBranchAddress("mdt_hit_z", &mdt_hit_z, &b_mdt_hit_z);
   //fChain->SetBranchAddress("mdt_hit_residual", &mdt_hit_residual, &b_mdt_hit_residual);
   //fChain->SetBranchAddress("segment_r", &segment_r, &b_segment_r);
   //fChain->SetBranchAddress("segment_z", &segment_z, &b_segment_z);
   //fChain->SetBranchAddress("segment_eta", &segment_eta, &b_segment_eta);
   //fChain->SetBranchAddress("segment_phi", &segment_phi, &b_segment_phi);
   //fChain->SetBranchAddress("segment_chamber", &segment_chamber, &b_segment_chamber);
   Notify();
}

Bool_t newValidationT::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void newValidationT::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t newValidationT::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef newValidationT_cxx
