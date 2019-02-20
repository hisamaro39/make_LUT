//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 10 09:49:11 2016 by ROOT version 6.02/08
// from TTree t_tap/TrigMuonTagAndProbe
// found on file: test2.root
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
#include "vector"

class check2muon {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        probe_offline_pt;
   Double_t        probe_offline_eta;
   Double_t        probe_offline_phi;
   Int_t           probe_offline_charge;
   Double_t        tag_offline_pt;
   Double_t        tag_offline_eta;
   Double_t        tag_offline_phi;
   Int_t           tag_offline_charge;
   Double_t        probe_ext_eta;
   Double_t        probe_ext_phi;
   Double_t        tag_ext_eta;
   Double_t        tag_ext_phi;
   Double_t        probe_roi_eta;
   Double_t        probe_roi_phi;
   Double_t        tag_roi_eta;
   Double_t        tag_roi_phi;
   Double_t        probe_sa_pt;
   Double_t        probe_sa_eta;
   Double_t        probe_sa_phi;
   Double_t        probe_sa_etaMS;
   Double_t        probe_sa_phiMS;
   Int_t           probe_sa_charge;
   Int_t           probe_saddress;
   Double_t        tag_sa_pt;
   Double_t        tag_sa_eta;
   Double_t        tag_sa_phi;
   Double_t        tag_sa_etaMS;
   Double_t        tag_sa_phiMS;
   Int_t           tag_sa_charge;
   Int_t           tag_saddress;
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
   Bool_t          tag_passL1;
   Bool_t          tag_passSA;
   Bool_t          tag_passComb;
   Bool_t          tag_passEF;
   Bool_t          tag_Combexist;
   Bool_t          tag_EFexist;
   Bool_t          probe_passL1;
   Bool_t          probe_passSA;
   Bool_t          probe_passComb;
   Bool_t          probe_passEF;
   Bool_t          probe_Combexist;
   Bool_t          probe_EFexist;
   vector<float>   *tag_sp_r;
   vector<float>   *tag_sp_z;
   vector<float>   *tag_road_aw;
   vector<float>   *tag_road_bw;
   vector<float>   *tag_mdt_hit_r;
   vector<float>   *tag_mdt_hit_z;
   vector<float>   *tag_mdt_hit_residual;
   vector<float>   *tag_mdt_hit_space;
   vector<int>   *tag_mdt_hit_offlineId;
   vector<int>   *tag_mdt_hit_chamber;
   vector<float>   *tag_tgc_hit_eta;
   vector<float>   *tag_tgc_hit_phi;
   vector<bool>    *tag_tgc_hit_isStrip;
   vector<int>     *tag_tgc_hit_stationNum;
   vector<float>   *tag_rpc_hit_x;
   vector<float>   *tag_rpc_hit_y;
   vector<float>   *tag_rpc_hit_z;
   vector<int>     *tag_rpc_hit_layer;
   vector<float>   *probe_sp_r;
   vector<float>   *probe_sp_z;
   vector<float>   *probe_road_aw;
   vector<float>   *probe_road_bw;
   vector<float>   *probe_mdt_hit_r;
   vector<float>   *probe_mdt_hit_z;
   vector<float>   *probe_mdt_hit_residual;
   vector<float>   *probe_mdt_hit_space;
   vector<int>   *probe_mdt_hit_offlineId;
   vector<int>   *probe_mdt_hit_chamber;
   vector<float>   *probe_tgc_hit_eta;
   vector<float>   *probe_tgc_hit_phi;
   vector<bool>    *probe_tgc_hit_isStrip;
   vector<int>     *probe_tgc_hit_stationNum;
   vector<float>   *probe_rpc_hit_x;
   vector<float>   *probe_rpc_hit_y;
   vector<float>   *probe_rpc_hit_z;
   vector<int>     *probe_rpc_hit_layer;

   // List of branches
   TBranch        *b_probe_offline_pt;   //!
   TBranch        *b_probe_offline_eta;   //!
   TBranch        *b_probe_offline_phi;   //!
   TBranch        *b_probe_offline_charge;   //!
   TBranch        *b_tag_offline_pt;   //!
   TBranch        *b_tag_offline_eta;   //!
   TBranch        *b_tag_offline_phi;   //!
   TBranch        *b_tag_offline_charge;   //!
   TBranch        *b_probe_ext_eta;   //!
   TBranch        *b_probe_ext_phi;   //!
   TBranch        *b_tag_ext_eta;   //!
   TBranch        *b_tag_ext_phi;   //!
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_tag_roi_eta;   //!
   TBranch        *b_tag_roi_phi;   //!
   TBranch        *b_probe_sa_pt;   //!
   TBranch        *b_probe_sa_eta;   //!
   TBranch        *b_probe_sa_phi;   //!
   TBranch        *b_probe_sa_etaMS;   //!
   TBranch        *b_probe_sa_phiMS;   //!
   TBranch        *b_probe_sa_charge;   //!
   TBranch        *b_probe_saddress;   //!
   TBranch        *b_tag_sa_pt;   //!
   TBranch        *b_tag_sa_eta;   //!
   TBranch        *b_tag_sa_phi;   //!
   TBranch        *b_tag_sa_etaMS;   //!
   TBranch        *b_tag_sa_phiMS;   //!
   TBranch        *b_tag_sa_charge;   //!
   TBranch        *b_tag_saddress;   //!
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
   TBranch        *b_tag_passL1;   //!
   TBranch        *b_tag_passSA;   //!
   TBranch        *b_tag_passComb;   //!
   TBranch        *b_tag_passEF;   //!
   TBranch        *b_tag_Combexist;   //!
   TBranch        *b_tag_EFexist;   //!
   TBranch        *b_probe_passL1;   //!
   TBranch        *b_probe_passSA;   //!
   TBranch        *b_probe_passComb;   //!
   TBranch        *b_probe_passEF;   //!
   TBranch        *b_probe_Combexist;   //!
   TBranch        *b_probe_EFexist;   //!
   TBranch        *b_tag_sp_r;   //!
   TBranch        *b_tag_sp_z;   //!
   TBranch        *b_tag_road_aw;   //!
   TBranch        *b_tag_road_bw;   //!
   TBranch        *b_tag_mdt_hit_r;   //!
   TBranch        *b_tag_mdt_hit_z;   //!
   TBranch        *b_tag_mdt_hit_residual;   //!
   TBranch        *b_tag_mdt_hit_space;   //!
   TBranch        *b_tag_mdt_hit_offlineId;   //!
   TBranch        *b_tag_mdt_hit_chamber;   //!
   TBranch        *b_tag_tgc_hit_eta;   //!
   TBranch        *b_tag_tgc_hit_phi;   //!
   TBranch        *b_tag_tgc_hit_isStrip;   //!
   TBranch        *b_tag_tgc_hit_stationNum;   //!
   TBranch        *b_tag_rpc_hit_x;   //!
   TBranch        *b_tag_rpc_hit_y;   //!
   TBranch        *b_tag_rpc_hit_z;   //!
   TBranch        *b_tag_rpc_hit_layer;   //!
   TBranch        *b_probe_sp_r;   //!
   TBranch        *b_probe_sp_z;   //!
   TBranch        *b_probe_road_aw;   //!
   TBranch        *b_probe_road_bw;   //!
   TBranch        *b_probe_mdt_hit_r;   //!
   TBranch        *b_probe_mdt_hit_z;   //!
   TBranch        *b_probe_mdt_hit_residual;   //!
   TBranch        *b_probe_mdt_hit_space;   //!
   TBranch        *b_probe_mdt_hit_offlineId;   //!
   TBranch        *b_probe_mdt_hit_chamber;   //!
   TBranch        *b_probe_tgc_hit_eta;   //!
   TBranch        *b_probe_tgc_hit_phi;   //!
   TBranch        *b_probe_tgc_hit_isStrip;   //!
   TBranch        *b_probe_tgc_hit_stationNum;   //!
   TBranch        *b_probe_rpc_hit_x;   //!
   TBranch        *b_probe_rpc_hit_y;   //!
   TBranch        *b_probe_rpc_hit_z;   //!
   TBranch        *b_probe_rpc_hit_layer;   //!

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
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLTmu4msonly.list";
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLT2mu4.list";
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLT2mu4Default.list";
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLT2mu4AddRoadSeparation.list";
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLT2mu4AddRoadSeparationNewestMuComb.list";
  const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLTmu4msonlyDefault.list";
  //const char* rec = "datalist/inputEfficiency2muonJpsimu4mu20HLTmu4msonlyAddRoadSeparation.list";

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

   // Set object pointer
   tag_sp_r = 0;
   tag_sp_z = 0;
   tag_road_aw = 0;
   tag_road_bw = 0;
   tag_mdt_hit_r = 0;
   tag_mdt_hit_z = 0;
   tag_mdt_hit_residual = 0;
   tag_mdt_hit_space = 0;
   tag_mdt_hit_offlineId = 0;
   tag_mdt_hit_chamber = 0;
   tag_tgc_hit_eta = 0;
   tag_tgc_hit_phi = 0;
   tag_tgc_hit_isStrip = 0;
   tag_tgc_hit_stationNum = 0;
   tag_rpc_hit_x = 0;
   tag_rpc_hit_y = 0;
   tag_rpc_hit_z = 0;
   tag_rpc_hit_layer = 0;
   probe_sp_r = 0;
   probe_sp_z = 0;
   probe_road_aw = 0;
   probe_road_bw = 0;
   probe_mdt_hit_r = 0;
   probe_mdt_hit_z = 0;
   probe_mdt_hit_residual = 0;
   probe_mdt_hit_space = 0;
   probe_mdt_hit_offlineId = 0;
   probe_mdt_hit_chamber = 0;
   probe_tgc_hit_eta = 0;
   probe_tgc_hit_phi = 0;
   probe_tgc_hit_isStrip = 0;
   probe_tgc_hit_stationNum = 0;
   probe_rpc_hit_x = 0;
   probe_rpc_hit_y = 0;
   probe_rpc_hit_z = 0;
   probe_rpc_hit_layer = 0;
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
   fChain->SetBranchAddress("probe_ext_eta", &probe_ext_eta, &b_probe_ext_eta);
   fChain->SetBranchAddress("probe_ext_phi", &probe_ext_phi, &b_probe_ext_phi);
   fChain->SetBranchAddress("tag_ext_eta", &tag_ext_eta, &b_tag_ext_eta);
   fChain->SetBranchAddress("tag_ext_phi", &tag_ext_phi, &b_tag_ext_phi);
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("tag_roi_eta", &tag_roi_eta, &b_tag_roi_eta);
   fChain->SetBranchAddress("tag_roi_phi", &tag_roi_phi, &b_tag_roi_phi);
   fChain->SetBranchAddress("probe_sa_pt", &probe_sa_pt, &b_probe_sa_pt);
   fChain->SetBranchAddress("probe_sa_eta", &probe_sa_eta, &b_probe_sa_eta);
   fChain->SetBranchAddress("probe_sa_phi", &probe_sa_phi, &b_probe_sa_phi);
   fChain->SetBranchAddress("probe_sa_etaMS", &probe_sa_etaMS, &b_probe_sa_etaMS);
   fChain->SetBranchAddress("probe_sa_phiMS", &probe_sa_phiMS, &b_probe_sa_phiMS);
   fChain->SetBranchAddress("probe_sa_charge", &probe_sa_charge, &b_probe_sa_charge);
   fChain->SetBranchAddress("probe_saddress", &probe_saddress, &b_probe_saddress);
   fChain->SetBranchAddress("tag_sa_pt", &tag_sa_pt, &b_tag_sa_pt);
   fChain->SetBranchAddress("tag_sa_eta", &tag_sa_eta, &b_tag_sa_eta);
   fChain->SetBranchAddress("tag_sa_phi", &tag_sa_phi, &b_tag_sa_phi);
   fChain->SetBranchAddress("tag_sa_etaMS", &tag_sa_etaMS, &b_tag_sa_etaMS);
   fChain->SetBranchAddress("tag_sa_phiMS", &tag_sa_phiMS, &b_tag_sa_phiMS);
   fChain->SetBranchAddress("tag_sa_charge", &tag_sa_charge, &b_tag_sa_charge);
   fChain->SetBranchAddress("tag_saddress", &tag_saddress, &b_tag_saddress);
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
   fChain->SetBranchAddress("tag_passL1", &tag_passL1, &b_tag_passL1);
   fChain->SetBranchAddress("tag_passSA", &tag_passSA, &b_tag_passSA);
   fChain->SetBranchAddress("tag_passComb", &tag_passComb, &b_tag_passComb);
   fChain->SetBranchAddress("tag_passEF", &tag_passEF, &b_tag_passEF);
   fChain->SetBranchAddress("tag_Combexist", &tag_Combexist, &b_tag_Combexist);
   fChain->SetBranchAddress("tag_EFexist", &tag_EFexist, &b_tag_EFexist);
   fChain->SetBranchAddress("probe_passL1", &probe_passL1, &b_probe_passL1);
   fChain->SetBranchAddress("probe_passSA", &probe_passSA, &b_probe_passSA);
   fChain->SetBranchAddress("probe_passComb", &probe_passComb, &b_probe_passComb);
   fChain->SetBranchAddress("probe_passEF", &probe_passEF, &b_probe_passEF);
   fChain->SetBranchAddress("probe_Combexist", &probe_Combexist, &b_probe_Combexist);
   fChain->SetBranchAddress("probe_EFexist", &probe_EFexist, &b_probe_EFexist);
   fChain->SetBranchAddress("tag_sp_r", &tag_sp_r, &b_tag_sp_r);
   fChain->SetBranchAddress("tag_sp_z", &tag_sp_z, &b_tag_sp_z);
   fChain->SetBranchAddress("tag_road_aw", &tag_road_aw, &b_tag_road_aw);
   fChain->SetBranchAddress("tag_road_bw", &tag_road_bw, &b_tag_road_bw);
   fChain->SetBranchAddress("tag_mdt_hit_r", &tag_mdt_hit_r, &b_tag_mdt_hit_r);
   fChain->SetBranchAddress("tag_mdt_hit_z", &tag_mdt_hit_z, &b_tag_mdt_hit_z);
   fChain->SetBranchAddress("tag_mdt_hit_residual", &tag_mdt_hit_residual, &b_tag_mdt_hit_residual);
   fChain->SetBranchAddress("tag_mdt_hit_space", &tag_mdt_hit_space, &b_tag_mdt_hit_space);
   fChain->SetBranchAddress("tag_mdt_hit_offlineId", &tag_mdt_hit_offlineId, &b_tag_mdt_hit_offlineId);
   fChain->SetBranchAddress("tag_mdt_hit_chamber", &tag_mdt_hit_chamber, &b_tag_mdt_hit_chamber);
   fChain->SetBranchAddress("tag_tgc_hit_eta", &tag_tgc_hit_eta, &b_tag_tgc_hit_eta);
   fChain->SetBranchAddress("tag_tgc_hit_phi", &tag_tgc_hit_phi, &b_tag_tgc_hit_phi);
   fChain->SetBranchAddress("tag_tgc_hit_isStrip", &tag_tgc_hit_isStrip, &b_tag_tgc_hit_isStrip);
   fChain->SetBranchAddress("tag_tgc_hit_stationNum", &tag_tgc_hit_stationNum, &b_tag_tgc_hit_stationNum);
   fChain->SetBranchAddress("tag_rpc_hit_x", &tag_rpc_hit_x, &b_tag_rpc_hit_x);
   fChain->SetBranchAddress("tag_rpc_hit_y", &tag_rpc_hit_y, &b_tag_rpc_hit_y);
   fChain->SetBranchAddress("tag_rpc_hit_z", &tag_rpc_hit_z, &b_tag_rpc_hit_z);
   fChain->SetBranchAddress("tag_rpc_hit_layer", &tag_rpc_hit_layer, &b_tag_rpc_hit_layer);
   fChain->SetBranchAddress("probe_sp_r", &probe_sp_r, &b_probe_sp_r);
   fChain->SetBranchAddress("probe_sp_z", &probe_sp_z, &b_probe_sp_z);
   fChain->SetBranchAddress("probe_road_aw", &probe_road_aw, &b_probe_road_aw);
   fChain->SetBranchAddress("probe_road_bw", &probe_road_bw, &b_probe_road_bw);
   fChain->SetBranchAddress("probe_mdt_hit_r", &probe_mdt_hit_r, &b_probe_mdt_hit_r);
   fChain->SetBranchAddress("probe_mdt_hit_z", &probe_mdt_hit_z, &b_probe_mdt_hit_z);
   fChain->SetBranchAddress("probe_mdt_hit_residual", &probe_mdt_hit_residual, &b_probe_mdt_hit_residual);
   fChain->SetBranchAddress("probe_mdt_hit_space", &probe_mdt_hit_space, &b_probe_mdt_hit_space);
   fChain->SetBranchAddress("probe_mdt_hit_offlineId", &probe_mdt_hit_offlineId, &b_probe_mdt_hit_offlineId);
   fChain->SetBranchAddress("probe_mdt_hit_chamber", &probe_mdt_hit_chamber, &b_probe_mdt_hit_chamber);
   fChain->SetBranchAddress("probe_tgc_hit_eta", &probe_tgc_hit_eta, &b_probe_tgc_hit_eta);
   fChain->SetBranchAddress("probe_tgc_hit_phi", &probe_tgc_hit_phi, &b_probe_tgc_hit_phi);
   fChain->SetBranchAddress("probe_tgc_hit_isStrip", &probe_tgc_hit_isStrip, &b_probe_tgc_hit_isStrip);
   fChain->SetBranchAddress("probe_tgc_hit_stationNum", &probe_tgc_hit_stationNum, &b_probe_tgc_hit_stationNum);
   fChain->SetBranchAddress("probe_rpc_hit_x", &probe_rpc_hit_x, &b_probe_rpc_hit_x);
   fChain->SetBranchAddress("probe_rpc_hit_y", &probe_rpc_hit_y, &b_probe_rpc_hit_y);
   fChain->SetBranchAddress("probe_rpc_hit_z", &probe_rpc_hit_z, &b_probe_rpc_hit_z);
   fChain->SetBranchAddress("probe_rpc_hit_layer", &probe_rpc_hit_layer, &b_probe_rpc_hit_layer);
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
