//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun  6 18:18:31 2016 by ROOT version 5.34/18
// from TTree validationT/validationT
// found on file: hist-DAOD_MUON0.08384499._000001.pool.root.1.root
//////////////////////////////////////////////////////////

#ifndef ScaleFactorTrack_h
#define ScaleFactorTrack_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
#include "macro/util.h"
#include "macro/makeEfficiencyPlot.h"
using namespace std;

// Header file for the classes stored in the TTree if any.
#include <string>
#include <vector>
#include "macro/util.h"

// Fixed size dimensions of array or collections stored in the TTree if any.

class ScaleFactorTrack {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   util m_util;
   makeEffPlot m_makeEffPlot;

   // Declaration of leaf types
   Float_t         tag_offline_pt;
   Float_t         tag_offline_eta;
   Float_t         tag_offline_phi;
   Float_t         tag_offline_charge;
   Float_t         tag_offline_iso;
   string          *tag_trigger_chain;
   Int_t          *tag_trigger_threshold;
   Float_t         probe_offline_pt;
   Float_t         probe_offline_eta;
   Float_t         probe_offline_phi;
   Float_t         probe_offline_charge;
   Float_t         probe_offline_iso;
   Float_t         probe_roi_eta;
   Float_t         probe_roi_phi;
   Float_t         probe_sa_pt;
   Float_t         probe_sa_eta;
   Float_t         probe_sa_phi;
   Int_t           probe_sa_charge;
   Int_t           probe_sa_saddress;
   Float_t         probe_comb_pt;
   Float_t         probe_comb_eta;
   Float_t         probe_comb_phi;
   Int_t           probe_comb_charge;
   Float_t         probe_ef_pt;
   Float_t         probe_ef_eta;
   Float_t         probe_ef_phi;
   Int_t           probe_ef_charge;
   vector<bool>    *probe_passL1;
   vector<bool>    *probe_passSA;
   vector<bool>    *probe_passComb;
   vector<bool>    *probe_passEF;
   vector<string>  *probe_trigger_chain;
   vector<int>  *probe_trigger_threshold;
   vector<float>    *track_pt;
   vector<float>    *track_eta;
   vector<float>    *track_phi;

   // List of branches
   TBranch        *b_tag_offline_pt;   //!
   TBranch        *b_tag_offline_eta;   //!
   TBranch        *b_tag_offline_phi;   //!
   TBranch        *b_tag_offline_charge;   //!
   TBranch        *b_tag_offline_iso;   //!
   TBranch        *b_tag_trigger_chain;   //!
   TBranch        *b_tag_trigger_threshold;   //!
   TBranch        *b_probe_offline_pt;   //!
   TBranch        *b_probe_offline_eta;   //!
   TBranch        *b_probe_offline_phi;   //!
   TBranch        *b_probe_offline_charge;   //!
   TBranch        *b_probe_offline_iso;   //!
   TBranch        *b_probe_roi_eta;   //!
   TBranch        *b_probe_roi_phi;   //!
   TBranch        *b_probe_sa_pt;   //!
   TBranch        *b_probe_sa_eta;   //!
   TBranch        *b_probe_sa_phi;   //!
   TBranch        *b_probe_sa_charge;   //!
   TBranch        *b_probe_sa_saddress;   //!
   TBranch        *b_probe_comb_pt;   //!
   TBranch        *b_probe_comb_eta;   //!
   TBranch        *b_probe_comb_phi;   //!
   TBranch        *b_probe_comb_charge;   //!
   TBranch        *b_probe_ef_pt;   //!
   TBranch        *b_probe_ef_eta;   //!
   TBranch        *b_probe_ef_phi;   //!
   TBranch        *b_probe_ef_charge;   //!
   TBranch        *b_probe_passL1;   //!
   TBranch        *b_probe_passSA;   //!
   TBranch        *b_probe_passComb;   //!
   TBranch        *b_probe_passEF;   //!
   TBranch        *b_probe_trigger_chain;   //!
   TBranch        *b_probe_trigger_threshold;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!

   ScaleFactorTrack();
   virtual ~ScaleFactorTrack();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const char* list, const char* output);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void getThreshold(const char *list,vector<int> &chain_list,vector<string> &chain_name_list);
   void     makeResolutionPlot(vector<string> list, string output);
   void     fitResolution(string input, string output);
};

#endif

#ifdef ScaleFactorTrack_cxx
ScaleFactorTrack::ScaleFactorTrack(){}

ScaleFactorTrack::~ScaleFactorTrack(){}

Int_t ScaleFactorTrack::GetEntry(Long64_t entry)
{
// Read contents of entry.
   //if (!fChain) return 0;
   //return fChain->GetEntry(entry);
}

Long64_t ScaleFactorTrack::LoadTree(Long64_t entry)
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

void ScaleFactorTrack::Init(TTree *tree)
{
   // Set object pointer
   tag_trigger_chain = 0;
   tag_trigger_threshold = 0;
   probe_passL1 = 0;
   probe_passSA = 0;
   probe_passComb = 0;
   probe_passEF = 0;
   probe_trigger_chain = 0;
   probe_trigger_threshold = 0;
   track_pt = 0;
   track_eta = 0;
   track_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("tag_offline_pt", &tag_offline_pt, &b_tag_offline_pt);
   fChain->SetBranchAddress("tag_offline_eta", &tag_offline_eta, &b_tag_offline_eta);
   fChain->SetBranchAddress("tag_offline_phi", &tag_offline_phi, &b_tag_offline_phi);
   fChain->SetBranchAddress("tag_offline_charge", &tag_offline_charge, &b_tag_offline_charge);
   fChain->SetBranchAddress("tag_offline_iso", &tag_offline_iso, &b_tag_offline_iso);
   fChain->SetBranchAddress("tag_trigger_chain", &tag_trigger_chain, &b_tag_trigger_chain);
   fChain->SetBranchAddress("tag_trigger_threshold", &tag_trigger_threshold, &b_tag_trigger_threshold);
   fChain->SetBranchAddress("probe_offline_pt", &probe_offline_pt, &b_probe_offline_pt);
   fChain->SetBranchAddress("probe_offline_eta", &probe_offline_eta, &b_probe_offline_eta);
   fChain->SetBranchAddress("probe_offline_phi", &probe_offline_phi, &b_probe_offline_phi);
   fChain->SetBranchAddress("probe_offline_charge", &probe_offline_charge, &b_probe_offline_charge);
   fChain->SetBranchAddress("probe_offline_iso", &probe_offline_iso, &b_probe_offline_iso);
   fChain->SetBranchAddress("probe_roi_eta", &probe_roi_eta, &b_probe_roi_eta);
   fChain->SetBranchAddress("probe_roi_phi", &probe_roi_phi, &b_probe_roi_phi);
   fChain->SetBranchAddress("probe_sa_pt", &probe_sa_pt, &b_probe_sa_pt);
   fChain->SetBranchAddress("probe_sa_eta", &probe_sa_eta, &b_probe_sa_eta);
   fChain->SetBranchAddress("probe_sa_phi", &probe_sa_phi, &b_probe_sa_phi);
   fChain->SetBranchAddress("probe_sa_charge", &probe_sa_charge, &b_probe_sa_charge);
   fChain->SetBranchAddress("probe_sa_saddress", &probe_sa_saddress, &b_probe_sa_saddress);
   fChain->SetBranchAddress("probe_comb_pt", &probe_comb_pt, &b_probe_comb_pt);
   fChain->SetBranchAddress("probe_comb_eta", &probe_comb_eta, &b_probe_comb_eta);
   fChain->SetBranchAddress("probe_comb_phi", &probe_comb_phi, &b_probe_comb_phi);
   fChain->SetBranchAddress("probe_comb_charge", &probe_comb_charge, &b_probe_comb_charge);
   fChain->SetBranchAddress("probe_ef_pt", &probe_ef_pt, &b_probe_ef_pt);
   fChain->SetBranchAddress("probe_ef_eta", &probe_ef_eta, &b_probe_ef_eta);
   fChain->SetBranchAddress("probe_ef_phi", &probe_ef_phi, &b_probe_ef_phi);
   fChain->SetBranchAddress("probe_ef_charge", &probe_ef_charge, &b_probe_ef_charge);
   fChain->SetBranchAddress("probe_passL1", &probe_passL1, &b_probe_passL1);
   fChain->SetBranchAddress("probe_passSA", &probe_passSA, &b_probe_passSA);
   fChain->SetBranchAddress("probe_passComb", &probe_passComb, &b_probe_passComb);
   fChain->SetBranchAddress("probe_passEF", &probe_passEF, &b_probe_passEF);
   fChain->SetBranchAddress("probe_trigger_chain", &probe_trigger_chain, &b_probe_trigger_chain);
   fChain->SetBranchAddress("probe_trigger_threshold", &probe_trigger_threshold, &b_probe_trigger_threshold);
   fChain->SetBranchAddress("track_pt", &track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_eta", &track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", &track_phi, &b_track_phi);
   Notify();
}

Bool_t ScaleFactorTrack::Notify()
{
   return kTRUE;
}

void ScaleFactorTrack::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ScaleFactorTrack::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ScaleFactorTrack_cxx
