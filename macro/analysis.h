//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 25 20:30:28 2015 by ROOT version 5.34/18
// from TTree analysis/analysis
// found on file: hist-valid3.147407.PowhegPythia8_AZNLO_Zmumu.merge.DAOD_MUON0.e3099_s2579_r6164_p2324_tid05261900_00.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <iostream>
using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Float_t         L2_pt;
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
   Int_t         L2_saddress;
   Float_t         tgcPt;
   Float_t         tgcMid1_r;
   Float_t         tgcMid1_z;
   Float_t         tgcMid1_eta;
   Float_t         tgcMid1_phi;
   Float_t         tgcMid2_r;
   Float_t         tgcMid2_z;
   Float_t         offline_pt;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_L2_pt;   //!
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
   TBranch        *b_L2_saddress;   //!
   TBranch        *b_tgcPt;   //!
   TBranch        *b_tgcMid1_r;   //!
   TBranch        *b_tgcMid1_z;   //!
   TBranch        *b_tgcMid1_eta;   //!
   TBranch        *b_tgcMid1_phi;   //!
   TBranch        *b_tgcMid2_r;   //!
   TBranch        *b_tgcMid2_z;   //!
   TBranch        *b_offline_pt;   //!

   analysis(TTree *tree=0);
   virtual ~analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_cxx
analysis::analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   /*if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("inputValidation/user.mtanaka.5578272._000001.hist-output.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("inputValidation/user.mtanaka.5578272._000001.hist-output.root");
      }
      f->GetObject("validationT",tree);

   }*/
  TChain *chain = new TChain("validationT");
  //chain->Add("inputValidation/data/user.mtanaka.5582965._000001.hist-output.root");
  //chain->Add("/data/maxi174/zp/mtanaka/myROOTAnalysisLUT/Run/user.mtanaka.test.201506141704.00267073.f594_m1435_p2361_hist-output.root.31081490/*root");
  //const char* rec = "datalist/inputValidationGrl.list";
  //const char* rec = "datalist/inputValidationNoVias.list";
  //const char* rec = "datalist/inputValidationTagAndProbe.list";
  //const char* rec = "datalist/inputValidationMinThr.list";
  //const char* rec = "datalist/inputValidationSingle.list";
  //const char* rec = "datalist/inputValidationSingleTGC.list";
  const char* rec = "datalist/inputValidationNew.list";
  //const char* rec = "datalist/inputValidationMCAll.list";
  //const char* rec = "datalist/inputValidationMCHLTmu4.list";
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

analysis::~analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
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

void analysis::Init(TTree *tree)
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("L2_pt", &L2_pt, &b_L2_pt);
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
   fChain->SetBranchAddress("L2_saddress", &L2_saddress, &b_L2_saddress);
   fChain->SetBranchAddress("tgcPt", &tgcPt, &b_tgcPt);
   fChain->SetBranchAddress("tgcMid1_r", &tgcMid1_r, &b_tgcMid1_r);
   fChain->SetBranchAddress("tgcMid1_z", &tgcMid1_z, &b_tgcMid1_z);
   fChain->SetBranchAddress("tgcMid1_eta", &tgcMid1_eta, &b_tgcMid1_eta);
   fChain->SetBranchAddress("tgcMid1_phi", &tgcMid1_phi, &b_tgcMid1_phi);
   fChain->SetBranchAddress("tgcMid2_r", &tgcMid2_r, &b_tgcMid2_r);
   fChain->SetBranchAddress("tgcMid2_z", &tgcMid2_z, &b_tgcMid2_z);
   fChain->SetBranchAddress("offline_pt", &offline_pt, &b_offline_pt);
   Notify();
}

Bool_t analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
