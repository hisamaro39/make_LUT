#define resolution_cxx
#include "resolution.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

void resolution::Loop()
{

  int nentries = fChain->GetEntries();

  //TFile *file = new TFile("outputEfficiency/resolution_online.root","recreate");
  TFile *file = new TFile("outputEfficiency/resolution_defaultLUT.root","recreate");
  //TFile *file = new TFile("outputEfficiency/resolution_newLUT.root","recreate");
  //TFile *file = new TFile("outputEfficiency/resolution_newLUT_final2.root","recreate");
  //TFile *file = new TFile("outputEfficiency/aho.root","recreate");

  TH1 *ptEndcapResolution = new TH1F("ptEndcapResolution",";pt resolution;Events",100,-1,1);
  TH1 *ptAlphaResolution = new TH1F("ptAlphaResolution",";pt resolution;Events",100,-1,1);
  TH1 *ptBetaResolution = new TH1F("ptBetaResolution",";pt resolution;Events",100,-1,1);
  TH1 *ptTgcResolution = new TH1F("ptTgcResolution",";pt resolution;Events",100,-1,1);

  cout << "Number of Events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    fChain->GetEntry(jentry);   
    float offline_pt = fabs(probe_offline_pt);
    float offline_eta = probe_offline_eta;
    float offline_phi = probe_offline_phi;
    float sa_pt = fabs(probe_sa_pt);
    float sa_eta = probe_sa_eta;
    float sa_pt_alpha = fabs(probe_sa_pt_alpha);
    float sa_pt_beta = fabs(probe_sa_pt_beta);
    float sa_pt_tgc = fabs(probe_sa_pt_tgc);
    int saddress = probe_sa_saddress;
    int etaBinEC = probe_sa_eta_bin;
    int phiBinEC = probe_sa_phi_bin;
    //cout << "etabin/phibin=" << etaBinEC << "/" << phiBinEC << endl;

    //if (offline_pt>6) continue;
    if (fabs(offline_eta)>1.05){//Endcap
      //cout << "pt offline/sa=" << offline_pt << "/" << sa_pt << endl;
      float pt_ec_reso = 1 - offline_pt/sa_pt;
      float pt_alpha_reso = 1 - offline_pt/sa_pt_alpha;
      float pt_beta_reso = 1 - offline_pt/sa_pt_beta;
      float pt_tgc_reso = 1 - offline_pt/sa_pt_tgc;
     
      if ( fabs(sa_eta) > 2.3 && fabs(sa_eta)<2.35){
        ptEndcapResolution->Fill(pt_ec_reso);
        ptAlphaResolution->Fill(pt_alpha_reso);
        ptBetaResolution->Fill(pt_beta_reso);
        ptTgcResolution->Fill(pt_tgc_reso);
      }
    }
  }

  file->Write();

}
