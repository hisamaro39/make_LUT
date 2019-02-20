#define validEfficiency_cxx
#include "validEfficiency.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <iostream>
using namespace std;

void validEfficiency::Loop()
{
  TFile *file = new TFile("outputEfficiency/no_dr.root","recreate");

  TH2 *drL1ExtVsPt = new TH2F("drL1ExtVsPt",";offline pt(GeV);dr",100,0,50,100,0,1);

  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;

  for (int jentry=0; jentry<nentries;jentry++) {
    fChain->GetEntry(jentry);   
    if(jentry%1000000==0) cout << "Events=" << jentry << endl;
    float deta = roi_eta - offline_ext_eta;
    float dphi = acos( cos( roi_phi - offline_ext_phi));
    float dr = sqrt(deta*deta+dphi*dphi);
    //cout << "dr/pt=" << dr << "/" << offline_pt << endl;
    drL1ExtVsPt->Fill(offline_pt,dr);
  }
  
  file->Write();

}
