#define analysis_cxx
#include "macro/analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#define PI 3.14159265258979

double const PI_OVER_4 = PI/4.0;
double const PI_OVER_8 = PI/8.0;
double const PHI_RANGE = 12.0/PI_OVER_8;
double const ZERO_LIMIT = 1e-5;

pair<int,int> GetBinNumber(float m_tgcMid1_phi, float m_tgcMid1_eta){

  int Octant = (int)(m_tgcMid1_phi/PI_OVER_4);
  double PhiInOctant = fabs(m_tgcMid1_phi - Octant*PI_OVER_4);
  if(PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  int phiBin = static_cast<int>(PhiInOctant*PHI_RANGE);
  int etaBin = static_cast<int>((fabs(m_tgcMid1_eta)-1.)/0.05);

  if (etaBin == -1) etaBin =  0;
  if (etaBin == 30) etaBin = 29;

  if (etaBin < -0.5 || etaBin > 29.5 || phiBin < -0.5 || phiBin > 11.5) return make_pair(-1,-1);

  return make_pair(etaBin,phiBin);
}
///////////////////////

pair<int,int> GetBinNumberEE(float phi){
  int SL=-1;//0:small 1:large
  if (fabs(phi)>0.196 && fabs(phi)<0.589) SL=0;
  if (fabs(phi)>0.982 && fabs(phi)<1.375) SL=0;
  if (fabs(phi)>1.767 && fabs(phi)<2.16) SL=0;
  if (fabs(phi)>2.553 && fabs(phi)<2.945) SL=0;
  if (fabs(phi)<0.196) SL=1;
  if (fabs(phi)>0.589 && fabs(phi)<0.982) SL=1;
  if (fabs(phi)>1.375 && fabs(phi)<1.767) SL=1;
  if (fabs(phi)>2.16 && fabs(phi)<2.553) SL=1;
  if (fabs(phi)>2.945) SL=1;

  int phiBin24=-1;
  int Octant = (int)(phi/PI_OVER_4);
  double PhiInOctant = fabs(phi - Octant * PI_OVER_4);
  if (PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  if ( SL==0 ){//Small
    int OctantSmall = Octant;
    double PhiInOctantSmall = PhiInOctant;
    if(phi<0) PhiInOctantSmall = fabs(phi - (OctantSmall-1)*PI_OVER_4);
    phiBin24 = PhiInOctantSmall * PHI_RANGE;
  }
  else {//Large
    phi = phi + PI_OVER_8;
    int OctantLarge = (int)(phi / PI_OVER_4);
    double PhiInOctantLarge = fabs(phi - OctantLarge * PI_OVER_4);
    if (phi<0) PhiInOctantLarge = fabs(phi - (OctantLarge-1)*PI_OVER_4);
    phiBin24 = PhiInOctantLarge * PHI_RANGE;
  }

  return make_pair(phiBin24,SL);
}

void analysis::Loop()
{
  TFile *fout = new TFile("resultAnalysis/data15_13TeV_final.root","recreate");

  //histgram
  ////////////////////////////////////
  TH2 *alphaVsL2Phi = TH2F("alphaVsL2Phi",";;")


  /////////////////////////////////////
  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000000==0) cout << "entry=" << jentry << endl;
    //if (jentry<9000000) continue;
    //cout << "entry=" << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry); 


  }

  fout->Write();
}
