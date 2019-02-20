#define CheckEff_cxx
#include "CheckEff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <iostream>
#include <sstream>
using namespace std;

void CheckEff::Loop()
{
  //TFile *file = new TFile("outputCheckEff/efficiency_zmumu_r7447.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_zmumu_r7463.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_zmumu_r7514.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_zmumu_r7534.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_noTandP_r7507.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_noTandP_r7540.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_zmumu_use_middle_eta_default.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_jpsimu4mu4_use_middle_eta_default.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_jpsimu2p5mu15.root","recreate");
  //TFile *file = new TFile("outputCheckEff/efficiency_jpsi_thonda_default2.root","recreate");
  TFile *file = new TFile("outputCheckEff/efficiency_jpsi_thonda_modified2.root","recreate");
  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;

  //////////////////////////initialization
  int nbinPt=30,divbinPt=2,ptMax=60,nChain=6;
  int nbinEta=50;
  int Thr[6]={4,6,10,14,18,20};
  int numPt[nChain][5][2][nbinPt];//chain,level,region,ptbin
  int numEta[nChain][5][nbinEta];//chain,level,region,ptbin
  float xbinPt[nbinPt],ybinPt[nbinPt],xbinerrPt[nbinPt],ybinerrPt[nbinPt];
  float xbinEta[nbinEta],ybinEta[nbinEta],xbinerrEta[nbinEta],ybinerrEta[nbinEta];
  float effPt[nChain][5][2][nbinPt],efferrPt[nChain][5][2][nbinPt]; 
  float effEta[nChain][5][nbinEta],efferrEta[nChain][5][nbinEta]; 
  for (int ipt=0;ipt<nbinPt;ipt++){
    ybinPt[ipt]=0;ybinerrPt[ipt]=0;
    xbinPt[ipt]=1.0*divbinPt/2+ipt*divbinPt;
    xbinerrPt[ipt]=1.0*divbinPt/2;
    for (int ch=0;ch<nChain;ch++){
      for (int lv=0;lv<5;lv++){
        for (int eb=0;eb<2;eb++){ 
          numPt[ch][lv][eb][ipt]=0;
          effPt[ch][lv][eb][ipt]=0;
          efferrPt[ch][lv][eb][ipt]=0;
        }
      }
    }
  }
  for (int ieta=0;ieta<nbinEta;ieta++){
    ybinEta[ieta]=0;ybinerrEta[ieta]=0;
    xbinEta[ieta]=-2.45+0.1*ieta;
    xbinerrEta[ieta]=0.05;
    for (int ch=0;ch<nChain;ch++){
      for (int lv=0;lv<5;lv++){
        numEta[ch][lv][ieta]=0;
        effEta[ch][lv][ieta]=0;
        efferrEta[ch][lv][ieta]=0;
      }
    }
  }
  /////////////////////////

  ////////////////////////Event loop
  for (int jentry=0; jentry<nentries;jentry++) {
    //if (jentry>0) break;
    fChain->GetEntry(jentry);  
    float offpt = probe_offline_pt;
    float offeta = probe_offline_eta;
    float offphi = probe_offline_phi;
    bool passL1[nChain],passSA[nChain],passComb[nChain],passEF[nChain];
    int threshold = probe_threshold;
    for (int ch=0;ch<nChain;ch++){
      if(threshold==Thr[ch]){
        passL1[ch] = probe_pass_L1;
        passSA[ch] = probe_pass_SA;
        passComb[ch] = probe_pass_Comb;
        passEF[ch] = probe_pass_EF;
        //cout << "threshold=" << threshold << endl;
        if (fabs(offeta)<1.05){//barrel
          for (int ipt=0;ipt<nbinPt;ipt++){
            if (offpt>ipt*divbinPt && offpt<(ipt+1)*divbinPt){
              numPt[ch][0][0][ipt]++;
              if (passL1[ch]) numPt[ch][1][0][ipt]++;
              if (passSA[ch]) numPt[ch][2][0][ipt]++;
              if (passComb[ch]) numPt[ch][3][0][ipt]++;
              if (passEF[ch]) numPt[ch][4][0][ipt]++;
            }
          }
        }
        else {//endcap
          for (int ipt=0;ipt<nbinPt;ipt++){
            if (offpt>ipt*divbinPt && offpt<(ipt+1)*divbinPt){
              numPt[ch][0][1][ipt]++;
              if (passL1[ch]) numPt[ch][1][1][ipt]++;
              if (passSA[ch]) numPt[ch][2][1][ipt]++;
              if (passComb[ch]) numPt[ch][3][1][ipt]++;
              if (passEF[ch]) numPt[ch][4][1][ipt]++;
            }
          }
        }

        if (offpt>Thr[ch]){
          for (int ieta=0;ieta<nbinEta;ieta++){
            //if (offeta>-2.5+0.1*ieta && offeta<-2.4+0.1*ieta){
            if (probe_roi_eta>-2.5+0.1*ieta && probe_roi_eta<-2.4+0.1*ieta){
              numEta[ch][0][ieta]++;
              if (passL1[ch]) numEta[ch][1][ieta]++;
              if (passSA[ch]) numEta[ch][2][ieta]++;
              if (passComb[ch]) numEta[ch][3][ieta]++;
              if (passEF[ch]) numEta[ch][4][ieta]++;
            }
          }
        }

      }
    }
  }
  /////////////////////////////

  /////////////////////////Calculate efficiency
  TGraphErrors* EffPt[nChain][4][2];
  string LEVEL[4] = {"L1","SA","Comb","EF"};
  string EB[2] = {"Barrel","Endcap"};
  stringstream effName; 
  for (int ch=0;ch<nChain;ch++){
    for (int level=0;level<4;level++){//Trigger level
      for (int eb=0;eb<2;eb++){//barrel or endcap
        for (int ipt=0;ipt<nbinPt;ipt++){
          //cout << "numPt level/eb/ipt=" << numPt[level][eb][ipt] << endl;
          if (numPt[ch][level][eb][ipt]>0){ 
            effPt[ch][level][eb][ipt] = 
              1.*numPt[ch][level+1][eb][ipt]/numPt[ch][level][eb][ipt];
            efferrPt[ch][level][eb][ipt] = 
              sqrt(effPt[ch][level][eb][ipt]*(1-effPt[ch][level][eb][ipt]) / numPt[ch][level][eb][ipt]);
          }
        }
        effName << "Efficiency" << LEVEL[level] << "Mu" << Thr[ch] << EB[eb] ;
        EffPt[ch][level][eb] = new TGraphErrors(nbinPt,xbinPt,effPt[ch][level][eb],xbinerrPt,efferrPt[ch][level][eb]);
        EffPt[ch][level][eb]->SetName(effName.str().c_str());
        EffPt[ch][level][eb]->SetTitle("");
        EffPt[ch][level][eb]->GetXaxis()->SetLabelSize(0.07);
        EffPt[ch][level][eb]->GetYaxis()->SetLabelSize(0.07);
        EffPt[ch][level][eb]->GetYaxis()->SetRangeUser(0,1.1);
        file->Add(EffPt[ch][level][eb]);
        effName.str("");
      }
    }
  }

  TGraphErrors* EffEta[nChain][4];
  for (int ch=0;ch<nChain;ch++){
    for (int level=0;level<4;level++){//Trigger level
      for (int ieta=0;ieta<nbinEta;ieta++){
        if (numEta[ch][level][ieta]>0){ 
          effEta[ch][level][ieta] = 
            1.*numEta[ch][level+1][ieta]/numEta[ch][level][ieta];
          efferrEta[ch][level][ieta] = 
            sqrt(effEta[ch][level][ieta]*(1-effEta[ch][level][ieta]) / numEta[ch][level][ieta]);
        }
      }
      effName << "EfficiencyEta" << LEVEL[level] << "Mu" << Thr[ch] ;
      EffEta[ch][level] = new TGraphErrors(nbinEta,xbinEta,effEta[ch][level],xbinerrEta,efferrEta[ch][level]);
      EffEta[ch][level]->SetName(effName.str().c_str());
      EffEta[ch][level]->SetTitle("");
      EffEta[ch][level]->GetXaxis()->SetLabelSize(0.07);
      EffEta[ch][level]->GetYaxis()->SetLabelSize(0.07);
      EffEta[ch][level]->GetYaxis()->SetRangeUser(0,1.1);
      file->Add(EffEta[ch][level]);
      effName.str("");
    }
  }

  file->Write();
}
