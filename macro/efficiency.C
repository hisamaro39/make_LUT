#define efficiency_cxx
#include "efficiency.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

void efficiency::Loop()
{

  int nentries = fChain->GetEntries();

  //TFile *file = new TFile("outputEfficiency/efficiency_online.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_eeRegion.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_newLUT.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_newLUT_final2.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_new_dr.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_new_dr2.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_plus_ee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_newLUT_new_dr2.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_new_dr3.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_new_dr4.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_defaultLUT_jpsi.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_newLUT_EENoSL2.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_combine_vesion3.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_combine_ee_small_large.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_combine_ee_small.root","recreate");
  TFile *file = new TFile("outputEfficiency/efficiency_combine_version4.root","recreate");

  int numProbe_ec_pt[60],numL1_ec_pt[2][60],numSA_ec_pt[2][60],numSA_ec_pt_new[2][60];
  int numProbe_br_pt[60],numL1_br_pt[2][60],numSA_br_pt[2][60],numSA_br_pt_new[2][60];
  int numL1_etaphi[100][64],numSA_etaphi[100][64];
  int numProbe_ec_phi[32],numL1_ec_phi[2][32],numSA_ec_phi[2][32];
  int numProbe_br_phi[32],numL1_br_phi[2][32],numSA_br_phi[2][32];
  int numProbe_eta[50],numL1_eta[2][50],numSA_eta[2][50];
  int numSA_ec_charge_pt[40],numSA_br_charge_pt[40];
  int numSA_ec_charge_correct_pt[40],numSA_br_charge_correct_pt[40];
  int numSA_ec_charge_correct_phi[32],numSA_br_charge_correct_phi[32];
  int numSA_charge_correct_eta[50];
  for (int i=0;i<100;i++){
    for (int j=0;j<64;j++){
      numL1_etaphi[i][j]=0;
      numSA_etaphi[i][j]=0;
    }
  }
  for (int k=0;k<60;k++){
    numProbe_ec_pt[k]=0;
    numProbe_br_pt[k]=0;
    numSA_ec_charge_pt[k]=0;
    numSA_br_charge_pt[k]=0;
    for (int o=0;o<2;o++){
      numL1_ec_pt[o][k]=0;
      numSA_ec_pt[o][k]=0;
      numL1_br_pt[o][k]=0;
      numSA_br_pt[o][k]=0;
      numSA_ec_pt_new[o][k]=0;
      numSA_br_pt_new[o][k]=0;
    }
  }
  for (int k=0;k<40;k++){
    numSA_ec_charge_correct_pt[k]=0;
    numSA_br_charge_correct_pt[k]=0;
  }

  for (int q=0;q<32;q++){
    numProbe_ec_phi[q]=0;
    numProbe_br_phi[q]=0;
    numSA_ec_charge_correct_phi[q]=0;
    numSA_br_charge_correct_phi[q]=0;
    for (int r=0;r<2;r++){
      numL1_br_phi[r][q]=0;
      numSA_br_phi[r][q]=0;
      numL1_ec_phi[r][q]=0;
      numSA_ec_phi[r][q]=0;
    }
  }

  for (int q=0;q<50;q++){
    numProbe_eta[q]=0;
    for (int r=0;r<2;r++){
      numL1_eta[r][q]=0;
      numSA_eta[r][q]=0;
    }
  }

  cout << "Number of Events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    fChain->GetEntry(jentry);   
    float offline_pt = fabs(probe_offline_pt);
    float offline_eta = probe_offline_eta;
    float offline_phi = probe_offline_phi;
    float sa_pt = fabs(probe_sa_pt);
    int saddress = probe_sa_saddress;
    bool pass_L1mu4 = probe_passL1_1;
    bool pass_SAmu4 = probe_passSA_1;
    bool pass_SAmu4_new = probe_passSA_1_new;
    bool pass_L1mu6 = probe_passL1_2;
    bool pass_SAmu6 = probe_passSA_2;
    bool pass_SAmu6_new = probe_passSA_2_new;
    bool charge_correct = (probe_offline_charge==probe_sa_charge)? true : false;
    if (pass_SAmu4){
      //cout << "charge offline/SA=" << probe_offline_charge << "/" << probe_sa_charge << endl;
      //cout << "offline & SA charge match is " << charge_correct << endl;
    }
    //cout << "offline_pt=" << offline_pt << endl;

    for (int i=0;i<50;i++){
      if (offline_eta>-2.5+0.1*i && offline_eta<-2.4+0.1*i){
        numProbe_eta[i]++;
        if (pass_L1mu4) numL1_eta[0][i]++;
        if (pass_SAmu4) {
          numSA_eta[0][i]++;
          if(charge_correct) numSA_charge_correct_eta[i]++;
        }
        if (pass_L1mu6) numL1_eta[1][i]++;
        if (pass_SAmu6) numSA_eta[1][i]++;
      }
    }
    
    //if (saddress<0){//Endcap
    //if (fabs(offline_eta)>1.05 && fabs(offline_eta)<1.35){//EE region
    //if (fabs(offline_eta)>1.05){//Endcap
    if (fabs(offline_eta)>1.05 && fabs(offline_eta)<1.35){//EE

      //cout << "endcap pt resolution is " << pt_ec_reso << endl;
      if (fabs(offline_eta)>1.05 && fabs(offline_eta)<1.35)
      for (int i=0;i<60;i++){
        if (offline_pt>0.25*i && offline_pt<0.25*(i+1)){
          numProbe_ec_pt[i]++;
          if (pass_L1mu4) numL1_ec_pt[0][i]++;
          if (pass_SAmu4) numSA_ec_pt[0][i]++;
          if (pass_SAmu4_new) numSA_ec_pt_new[0][i]++;
          if (pass_L1mu6) numL1_ec_pt[1][i]++;
          if (pass_SAmu6) numSA_ec_pt[1][i]++;
          if (pass_SAmu6_new) numSA_ec_pt_new[1][i]++;
        }
      }
      for (int i=0;i<100;i++){
        for (int j=0;j<64;j++){
          if (offline_eta>0.05*i-2.5 && offline_eta<0.05*i-2.45){
            if (offline_phi>0.1*i-3.2 && offline_phi<0.1*i-3.1){
              if (pass_L1mu6) numL1_etaphi[i][j]++;
              if (pass_SAmu6) numSA_etaphi[i][j]++;
            }
          }
        }
      }
      for (int i=0;i<40;i++){
        if (offline_pt>i*0.5 && offline_pt<0.5*(i+1)){
          if (pass_SAmu4) {
            numSA_ec_charge_pt[i]++;
            if (charge_correct) numSA_ec_charge_correct_pt[i]++;
          }
        }
      }
      for (int s=0;s<32;s++){
        if (offline_phi>-3.2+0.2*s && offline_phi<-3+0.2*s){
          numProbe_ec_phi[s]++;
          if (pass_L1mu4) numL1_ec_phi[0][s]++;
          if (pass_SAmu4) {
            numSA_ec_phi[0][s]++;
            if (charge_correct) numSA_ec_charge_correct_phi[s]++;
          }
          if (pass_L1mu6) numL1_ec_phi[1][s]++;
          if (pass_SAmu6) numSA_ec_phi[1][s]++;
        }
      }
    }
    else  {//Barrel
      for (int i=0;i<60;i++){
        if (offline_pt>0.25*i && offline_pt<0.25*(i+1)){
          numProbe_br_pt[i]++;
          if (pass_L1mu4) numL1_br_pt[0][i]++;
          if (pass_SAmu4) numSA_br_pt[0][i]++;
          if (pass_SAmu4_new) numSA_br_pt_new[0][i]++;
          if (pass_L1mu6) numL1_br_pt[1][i]++;
          if (pass_SAmu6) numSA_br_pt[1][i]++;
          if (pass_SAmu6_new) numSA_br_pt_new[1][i]++;
        }
      }
      for (int i=0;i<40;i++){
        if (offline_pt>i && offline_pt<i+1){
          if (pass_SAmu4) {
            numSA_br_charge_pt[i]++;
            if (charge_correct) numSA_br_charge_correct_pt[i]++;
          }
        }
      }
      for (int s=0;s<32;s++){
        if (offline_phi>-3.2+0.2*s && offline_phi<-3+0.2*s){
          numProbe_br_phi[s]++;
          if (pass_L1mu4) numL1_br_phi[0][s]++;
          if (pass_SAmu4) {
            numSA_br_phi[0][s]++;
            if (charge_correct) numSA_br_charge_correct_phi[s]++;
          }
          if (pass_L1mu6) numL1_br_phi[1][s]++;
          if (pass_SAmu6) numSA_br_phi[1][s]++;
        }
      }

    }
  }

  float eff_ec_pt_L1[2][60],efferr_ec_pt_L1[2][60],eff_ec_pt_SA[2][60],efferr_ec_pt_SA[2][60],eff_ec_pt_SA_new[2][60],efferr_ec_pt_SA_new[2][60];
  float eff_br_pt_L1[2][60],efferr_br_pt_L1[2][60],eff_br_pt_SA[2][60],efferr_br_pt_SA[2][60],eff_br_pt_SA_new[2][60],efferr_br_pt_SA_new[2][60];
  float ptbin[2][60],ptbinerr[2][60];
  for (int p=0;p<2;p++){
    cout << "p=" << p << endl;
    for (int m=0;m<60;m++){
      cout << "m=" << m << endl;
      //endcap
      //default LUT L1 SA
      eff_ec_pt_L1[p][m] = (numProbe_ec_pt[m]) ? 
        1.*numL1_ec_pt[p][m]/numProbe_ec_pt[m] : 0 ;
      efferr_ec_pt_L1[p][m]=sqrt( eff_ec_pt_L1[p][m]*(1-eff_ec_pt_L1[p][m]) )/numProbe_ec_pt[m] ;
      eff_ec_pt_SA[p][m]= (numL1_ec_pt[p][m]) ?
        1.*numSA_ec_pt[p][m]/numL1_ec_pt[p][m] : 0;
      efferr_ec_pt_SA[p][m]=sqrt( eff_ec_pt_SA[p][m]*(1-eff_ec_pt_SA[p][m]) )/numL1_ec_pt[p][m] ;
      //new LUT SA
      eff_ec_pt_SA_new[p][m]= (numL1_ec_pt[p][m]) ?
        1.*numSA_ec_pt_new[p][m]/numL1_ec_pt[p][m] : 0;
      efferr_ec_pt_SA_new[p][m]=sqrt( eff_ec_pt_SA_new[p][m]*(1-eff_ec_pt_SA_new[p][m]) )/numL1_ec_pt[p][m] ;
      //barrel
      //default LUT L1 SA
      eff_br_pt_L1[p][m] = (numProbe_br_pt[m]) ? 
        1.*numL1_br_pt[p][m]/numProbe_br_pt[m] : 0 ;
      efferr_br_pt_L1[p][m]=sqrt( eff_br_pt_L1[p][m]*(1-eff_br_pt_L1[p][m]) )/numProbe_br_pt[m] ;
      eff_br_pt_SA[p][m]= (numL1_br_pt[p][m]) ?
        1.*numSA_br_pt[p][m]/numL1_br_pt[p][m] : 0;
      efferr_br_pt_SA[p][m]=sqrt( eff_br_pt_SA[p][m]*(1-eff_br_pt_SA[p][m]) )/numL1_br_pt[p][m] ;
      //new LUT SA
      eff_br_pt_SA_new[p][m]= (numL1_br_pt[p][m]) ?
        1.*numSA_br_pt_new[p][m]/numL1_br_pt[p][m] : 0;
      efferr_br_pt_SA_new[p][m]=sqrt( eff_br_pt_SA_new[p][m]*(1-eff_br_pt_SA_new[p][m]) )/numL1_br_pt[p][m] ;
      ptbin[p][m]=0.125+0.25*m;
      ptbinerr[p][m] = 0.125;
      if (p==0) cout << "HLT_mu4" << endl;
      else cout << "HLT_mu6" << endl;
      cout << "bin" << m << ": Probe/L1/SA=" 
        << numProbe_ec_pt[m] << "/" << numL1_ec_pt[p][m] << "/" << numSA_ec_pt[p][m] << endl;
      cout << "L1 efficiency is " << eff_ec_pt_L1[p][m] << " +- " << efferr_ec_pt_L1[p][m] << endl;
      cout << "SA efficiency is " << eff_ec_pt_SA[p][m] << " +- " << efferr_ec_pt_SA[p][m] << endl;
      cout << "SA efficiency new is " << eff_ec_pt_SA_new[p][m] << " +- " << efferr_ec_pt_SA_new[p][m] << endl;
    }
  }

  float eff_ec_phi_L1[2][32],efferr_ec_phi_L1[2][32],eff_ec_phi_SA[2][32],efferr_ec_phi_SA[2][32];
  float eff_br_phi_L1[2][32],efferr_br_phi_L1[2][32],eff_br_phi_SA[2][32],efferr_br_phi_SA[2][32];
  float phibin[2][32],phibinerr[2][32];
  for (int p=0;p<2;p++){
    for (int m=0;m<32;m++){
      //endcap
      eff_ec_phi_L1[p][m] = (numProbe_ec_phi[m]) ? 
        1.*numL1_ec_phi[p][m]/numProbe_ec_phi[m] : 0 ;
      efferr_ec_phi_L1[p][m]=sqrt( eff_ec_phi_L1[p][m]*(1-eff_ec_phi_L1[p][m]) )/numProbe_ec_phi[m] ;
      eff_ec_phi_SA[p][m]= (numL1_ec_phi[p][m]) ?
        1.*numSA_ec_phi[p][m]/numL1_ec_phi[p][m] : 0;
      efferr_ec_phi_SA[p][m]=sqrt( eff_ec_phi_SA[p][m]*(1-eff_ec_phi_SA[p][m]) )/numL1_ec_phi[p][m] ;
      phibin[p][m]=-3.1+0.2*m;
      phibinerr[p][m] = 0.1;
      //barrel
      eff_br_phi_L1[p][m] = (numProbe_br_phi[m]) ? 
        1.*numL1_br_phi[p][m]/numProbe_br_phi[m] : 0 ;
      efferr_br_phi_L1[p][m]=sqrt( eff_br_phi_L1[p][m]*(1-eff_br_phi_L1[p][m]) )/numProbe_br_phi[m] ;
      eff_br_phi_SA[p][m]= (numL1_br_phi[p][m]) ?
        1.*numSA_br_phi[p][m]/numL1_br_phi[p][m] : 0;
      efferr_br_phi_SA[p][m]=sqrt( eff_br_phi_SA[p][m]*(1-eff_br_phi_SA[p][m]) )/numL1_br_phi[p][m] ;
      phibin[p][m]=-3.1+0.2*m;
      phibinerr[p][m] = 0.1;
    }
  }

  float eff_eta_L1[2][50],efferr_eta_L1[2][50],eff_eta_SA[2][50],efferr_eta_SA[2][50];
  float etabin[2][50],etabinerr[2][50];
  for (int p=0;p<2;p++){
    for (int m=0;m<50;m++){
      //endcap
      eff_eta_L1[p][m] = (numProbe_eta[m]) ? 
        1.*numL1_eta[p][m]/numProbe_eta[m] : 0 ;
      efferr_eta_L1[p][m]=sqrt( eff_eta_L1[p][m]*(1-eff_eta_L1[p][m]) )/numProbe_eta[m] ;
      eff_eta_SA[p][m]= (numL1_eta[p][m]) ?
        1.*numSA_eta[p][m]/numL1_eta[p][m] : 0;
      efferr_eta_SA[p][m]=sqrt( eff_eta_SA[p][m]*(1-eff_eta_SA[p][m]) )/numL1_eta[p][m] ;
      etabin[p][m]=-2.45+0.1*m;
      etabinerr[p][m] = 0.05;
    }
  }

  float eff_ec_charge_correct_pt[40],eff_br_charge_correct_pt[40],eff_ec_charge_correct_phi[32],eff_br_charge_correct_phi[32],eff_charge_correct_eta[32];
  float efferr_ec_charge_correct_pt[40],efferr_br_charge_correct_pt[40],efferr_ec_charge_correct_phi[32],efferr_br_charge_correct_phi[32],efferr_charge_correct_eta[32];
  float ptbin2[2][40],ptbinerr2[2][40];
  for (int m=0;m<40;m++){
    cout << "m=" << m << endl;
    //endcap
    eff_ec_charge_correct_pt[m]= (numSA_ec_charge_pt[m]) ?
      1.*numSA_ec_charge_correct_pt[m]/numSA_ec_charge_pt[m] : 0;
    efferr_ec_charge_correct_pt[m]=sqrt( eff_ec_charge_correct_pt[m]*(1-eff_ec_charge_correct_pt[m]) )/numSA_ec_charge_pt[m] ;
    //barrel
    eff_br_charge_correct_pt[m]= (numSA_br_charge_pt[m]) ?
      1.*numSA_br_charge_correct_pt[m]/numSA_br_charge_pt[m] : 0;
    efferr_br_charge_correct_pt[m]=sqrt( eff_br_charge_correct_pt[m]*(1-eff_br_charge_correct_pt[m]) )/numSA_br_charge_pt[m] ;
    ptbin2[0][m]=0.5+m;
    ptbinerr2[0][m] = 0.5;
    cout << "bin" << m << ": Endcap Probe/SA/charge_correct=" 
      << numProbe_ec_pt[m] << "/" << numSA_ec_pt[0][m] << "/" << numSA_ec_charge_correct_pt[m] << endl;
    cout << "SA Endcap charge efficiency is " << eff_ec_charge_correct_pt[m] << " +- " << efferr_ec_charge_correct_pt[m] << endl;
    cout << "bin" << m << ": Barrel Probe/SA/charge_correct=" 
      << numProbe_br_pt[m] << "/" << numSA_br_pt[0][m] << "/" << numSA_br_charge_correct_pt[m] << endl;
    cout << "SA Barrel charge efficiency is " << eff_br_charge_correct_pt[m] << " +- " << efferr_br_charge_correct_pt[m] << endl;
  }
  
  //endcap mu4 pt
  TGraphErrors *effEndcapL1mu4 = new TGraphErrors(60,ptbin[0],eff_ec_pt_L1[0],ptbinerr[0],efferr_ec_pt_L1[0]);
  effEndcapL1mu4->SetName("effEndcapL1mu4");
  effEndcapL1mu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapL1mu4);
  TGraphErrors *effEndcapSAmu4 = new TGraphErrors(60,ptbin[0],eff_ec_pt_SA[0],ptbinerr[0],efferr_ec_pt_SA[0]);
  effEndcapSAmu4->SetName("effEndcapSAmu4");
  effEndcapSAmu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu4);
  TGraphErrors *effEndcapSAmu4New = new TGraphErrors(60,ptbin[0],eff_ec_pt_SA_new[0],ptbinerr[0],efferr_ec_pt_SA_new[0]);
  effEndcapSAmu4New->SetName("effEndcapSAmu4New");
  effEndcapSAmu4New->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu4New);
  //endcap mu6 pt
  TGraphErrors *effEndcapL1mu6 = new TGraphErrors(60,ptbin[1],eff_ec_pt_L1[1],ptbinerr[1],efferr_ec_pt_L1[1]);
  effEndcapL1mu6->SetName("effEndcapL1mu6");
  effEndcapL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapL1mu6);
  TGraphErrors *effEndcapSAmu6 = new TGraphErrors(60,ptbin[1],eff_ec_pt_SA[1],ptbinerr[1],efferr_ec_pt_SA[1]);
  effEndcapSAmu6->SetName("effEndcapSAmu6");
  effEndcapSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu6);
  TGraphErrors *effEndcapSAmu6New = new TGraphErrors(60,ptbin[1],eff_ec_pt_SA_new[1],ptbinerr[1],efferr_ec_pt_SA_new[1]);
  effEndcapSAmu6New->SetName("effEndcapSAmu6New");
  effEndcapSAmu6New->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu6New);
  //barrel mu4 pt
  TGraphErrors *effBarrelL1mu4 = new TGraphErrors(60,ptbin[0],eff_br_pt_L1[0],ptbinerr[0],efferr_br_pt_L1[0]);
  effBarrelL1mu4->SetName("effBarrelL1mu4");
  effBarrelL1mu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelL1mu4);
  TGraphErrors *effBarrelSAmu4 = new TGraphErrors(60,ptbin[0],eff_br_pt_SA[0],ptbinerr[0],efferr_br_pt_SA[0]);
  effBarrelSAmu4->SetName("effBarrelSAmu4");
  effBarrelSAmu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu4);
  TGraphErrors *effBarrelSAmu4New = new TGraphErrors(60,ptbin[0],eff_br_pt_SA_new[0],ptbinerr[0],efferr_br_pt_SA_new[0]);
  effBarrelSAmu4New->SetName("effBarrelSAmu4New");
  effBarrelSAmu4New->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu4New);
  //barrel mu6 pt
  TGraphErrors *effBarrelL1mu6 = new TGraphErrors(60,ptbin[1],eff_br_pt_L1[1],ptbinerr[1],efferr_br_pt_L1[1]);
  effBarrelL1mu6->SetName("effBarrelL1mu6");
  effBarrelL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelL1mu6);
  TGraphErrors *effBarrelSAmu6 = new TGraphErrors(60,ptbin[1],eff_br_pt_SA[1],ptbinerr[1],efferr_br_pt_SA[1]);
  effBarrelSAmu6->SetName("effBarrelSAmu6");
  effBarrelSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu6);
  TGraphErrors *effBarrelSAmu6New = new TGraphErrors(60,ptbin[1],eff_br_pt_SA_new[1],ptbinerr[1],efferr_br_pt_SA_new[1]);
  effBarrelSAmu6New->SetName("effBarrelSAmu6New");
  effBarrelSAmu6New->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu6New);
  //endcap mu4 phi
  TGraphErrors *effPhiEndcapL1mu4 = new TGraphErrors(32,phibin[0],eff_ec_phi_L1[0],phibinerr[0],efferr_ec_phi_L1[0]);
  effPhiEndcapL1mu4->SetName("effPhiEndcapL1mu4");
  effPhiEndcapL1mu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapL1mu4);
  TGraphErrors *effPhiEndcapSAmu4 = new TGraphErrors(32,phibin[0],eff_ec_phi_SA[0],phibinerr[0],efferr_ec_phi_SA[0]);
  effPhiEndcapSAmu4->SetName("effPhiEndcapSAmu4");
  effPhiEndcapSAmu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapSAmu4);
  //endcap mu6 phi
  TGraphErrors *effPhiEndcapL1mu6 = new TGraphErrors(32,phibin[1],eff_ec_phi_L1[1],phibinerr[1],efferr_ec_phi_L1[1]);
  effPhiEndcapL1mu6->SetName("effPhiEndcapL1mu6");
  effPhiEndcapL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapL1mu6);
  TGraphErrors *effPhiEndcapSAmu6 = new TGraphErrors(32,phibin[1],eff_ec_phi_SA[1],phibinerr[1],efferr_ec_phi_SA[1]);
  effPhiEndcapSAmu6->SetName("effPhiEndcapSAmu6");
  effPhiEndcapSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapSAmu6);
  //barrel mu4 phi
  TGraphErrors *effPhiBarrelL1mu4 = new TGraphErrors(32,phibin[0],eff_br_phi_L1[0],phibinerr[0],efferr_br_phi_L1[0]);
  effPhiBarrelL1mu4->SetName("effPhiBarrelL1mu4");
  effPhiBarrelL1mu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelL1mu4);
  TGraphErrors *effPhiBarrelSAmu4 = new TGraphErrors(32,phibin[0],eff_br_phi_SA[0],phibinerr[0],efferr_br_phi_SA[0]);
  effPhiBarrelSAmu4->SetName("effPhiBarrelSAmu4");
  effPhiBarrelSAmu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelSAmu4);
  //barrel mu6 phi
  TGraphErrors *effPhiBarrelL1mu6 = new TGraphErrors(32,phibin[1],eff_br_phi_L1[1],phibinerr[1],efferr_br_phi_L1[1]);
  effPhiBarrelL1mu6->SetName("effPhiBarrelL1mu6");
  effPhiBarrelL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelL1mu6);
  TGraphErrors *effPhiBarrelSAmu6 = new TGraphErrors(32,phibin[1],eff_br_phi_SA[1],phibinerr[1],efferr_br_phi_SA[1]);
  effPhiBarrelSAmu6->SetName("effPhiBarrelSAmu6");
  effPhiBarrelSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelSAmu6);
  //mu4 eta
  TGraphErrors *effEtaL1mu4 = new TGraphErrors(50,etabin[0],eff_eta_L1[0],etabinerr[0],efferr_eta_L1[0]);
  effEtaL1mu4->SetName("effEtaL1mu4");
  effEtaL1mu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEtaL1mu4);
  TGraphErrors *effEtaSAmu4 = new TGraphErrors(50,etabin[0],eff_eta_SA[0],etabinerr[0],efferr_eta_SA[0]);
  effEtaSAmu4->SetName("effEtaSAmu4");
  effEtaSAmu4->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEtaSAmu4);
  //mu6 eta
  TGraphErrors *effEtaL1mu6 = new TGraphErrors(50,etabin[1],eff_eta_L1[1],etabinerr[1],efferr_eta_L1[1]);
  effEtaL1mu6->SetName("effEtaL1mu6");
  effEtaL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEtaL1mu6);
  TGraphErrors *effEtaSAmu6 = new TGraphErrors(50,etabin[1],eff_eta_SA[1],etabinerr[1],efferr_eta_SA[1]);
  effEtaSAmu6->SetName("effEtaSAmu6");
  effEtaSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEtaSAmu6);
  //mu4 charge endcap pt
  TGraphErrors *effEndcapSAmu4Charge = new TGraphErrors(40,ptbin2[0],eff_ec_charge_correct_pt,ptbinerr2[0],efferr_ec_charge_correct_pt);
  effEndcapSAmu4Charge->SetName("effEndcapSAmu4Charge");
  effEndcapSAmu4Charge->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu4Charge);
  //mu4 charge barrel pt
  TGraphErrors *effBarrelSAmu4Charge = new TGraphErrors(40,ptbin2[0],eff_br_charge_correct_pt,ptbinerr2[0],efferr_br_charge_correct_pt);
  effBarrelSAmu4Charge->SetName("effBarrelSAmu4Charge");
  effBarrelSAmu4Charge->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu4Charge);
  
  file->Write();

  }
