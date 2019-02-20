#define efficiency_mc_cxx
#include "efficiency_mc.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include <sstream>

void efficiency_mc::Loop()
{

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetOptStat(0);

  int nentries = fChain->GetEntries();

  //TFile *file = new TFile("outputEfficiency/efficiency_mc_zmumu_default2.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_mc_zmumu_default_pt0.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_mc_zmumu_default_charge.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_data_zmumu.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_data_zmumu_mu6.root","recreate");
  //TFile *file = new TFile("outputEfficiency/efficiency_data_zmumu_deltaPt.root","recreate");
  TFile *file = new TFile("outputEfficiency/efficiency_data_zmumu_version4.root","recreate");
  //TFile *file = new TFile("outputEfficiency/test.root","recreate");

  TH1* saPt0barrel = new TH1F("saPt0barrel",";offline phi(GeV);Events",64,-3.2,3.2);
  TH2* offPtVsInvQptResidualEndcap[3];
  TH2* offPtVsInvQptResidualBarrel[3];
  TH2* offPtVsInvPtResidualEndcap[3];
  TH2* offPtVsInvPtResidualBarrel[3];
  string histName1[3] = {"offPtVsInvQptResidualEndcap","offPtVsInvQptResidualEndcapChargeCorrect","offPtVsInvQptResidualEndcapChargeIncorrect"};
  string histName2[3] = {"offPtVsInvQptResidualBarrel","offPtVsInvQptResidualBarrelChargeCorrect","offPtVsInvQptResidualBarrelChargeIncorrect"};
  string histName3[3] = {"offPtVsInvPtResidualEndcap","offPtVsInvPtResidualEndcapChargeCorrect","offPtVsInvPtResidualEndcapChargeIncorrect"};
  string histName4[3] = {"offPtVsInvPtResidualBarrel","offPtVsInvPtResidualBarrelChargeCorrect","offPtVsInvPtResidualBarrelChargeIncorrect"};
  for (int i=0;i<3;i++){
    offPtVsInvQptResidualEndcap[i] = new TH2F(histName1[i].c_str(),";offline pt(GeV);(1/Q*pt)_{off}-(1/Q*pt)_{SA}",100,0,100,300,-0.3,0.3);
    offPtVsInvQptResidualBarrel[i] = new TH2F(histName2[i].c_str(),";offline pt(GeV);(1/Q*pt)_{off}-(1/Q*pt)_{SA}",100,0,100,300,-0.3,0.3);
    offPtVsInvPtResidualEndcap[i] = new TH2F(histName3[i].c_str(),";offline pt(GeV);(1/pt)_{off}-(1/pt)_{SA}",100,0,100,300,-0.3,0.3);
    offPtVsInvPtResidualBarrel[i] = new TH2F(histName4[i].c_str(),";offline pt(GeV);(1/pt)_{off}-(1/pt)_{SA}",100,0,100,300,-0.3,0.3);
  }
  TH2* PtVsDeltaPtEndcap = new TH2F("PtVsDeltaPtEndcap",";offline pt(GeV);deltaPt",100,0,100,200,-10,10); 
  TH2* offPtVsPtResidualEndcap = new TH2F("offPtVsPtResidualEndcap",";offline pt(GeV);pt_{off}-pt_{SA}",100,0,100,200,-10,10);
  TH2* PtVsDeltaPtBarrel[5];
  TH2* PtVsNewDeltaPtBarrel[5];
  TH2* offPtVsPtResidualBarrel[5];
  string histName5[5] = {"PtVsDeltaPtBarrelAll","PtVsDeltaPtBarrelSmall","PtVsDeltaPtBarrelSmallSpecial","PtVsDeltaPtBarrelLarge","PtVsDeltaPtBarrelLargeSpecial"};
  string histName6[5] = {"offPtVsPtResidualBarrelAll","offPtVsPtResidualBarrelSmall","offPtVsPtResidualBarrelSmallSpecial","offPtVsPtResidualBarrelLarge","offPtVsPtResidualBarrelLargeSpecial"};
  string histName7[5] = {"PtVsNewDeltaPtBarrelAll","PtVsNewDeltaPtBarrelSmall","PtVsNewDeltaPtBarrelSmallSpecial","PtVsNewDeltaPtBarrelLarge","PtVsNewDeltaPtBarrelLargeSpecial"};
  for (int i=0;i<5;i++){
    PtVsDeltaPtBarrel[i] = new TH2F(histName5[i].c_str(),";offline pt(GeV);deltaPt",100,0,100,200,-10,10); 
    PtVsNewDeltaPtBarrel[i] = new TH2F(histName7[i].c_str(),";offline pt(GeV);new deltaPt",100,0,100,200,-10,10); 
    offPtVsPtResidualBarrel[i] = new TH2F(histName6[i].c_str(),";offline pt(GeV);pt_{off}-pt_{SA}",100,0,100,200,-10,10); 
  }
  TH2* saPtVsDeltaPtBarrel = new TH2F("saPtVsDeltaPtBarrel",";pt_{SA}(GeV);deltaPt",100,0,100,200,-10,10); 
  TH2* saPtVsDeltaPtEndcap = new TH2F("saPtVsDeltaPtEndcap",";pt_{SA}(GeV);deltaPt",100,0,100,200,-50,50); 
  TH2* saPtVsPtDifferenceEndcap = new TH2F("saPtVsPtDifferenceEndcap",";pt_{SA}(GeV);pt_{off}-pt_{SA}",100,0,100,200,-50,50);
  TH2* saPtVsPtDifferenceBarrel = new TH2F("saPtVsPtDifferenceBarrel",";pt_{SA}(GeV);pt_{off}-pt_{SA}",100,0,100,200,-50,50);
  TH2* saPtVsNewDeltaPtBarrel = new TH2F("saPtVsNewDeltaPtBarrel",";pt_{SA}(GeV);new deltaPt",100,0,100,200,-50,50); 
  TH2* saPtVsDefDeltaPtBarrel = new TH2F("saPtVsDefDeltaPtBarrel",";pt_{SA}(GeV);default deltaPt",100,0,100,200,-50,50); 
  TH2* offPtVsCalcInvPtResidualEndcap = new TH2F("offPtVsCalcInvPtResidualEndcap",";offline pt(GeV);deltaPt/(pt_{SA}*pt_{SA})",100,0,100,200,-0.01,0.01);
  TH2* offPtVsCalcInvPtResidualBarrel = new TH2F("offPtVsCalcInvPtResidualBarrel",";offline pt(GeV);deltaPt/(pt_{SA}*pt_{SA})",100,0,100,200,-0.01,0.01);
  TH2 *invSAPtVsPtDifferenceEndcap = new TH2F("invSAPtVsPtDifferenceEndcap",";1/pt_{SA};pt_{off}-pt_{SA}",100,0,0.1,100,-10,10);
  TH2 *invSAPtVsPtDifferenceBarrel = new TH2F("invSAPtVsPtDifferenceBarrel",";1/pt_{SA};pt_{off}-pt_{SA}",100,0,0.1,200,-10,10);
  TH2 *invoffPtVsPtDifferenceEndcap = new TH2F("invoffPtVsPtDifferenceEndcap",";1/pt_{off};pt_{off}-pt_{SA}",200,0,0.1,200,-10,10);
  TH2 *invoffPtVsPtDifferenceBarrel = new TH2F("invoffPtVsPtDifferenceBarrel",";1/pt_{off};pt_{off}-pt_{SA}",200,0,0.1,200,-10,10);

  int numProbe_ec_pt[40],numL1_ec_pt[40],numSA_ec_pt[40];
  int numProbe_br_pt[40],numL1_br_pt[40],numSA_br_pt[40];
  int numSA_ec_charge_pt[80],numSA_br_charge_pt[80];
  int numSA_ec_charge_correct_pt[80],numSA_br_charge_correct_pt[80];
  int numSA_ec_charge_phi[64],numSA_br_charge_phi[64];
  int numSA_ec_charge_correct_phi[64],numSA_br_charge_correct_phi[64];
  int numProbe_ec_phi[32],numL1_ec_phi[32],numSA_ec_phi[32];
  int numProbe_br_phi[32],numL1_br_phi[32],numSA_br_phi[32];
  int numProbe_eta[50],numL1_eta[50],numSA_eta[50];
  for (int k=0;k<40;k++){
    numProbe_ec_pt[k]=0;
    numProbe_br_pt[k]=0;
    numL1_ec_pt[k]=0;
    numSA_ec_pt[k]=0;
    numL1_br_pt[k]=0;
    numSA_br_pt[k]=0;
  }
  for (int k=0;k<80;k++){
    numSA_ec_charge_pt[k]=0;
    numSA_br_charge_pt[k]=0;
    numSA_ec_charge_correct_pt[k]=0;
    numSA_br_charge_correct_pt[k]=0;
  }
  for (int k=0;k<64;k++){
    numSA_ec_charge_phi[k]=0;
    numSA_br_charge_phi[k]=0;
    numSA_ec_charge_correct_phi[k]=0;
    numSA_br_charge_correct_phi[k]=0;
  }

  for (int q=0;q<32;q++){
    numProbe_ec_phi[q]=0;
    numProbe_br_phi[q]=0;
    numL1_br_phi[q]=0;
    numSA_br_phi[q]=0;
    numL1_ec_phi[q]=0;
    numSA_ec_phi[q]=0;
  }

  for (int q=0;q<50;q++){
    numProbe_eta[q]=0;
    numL1_eta[q]=0;
    numSA_eta[q]=0;
  }

  float vec[6] = {0.096348,-0.00788948,-0.726915,0.347221,0.0087001,0.001};
  float defvec[6] = {0.042, -0.00046, 3.5, -1.8, 0.35, -0.017};

  cout << "Number of Events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    fChain->GetEntry(jentry);   
    float offline_pt = fabs(probe_offline_pt);
    float offline_eta = probe_offline_eta;
    float offline_phi = probe_offline_phi;
    float sa_pt = fabs(probe_sa_pt);
    int saddress = probe_sa_saddress;
    bool pass_L1mu6 = probe_passL1mu6;
    bool pass_SAmu6 = probe_passSAmu6;
    bool charge_correct = (probe_offline_charge==probe_sa_charge)? true : false;
    float invQptResidual = 1000.;
    float invPtResidual = 1000.;
    float calcInvPtResidual = 1000.;
    float PtResidual = 1000.;
    float deltaPt = probe_sa_delta_pt;
    float newDeltaPt = 1000.;
    float defDeltaPt = 1000.;
    float invSApt = 1000.;
    float invSApt2 = 1000.;
    float invSApt3 = 1000.;
    //cout << "deltaPt=" << probe_sa_delta_pt << endl;
    if(sa_pt>1e-5) {
      invQptResidual = probe_offline_charge/offline_pt - probe_sa_charge/sa_pt;
      invPtResidual = 1./offline_pt - 1./sa_pt;
      calcInvPtResidual = deltaPt/(sa_pt*sa_pt);
      PtResidual = offline_pt-sa_pt;
      invSApt = 1./sa_pt;
      invSApt2 = invSApt*invSApt;
      invSApt3 = invSApt*invSApt*invSApt;
    }
    //cout << "invQptResidual=" << invQptResidual << endl;

    for (int i=0;i<50;i++){
      if (offline_eta>-2.5+0.1*i && offline_eta<-2.4+0.1*i){
        numProbe_eta[i]++;
        if (pass_L1mu6) numL1_eta[i]++;
        //if (pass_SAmu6) numSA_eta[i]++;
        if (sa_pt>1e-5) {
          numSA_eta[i]++;
          //cout << "sa pt=" << sa_pt << endl;
        }
        else if (pass_L1mu6 && fabs(offline_eta)<1.05) saPt0barrel->Fill(probe_offline_phi);
        
      }
    }

    if (saddress==-1){//Endcap
    //if (fabs(offline_eta)>1.05 && fabs(offline_eta)<1.35){//EE region
    //if (fabs(offline_eta)>1.05){//Endcap
      if (sa_pt>1e-5) {
        saPtVsPtDifferenceEndcap->Fill(sa_pt,PtResidual);
        invSAPtVsPtDifferenceEndcap->Fill(1./sa_pt,PtResidual);
        invoffPtVsPtDifferenceEndcap->Fill(1./offline_pt,PtResidual);
        saPtVsDeltaPtEndcap->Fill(sa_pt,deltaPt);
        //cout << "endcap default sa pt=" << sa_pt << " +- " << deltaPt << endl;
      }
      if (fabs(deltaPt)>1e-5) {
        //cout << "Endcap sa pt=" << sa_pt << " +- " << deltaPt << endl; 
        PtVsDeltaPtEndcap->Fill(offline_pt,deltaPt);
        offPtVsCalcInvPtResidualEndcap->Fill(offline_pt,calcInvPtResidual);
      }
      offPtVsPtResidualEndcap->Fill(offline_pt,PtResidual);
      if (pass_SAmu6){
        offPtVsInvQptResidualEndcap[0]->Fill(offline_pt,invQptResidual);
        offPtVsInvPtResidualEndcap[0]->Fill(offline_pt,invPtResidual);
        if (charge_correct) {
          offPtVsInvQptResidualEndcap[1]->Fill(offline_pt,invQptResidual);
          offPtVsInvPtResidualEndcap[1]->Fill(offline_pt,invPtResidual);
        }
        else {
          offPtVsInvQptResidualEndcap[2]->Fill(offline_pt,invQptResidual);
          offPtVsInvPtResidualEndcap[2]->Fill(offline_pt,invPtResidual);
        }
      }

      //cout << "endcap pt resolution is " << pt_ec_reso << endl;
      for (int i=0;i<40;i++){
        if (offline_pt>i && offline_pt<i+1){
          numProbe_ec_pt[i]++;
          if (pass_L1mu6) numL1_ec_pt[i]++;
          if (pass_SAmu6) numSA_ec_pt[i]++;
        }
      }
      for (int i=0;i<80;i++){
        if (offline_pt>i && offline_pt<i+1){
          if (pass_SAmu6) {
            numSA_ec_charge_pt[i]++;
            if (charge_correct) numSA_ec_charge_correct_pt[i]++;
          }
        }
      }
      for (int s=0;s<32;s++){
        if (offline_phi>-3.2+0.2*s && offline_phi<-3+0.2*s){
          numProbe_ec_phi[s]++;
          if (pass_L1mu6) numL1_ec_phi[s]++;
          if (pass_SAmu6) numSA_ec_phi[s]++;
        }
      }
      for (int i=0;i<64;i++){
        if (offline_phi>0.1*i-3.2 && offline_phi<0.1*i-3.1){
          if (pass_SAmu6) {
            numSA_ec_charge_phi[i]++;
            if (charge_correct) numSA_ec_charge_correct_phi[i]++;
          }
        }
      }
    }
    else  {//Barrel
      if (sa_pt>1e-5) {
        saPtVsPtDifferenceBarrel->Fill(sa_pt,PtResidual);
        invSAPtVsPtDifferenceBarrel->Fill(1./sa_pt,PtResidual);
        invoffPtVsPtDifferenceBarrel->Fill(1./offline_pt,PtResidual);
        float res=0,defres=0;
        if (invSApt<0.184) res = invSApt3*vec[2] + invSApt2*vec[3] + invSApt*vec[4] + vec[5];
        else res = invSApt*vec[0] + vec[1];
        if (invSApt<0.186) defres = invSApt*defvec[0] + defvec[1];
        else defres = invSApt3*defvec[2] + invSApt2*defvec[3] + invSApt*defvec[4] + defvec[5];
        //cout << "res=" << res << endl; 
        newDeltaPt = res * sa_pt * sa_pt ;
        defDeltaPt = defres * sa_pt * sa_pt ;
        //cout << "new definition sa pt=" << sa_pt << " +- " << newDeltaPt << endl;
        PtVsNewDeltaPtBarrel[0]->Fill(offline_pt,newDeltaPt);
        if (fabs(deltaPt-defDeltaPt)>1){
          //cout << "saddress=" << saddress << endl;
          //cout << "deltapt online/calc=" << deltaPt << "/" << defDeltaPt << endl;
        }
        //cout << "barrel default sa pt=" << sa_pt << " +- " << deltaPt << endl;
        saPtVsDeltaPtBarrel->Fill(sa_pt,deltaPt);
        saPtVsNewDeltaPtBarrel->Fill(sa_pt,newDeltaPt);
        saPtVsDefDeltaPtBarrel->Fill(sa_pt,defDeltaPt);
        if (saddress==0) PtVsNewDeltaPtBarrel[1]->Fill(offline_pt,newDeltaPt);
        else if (saddress==1) PtVsNewDeltaPtBarrel[2]->Fill(offline_pt,newDeltaPt);
        else if (saddress==2) PtVsNewDeltaPtBarrel[3]->Fill(offline_pt,newDeltaPt);
        else if (saddress==3) PtVsNewDeltaPtBarrel[4]->Fill(offline_pt,newDeltaPt);
      }
      if(fabs(deltaPt)>1e-5) {
        //cout << "Barrel sa pt=" << sa_pt << " +- " << deltaPt << endl; 
        //deltaPt = deltaPt*10;
        PtVsDeltaPtBarrel[0]->Fill(offline_pt,deltaPt);
        if (saddress==0) PtVsDeltaPtBarrel[1]->Fill(offline_pt,deltaPt);
        else if (saddress==1) PtVsDeltaPtBarrel[2]->Fill(offline_pt,deltaPt);
        else if (saddress==2) PtVsDeltaPtBarrel[3]->Fill(offline_pt,deltaPt);
        else if (saddress==3) PtVsDeltaPtBarrel[4]->Fill(offline_pt,deltaPt);
        offPtVsCalcInvPtResidualBarrel->Fill(offline_pt,calcInvPtResidual);
      }
      offPtVsPtResidualBarrel[0]->Fill(offline_pt,PtResidual);
      if (saddress==0) offPtVsPtResidualBarrel[1]->Fill(offline_pt,PtResidual);
      else if (saddress==1) offPtVsPtResidualBarrel[2]->Fill(offline_pt,PtResidual);
      else if (saddress==2) offPtVsPtResidualBarrel[3]->Fill(offline_pt,PtResidual);
      else if (saddress==3) offPtVsPtResidualBarrel[4]->Fill(offline_pt,PtResidual);
      if (pass_SAmu6){
        offPtVsInvQptResidualBarrel[0]->Fill(offline_pt,invQptResidual);
        offPtVsInvPtResidualBarrel[0]->Fill(offline_pt,invPtResidual);
        if (charge_correct) {
          offPtVsInvQptResidualBarrel[1]->Fill(offline_pt,invQptResidual);
          offPtVsInvPtResidualBarrel[1]->Fill(offline_pt,invPtResidual);
        }
        else {
          offPtVsInvQptResidualBarrel[2]->Fill(offline_pt,invQptResidual);
          offPtVsInvPtResidualBarrel[2]->Fill(offline_pt,invPtResidual);
        }
      }
      for (int i=0;i<40;i++){
        if (offline_pt>i && offline_pt<i+1){
          numProbe_br_pt[i]++;
          if (pass_L1mu6) numL1_br_pt[i]++;
          if (pass_SAmu6) numSA_br_pt[i]++;
        }
      }
      for (int i=0;i<80;i++){
        if (offline_pt>i && offline_pt<i+1){
          if (pass_SAmu6) {
            numSA_br_charge_pt[i]++;
            if (charge_correct) numSA_br_charge_correct_pt[i]++;
          }
        }
      }
      for (int s=0;s<32;s++){
        if (offline_phi>-3.2+0.2*s && offline_phi<-3+0.2*s){
          numProbe_br_phi[s]++;
          if (pass_L1mu6) numL1_br_phi[s]++;
          if (pass_SAmu6) numSA_br_phi[s]++;
        }
      }
      for (int i=0;i<64;i++){
        if (offline_phi>0.1*i-3.2 && offline_phi<0.1*i-3.1){
          if (pass_SAmu6) {
            numSA_br_charge_phi[i]++;
            if (charge_correct) numSA_br_charge_correct_phi[i]++;
          }
        }
      }
    }
  }

  float eff_ec_pt_L1[40],efferr_ec_pt_L1[40],eff_ec_pt_SA[40],efferr_ec_pt_SA[40];
  float eff_br_pt_L1[40],efferr_br_pt_L1[40],eff_br_pt_SA[40],efferr_br_pt_SA[40];
  float ptbin[40],ptbinerr[40];
  for (int m=0;m<40;m++){
    //cout << "m=" << m << endl;
    //endcap
    eff_ec_pt_L1[m] = (numProbe_ec_pt[m]) ? 
      1.*numL1_ec_pt[m]/numProbe_ec_pt[m] : 0 ;
    efferr_ec_pt_L1[m]=sqrt( eff_ec_pt_L1[m]*(1-eff_ec_pt_L1[m]) )/numProbe_ec_pt[m] ;
    eff_ec_pt_SA[m]= (numL1_ec_pt[m]) ?
      1.*numSA_ec_pt[m]/numL1_ec_pt[m] : 0;
    efferr_ec_pt_SA[m]=sqrt( eff_ec_pt_SA[m]*(1-eff_ec_pt_SA[m]) )/numL1_ec_pt[m] ;
    //barrel
    eff_br_pt_L1[m] = (numProbe_br_pt[m]) ? 
      1.*numL1_br_pt[m]/numProbe_br_pt[m] : 0 ;
    efferr_br_pt_L1[m]=sqrt( eff_br_pt_L1[m]*(1-eff_br_pt_L1[m]) )/numProbe_br_pt[m] ;
    eff_br_pt_SA[m]= (numL1_br_pt[m]) ?
      1.*numSA_br_pt[m]/numL1_br_pt[m] : 0;
    efferr_br_pt_SA[m]=sqrt( eff_br_pt_SA[m]*(1-eff_br_pt_SA[m]) )/numL1_br_pt[m] ;
    ptbin[m]=0.5+m;
    ptbinerr[m] = 0.5;
    //cout << "bin" << m << ": Probe/L1/SA=" 
      //<< numProbe_ec_pt[m] << "/" << numL1_ec_pt[m] << "/" << numSA_ec_pt[m] << endl;
    //cout << "L1 efficiency is " << eff_ec_pt_L1[m] << " +- " << efferr_ec_pt_L1[m] << endl;
    //cout << "SA efficiency is " << eff_ec_pt_SA[m] << " +- " << efferr_ec_pt_SA[m] << endl;
  }

  float eff_ec_phi_L1[32],efferr_ec_phi_L1[32],eff_ec_phi_SA[32],efferr_ec_phi_SA[32];
  float eff_br_phi_L1[32],efferr_br_phi_L1[32],eff_br_phi_SA[32],efferr_br_phi_SA[32];
  float phibin[32],phibinerr[32];
  for (int m=0;m<32;m++){
    //endcap
    eff_ec_phi_L1[m] = (numProbe_ec_phi[m]) ? 
      1.*numL1_ec_phi[m]/numProbe_ec_phi[m] : 0 ;
    efferr_ec_phi_L1[m]=sqrt( eff_ec_phi_L1[m]*(1-eff_ec_phi_L1[m]) )/numProbe_ec_phi[m] ;
    eff_ec_phi_SA[m]= (numL1_ec_phi[m]) ?
      1.*numSA_ec_phi[m]/numL1_ec_phi[m] : 0;
    efferr_ec_phi_SA[m]=sqrt( eff_ec_phi_SA[m]*(1-eff_ec_phi_SA[m]) )/numL1_ec_phi[m] ;
    phibin[m]=-3.1+0.2*m;
    phibinerr[m] = 0.1;
    //barrel
    eff_br_phi_L1[m] = (numProbe_br_phi[m]) ? 
      1.*numL1_br_phi[m]/numProbe_br_phi[m] : 0 ;
    efferr_br_phi_L1[m]=sqrt( eff_br_phi_L1[m]*(1-eff_br_phi_L1[m]) )/numProbe_br_phi[m] ;
    eff_br_phi_SA[m]= (numL1_br_phi[m]) ?
      1.*numSA_br_phi[m]/numL1_br_phi[m] : 0;
    efferr_br_phi_SA[m]=sqrt( eff_br_phi_SA[m]*(1-eff_br_phi_SA[m]) )/numL1_br_phi[m] ;
    phibin[m]=-3.1+0.2*m;
    phibinerr[m] = 0.1;
  }

  float eff_eta_L1[50],efferr_eta_L1[50],eff_eta_SA[50],efferr_eta_SA[50];
  float etabin[50],etabinerr[50];
  for (int m=0;m<50;m++){
    //endcap
    eff_eta_L1[m] = (numProbe_eta[m]) ? 
      1.*numL1_eta[m]/numProbe_eta[m] : 0 ;
    efferr_eta_L1[m]=sqrt( eff_eta_L1[m]*(1-eff_eta_L1[m]) )/numProbe_eta[m] ;
    eff_eta_SA[m]= (numL1_eta[m]) ?
      1.*numSA_eta[m]/numL1_eta[m] : 0;
    efferr_eta_SA[m]=sqrt( eff_eta_SA[m]*(1-eff_eta_SA[m]) )/numL1_eta[m] ;
    etabin[m]=-2.45+0.1*m;
    etabinerr[m] = 0.05;
  }

  float eff_ec_charge_correct_pt[80],eff_br_charge_correct_pt[80],eff_ec_charge_correct_phi[64],eff_br_charge_correct_phi[64],eff_charge_correct_eta[32];
  float efferr_ec_charge_correct_pt[80],efferr_br_charge_correct_pt[80],efferr_ec_charge_correct_phi[64],efferr_br_charge_correct_phi[64],efferr_charge_correct_eta[32];
  float ptbin2[80],ptbinerr2[80];
  float phibin2[80],phibinerr2[80];
  
  for (int m=0;m<80;m++){
    //cout << "m=" << m << endl;
    //endcap
    eff_ec_charge_correct_pt[m]= (numSA_ec_charge_pt[m]) ?
      1.*numSA_ec_charge_correct_pt[m]/numSA_ec_charge_pt[m] : 0;
    efferr_ec_charge_correct_pt[m]=sqrt( eff_ec_charge_correct_pt[m]*(1-eff_ec_charge_correct_pt[m]) /numSA_ec_charge_pt[m] );
    //barrel
    eff_br_charge_correct_pt[m]= (numSA_br_charge_pt[m]) ?
      1.*numSA_br_charge_correct_pt[m]/numSA_br_charge_pt[m] : 0;
    efferr_br_charge_correct_pt[m]=sqrt( eff_br_charge_correct_pt[m]*(1-eff_br_charge_correct_pt[m]) /numSA_br_charge_pt[m] );
    ptbin2[m]=0.5+m;
    ptbinerr2[m] = 0.5;
    //cout << "bin" << m << ": Endcap Probe/SA/charge_correct=" 
      //<< numProbe_ec_pt[m] << "/" << numSA_ec_charge_pt[m] << "/" << numSA_ec_charge_correct_pt[m] << endl;
    //cout << "SA Endcap charge efficiency is " << eff_ec_charge_correct_pt[m] << " +- " << efferr_ec_charge_correct_pt[m] << endl;
    //cout << "bin" << m << ": Barrel Probe/SA/charge_correct=" 
      //<< numProbe_br_pt[m] << "/" << numSA_br_charge_pt[m] << "/" << numSA_br_charge_correct_pt[m] << endl;
    //cout << "SA Barrel charge efficiency is " << eff_br_charge_correct_pt[m] << " +- " << efferr_br_charge_correct_pt[m] << endl;
  }

  for (int m=0;m<64;m++){
    //cout << "m=" << m << endl;
    //endcap
    eff_ec_charge_correct_phi[m]= (numSA_ec_charge_phi[m]) ?
      1.*numSA_ec_charge_correct_phi[m]/numSA_ec_charge_phi[m] : 0;
    efferr_ec_charge_correct_phi[m]=sqrt( eff_ec_charge_correct_phi[m]*(1-eff_ec_charge_correct_phi[m]) /numSA_ec_charge_phi[m] );
    //barrel
    eff_br_charge_correct_phi[m]= (numSA_br_charge_phi[m]) ?
      1.*numSA_br_charge_correct_phi[m]/numSA_br_charge_phi[m] : 0;
    efferr_br_charge_correct_phi[m]=sqrt( eff_br_charge_correct_phi[m]*(1-eff_br_charge_correct_phi[m]) /numSA_br_charge_phi[m] );
    phibin2[m]=0.1*m-3.15;
    phibinerr2[m] = 0.05;
    //cout << "bin" << m << ": Endcap Probe/SA/charge_correct=" 
      //<< numProbe_ec_phi[m] << "/" << numSA_ec_charge_phi[m] << "/" << numSA_ec_charge_correct_phi[m] << endl;
    //cout << "SA Endcap charge efficiency is " << eff_ec_charge_correct_phi[m] << " +- " << efferr_ec_charge_correct_phi[m] << endl;
    //cout << "bin" << m << ": Barrel Probe/SA/charge_correct=" 
      //<< numProbe_br_phi[m] << "/" << numSA_br_charge_phi[m] << "/" << numSA_br_charge_correct_phi[m] << endl;
    //cout << "SA Barrel charge efficiency is " << eff_br_charge_correct_phi[m] << " +- " << efferr_br_charge_correct_phi[m] << endl;
  }

  TGraphErrors *effEndcapL1mu6 = new TGraphErrors(40,ptbin,eff_ec_pt_L1,ptbinerr,efferr_ec_pt_L1);
  effEndcapL1mu6->SetName("effEndcapL1mu6");
  effEndcapL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapL1mu6);
  TGraphErrors *effEndcapSAmu6 = new TGraphErrors(40,ptbin,eff_ec_pt_SA,ptbinerr,efferr_ec_pt_SA);
  effEndcapSAmu6->SetName("effEndcapSAmu6");
  effEndcapSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu6);
  TGraphErrors *effBarrelL1mu6 = new TGraphErrors(40,ptbin,eff_br_pt_L1,ptbinerr,efferr_br_pt_L1);
  effBarrelL1mu6->SetName("effBarrelL1mu6");
  effBarrelL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelL1mu6);
  TGraphErrors *effBarrelSAmu6 = new TGraphErrors(40,ptbin,eff_br_pt_SA,ptbinerr,efferr_br_pt_SA);
  effBarrelSAmu6->SetName("effBarrelSAmu6");
  effBarrelSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu6);
  TGraphErrors *effPhiEndcapL1mu6 = new TGraphErrors(32,phibin,eff_ec_phi_L1,phibinerr,efferr_ec_phi_L1);
  effPhiEndcapL1mu6->SetName("effPhiEndcapL1mu6");
  effPhiEndcapL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapL1mu6);
  TGraphErrors *effPhiEndcapSAmu6 = new TGraphErrors(32,phibin,eff_ec_phi_SA,phibinerr,efferr_ec_phi_SA);
  effPhiEndcapSAmu6->SetName("effPhiEndcapSAmu6");
  effPhiEndcapSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapSAmu6);
  TGraphErrors *effPhiBarrelL1mu6 = new TGraphErrors(32,phibin,eff_br_phi_L1,phibinerr,efferr_br_phi_L1);
  effPhiBarrelL1mu6->SetName("effPhiBarrelL1mu6");
  effPhiBarrelL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelL1mu6);
  TGraphErrors *effPhiBarrelSAmu6 = new TGraphErrors(32,phibin,eff_br_phi_SA,phibinerr,efferr_br_phi_SA);
  effPhiBarrelSAmu6->SetName("effPhiBarrelSAmu6");
  effPhiBarrelSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelSAmu6);
  TGraphErrors *effEtaL1mu6 = new TGraphErrors(50,etabin,eff_eta_L1,etabinerr,efferr_eta_L1);
  effEtaL1mu6->SetName("effEtaL1mu6");
  effEtaL1mu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEtaL1mu6);
  TGraphErrors *effEtaSAmu6 = new TGraphErrors(50,etabin,eff_eta_SA,etabinerr,efferr_eta_SA);
  effEtaSAmu6->SetName("effEtaSAmu6");
  effEtaSAmu6->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEtaSAmu6);
  //mu4 charge endcap pt
  TGraphErrors *effEndcapSAmu14Charge = new TGraphErrors(80,ptbin2,eff_ec_charge_correct_pt,ptbinerr2,efferr_ec_charge_correct_pt);
  effEndcapSAmu14Charge->SetName("effEndcapSAmu14Charge");
  effEndcapSAmu14Charge->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effEndcapSAmu14Charge);
  //mu4 charge barrel pt
  TGraphErrors *effBarrelSAmu14Charge = new TGraphErrors(80,ptbin2,eff_br_charge_correct_pt,ptbinerr2,efferr_br_charge_correct_pt);
  effBarrelSAmu14Charge->SetName("effBarrelSAmu14Charge");
  effBarrelSAmu14Charge->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effBarrelSAmu14Charge);
  //mu4 charge endcap phi
  TGraphErrors *effPhiEndcapSAmu14Charge = new TGraphErrors(64,phibin2,eff_ec_charge_correct_phi,phibinerr2,efferr_ec_charge_correct_phi);
  effPhiEndcapSAmu14Charge->SetName("effPhiEndcapSAmu14Charge");
  effPhiEndcapSAmu14Charge->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiEndcapSAmu14Charge);
  //mu4 charge barrel phi
  TGraphErrors *effPhiBarrelSAmu14Charge = new TGraphErrors(64,phibin2,eff_br_charge_correct_phi,phibinerr2,efferr_br_charge_correct_phi);
  effPhiBarrelSAmu14Charge->SetName("effPhiBarrelSAmu14Charge");
  effPhiBarrelSAmu14Charge->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(effPhiBarrelSAmu14Charge);

  file->Write();

  }
