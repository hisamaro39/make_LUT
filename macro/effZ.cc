#define effZ_cxx
#include "effZ.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <sstream>
#include <TLorentzVector.h>

void effZ::Loop()
{
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetOptStat(0);
  int nentries = fChain->GetEntries();
  cout << "Nubmer of events is " << nentries << endl;
  
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_ee_sym.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_separate.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_distance_separate_eept.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_final_separate_eept.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_selected.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_allregion.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_change_ee_default.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_change_ee_after.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_ee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_onlyee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_version9_nowindow.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_sampleA_r7447.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_z_sampleA_r7463.root","recreate");
  TFile *file = new TFile("outputEfficiency/result/efficiency_z_jira_12425.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/aho.root","recreate");

  int numProbePt[2][80],numPt[2][2][2][3][80];
  int numProbePhi[2][64],numPhi[2][2][2][3][64];
  float eff_pt[2][2][2][3][80],efferr_pt[2][2][2][3][80],xbin_pt[80],xbinerr_pt[80];
  float eff_phi[2][2][2][3][64],efferr_phi[2][2][2][3][64],xbin_phi[64],xbinerr_phi[64];
  float eff_eta_phi[2][2][3][50][64];
  int numProbeEtaPhi[50][64],numEtaPhi[2][2][3][50][64];
  for (int pt=0;pt<80;pt++){//pt bin
    xbin_pt[pt]=0.5+pt;
    xbinerr_pt[pt]=0.5;
    for (int side=0;side<2;side++){//barrel or endcap
      numProbePt[side][pt]=0;
      for (int level=0;level<2;level++){//L1 or SA
        for (int chain=0;chain<2;chain++){//mu4 or mu6
          for (int type=0;type<3;type++){//default, new or ee LUT
            numPt[level][chain][side][type][pt]=0;
            eff_pt[level][chain][side][type][pt]=0;
            efferr_pt[level][chain][side][type][pt]=0;
          }
        }
      }
    }
  }

  for (int phi=0;phi<64;phi++){//phi bin
    xbin_phi[phi]=-3.15+0.1*phi;
    xbinerr_phi[phi]=0.05;
    for (int side=0;side<2;side++){//barrel or endcap
      numProbePhi[side][phi]=0;
      for (int level=0;level<2;level++){//L1 or SA
        for (int chain=0;chain<2;chain++){//mu4 or mu6
          for (int type=0;type<3;type++){//default, new or ee LUT
            numPhi[level][chain][side][type][phi]=0;
            eff_phi[level][chain][side][type][phi]=0;
            efferr_phi[level][chain][side][type][phi]=0;
          }
        }
      }
    }
  }
  TH1 *SAPtResolutionEndcap[12],*SAPtEEResolutionEndcap[12],*SAPtNewResolutionEndcap[12];
  TH1 *SAPhiResolutionEndcap[32],*SAPhiEEResolutionEndcap[32];
  TH1 *SAEtaResolutionEndcap[100],*SAEtaEEResolutionEndcap[100];
  stringstream saResoDef,saResoEE,saResoNew;
  for (int i=0;i<12;i++){
    saResoDef << "SAPtResolutionEndcap";
    saResoEE << "SAPtEEResolutionEndcap";
    saResoNew << "SAPtNewResolutionEndcap";
    if (i<9){
      saResoDef << 15+i*5 ;saResoEE << 15+i*5; saResoNew << 15+i*5 ;
    }
    else {
      saResoDef << 60+(i-9)*10; saResoEE << 60+(i-9)*10; saResoNew << 60+(i-9)*10;
    }
    SAPtResolutionEndcap[i] = new TH1F(saResoDef.str().c_str(),";pt resolution;Events",100,-1,1);
    SAPtEEResolutionEndcap[i] = new TH1F(saResoEE.str().c_str(),";pt resolution;Events",100,-1,1);
    SAPtNewResolutionEndcap[i] = new TH1F(saResoNew.str().c_str(),";pt resolution;Events",100,-1,1);
    saResoDef.str("");saResoEE.str(""),saResoNew.str("");
  }
  for (int i=0;i<32;i++){
    saResoDef << "SAPhiResolutionEndcap";
    saResoEE << "SAPhiEEResolutionEndcap";
    saResoDef << i ;saResoEE << i ;
    SAPhiResolutionEndcap[i] = new TH1F(saResoDef.str().c_str(),";pt resolution;Events",100,-1,1);
    SAPhiEEResolutionEndcap[i] = new TH1F(saResoEE.str().c_str(),";pt resolution;Events",100,-1,1);
    saResoDef.str("");saResoEE.str("");
  }
  for (int i=0;i<100;i++){
    saResoDef << "SAEtaResolution";
    saResoEE << "SAEtaEEResolution";
    saResoDef << i ;saResoEE << i ;
    SAEtaResolutionEndcap[i] = new TH1F(saResoDef.str().c_str(),";pt resolution;Events",100,-1,1);
    SAEtaEEResolutionEndcap[i] = new TH1F(saResoEE.str().c_str(),";pt resolution;Events",100,-1,1);
    saResoDef.str("");saResoEE.str("");
  }
  TH2 *badEtaPhi = new TH2F("badEtaPhi",";SA #eta;SA #phi",100,-2.5,2.5,100,-3.2,3.2);
  TH1 *probePt = new TH1F("probePt",";p_{T,off};Number of events",200,0,100);
  TH1 *DiMuonMass = new TH1F("DiMuonMass",";M_{#mu#mu}(GeV);Number of events",100,60,120);
  TH2 *RoIEtaPhi = new TH2F("RoIEtaPhi",";RoI #eta;RoI #phi",100,-2.5,2.5,100,-3.2,3.2);
  TH2 *RoIEtaPhiAll = new TH2F("RoIEtaPhiAll",";RoI #eta;RoI #phi",100,-2.5,2.5,100,-3.2,3.2);
  TH2 *RoIEtaPhiSAFail = new TH2F("RoIEtaPhiSAFail",";RoI #eta;RoI #phi",100,-2.5,2.5,100,-3.2,3.2);
  TH1 *RoIEtaSAFail = new TH1F("RoIEtaSAFail",";RoI #eta;",100,-2.5,2.5);
  TH1 *RoIPhiSAFail = new TH1F("RoIPhiSAFail",";RoI #phi;",100,-3.2,3.2);

  for (int eta=0;eta<50;eta++){//eta bin
    for (int phi=0;phi<64;phi++){//phi bin
      numProbeEtaPhi[eta][phi]=0;
      for (int level=0;level<2;level++){//L1 or SA
        for (int chain=0;chain<2;chain++){//mu4 or mu6
          for (int type=0;type<3;type++){//default, new or ee LUT
            numEtaPhi[level][chain][type][eta][phi]=0;
            eff_eta_phi[level][chain][type][eta][phi]=0;
          }
        }
      }
    }
  }


  for (int jentry=0; jentry<nentries;jentry++) {
    int ientry = LoadTree(jentry);
    fChain->GetEntry(jentry); 

    float offline_pt = probe_offline_pt; 
    float offline_eta = probe_offline_eta; 
    float offline_phi = probe_offline_phi; 
    probePt->Fill(offline_pt);;
    float offline_pt2 = tag_offline_pt; 
    float offline_eta2 = tag_offline_eta; 
    float offline_phi2 = tag_offline_phi; 
    float muon_mass = 0.1056; 
    TLorentzVector muon1,muon2;
    muon1.SetPtEtaPhiM(offline_pt,offline_eta,offline_phi,muon_mass);
    muon2.SetPtEtaPhiM(offline_pt2,offline_eta2,offline_phi2,muon_mass);
    float mass = (muon1+muon2).M();
    DiMuonMass->Fill(mass);
    //cout << "probe pt/eta/phi=" << offline_pt << "/" << offline_eta << "/" << offline_phi << endl; 
    //cout << "tag pt/eta/phi=" << tag_offline_pt << "/" << tag_offline_eta << "/" << tag_offline_phi << endl; 
    //cout << "dimuon mass=" << mass << endl;
    float sa_pt = fabs(probe_sa_pt);
    float sa_pt_ee = fabs(probe_sa_pt_new_ee);
    float sa_pt_new = fabs(probe_sa_pt_new);
    bool passL1mu4 = probe_passL1mu4;
    bool passSAmu4 = probe_passSAmu4;
    bool passSAmu4_new = probe_passSAmu4new;
    bool passSAmu4_ee = probe_passSAmu4newee;
    bool passL1mu6 = probe_passL1mu6;
    bool passSAmu6 = probe_passSAmu6;
    bool passSAmu6_new = probe_passSAmu6new;
    bool passSAmu6_ee = probe_passSAmu6newee;
    int saddress = probe_sa_saddress;
    float saPtResolution = (sa_pt>1e-5)? 1-offline_pt/sa_pt : 1000;
    float saPtEEResolution = (sa_pt_ee>1e-5)? 1-offline_pt/sa_pt_ee : 1000;
    float saPtNewResolution = (sa_pt_new>1e-5)? 1-offline_pt/sa_pt_new : 1000;
    if (sa_pt<1e-5 && (passL1mu4 || passL1mu6)) {
      RoIEtaPhiSAFail->Fill(probe_roi_eta,probe_roi_phi);
      RoIEtaSAFail->Fill(probe_roi_eta);
      RoIPhiSAFail->Fill(probe_roi_phi);
      if(fabs(probe_roi_eta)>2){
        cout << "RoI eta/phi=" << probe_roi_eta << "/" << probe_roi_phi << endl;
        cout << "sa alpha/beta=" << probe_sa_alpha << "/" << probe_sa_beta << endl;
      }
    }
    if (fabs(offline_eta)<1.0){//barrel
    //if (saddress!=-1){//barrel
      for (int ipt=0;ipt<80;ipt++){
        if (offline_pt>ipt && offline_pt<ipt+1){
          numProbePt[0][ipt]++;
          if (passL1mu4) {
            numPt[0][0][0][0][ipt]++; 
            numPt[0][0][0][1][ipt]++; 
            numPt[0][0][0][2][ipt]++; 
          }
          if (passL1mu6) {
            numPt[0][1][0][0][ipt]++; 
            numPt[0][1][0][1][ipt]++; 
            numPt[0][1][0][2][ipt]++; 
          }
          if (passSAmu4) numPt[1][0][0][0][ipt]++; 
          if (passSAmu4_new) numPt[1][0][0][1][ipt]++; 
          if (passSAmu4_ee) numPt[1][0][0][2][ipt]++; 
          if (passSAmu6) numPt[1][1][0][0][ipt]++; 
          if (passSAmu6_new) numPt[1][1][0][1][ipt]++; 
          if (passSAmu6_ee) numPt[1][1][0][2][ipt]++; 
        }
      }
      for (int iphi=0;iphi<64;iphi++){
        if (offline_phi>-3.2+iphi*0.1 && offline_phi<-3.1+iphi*0.1){
          numProbePhi[0][iphi]++;
          if (passL1mu4) {
            numPhi[0][0][0][0][iphi]++; 
            numPhi[0][0][0][1][iphi]++; 
            numPhi[0][0][0][2][iphi]++; 
          }
          if (passL1mu6) {
            numPhi[0][1][0][0][iphi]++; 
            numPhi[0][1][0][1][iphi]++; 
            numPhi[0][1][0][2][iphi]++; 
          }
          if (passSAmu4) numPhi[1][0][0][0][iphi]++; 
          if (passSAmu4_new) numPhi[1][0][0][1][iphi]++; 
          if (passSAmu4_ee) numPhi[1][0][0][2][iphi]++; 
          if (passSAmu6) numPhi[1][1][0][0][iphi]++; 
          if (passSAmu6_new) numPhi[1][1][0][1][iphi]++; 
          if (passSAmu6_ee) numPhi[1][1][0][2][iphi]++; 
        }
      }//phi loop end
    }
    else {//endcap
      //if (fabs(offline_eta)>1.4) continue;//EE region
      //if (fabs(sa_pt-sa_pt_ee)<1e-5) continue;
      if (sa_pt_ee-sa_pt>50){
        //cout << "****************" << endl;
        //cout << "pt offline/default/ee=" << offline_pt << "/" << sa_pt << "/" << sa_pt_ee << endl;
        //cout << "sa eta/phi=" << probe_sa_eta << "/" << probe_sa_phi << endl;
        //cout << "resolution default/ee=" << saPtResolution << "/" << saPtEEResolution << endl;
        badEtaPhi->Fill(probe_sa_eta,probe_sa_phi);
      }
      if (saPtResolution<999){
        for (int i=0;i<9;i++)
          if (offline_pt>15+i*5 && offline_pt<20+i*5)
            SAPtResolutionEndcap[i]->Fill(saPtResolution);
        for (int i=9;i<12;i++)
          if (offline_pt>60+(i-9)*10 && offline_pt<60+(i-8)*10)
            SAPtResolutionEndcap[i]->Fill(saPtResolution);
        for (int j=0;j<32;j++)
          if (offline_phi>-3.2+0.2*j && offline_phi<-3+0.2*j)
            SAPhiResolutionEndcap[j]->Fill(saPtResolution);
      }
      if (saPtEEResolution<999){
        for (int i=0;i<9;i++)
          if (offline_pt>15+i*5 && offline_pt<20+i*5)
            SAPtEEResolutionEndcap[i]->Fill(saPtEEResolution);
        for (int i=9;i<12;i++)
          if (offline_pt>60+(i-9)*10 && offline_pt<60+(i-8)*10)
            SAPtEEResolutionEndcap[i]->Fill(saPtEEResolution);
        for (int j=0;j<32;j++)
          if (offline_phi>-3.2+0.2*j && offline_phi<-3+0.2*j)
            SAPhiEEResolutionEndcap[j]->Fill(saPtEEResolution);
        for (int k=0;k<100;k++)
          if (offline_eta>-2.5+0.05*k && offline_eta<-2.45+0.05*k)
            SAEtaEEResolutionEndcap[k]->Fill(saPtEEResolution);
      }
      if (saPtNewResolution<999){
        for (int i=0;i<9;i++)
          if (offline_pt>15+i*5 && offline_pt<20+i*5)
            SAPtNewResolutionEndcap[i]->Fill(saPtNewResolution);
        for (int i=9;i<12;i++)
          if (offline_pt>60+(i-9)*10 && offline_pt<60+(i-8)*10)
            SAPtNewResolutionEndcap[i]->Fill(saPtNewResolution);
      }

      for (int ipt=0;ipt<80;ipt++){
        if (offline_pt>ipt && offline_pt<ipt+1){
          numProbePt[1][ipt]++;
          if (passL1mu4) {
            numPt[0][0][1][0][ipt]++; 
            numPt[0][0][1][1][ipt]++; 
            numPt[0][0][1][2][ipt]++; 
          }
          if (passL1mu6) {
            numPt[0][1][1][0][ipt]++; 
            numPt[0][1][1][1][ipt]++; 
            numPt[0][1][1][2][ipt]++; 
          }
          if (passSAmu4) numPt[1][0][1][0][ipt]++; 
          if (passSAmu4_new) numPt[1][0][1][1][ipt]++; 
          if (passSAmu4_ee) numPt[1][0][1][2][ipt]++; 
          if (passSAmu6) numPt[1][1][1][0][ipt]++; 
          if (passSAmu6_new) numPt[1][1][1][1][ipt]++; 
          if (passSAmu6_ee) numPt[1][1][1][2][ipt]++; 
        }
      }
      for (int iphi=0;iphi<64;iphi++){
        if (offline_phi>-3.2+iphi*0.1 && offline_phi<-3.1+iphi*0.1){
          numProbePhi[1][iphi]++;
          if (passL1mu4) {
            numPhi[0][0][1][0][iphi]++; 
            numPhi[0][0][1][1][iphi]++; 
            numPhi[0][0][1][2][iphi]++; 
          }
          if (passL1mu6) {
            numPhi[0][1][1][0][iphi]++; 
            numPhi[0][1][1][1][iphi]++; 
            numPhi[0][1][1][2][iphi]++; 
          }
          if (passSAmu4) numPhi[1][0][1][0][iphi]++; 
          if (passSAmu4_new) numPhi[1][0][1][1][iphi]++; 
          if (passSAmu4_ee) numPhi[1][0][1][2][iphi]++; 
          if (passSAmu6) numPhi[1][1][1][0][iphi]++; 
          if (passSAmu6_new) numPhi[1][1][1][1][iphi]++; 
          if (passSAmu6_ee) numPhi[1][1][1][2][iphi]++; 
        }
      }//phi loop end
    }//endcap end

    if (saPtResolution<999){
      for (int k=0;k<100;k++)
        if (offline_eta>-2.5+0.05*k && offline_eta<-2.45+0.05*k)
          SAEtaResolutionEndcap[k]->Fill(saPtResolution);
    }
    if (saPtEEResolution<999){
      for (int k=0;k<100;k++)
        if (offline_eta>-2.5+0.05*k && offline_eta<-2.45+0.05*k)
          SAEtaEEResolutionEndcap[k]->Fill(saPtEEResolution);
    }
    
    RoIEtaPhiAll->Fill(probe_roi_eta,probe_roi_phi);
    for (int ieta=0;ieta<50;ieta++){
      for (int iphi=0;iphi<64;iphi++){
        if (offline_eta>-2.5+ieta*0.1 && offline_eta<-2.4+ieta*0.1){
          if (offline_phi>-3.2+iphi*0.1 && offline_phi<-3.1+iphi*0.1){
            if(offline_pt>4) {
              numProbeEtaPhi[ieta][iphi]++;
              if(fabs(probe_roi_eta)>1e-5) RoIEtaPhi->Fill(probe_roi_eta,probe_roi_phi);
            }
            if (passL1mu4 && offline_pt>4) {
              numEtaPhi[0][0][0][ieta][iphi]++; 
              numEtaPhi[0][0][1][ieta][iphi]++; 
              numEtaPhi[0][0][2][ieta][iphi]++; 
            }
            if (passL1mu6 && offline_pt>6) {
              numEtaPhi[0][1][0][ieta][iphi]++; 
              numEtaPhi[0][1][1][ieta][iphi]++; 
              numEtaPhi[0][1][2][ieta][iphi]++; 
            }
            if (passSAmu4 && offline_pt>4) numEtaPhi[1][0][0][ieta][iphi]++; 
            if (passSAmu4_new) numEtaPhi[1][0][1][ieta][iphi]++; 
            if (passSAmu4_ee) numEtaPhi[1][0][2][ieta][iphi]++; 
            if (passSAmu6 && offline_pt>6) numEtaPhi[1][1][0][ieta][iphi]++; 
            if (passSAmu6_new) numEtaPhi[1][1][1][ieta][iphi]++; 
            if (passSAmu6_ee) numEtaPhi[1][1][2][ieta][iphi]++; 
          }
        }
      }
    }//eta-phi loop end
  }

  TGraphErrors* EffPt[2][2][2][3];
  string LEVEL[2] = {"L1","SA"};
  string CHAIN[2] = {"Mu4","Mu6"};
  string SIDE[2] = {"Barrel","Endcap"};
  string TYPE[3] = {"Default","New","EE"};
  stringstream effName; 
  for (int level=0;level<2;level++){//L1 or SA
    for (int chain=0;chain<2;chain++){//mu4 or mu6
      for (int side=0;side<2;side++){//barrel or endcap
        for (int type=0;type<3;type++){//default, new or ee LUT
          for (int pt=0;pt<80;pt++){
            if (type==1){
              //cout << "level/chain/side/type/pt=" << level << "/" << chain << "/" << side << "/" << type << "/" << pt << endl;
              //cout << numPt[level][chain][side][type][pt] << endl;
            }
            if (level==0){//L1 efficiency
              if (numProbePt[side][pt]>0){ 
                eff_pt[level][chain][side][type][pt] = 
                  1.*numPt[level][chain][side][type][pt]/numProbePt[side][pt];
                efferr_pt[level][chain][side][type][pt] = 
                  sqrt(eff_pt[level][chain][side][type][pt]*(1-eff_pt[level][chain][side][type][pt]) / numProbePt[side][pt]);
              }
            }
            if (level==1){//SA efficiency
              if (numPt[0][chain][side][type][pt]>0){ 
                eff_pt[level][chain][side][type][pt] = 
                  1.*numPt[level][chain][side][type][pt]/numPt[0][chain][side][type][pt];
                efferr_pt[level][chain][side][type][pt] = 
                  sqrt(eff_pt[level][chain][side][type][pt]*(1-eff_pt[level][chain][side][type][pt]) / numPt[0][chain][side][type][pt]);
              }
            }
          }
          effName << "Efficiency" << LEVEL[level] << CHAIN[chain] << SIDE[side] << TYPE[type];
          EffPt[level][chain][side][type] = new TGraphErrors(80,xbin_pt,eff_pt[level][chain][side][type],xbinerr_pt,efferr_pt[level][chain][side][type]);
          EffPt[level][chain][side][type]->SetName(effName.str().c_str());
          EffPt[level][chain][side][type]->SetTitle("");
          EffPt[level][chain][side][type]->GetXaxis()->SetLabelSize(0.07);
          EffPt[level][chain][side][type]->GetYaxis()->SetLabelSize(0.07);
          file->Add(EffPt[level][chain][side][type]);
          effName.str("");
        }
      }
    }
  }

  TGraphErrors* EffPhi[2][2][2][3];
  for (int level=0;level<2;level++){//L1 or SA
    for (int chain=0;chain<2;chain++){//mu4 or mu6
      for (int side=0;side<2;side++){//barrel or endcap
        for (int type=0;type<3;type++){//default, new or ee LUT
          for (int phi=0;phi<64;phi++){
            if (level==0){//L1 efficiency
              if (numProbePhi[side][phi]>0){ 
                eff_phi[level][chain][side][type][phi] = 
                  1.*numPhi[level][chain][side][type][phi]/numProbePhi[side][phi];
                efferr_phi[level][chain][side][type][phi] = 
                  sqrt(eff_phi[level][chain][side][type][phi]*(1-eff_phi[level][chain][side][type][phi]) / numProbePhi[side][phi]);
              }
            }
            if (level==1){//SA efficiency
              if (numPhi[0][chain][side][type][phi]>0){ 
                eff_phi[level][chain][side][type][phi] = 
                  1.*numPhi[level][chain][side][type][phi]/numPhi[0][chain][side][type][phi];
                efferr_phi[level][chain][side][type][phi] = 
                  sqrt(eff_phi[level][chain][side][type][phi]*(1-eff_phi[level][chain][side][type][phi]) / numPhi[0][chain][side][type][phi]);
              }
            }
          }
          effName << "EfficiencyPhi" << LEVEL[level] << CHAIN[chain] << SIDE[side] << TYPE[type];
          EffPhi[level][chain][side][type] = new TGraphErrors(64,xbin_phi,eff_phi[level][chain][side][type],xbinerr_phi,efferr_phi[level][chain][side][type]);
          EffPhi[level][chain][side][type]->SetName(effName.str().c_str());
          EffPhi[level][chain][side][type]->SetTitle("");
          EffPhi[level][chain][side][type]->GetXaxis()->SetLabelSize(0.07);
          EffPhi[level][chain][side][type]->GetYaxis()->SetLabelSize(0.07);
          EffPhi[level][chain][side][type]->GetYaxis()->SetRangeUser(0,1.2);
          file->Add(EffPhi[level][chain][side][type]);
          effName.str("");
        }
      }
    }
  }
  
  TH2 *EffEtaPhi[2][2][3];
  for (int level=0;level<2;level++){//L1 or SA
    for (int chain=0;chain<2;chain++){//mu4 or mu6
      for (int type=0;type<3;type++){//default, new or ee LUT
        effName << "EfficiencyEtaPhi" << LEVEL[level] << CHAIN[chain] << TYPE[type];
        EffEtaPhi[level][chain][type] = new TH2F(effName.str().c_str(),";offline #eta;offline #phi",50,-2.5,2.5,64,-3.2,3.2);
        EffEtaPhi[level][chain][type]->GetXaxis()->SetLabelSize(0.07);
        EffEtaPhi[level][chain][type]->GetYaxis()->SetLabelSize(0.07);
        for (int eta=0;eta<50;eta++){
          for (int phi=0;phi<64;phi++){
            if (level==0){//L1 efficiency
              if (numProbeEtaPhi[eta][phi]>0){ 
                eff_eta_phi[level][chain][type][eta][phi] = 
                  1.*numEtaPhi[level][chain][type][eta][phi]/numProbeEtaPhi[eta][phi];
              }
            }
            if (level==1){//SA efficiency
              if (numEtaPhi[0][chain][type][eta][phi]>0){ 
                eff_eta_phi[level][chain][type][eta][phi] = 
                  1.*numEtaPhi[level][chain][type][eta][phi]/numEtaPhi[0][chain][type][eta][phi];
              }
            }
            //cout << CHAIN[chain] << "/" << TYPE[type] << endl;
            //cout << "eta/phi=" << eta << "/" << phi << "/" << endl; 
            //cout << "probe/L1/SA=" << numProbeEtaPhi[eta][phi] << "/" << numEtaPhi[0][chain][type][eta][phi] << "/" 
              //<< numEtaPhi[1][chain][type][eta][phi] << endl; 
            //cout << "efficiency L1/Probe is "  << eff_eta_phi[0][chain][type][eta][phi] << endl;
            //cout << "efficiency SA/L1 is "  << eff_eta_phi[1][chain][type][eta][phi] << endl;
            EffEtaPhi[level][chain][type]->Fill(-2.45+0.1*eta,-3.15+0.1*phi,eff_eta_phi[level][chain][type][eta][phi]);
          }
        }
        effName.str("");
      }
    }
  }
  
  file->Write();

}
