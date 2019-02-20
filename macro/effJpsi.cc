#define effJpsi_cxx
#include "effJpsi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraphErrors.h"
#include <iostream>
#include <sstream>
#include <TLorentzVector.h>
using namespace std;

#define PI 3.14159265258979

double const PI_OVER_4 = PI/4.0;
double const PI_OVER_8 = PI/8.0;
double const PI_OVER_16 = PI/16.0;
double const PHI_RANGE = 12.0/PI_OVER_8;
double const ZERO_LIMIT = 1e-5;

int GetPhiBinNumber(float phi){
  int Octant = (int)(phi/PI_OVER_4);
  double PhiInOctant = fabs(phi - Octant*PI_OVER_4);
  if(PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  int phiBin = static_cast<int>(PhiInOctant*PHI_RANGE);
  if (phiBin < -0.5 || phiBin > 11.5) return -1;
  return phiBin;
}

int GetEtaBinNumber(float eta){
  int etaBin = static_cast<int>((fabs(eta)-1.)/0.05);
  if (etaBin == -1) etaBin =  0;
  if (etaBin == 30) etaBin = 29;
  if (etaBin < -0.5 || etaBin > 29.5) return -1;
  return etaBin;
}

int GetPhiBinAllNumber(float phi){
  if (phi<0) phi = phi + 2*PI;
  float phibin = (int) (phi * 96/PI);

  return phibin;
}

float calcDistance(float x1,float y1,float x2,float y2,float x3,float y3){
  float xm1=(x1+x2)/2;
  float xm2=(x2+x3)/2;
  float ym1=(y1+y2)/2;
  float ym2=(y2+y3)/2;
  float a1=(x1-x2)/(y2-y1);
  float a2=(x2-x3)/(y3-y2);
  float x0=(a2*xm2-a1*xm1-ym2+ym1)/(a2-a1);//center of circle
  float y0=a1*(x0-xm1)+ym1;//center of circle
  float a = (x0-x1)/(y1-y0);//slope of sessen
  float b = y1+x1*(x1-x0)/(y1-y0);//intercept of sessen
  float d=fabs(b)/sqrt(a*a+1);
  return d;
}

bool isSmall(float eez){
  bool small = false;
  if (fabs(eez)>10000 && fabs(eez)<10600) small =true;
  else if(fabs(eez)>10600 && fabs(eez)<12000) small = false;
  return small;
}

double computeRadius3Points(double InnerZ, double InnerR, 
    double EEZ, double EER,
    double MiddleZ, double MiddleR)
{
  double radius_EE;
  double a3;
  double m = 0.;
  double cost = 0.;
  double x0 = 0., y0 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0.;
  double tm = 0.;
  a3 = ( MiddleZ - InnerZ ) / ( MiddleR - InnerR );
  m = a3;
  cost = cos(atan(m));
  x2 = EER - InnerR;
  y2 = EEZ - InnerZ;
  x3 = MiddleR - InnerR;
  y3 = MiddleZ - InnerZ;
  tm = x2;
  x2 = ( x2   + y2*m)*cost;
  y2 = (-tm*m + y2  )*cost;
  tm = x3;
  x3 = ( x3   + y3*m)*cost;
  y3 = (-tm*m + y3  )*cost;
  x0 = x3/2.;
  y0 = (y2*y2 + x2*x2 -x2*x3)/(2*y2);
  radius_EE = sqrt(x0*x0 + y0*y0);
  return radius_EE;
}

float calcDis(float r1,float z1,float r2,float z2){
  float dr = sqrt((r1-r2)*(r1-r2)+(z1-z2)*(z1-z2));
  return dr;
}

void effJpsi::Loop()
{
  /*gStyle->SetLabelSize(0.07,"x");
  gStyle->SetLabelSize(0.07,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetOptStat(0);
  */

  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee_distance.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee_distance_separate_eept.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee_distance_separate_eept_noee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee_distance_all_separate_eept_noee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee_distance_all_allregion.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee_all.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_final.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_selected.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_separate.root","recreate");
  TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_separate_sa_offline.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_ee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_onlyee.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_no_window.root","recreate");
  //TFile *file = new TFile("outputEfficiency/result/efficiency_jpsi_version9_offline_charge.root","recreate");
  //TFile *file = new TFile("outputEfficiency/aho.root","recreate");

  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;

  int numProbePt[2][60],numPt[2][2][2][3][60];
  int numProbeEta[50],numEta[2][2][3][50];
  int numProbePhi[2][64],numPhi[2][2][2][3][64];
  int numProbeEtaPhi[50][64],numEtaPhi[2][2][3][50][64];
  int numEERegionPt[60],numUseEEPt[60];
  int numEERegionEta[60],numUseEEEta[60];
  int numEERegionPhi[64],numUseEEPhi[64];
  int numPassSA[2][15],numChargeCorrect[2][15];
  float eff_pt[2][2][2][3][60],efferr_pt[2][2][2][3][60],xbin_pt[60],xbinerr_pt[60];
  float eff_eta[2][2][3][50],efferr_eta[2][2][3][50],xbin_eta[50],xbinerr_eta[50];
  float eff_phi[2][2][2][3][64],efferr_phi[2][2][2][3][64],xbin_phi[64],xbinerr_phi[64];
  float eff_use_ee_pt[60],efferr_use_ee_pt[60],xbin_use_ee_pt[60],xbinerr_use_ee_pt[60];
  float eff_use_ee_eta[60],efferr_use_ee_eta[60],xbin_use_ee_eta[60],xbinerr_use_ee_eta[60];
  float eff_use_ee_phi[64],efferr_use_ee_phi[64],xbin_use_ee_phi[64],xbinerr_use_ee_phi[64];
  float eff_eta_phi[2][2][3][50][64],xbin_charge[15],xbinerr_charge[15];
  float eff_charge_correct[2][15],efferr_charge_correct[2][15];

  for (int i=0;i<2;i++){
    for (int j=0;j<15;j++){
      numPassSA[i][j]=0;
      numChargeCorrect[i][j]=0;
      eff_charge_correct[i][j]=0;
      xbin_charge[j]=0.5+j;
      xbinerr_charge[j]=0.5;
    }
  }

  for (int ee=0;ee<60;ee++){
    eff_use_ee_pt[ee]=0;
    xbin_use_ee_pt[ee]=0.125+0.25*ee;
    xbinerr_use_ee_pt[ee]=0.125;
  }

  for (int ee=0;ee<60;ee++){
    eff_use_ee_eta[ee]=0;
    xbin_use_ee_eta[ee]=-1.475+0.05*ee;
    xbinerr_use_ee_eta[ee]=0.025;
  }

  for (int ee=0;ee<64;ee++){
    eff_use_ee_phi[ee]=0;
    xbin_use_ee_phi[ee]=-3.15+0.1*ee;
    xbinerr_use_ee_phi[ee]=0.05;
  }

  for (int pt=0;pt<60;pt++){//pt bin
    xbin_pt[pt]=0.125+0.25*pt;
    xbinerr_pt[pt]=0.125;
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

  for (int eta=0;eta<50;eta++){//eta bin
    xbin_eta[eta]=-2.45+0.1*eta;
    xbinerr_eta[eta]=0.05;
    numProbeEta[eta]=0;
    for (int level=0;level<2;level++){//L1 or SA
      for (int chain=0;chain<2;chain++){//mu4 or mu6
        for (int type=0;type<3;type++){//default, new or ee LUT
          numEta[level][chain][type][eta]=0;
          eff_eta[level][chain][type][eta]=0;
          efferr_eta[level][chain][type][eta]=0;
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

  TH1 *SAPtResolutionEndcap[20],*SAPtNewResolutionEndcap[20],*SAPtEEResolutionEndcap[20];
  TH1 *SAPtResolutionBarrel[20];
  TH1 *SAEtaResolutionEndcap[100],*SAEtaNewResolutionEndcap[100],*SAEtaEEResolutionEndcap[100];
  stringstream saResoDef,saResoNew,saResoEE,saResoBarrel;
  for (int i=0;i<20;i++){
    saResoDef << "SAPtResolutionEndcap";
    saResoBarrel << "SAPtResolutionBarrel";
    saResoNew << "SAPtNewResolutionEndcap";
    saResoEE << "SAPtEEResolutionEndcap";
    saResoDef << i; saResoNew << i; saResoEE << i ; saResoBarrel << i;
    SAPtResolutionEndcap[i] = new TH1F(saResoDef.str().c_str(),";;",100,-1,1);
    SAPtResolutionBarrel[i] = new TH1F(saResoBarrel.str().c_str(),";;",100,-1,1);
    SAPtNewResolutionEndcap[i] = new TH1F(saResoNew.str().c_str(),";;",100,-1,1);
    SAPtEEResolutionEndcap[i] = new TH1F(saResoEE.str().c_str(),";;",100,-1,1);
    saResoDef.str("");saResoNew.str(""),saResoEE.str(""),saResoBarrel.str("");
  }
  for (int i=0;i<100;i++){
    saResoDef << "SAEtaResolution";
    saResoNew << "SAEtaNewResolution";
    saResoEE << "SAEtaEEResolution";
    saResoDef << i ;saResoNew << i ;saResoEE << i;
    SAEtaResolutionEndcap[i] = new TH1F(saResoDef.str().c_str(),";pt resolution;Events",100,-1,1);
    SAEtaNewResolutionEndcap[i] = new TH1F(saResoNew.str().c_str(),";pt resolution;Events",100,-1,1);
    SAEtaEEResolutionEndcap[i] = new TH1F(saResoEE.str().c_str(),";pt resolution;Events",100,-1,1);
    saResoDef.str("");saResoNew.str("");saResoEE.str("");
  }
  TH2 *SAPtVsPtDifferenceBarrel = new TH2F("SAPtVsPtDifferenceBarrel",";pt_{SA}(GeV);pt_{off}-pt_{SA}(GeV)",100,0,20,200,-10,10);
  TH2 *SAPtVsDeltaPtBarrel = new TH2F("SAPtVsDeltaPtBarrel",";pt_{SA}(GeV);deltaPt(GeV)",100,0,20,200,-10,10);
  TH2 *SAPtVsDeltaPtEndcap = new TH2F("SAPtVsDeltaPtEndcap",";pt_{SA}(GeV);deltaPt(GeV)",100,0,20,200,-10,10);
  TH2 *invSAPtVsInvPtDifferenceBarrel = new TH2F("invSAPtVsInvPtDifferenceBarrel",";1/pt_{SA};1/pt_{off}-1/pt_{SA}",250,0,0.25,400,-0.05,0.05);
  TH2 *badBinEtaPhiAll = new TH2F("badBinEtaPhiAll",";etaBin;phiBin",30,0,30,192,0,192);
  TH2 *badBinEtaPhi = new TH2F("badBinEtaPhi",";etaBin;phiBin",30,0,30,12,0,12);
  TH2 *goodBinEtaPhiAll = new TH2F("goodBinEtaPhiAll",";etaBin;phiBin",30,0,30,192,0,192);
  TH2 *goodBinEtaPhi = new TH2F("goodBinEtaPhi",";etaBin;phiBin",30,0,30,12,0,12);
  TH1 *ptRatioBad = new TH1F("ptRatioBad",";pt_{SA}/pt_{EE};Events",100,0,3);
  TH1 *ptRatioGood = new TH1F("ptRatioGood",";pt_{SA}/pt_{EE};Events",100,0,3);
  TH1 *ptRatioCut = new TH1F("ptRatioCut",";pt_{SA}/pt_{EE};Events",100,0,3);
  TH1 *ptBad = new TH1F("ptBad",";offline pt(GeV);Events",20,0,20);
  TH1 *ptGood = new TH1F("ptGood",";offline pt(GeV);Events",20,0,20);
  TH2 *ptRatioVsPtResolution = new TH2F("ptRatioVsPtResolution",";pt_{SA}/pt_{EE};pt_{EE} resolution",100,0,3,200,-5,5);
  TH1 *SAEEresoL6DefEENo = new TH1F("SAEEresoL6DefEENo",";reso;Events",200,-10,10);
  TH1 *SAEEresoS6DefEENo = new TH1F("SAEEresoS6DefEENo",";reso;Events",200,-10,10);
  TH1 *SAEEresoL4DefEENo = new TH1F("SAEEresoL4DefEENo",";reso;Events",200,-10,10);
  TH1 *SAEEresoS4DefEENo = new TH1F("SAEEresoS4DefEENo",";reso;Events",200,-10,10);
  TH2* EESAvsEEOffReso = new TH2F("EESAvsEEOffReso",";1-pt_{SA}/pt_{EE};1-pt{EE}/pt_{off}",400,-2,2,400,-2,2);
  TH2* invSAPtVsAlpha00270 = new TH2F("invSAPtVsAlpha00270",";1/pt_{SA};alpha",250,0,0.25,800,0,0.4);
  TH2* invPtVsAlpha00270 = new TH2F("invPtVsAlpha00270",";1/pt_{SA};alpha",250,0,0.25,800,0,0.4);
  TH1* ptResolutionDefault[12][2];
  TH1* ptResolutionNew[12][2];
  TH1* Distance = new TH1F("Distance",";distance(mm);",100,0,5000);
  TH1* DistanceBad = new TH1F("DistanceBad",";distance(mm);",100,0,5000);
  TH1* DistanceGood = new TH1F("DistanceGood",";distance(mm);",100,0,5000);
  TH1* DistanceCut = new TH1F("DistanceCut",";distance(mm);",100,0,5000);
  TH2* DistanceVsPtResolution = new TH2F("DistanceVsPtResolution",";distance(mm);pt_{EE} resolution",500,0,5000,200,-5,5);
  TH2* SegmentChamberVsR = new TH2F("SegmentChamberVsR",";segment chamber;R(mm)",16,0,16,100,0,12000);
  TH2* SegmentChamberVsZ = new TH2F("SegmentChamberVsZ",";segment chamber;Z(mm)",16,0,16,100,0,22000);
  TH1 *Rdiff = new TH1F("Rdiff",";(R_{off}-R_{SA})/R_{off};",200,-1,1);
  TH1 *RdiffBad = new TH1F("RdiffBad",";(R_{off}-R_{SA})/R_{off};",200,-1,1);
  TH2 *RdiffVsReso = new TH2F("RdiffVsReso",";(R_{off}-R_{SA})/R_{off};EE pt residual",200,-1,1,200,-1,1);
  TH1 *EEDistance = new TH1F("EEDistance",";distance(mm);",100,0,1000);
  TH1 *EEDistanceBad = new TH1F("EEDistanceBad",";distance(mm);",100,0,1000);
  TH1 *EIDistance = new TH1F("EIDistance",";distance(mm);",100,0,1000);
  TH1 *EIDistanceSmallBad = new TH1F("EIDistanceSmallBad",";distance(mm);",100,0,1000);
  TH1 *EIDistanceLargeBad = new TH1F("EIDistanceLargeBad",";distance(mm);",100,0,1000);
  TH1 *EMDistance = new TH1F("EMDistance",";distance(mm);",100,0,1000);
  TH1 *EMDistanceBad = new TH1F("EMDistanceBad",";distance(mm);",100,0,1000);
  TH1 *BIDistance = new TH1F("BIDistance",";distance(mm);",100,0,1000);
  TH1 *BIDistanceBad = new TH1F("BIDistanceBad",";distance(mm);",100,0,1000);
  TH1 *probePt = new TH1F("probePt",";p_{T,off}(GeV);Number of events",100,0,50);
  TH1 *DiMuonMass = new TH1F("DiMuonMass",";M_{#mu#mu}(GeV);Number of events",100,2.5,3.5);
  TH2 *EtaPhiMu4UnderThr = new TH2F("EtaPhiMu4UnderThr",";#eta;#phi",100,-2.5,2.5,100,-3.2,3.2);
  TH2 *EtaPhiMu6UnderThr = new TH2F("EtaPhiMu6UnderThr",";#eta;#phi",100,-2.5,2.5,100,-3.2,3.2);
  TH1 *SAPtResolutionMu4UnderThr = new TH1F("SAPtResolutionMu4UnderThr",";p_{T} residual;",100,-1,1);
  TH1 *SAPtResolutionMu6UnderThr = new TH1F("SAPtResolutionMu6UnderThr",";p_{T} residual;",100,-1,1);
  TH1 *SAPtResolutionMu4AboveThr = new TH1F("SAPtResolutionMu4AboveThr",";p_{T} residual;",100,-1,1);
  TH1 *SAPtResolutionMu6AboveThr = new TH1F("SAPtResolutionMu6AboveThr",";p_{T} residual;",100,-1,1);
  stringstream test1,test2;
  for (int i=0;i<12;i++){
    for (int j=0;j<2;j++){
      test1 << "ptResolutionDefaultPhi" << i << "Eta27Qeta" << j ;
      test2 << "ptResolutionNewPhi" << i << "Eta27Qeta" << j ;
      ptResolutionDefault[i][j] = new TH1F(test1.str().c_str(),";pt resolution;Events",100,-1,1);
      ptResolutionNew[i][j] = new TH1F(test2.str().c_str(),";pt resolution;Events",100,-1,1);
      test1.str("");test2.str("");
    }
  }

  int aho1=0,aho2=0,aho3,aho4;
  int num2=0,num2_false=0;
  int num3=0,num3_false=0;
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000000==0) cout << "jetnry=" << jentry << endl;
    fChain->GetEntry(jentry); 
    //float tag_pt = tag_offline_pt;
    float offline_pt = probe_offline_pt; 
    float offline_eta = probe_offline_eta; 
    float offline_phi = probe_offline_phi; 
    float offline_charge = probe_offline_charge; 
    bool charge_correct=false;
    if (offline_charge==probe_sa_charge) charge_correct=true;
    probePt->Fill(offline_pt);
    //float offline_pt2 = tag_offline_pt; 
    //float offline_eta2 = tag_offline_eta; 
    //float offline_phi2 = tag_offline_phi; 
    //float offline_charge2 = tag_offline_charge; 
    float muon_mass = 0.1056; 
    TLorentzVector muon1,muon2;
    muon1.SetPtEtaPhiM(offline_pt,offline_eta,offline_phi,muon_mass);
    //muon2.SetPtEtaPhiM(offline_pt2,offline_eta2,offline_phi2,muon_mass);
    //float mass = (muon1+muon2).M();
    //DiMuonMass->Fill(mass);
    float sa_pt = fabs(probe_sa_pt);
    float sa_pt_alpha = fabs(probe_sa_pt_alpha);
    float sa_pt_beta = fabs(probe_sa_pt_beta);
    float sa_pt_tgc = fabs(probe_sa_pt_tgc);
    float sa_pt_new = fabs(probe_sa_pt_new);
    float sa_eta = probe_sa_eta;
    float sa_phi = probe_sa_phi;
    float sa_pt_ee = fabs(probe_sa_pt_ee);
    float deltaPt = probe_sa_deltaPt;
    float PtDifference = (sa_pt>1e-5)? offline_pt-sa_pt : 1000;
    float invPtDifference = (sa_pt>1e-5)? 1/offline_pt-1/sa_pt : 1000;
    float saPtResolution = (sa_pt>1e-5)? 1-offline_pt/sa_pt : 1000;
    float saPtNewResolution = (sa_pt_new>1e-5)? 1-offline_pt/sa_pt_new : 1000;
    float saPtEEResolution = (sa_pt_ee>1e-5)? 1-offline_pt/sa_pt_ee : 1000;
    bool passL1mu4 = probe_pass_L1mu4;
    bool passSAmu4 = probe_pass_SAmu4;
    bool passSAmu4_new = probe_pass_SAmu4_new;
    bool passSAmu4_ee = probe_pass_SAmu4_ee;
    bool passL1mu6 = probe_pass_L1mu6;
    bool passSAmu6 = probe_pass_SAmu6;
    bool passSAmu6_new = probe_pass_SAmu6_new;
    bool passSAmu6_ee = probe_pass_SAmu6_ee;
    int saddress = probe_sa_saddress;
    int phibin = GetPhiBinNumber(sa_phi);
    int phibinall = GetPhiBinAllNumber(sa_phi);
    int etabin = GetEtaBinNumber(sa_eta);

    float SP_ec_inner_z=0,SP_ec_inner_r=0,SP_ec_ee_z=0,SP_ec_ee_r=0;
    float SP_br_inner_z=0,SP_br_inner_r=0,SP_ec_middle_z=0,SP_ec_middle_r=0;
    /*if (sp_z->size()){
      SP_ec_inner_z = sp_z->at(3);
      SP_br_inner_z = sp_z->at(0);
      SP_ec_ee_z = sp_z->at(6);
      SP_ec_middle_z = sp_z->at(4);
      SP_ec_inner_r = sp_r->at(3);
      SP_br_inner_r = sp_r->at(0);
      SP_ec_ee_r = sp_r->at(6);
      SP_ec_middle_r = sp_r->at(4);
    }*/
    bool small = isSmall(SP_ec_ee_z);

    if (fabs(offline_eta)<1.0){//barrel
      if (passSAmu6){
        for (int ipt=0;ipt<15;ipt++){
          if (offline_pt>ipt && offline_pt<ipt+1){
            numPassSA[0][ipt]++;
            if (charge_correct) numChargeCorrect[0][ipt]++;
          }
        }
      }
      if (PtDifference<999) SAPtVsPtDifferenceBarrel->Fill(sa_pt,PtDifference);
      if (sa_pt>1e-5) {
        SAPtVsDeltaPtBarrel->Fill(sa_pt,deltaPt);
        invSAPtVsInvPtDifferenceBarrel->Fill(1/sa_pt,invPtDifference);
      }
      for (int ipt=0;ipt<60;ipt++){
        if (offline_pt>ipt*0.25 && offline_pt<(ipt+1)*0.25){
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
      }//pt loop end
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

      if (saPtResolution<999){
        for (int i=0;i<20;i++)
          if (offline_pt>i && offline_pt<i+1)
            SAPtResolutionBarrel[i]->Fill(saPtResolution);
      }
    }//Barrel end
    else {//endcap
      if (passSAmu6){
        for (int ipt=0;ipt<15;ipt++){
          if (offline_pt>ipt && offline_pt<ipt+1){
            numPassSA[1][ipt]++;
            if (charge_correct) numChargeCorrect[1][ipt]++;
          }
        }
      }
      //if (fabs(offline_eta)>1.4) continue;//EE region
      //if (offline_eta>0) continue;
      //if (fabs(sa_pt-sa_pt_ee)<1e-5) continue;
      if (passL1mu4 && passSAmu4 && offline_pt<3) {
        if(saPtResolution>0.8 && saPtResolution<999) EtaPhiMu4UnderThr->Fill(offline_eta,offline_phi);
        if(saPtResolution<999) SAPtResolutionMu4UnderThr->Fill(saPtResolution);
        //cout << "passL1mu4! passL1mu6=" << passL1mu6 << endl;
      }
      else if (passL1mu4 && passSAmu4 && offline_pt>3)
        if(saPtResolution<999) SAPtResolutionMu4AboveThr->Fill(saPtResolution);
      
      if (passL1mu6 && passSAmu6 && offline_pt<4) {
        if(saPtResolution>0.8 && saPtResolution<999) EtaPhiMu6UnderThr->Fill(offline_eta,offline_phi);
        if(saPtResolution<999) SAPtResolutionMu6UnderThr->Fill(saPtResolution);
        //cout << "passL1mu6! passL1mu4=" << passL1mu4 << endl;
        if (saPtResolution>0.8){
          //cout << "**********************" << endl;
          //cout << "offline pt=" << offline_pt << endl;
          //cout << "pt sa/alpha/beta/tgc=" << sa_pt << "/" << sa_pt_alpha << "/" << sa_pt_beta << "/" << sa_pt_tgc << endl; 
          //cout << "charge offline/sa=" << offline_charge << "/" << probe_sa_charge << endl;
          num2++;
          if(offline_charge!=probe_sa_charge) num2_false++;
        }
      }
      else if (passL1mu6 && passSAmu6 && offline_pt>4){ 
        if(saPtResolution<999) SAPtResolutionMu6AboveThr->Fill(saPtResolution);
        num3++;
        if(offline_charge!=probe_sa_charge) num3_false++;

      }
      
      //int nSegments = segment_r->size();
      float EIS_segr=0,EIL_segr=0,BIS_segr=0,BIL_segr=0;
      float EIS_segz=0,EIL_segz=0,BIS_segz=0,BIL_segz=0;
      float EES_segr=0,EEL_segr=0,EMS_segr=0,EML_segr=0;
      float EES_segz=0,EEL_segz=0,EMS_segz=0,EML_segz=0;
      //cout << "****************" << endl;
      /*for (int seg=0;seg<nSegments;seg++){
        float segr = segment_r->at(seg);
        float segz = segment_z->at(seg);
        int segch = segment_chamber->at(seg);
        SegmentChamberVsR->Fill(segch,segr);
        SegmentChamberVsZ->Fill(segch,segz);
        switch(segch){
          case 0:
            BIS_segr=segr;
            BIS_segz=segz;
            break;
          case 1:
            BIL_segr=segr;
            BIL_segz=segz;
            break;
          case 7:
            EIS_segr=segr;
            EIS_segz=segz;
            break;
          case 8:
            EIL_segr=segr;
            EIL_segz=segz;
            break;
          case 13:
            EES_segr=segr;
            EES_segz=segz;
            break;
          case 14:
            EEL_segr=segr;
            EEL_segz=segz;
            break;
          case 9:
            EMS_segr=segr;
            EMS_segz=segz;
            break;
          case 10:
            EML_segr=segr;
            EML_segz=segz;
            break;
        }
      }*/
      float calc_ec_radius=0,distance=0,offline_radius=0;
      if (BIS_segr>1e-5 && EES_segr>1e-5 && EMS_segr>1e-5) 
        offline_radius = computeRadius3Points(BIS_segz,BIS_segr,EES_segz,EES_segr,EMS_segz,EMS_segr);
      if (BIL_segr>1e-5 && EEL_segr>1e-5 && EML_segr>1e-5) 
        offline_radius = computeRadius3Points(BIL_segz,BIL_segr,EEL_segz,EEL_segr,EML_segz,EML_segr);
      if (EIS_segr>1e-5 && EES_segr>1e-5 && EMS_segr>1e-5) 
        offline_radius = computeRadius3Points(EIS_segz,EIS_segr,EES_segz,EES_segr,EMS_segz,EMS_segr);
      if (EIL_segr>1e-5 && EEL_segr>1e-5 && EML_segr>1e-5) 
        offline_radius = computeRadius3Points(EIL_segz,EIL_segr,EEL_segz,EEL_segr,EML_segz,EML_segr);
      if(small){
        if (SP_br_inner_r>1e-5 && SP_ec_ee_r>1e-5 && SP_ec_middle_r>1e-5){
          calc_ec_radius = computeRadius3Points(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          distance = calcDistance(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
        }
      }
      else{
        if (SP_ec_inner_r>1e-5 && SP_ec_ee_r>1e-5 && SP_ec_middle_r>1e-5){
          calc_ec_radius = computeRadius3Points(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          distance = calcDistance(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
        }
      }
      if (fabs(sa_pt_ee-sa_pt)<1e-5 && calc_ec_radius>1e-5 ) 
        badBinEtaPhiAll->Fill(probe_sa_ec_eta_bin,phibinall);
      if (calc_ec_radius>1e-5 && offline_radius>1e-5){
        float r_diff = (offline_radius-calc_ec_radius)/offline_radius;
        Rdiff->Fill(r_diff);
        if(saPtEEResolution<999) RdiffVsReso->Fill(r_diff,saPtEEResolution);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) RdiffBad->Fill(r_diff);
      }
      if (EES_segr>1e-5 && SP_ec_ee_r>1e-5){
        float ee_dis = calcDis(EES_segr,EES_segz,SP_ec_ee_r,SP_ec_ee_z);
        EEDistance->Fill(ee_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) EEDistanceBad->Fill(ee_dis);
      }
      if (EEL_segr>1e-5 && SP_ec_ee_r>1e-5){
        float ee_dis = calcDis(EEL_segr,EEL_segz,SP_ec_ee_r,SP_ec_ee_z);
        EEDistance->Fill(ee_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) EEDistanceBad->Fill(ee_dis);
      }
      if (EIS_segr>1e-5 && SP_ec_inner_r>1e-5){
        float inner_dis = calcDis(EIS_segr,EIS_segz,SP_ec_inner_r,SP_ec_inner_z);
        EIDistance->Fill(inner_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) EIDistanceSmallBad->Fill(inner_dis);
      }
      if (EIL_segr>1e-5 && SP_ec_inner_r>1e-5){
        float inner_dis = calcDis(EIL_segr,EIL_segz,SP_ec_inner_r,SP_ec_inner_z);
        EIDistance->Fill(inner_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) {
          EIDistanceLargeBad->Fill(inner_dis);
          //cout << "********************" << endl;
          //cout << "EI SA r/z=" << SP_ec_inner_r << "/" << SP_ec_inner_z << endl;
          //cout << "EI offline r/z=" << EIL_segr << "/" << EIL_segz << endl;
          //cout << "EI distance=" << inner_dis << endl;
          //cout << "alpha/beta=" << probe_sa_alpha << "/" << probe_sa_beta << endl; 
          //cout << "pt alpha/beta/tgc=" << probe_sa_pt_alpha << "/" << probe_sa_pt_beta << "/" << probe_sa_pt_tgc << endl;
        }
      }
      if (EMS_segr>1e-5 && SP_ec_middle_r>1e-5){
        float middle_dis = calcDis(EMS_segr,EMS_segz,SP_ec_middle_r,SP_ec_middle_z);
        EMDistance->Fill(middle_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) EMDistanceBad->Fill(middle_dis);
      }
      if (EML_segr>1e-5 && SP_ec_middle_r>1e-5){
        float middle_dis = calcDis(EML_segr,EML_segz,SP_ec_middle_r,SP_ec_middle_z);
        EMDistance->Fill(middle_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) EMDistanceBad->Fill(middle_dis);
      }
      if (BIS_segr>1e-5 && SP_br_inner_r>1e-5){
        float br_inner_dis = calcDis(BIS_segr,BIS_segz,SP_br_inner_r,SP_br_inner_z);
        BIDistance->Fill(br_inner_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) BIDistanceBad->Fill(br_inner_dis);
      }
      if (BIL_segr>1e-5 && SP_br_inner_r>1e-5){
        float br_inner_dis = calcDis(BIL_segr,BIL_segz,SP_br_inner_r,SP_br_inner_z);
        BIDistance->Fill(br_inner_dis);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) BIDistanceBad->Fill(br_inner_dis);
      }

      if (distance>1e-5){
        DistanceVsPtResolution->Fill(distance,saPtEEResolution);
        Distance->Fill(distance);
        if (passSAmu6 && !passSAmu6_ee && offline_pt<6) DistanceCut->Fill(distance);
        if (!passSAmu6 && passSAmu6_ee && offline_pt>6) DistanceGood->Fill(distance);
        if (passSAmu6 && !passSAmu6_ee && offline_pt>6) {
          DistanceBad->Fill(distance);
        }
      }
      if (!passSAmu6 && passSAmu6_new && offline_pt<6){//check endcap LUT
        if (fabs(sa_pt_new-probe_sa_pt_alpha_new)<1e-5){
          if (probe_sa_ec_eta_bin==13 && probe_sa_ec_phi_bin==11 /*&& probe_sa_charge*sa_eta<0*/){
            /*cout << "*****************" << endl;
              cout << "pt offline/sa=" << offline_pt << endl;
              cout << "offline eta/phi/charge=" << offline_eta << "/" << offline_phi << "/" << offline_charge << endl;
              cout << "sa eta/phi/charge=" << sa_eta << "/" << sa_phi << "/" << probe_sa_charge << endl;
              cout << "etabin/phibin/invradius=" << probe_sa_ec_eta_bin << "/" << probe_sa_ec_phi_bin << "/" << 1./probe_sa_ec_radius << endl; 
              cout << "alpha/beta=" << probe_sa_alpha << "/" << probe_sa_beta << endl;
              cout << "sa pt default/new=" << sa_pt << "/" << sa_pt_new << endl;
              cout << "pt alpha/beta/tgc=" << probe_sa_pt_alpha << "/" << probe_sa_pt_beta << "/" << probe_sa_pt_tgc << endl;
              cout << "new pt alpha/beta/tgc=" << probe_sa_pt_alpha_new << "/" << probe_sa_pt_beta_new << "/" << probe_sa_pt_tgc_new << endl;
              */
            //badBinEtaPhi->Fill(probe_sa_ec_eta_bin,probe_sa_ec_phi_bin);
            invSAPtVsAlpha00270->Fill(1./sa_pt,probe_sa_alpha);
            invPtVsAlpha00270->Fill(1./offline_pt,probe_sa_alpha);
          }
        }
      }
      for (int bin=0;bin<12;bin++){
        if (probe_sa_ec_eta_bin==27 && probe_sa_ec_phi_bin==bin){
          int qeta=-1;
          if (sa_eta*probe_sa_charge<0) qeta=0;
          else qeta=1;
          if (saPtResolution<999) ptResolutionDefault[bin][qeta]->Fill(saPtResolution);
          if (saPtNewResolution<999) ptResolutionNew[bin][qeta]->Fill(saPtNewResolution);
        }
      }
      if (passSAmu6 && !passSAmu6_new && offline_pt<6) {
        if (fabs(sa_pt_new-probe_sa_pt_alpha_new)<1e-5){
          //goodBinEtaPhi->Fill(probe_sa_ec_eta_bin,probe_sa_ec_phi_bin);
        }
      }
      float par=1-sa_pt/sa_pt_ee;
      float par2=1-offline_pt/sa_pt_ee;
      if (offline_pt<10 && sa_pt_ee>1e-5 && fabs(sa_pt-sa_pt_ee)>1e-5) EESAvsEEOffReso->Fill(par,par2);
      if (passSAmu6 && !passSAmu6_ee ) {
        if (probe_sa_ec_eta_bin==2 && phibinall==141){
          //cout << "etabin/phibin=" << probe_sa_ec_eta_bin << "/" << phibinall << endl;
          //cout << "SP EE R/Z=" << SP_ec_ee_r << "/" << SP_ec_ee_z << endl; 
        }
      }
      if (passSAmu6 && !passSAmu6_ee && offline_pt>6){//check EE LUT
        badBinEtaPhi->Fill(probe_sa_ec_eta_bin,probe_sa_ec_phi_bin);
        ptRatioBad->Fill(sa_pt/sa_pt_ee);
        SAEEresoL6DefEENo->Fill(par);
        if (sa_pt/sa_pt_ee<1.5) ptBad->Fill(offline_pt);
        //cout << "*************" << endl;
        //cout << "offline pt/eta/phi=" << offline_pt << "/" << offline_eta << "/" << offline_phi << endl;
        //cout << "SP BI R/Z=" << SP_br_inner_r << "/" << SP_br_inner_z << endl; 
        //cout << "SP EI R/Z=" << SP_ec_inner_r << "/" << SP_ec_inner_z << endl; 
        //cout << "SP EE R/Z=" << SP_ec_ee_r << "/" << SP_ec_ee_z << endl; 
        //cout << "SP EM R/Z=" << SP_ec_middle_r << "/" << SP_ec_middle_z << endl; 
      }
      if (!passSAmu6 && passSAmu6_ee && offline_pt>6) {
        ptRatioGood->Fill(sa_pt/sa_pt_ee);
        if (sa_pt/sa_pt_ee<1.5) ptGood->Fill(offline_pt);
      }
      if (passSAmu6 && !passSAmu6_ee && offline_pt<6) {
        SAEEresoS6DefEENo->Fill(par);
        ptRatioCut->Fill(sa_pt/sa_pt_ee);
        goodBinEtaPhi->Fill(probe_sa_ec_eta_bin,probe_sa_ec_phi_bin);
        goodBinEtaPhiAll->Fill(probe_sa_ec_eta_bin,phibinall);
      }
      if (passSAmu4 && !passSAmu4_ee && offline_pt>4) SAEEresoL4DefEENo->Fill(par);
      if (passSAmu4 && !passSAmu4_ee && offline_pt<4) SAEEresoS4DefEENo->Fill(par);
      if (fabs(sa_eta)>1e-5){//MuonSA exist
        for (int ieta=0;ieta<60;ieta++){
          if (offline_eta>-1.5+0.05*ieta && offline_eta<-1.45+0.05*ieta){
            numEERegionEta[ieta]++;
            if (fabs(sa_pt-sa_pt_ee)>1e-5) numUseEEEta[ieta]++;
          }
        }
        for (int iphi=0;iphi<64;iphi++){
          if (fabs(offline_eta)>1 && fabs(offline_eta)<1.35){
            if (offline_phi>-3.2+0.1*iphi && offline_phi<-3.1+0.1*iphi){
              numEERegionPhi[iphi]++;
              if (fabs(sa_pt-sa_pt_ee)>1e-5) numUseEEPhi[iphi]++;
            }
          }
        }
        for (int ipt=0;ipt<60;ipt++){
          if (offline_pt>0.25*ipt && offline_pt<0.25*(ipt+1)){
            numEERegionPt[ipt]++;
            if (fabs(sa_pt-sa_pt_ee)>1e-5) numUseEEPt[ipt]++;
          }
        }
      }

      if (saPtResolution<999){
        for (int i=0;i<20;i++)
          if (offline_pt>i && offline_pt<i+1)
            SAPtResolutionEndcap[i]->Fill(saPtResolution);
      }
      if (saPtEEResolution<999){
        for (int i=0;i<20;i++)
          if (offline_pt>i && offline_pt<i+1)
            SAPtEEResolutionEndcap[i]->Fill(saPtEEResolution);
        if (fabs(sa_pt-sa_pt_ee)>1e-5) ptRatioVsPtResolution->Fill(sa_pt/sa_pt_ee,saPtEEResolution);
      }

      if (saPtNewResolution<999){
        for (int i=0;i<20;i++)
          if (offline_pt>i && offline_pt<i+1)
            SAPtNewResolutionEndcap[i]->Fill(saPtNewResolution);
      }
      if (sa_pt>1e-5) SAPtVsDeltaPtEndcap->Fill(sa_pt,deltaPt);

      for (int ipt=0;ipt<60;ipt++){
        if (offline_pt>ipt*0.25 && offline_pt<(ipt+1)*0.25){
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
      }//pt loop end

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

    }//Endcap end

    if (saPtResolution<999){
      for (int i=0;i<100;i++)
        if (offline_eta>-2.5+0.05*i && offline_eta<-2.45+0.05*i)
          SAEtaResolutionEndcap[i]->Fill(saPtResolution);
    }
    if (saPtNewResolution<999){
      for (int i=0;i<100;i++)
        if (offline_eta>-2.5+0.05*i && offline_eta<-2.45+0.05*i)
          SAEtaNewResolutionEndcap[i]->Fill(saPtNewResolution);
    }
    if (saPtEEResolution<999){
      for (int i=0;i<100;i++)
        if (offline_eta>-2.5+0.05*i && offline_eta<-2.45+0.05*i)
          SAEtaEEResolutionEndcap[i]->Fill(saPtEEResolution);
    }


    if (fabs(offline_pt)>6){
      for (int ieta=0;ieta<50;ieta++){
        if (offline_eta>-2.5+ieta*0.1 && offline_eta<-2.4+ieta*0.1){
          numProbeEta[ieta]++;
          if (passL1mu4) {
            numEta[0][0][0][ieta]++; 
            numEta[0][0][1][ieta]++; 
            numEta[0][0][2][ieta]++; 
          }
          if (passL1mu6) {
            numEta[0][1][0][ieta]++; 
            numEta[0][1][1][ieta]++; 
            numEta[0][1][2][ieta]++; 
          }
          if (passSAmu4) numEta[1][0][0][ieta]++; 
          if (passSAmu4_new) numEta[1][0][1][ieta]++; 
          if (passSAmu4_ee) numEta[1][0][2][ieta]++; 
          if (passSAmu6) numEta[1][1][0][ieta]++; 
          if (passSAmu6_new) numEta[1][1][1][ieta]++; 
          if (passSAmu6_ee) numEta[1][1][2][ieta]++; 
        }
      }//eta loop end
    }

    for (int ieta=0;ieta<50;ieta++){
      for (int iphi=0;iphi<64;iphi++){
        if (offline_eta>-2.5+ieta*0.1 && offline_eta<-2.4+ieta*0.1){
          if (offline_phi>-3.2+iphi*0.1 && offline_phi<-3.1+iphi*0.1){
            numProbeEtaPhi[ieta][iphi]++;
            if (passL1mu4) {
              numEtaPhi[0][0][0][ieta][iphi]++; 
              numEtaPhi[0][0][1][ieta][iphi]++; 
              numEtaPhi[0][0][2][ieta][iphi]++; 
            }
            if (passL1mu6) {
              numEtaPhi[0][1][0][ieta][iphi]++; 
              numEtaPhi[0][1][1][ieta][iphi]++; 
              numEtaPhi[0][1][2][ieta][iphi]++; 
            }
            if (passSAmu4) numEtaPhi[1][0][0][ieta][iphi]++; 
            if (passSAmu4_new) numEtaPhi[1][0][1][ieta][iphi]++; 
            if (passSAmu4_ee) numEtaPhi[1][0][2][ieta][iphi]++; 
            if (passSAmu6) numEtaPhi[1][1][0][ieta][iphi]++; 
            if (passSAmu6_new) numEtaPhi[1][1][1][ieta][iphi]++; 
            if (passSAmu6_ee) numEtaPhi[1][1][2][ieta][iphi]++; 
          }
        }
      }
    }//eta-phi loop end

  }//Event loop end

  string LEVEL[2] = {"L1","SA"};
  string CHAIN[2] = {"Mu4","Mu6"};
  string SIDE[2] = {"Barrel","Endcap"};
  string TYPE[3] = {"Default","New","EE"};
  stringstream effName; 

  TGraphErrors* EffPt[2][2][2][3];
  for (int level=0;level<2;level++){//L1 or SA
    for (int chain=0;chain<2;chain++){//mu4 or mu6
      for (int side=0;side<2;side++){//barrel or endcap
        for (int type=0;type<3;type++){//default, new or ee LUT
          for (int pt=0;pt<60;pt++){
            if (level==0){//L1 efficiency
              if (numProbePt[side][pt]>0){ 
                eff_pt[level][chain][side][type][pt] = 
                  1.*numPt[level][chain][side][type][pt]/numProbePt[side][pt];
                efferr_pt[level][chain][side][type][pt] = 
                  sqrt(eff_pt[level][chain][side][type][pt]*(1-eff_pt[level][chain][side][type][pt]) / numProbePt[side][pt]);
              }
            }
            if (level==1){//SA efficiency
              /*if (numPt[0][chain][side][type][pt]>0){ 
                eff_pt[level][chain][side][type][pt] = 
                  1.*numPt[level][chain][side][type][pt]/numPt[0][chain][side][type][pt];
                efferr_pt[level][chain][side][type][pt] = 
                  sqrt(eff_pt[level][chain][side][type][pt]*(1-eff_pt[level][chain][side][type][pt]) / numPt[0][chain][side][type][pt]);
              }*/
              if (numProbePt[side][pt]>0){ 
                eff_pt[level][chain][side][type][pt] = 
                  1.*numPt[level][chain][side][type][pt]/numProbePt[side][pt];
                efferr_pt[level][chain][side][type][pt] = 
                  sqrt(eff_pt[level][chain][side][type][pt]*(1-eff_pt[level][chain][side][type][pt]) / numProbePt[side][pt]);
              }
            }
          }
          effName << "EfficiencyPt" << LEVEL[level] << CHAIN[chain] << SIDE[side] << TYPE[type];
          EffPt[level][chain][side][type] = new TGraphErrors(60,xbin_pt,eff_pt[level][chain][side][type],xbinerr_pt,efferr_pt[level][chain][side][type]);
          EffPt[level][chain][side][type]->SetName(effName.str().c_str());
          EffPt[level][chain][side][type]->SetTitle("");
          EffPt[level][chain][side][type]->GetXaxis()->SetLabelSize(0.07);
          EffPt[level][chain][side][type]->GetYaxis()->SetLabelSize(0.07);
          EffPt[level][chain][side][type]->GetYaxis()->SetRangeUser(0,1.2);
          file->Add(EffPt[level][chain][side][type]);
          effName.str("");
        }
      }
    }
  }

  TGraphErrors* EffEta[2][2][3];
  for (int level=0;level<2;level++){//L1 or SA
    for (int chain=0;chain<2;chain++){//mu4 or mu6
      for (int type=0;type<3;type++){//default, new or ee LUT
        for (int eta=0;eta<50;eta++){
          if (level==0){//L1 efficiency
            if (numProbeEta[eta]>0){ 
              eff_eta[level][chain][type][eta] = 
                1.*numEta[level][chain][type][eta]/numProbeEta[eta];
              efferr_eta[level][chain][type][eta] = 
                sqrt(eff_eta[level][chain][type][eta]*(1-eff_eta[level][chain][type][eta]) / numProbeEta[eta]);
            }
          }
          if (level==1){//SA efficiency
            if (numEta[0][chain][type][eta]>0){ 
              eff_eta[level][chain][type][eta] = 
                1.*numEta[level][chain][type][eta]/numEta[0][chain][type][eta];
              efferr_eta[level][chain][type][eta] = 
                sqrt(eff_eta[level][chain][type][eta]*(1-eff_eta[level][chain][type][eta]) / numEta[0][chain][type][eta]);
            }
          }
        }
        effName << "EfficiencyEta" << LEVEL[level] << CHAIN[chain] << TYPE[type];
        EffEta[level][chain][type] = new TGraphErrors(50,xbin_eta,eff_eta[level][chain][type],xbinerr_eta,efferr_eta[level][chain][type]);
        EffEta[level][chain][type]->SetName(effName.str().c_str());
        EffEta[level][chain][type]->SetTitle("");
        EffEta[level][chain][type]->GetXaxis()->SetLabelSize(0.07);
        EffEta[level][chain][type]->GetYaxis()->SetLabelSize(0.07);
        EffEta[level][chain][type]->GetYaxis()->SetRangeUser(0,1.2);
        file->Add(EffEta[level][chain][type]);
        effName.str("");
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
            /*cout << CHAIN[chain] << "/" << TYPE[type] << endl;
              cout << "eta/phi=" << eta << "/" << phi << "/" << endl; 
              cout << "probe/L1/SA=" << numProbeEtaPhi[eta][phi] << "/" << numEtaPhi[0][chain][type][eta][phi] << "/" 
              << numEtaPhi[1][chain][type][eta][phi] << endl; 
              cout << "efficiency L1/Probe is "  << eff_eta_phi[0][chain][type][eta][phi] << endl;
              cout << "efficiency SA/L1 is "  << eff_eta_phi[1][chain][type][eta][phi] << endl;*/
            EffEtaPhi[level][chain][type]->Fill(-2.45+0.1*eta,-3.15+0.1*phi,eff_eta_phi[level][chain][type][eta][phi]);
            //EffEtaPhi[level][chain][type]->Fill(-2.45+0.1*eta,-3.15+0.1*phi);
          }
        }
        //file->Add(EffEtaPhi[level][chain][type]);
        effName.str("");
      }
    }
  }

  TGraphErrors* EffUseEEEta;
  for (int eta=0;eta<60;eta++){
    if (numEERegionEta[eta]>0){ 
      eff_use_ee_eta[eta] = 
        1.*numUseEEEta[eta]/numEERegionEta[eta];
      efferr_use_ee_eta[eta] = 
        sqrt(eff_use_ee_eta[eta]*(1-eff_use_ee_eta[eta]) / numEERegionEta[eta]);
      //cout << "ieta=" << eta << endl;
      //cout << "num eeregion/useee=" << numEERegionEta[eta] << "/" << numUseEEEta[eta] << endl;
      //cout << "efficiency=" << eff_use_ee_eta[eta] << " +- " << efferr_use_ee_eta[eta] << endl;
    }
  }
  EffUseEEEta = new TGraphErrors(60,xbin_use_ee_eta,eff_use_ee_eta,xbinerr_use_ee_eta,efferr_use_ee_eta);
  EffUseEEEta->SetName("EffUseEEEta");
  EffUseEEEta->SetTitle("");
  EffUseEEEta->GetXaxis()->SetLabelSize(0.07);
  EffUseEEEta->GetYaxis()->SetLabelSize(0.07);
  EffUseEEEta->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(EffUseEEEta);

  TGraphErrors* EffUseEEPhi;
  for (int phi=0;phi<64;phi++){
    if (numEERegionPhi[phi]>0){ 
      eff_use_ee_phi[phi] = 
        1.*numUseEEPhi[phi]/numEERegionPhi[phi];
      efferr_use_ee_phi[phi] = 
        sqrt(eff_use_ee_phi[phi]*(1-eff_use_ee_phi[phi]) / numEERegionPhi[phi]);
      //cout << "iphi=" << phi << endl;
      //cout << "num eeregion/useee=" << numEERegionPhi[phi] << "/" << numUseEEPhi[phi] << endl;
      //cout << "efficiency=" << eff_use_ee_phi[phi] << " +- " << efferr_use_ee_phi[phi] << endl;
    }
  }
  EffUseEEPhi = new TGraphErrors(64,xbin_use_ee_phi,eff_use_ee_phi,xbinerr_use_ee_phi,efferr_use_ee_phi);
  EffUseEEPhi->SetName("EffUseEEPhi");
  EffUseEEPhi->SetTitle("");
  EffUseEEPhi->GetXaxis()->SetLabelSize(0.07);
  EffUseEEPhi->GetYaxis()->SetLabelSize(0.07);
  EffUseEEPhi->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(EffUseEEPhi);

  TGraphErrors* EffUseEEPt;
  for (int pt=0;pt<60;pt++){
    if (numEERegionPt[pt]>0){ 
      eff_use_ee_pt[pt] = 
        1.*numUseEEPt[pt]/numEERegionPt[pt];
      efferr_use_ee_pt[pt] = 
        sqrt(eff_use_ee_pt[pt]*(1-eff_use_ee_pt[pt]) / numEERegionPt[pt]);
      //cout << "ipt=" << pt << endl;
      //cout << "num eeregion/useee=" << numEERegionPt[pt] << "/" << numUseEEPt[pt] << endl;
      //cout << "efficiency=" << eff_use_ee_pt[pt] << " +- " << efferr_use_ee_pt[pt] << endl;
    }
  }
  EffUseEEPt = new TGraphErrors(60,xbin_use_ee_pt,eff_use_ee_pt,xbinerr_use_ee_pt,efferr_use_ee_pt);
  EffUseEEPt->SetName("EffUseEEPt");
  EffUseEEPt->SetTitle("");
  EffUseEEPt->GetXaxis()->SetLabelSize(0.07);
  EffUseEEPt->GetYaxis()->SetLabelSize(0.07);
  EffUseEEPt->GetYaxis()->SetRangeUser(0,1.2);
  file->Add(EffUseEEPt);

  TGraphErrors* EffChargeCorrect[2];
  for (int side=0;side<2;side++){//barrel or endcap
    for (int pt=0;pt<15;pt++){
      if (numPassSA[side][pt]>0){ 
        eff_charge_correct[side][pt] = 
          1.*numChargeCorrect[side][pt]/numPassSA[side][pt];
        efferr_charge_correct[side][pt] = 
          sqrt(eff_charge_correct[side][pt]*(1-eff_charge_correct[side][pt]) / numPassSA[side][pt]);
      }
    }
    effName << "EfficiencyChargeCorrect" << SIDE[side];
    EffChargeCorrect[side] = new TGraphErrors(15,xbin_charge,eff_charge_correct[side],xbinerr_charge,efferr_charge_correct[side]);
    EffChargeCorrect[side]->SetName(effName.str().c_str());
    EffChargeCorrect[side]->SetTitle("");
    EffChargeCorrect[side]->GetXaxis()->SetLabelSize(0.07);
    EffChargeCorrect[side]->GetYaxis()->SetLabelSize(0.07);
    EffChargeCorrect[side]->GetYaxis()->SetRangeUser(0,1.2);
    file->Add(EffChargeCorrect[side]);
    effName.str("");
  }

  cout << "num2/false=" << num2 << "/" << num2_false << endl;
  cout << "num3/false=" << num3 << "/" << num3_false << endl;

  file->Write();
}

