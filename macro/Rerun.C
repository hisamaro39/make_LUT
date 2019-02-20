#define Rerun_cxx
#include "macro/Rerun.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include <sstream>

#define PI 3.14159265258979

double const PI_OVER_4 = PI/4.0;
double const PI_OVER_8 = PI/8.0;
double const PHI_RANGE = 12.0/PI_OVER_8;
double const ZERO_LIMIT = 10e-5;

pair<int,int> GetBinNumber(float eta, float phi){
  int Octant = (int)(phi/PI_OVER_4);
  double PhiInOctant = fabs(phi - Octant*PI_OVER_4);
  if(PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  int phiBin = static_cast<int>(PhiInOctant*PHI_RANGE);
  int etaBin = static_cast<int>((fabs(eta)-1.)/0.05);
  if (etaBin == -1) etaBin =  0;
  if (etaBin == 30) etaBin = 29;
  if (etaBin < -0.5 || etaBin > 29.5 || phiBin < -0.5 || phiBin > 11.5) return make_pair(-1,-1);
  return make_pair(etaBin,phiBin);
}

pair<int,int> GetBinNumberBarrel(float eta, float phi, int chamber){
  double EtaMin[4] = {-1.145, -1.150, -1.050, -1.050};
  double PhiMin[4] = {-0.230, -0.230, -0.181, -0.181};
  double EtaMax[4] = {1.145, 1.150, 1.050, 1.050};
  double PhiMax[4] = {0.230, 0.230, 0.181, 0.181};
  double EtaStep[4],PhiStep[4];
  for (int i=0;i<4;i++){
    EtaStep[i] = (EtaMax[i]-EtaMin[i])/30;
    PhiStep[i] = (PhiMax[i]-PhiMin[i])/30;

  }
  int etabin = (int)((eta-EtaMin[chamber])/EtaStep[chamber]);
  int phibin = (int)((phi-PhiMin[chamber])/PhiStep[chamber]);
  if(etabin<=-1) etabin = 0;
  if(etabin>=30) etabin = 29;
  if(phibin<=-1) phibin = 0;
  if(phibin>=30) phibin = 29;
  return make_pair(etabin,phibin);
}

int GetBinNumberAllPhi(float phi){
  if (phi<0) phi = phi + 2*PI;
  float phibin = (int) (phi * 96/PI);

  return phibin;
}

bool isBadPhi(int phibinall){
  //int badNum[18] = {67,68,69,70,76,77,78,113,114,115,118,119,120,167,168,172,173,174};//v1
  int badNum[19] = {68,69,70,71,75,76,77,78,113,114,115,118,119,120,167,168,169,173,174};//v2
  bool ans=false;
  for (int i=0;i<19;i++){
    if (phibinall==badNum[i])
      ans = true;
  }
  return ans;
}

float whichSAPT(float alphapt,float betapt,float tgcpt,float innerspz,float middlespz,float outerspz, float tgcmid1z, float beta){

  bool useMDT=false;
  if (fabs(middlespz)>1e-5 && fabs(outerspz)>1e-5) useMDT=true;
  else{
    if (tgcpt>=8 || fabs(tgcmid1z)<1e-5){
      if(fabs(middlespz)>1e-5) useMDT=true;
    }
  }

  if (!useMDT) return tgcpt;
  
  float l2pt=alphapt;
  if (alphapt>10 && beta>1e-5){
    float ratio1 = fabs(betapt-alphapt)/alphapt;
    float ratio2 = fabs(tgcpt-alphapt)/alphapt;
    float ratio3 = fabs(tgcpt-betapt)/betapt;
    if (ratio1<0.5) l2pt=betapt;
    else if (fabs(outerspz)<1e-5){
      if (betapt>alphapt || ratio2>ratio3) l2pt=betapt;
    }
  }

  if (fabs(l2pt)<1e-5) l2pt = tgcpt;

  return l2pt;
}

float calcIntercept(float r1,float z1,float r2,float z2){
  float slope = (r1-r2)/(z1-z2);
  float intercept = r1 - slope*z1;
  return intercept;
}

float calcAlpha(float r1,float z1,float r2,float z2){
  float slope1 = r1/z1;
  float slope2 = (r2 - r1)/(z2 - z1);
  float alpha = 0;
  alpha = fabs(atan(slope1) - atan(slope2));
  return alpha;
}

void Rerun::Loop()
{
  //version8
  //TFile *file = new TFile("outputRerun/default/data/defaultLUT_version8.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version8.root","recreate");
  //TFile *file = new TFile("outputRerun/default/data/defaultLUT_version8_noee.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version8_ee.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version8_ee_allphi.root","recreate");
  //TFile *file = new TFile("outputRerun/default/data/defaultLUT_version8_alpha_combine.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version8_alpha_combine.root","recreate");
  //TFile *file = new TFile("outputRerun/default/data/defaultLUT_version8_alpha_combine_again.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version8_alpha_combine_again.root","recreate");
  //TFile *file = new TFile("outputRerun/default/data/defaultLUT_version9.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version9.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version9_ratio.root","recreate");
  TFile *file = new TFile("outputRerun/default/data/defaultLUT_version9_separate.root","recreate");
  //TFile *file = new TFile("outputRerun/new/data/newLUT_version9_separate.root","recreate");
  //TFile *file = new TFile("outputRerun/aho.root","recreate");
  
  bool seeResolution = true;

  bool useDefaultLUT=true;
  bool useEE=false;
  int phinumber=12,etanumber=30;
  //endcap
  TH1 *ptEndcapRadiusResolution = new TH1F("ptEndcapRadiusResolution",";ptEndcapRadiusResolution;Events",100,-1,1);
  TH1 *ptEndcapRadiusRerunResolution = new TH1F("ptEndcapRadiusRerunResolution",";ptEndcapRadiusResolution;Events",100,-1,1);
  TH2 *ptAlphaCheck = new TH2F("ptAlphaCheck",";calculated pt(GeV);raw pt(GeV)",500,0,500,500,0,500);
  TH2 *ptTgcAlphaCheck = new TH2F("ptTgcAlphaCheck",";calculated pt(GeV);raw pt(GeV)",500,0,500,500,0,500);
  TH2 *ptBetaCheck = new TH2F("ptBetaCheck",";calculated pt(GeV);raw pt(GeV)",500,0,500,500,0,500);
  TH2 *ptBarrelCheck = new TH2F("ptBarrelCheck",";calculated pt(GeV);raw pt(GeV)",500,0,500,500,0,500);
  TH2 *ptECRadiusCheck = new TH2F("ptECRadiusCheck",";calculated pt(GeV);raw pt(GeV)",500,0,500,500,0,500);
  TH1 *pt08041 = new TH1F("pt08041",";pt;Events",60,0,60);
  TH1F *ptSAResolutionLow[phinumber][etanumber][2];
  TH1F *ptSAResolutionMiddle[phinumber][etanumber][2];
  TH1F *ptSAResolutionHigh[phinumber][etanumber][2];
  TH1F *ptAlphaResolutionLow[phinumber][etanumber][2];
  TH1F *ptAlphaResolutionMiddle[phinumber][etanumber][2];
  TH1F *ptAlphaResolutionHigh[phinumber][etanumber][2];
  TH1F *ptBetaResolutionLow[phinumber][etanumber][2];
  TH1F *ptBetaResolutionMiddle[phinumber][etanumber][2];
  TH1F *ptBetaResolutionHigh[phinumber][etanumber][2];
  TH1F *ptTgcResolutionLow[phinumber][etanumber][2];
  TH1F *ptTgcResolutionMiddle[phinumber][etanumber][2];
  TH1F *ptTgcResolutionHigh[phinumber][etanumber][2];
  TH1F *ptBarrelResolutionLow[2][4][30][30];//charge/chamber/etanum/phinum
  TH1F *ptBarrelResolutionMiddle[2][4][30][30];
  TH1F *ptBarrelResolutionHigh[2][4][30][30];
  //TH1F *ptEEResolutionLow[8][24][2][2];//etabin/phibin/qeta/sl
  //TH1F *ptEEResolutionMiddle[8][24][2][2];
  //TH1F *ptEEResolutionHigh[8][24][2][2];
  TH1F *ptEENoSLResolutionLow[phinumber][8][2];
  TH1F *ptEENoSLResolutionMiddle[phinumber][8][2];
  TH1F *ptEENoSLResolutionHigh[phinumber][8][2];
  TH1F *ptEEAllPhiNoQetaResolutionLow[192][8][2][2];
  TH1F *ptEEAllPhiNoQetaResolutionMiddle[192][8][2][2];
  TH1F *ptEEAllPhiNoQetaResolutionHigh[192][8][2][2];
  TH1F *ptSAAllPhiNoQetaResolutionLow[192][8][2][2];
  TH1F *ptSAAllPhiNoQetaResolutionMiddle[192][8][2][2];
  TH1F *ptSAAllPhiNoQetaResolutionHigh[192][8][2][2];
  stringstream SAdrawnameLow,SAdrawnameMiddle,SAdrawnameHigh;
  stringstream AlphadrawnameLow,AlphadrawnameMiddle,AlphadrawnameHigh;
  stringstream BetadrawnameLow,BetadrawnameMiddle,BetadrawnameHigh;
  stringstream TgcdrawnameLow,TgcdrawnameMiddle,TgcdrawnameHigh;
  stringstream BarreldrawnameLow,BarreldrawnameMiddle,BarreldrawnameHigh;
  //stringstream EEdrawnameLow,EEdrawnameMiddle,EEdrawnameHigh;
  stringstream EENoSLdrawnameLow,EENoSLdrawnameMiddle,EENoSLdrawnameHigh;
  stringstream EEAllPhiNoQetadrawnameLow,EEAllPhiNoQetadrawnameMiddle,EEAllPhiNoQetadrawnameHigh;
  stringstream SAAllPhiNoQetadrawnameLow,SAAllPhiNoQetadrawnameMiddle,SAAllPhiNoQetadrawnameHigh;
  stringstream iphi_ec,ieta_ec;
  stringstream iphi_ee,ieta_ee;
  stringstream iphi_br,ieta_br;

  if (seeResolution){
    //endcap
    for(int i = 0;i<phinumber;i++){
      for(int j = 0;j<etanumber;j++){
        for(int k = 0;k<2;k++){//Qeta
          if(i < 10) iphi_ec << "0" << i;
          else iphi_ec << i;
          if(j < 9) ieta_ec << "0" << j+1;
          else ieta_ec << j+1;
          SAdrawnameLow << "ptSAResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          SAdrawnameMiddle << "ptSAResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          SAdrawnameHigh << "ptSAResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          AlphadrawnameLow << "ptAlphaResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          AlphadrawnameMiddle << "ptAlphaResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          AlphadrawnameHigh << "ptAlphaResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          BetadrawnameLow << "ptBetaResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          BetadrawnameMiddle << "ptBetaResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          BetadrawnameHigh << "ptBetaResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          TgcdrawnameLow << "ptTgcResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          TgcdrawnameMiddle << "ptTgcResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          TgcdrawnameHigh << "ptTgcResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          EENoSLdrawnameLow << "ptEENoSLResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          EENoSLdrawnameMiddle << "ptEENoSLResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          EENoSLdrawnameHigh << "ptEENoSLResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          ptSAResolutionLow[i][j][k] = new TH1F(SAdrawnameLow.str().c_str(),"",200,-1,1);
          ptSAResolutionMiddle[i][j][k] = new TH1F(SAdrawnameMiddle.str().c_str(),"",200,-1,1);
          ptSAResolutionHigh[i][j][k] = new TH1F(SAdrawnameHigh.str().c_str(),"",200,-1,1);
          ptAlphaResolutionLow[i][j][k] = new TH1F(AlphadrawnameLow.str().c_str(),"",200,-1,1);
          ptAlphaResolutionMiddle[i][j][k] = new TH1F(AlphadrawnameMiddle.str().c_str(),"",200,-1,1);
          ptAlphaResolutionHigh[i][j][k] = new TH1F(AlphadrawnameHigh.str().c_str(),"",200,-1,1);
          ptBetaResolutionLow[i][j][k] = new TH1F(BetadrawnameLow.str().c_str(),"",200,-1,1);
          ptBetaResolutionMiddle[i][j][k] = new TH1F(BetadrawnameMiddle.str().c_str(),"",200,-1,1);
          ptBetaResolutionHigh[i][j][k] = new TH1F(BetadrawnameHigh.str().c_str(),"",200,-1,1);
          ptTgcResolutionLow[i][j][k] = new TH1F(TgcdrawnameLow.str().c_str(),"",200,-1,1);
          ptTgcResolutionMiddle[i][j][k] = new TH1F(TgcdrawnameMiddle.str().c_str(),"",200,-1,1);
          ptTgcResolutionHigh[i][j][k] = new TH1F(TgcdrawnameHigh.str().c_str(),"",200,-1,1);
          if (j<8){
            ptEENoSLResolutionLow[i][j][k] = new TH1F(EENoSLdrawnameLow.str().c_str(),"",100,-1,1);
            ptEENoSLResolutionMiddle[i][j][k] = new TH1F(EENoSLdrawnameMiddle.str().c_str(),"",100,-1,1);
            ptEENoSLResolutionHigh[i][j][k] = new TH1F(EENoSLdrawnameHigh.str().c_str(),"",100,-1,1);
          }
          SAdrawnameLow.str("");
          SAdrawnameMiddle.str("");
          SAdrawnameHigh.str("");
          AlphadrawnameLow.str("");
          AlphadrawnameMiddle.str("");
          AlphadrawnameHigh.str("");
          BetadrawnameLow.str("");
          BetadrawnameMiddle.str("");
          BetadrawnameHigh.str("");
          TgcdrawnameLow.str("");
          TgcdrawnameMiddle.str("");
          TgcdrawnameHigh.str("");
          EENoSLdrawnameLow.str("");
          EENoSLdrawnameMiddle.str("");
          EENoSLdrawnameHigh.str("");
          iphi_ec.str("");
          ieta_ec.str("");
        }
      }
    }

    //EE all phi no qeta
    for(int i = 0;i<192;i++){
      for(int j = 0;j<8;j++){
        for(int k = 0;k<2;k++){//charge
          for(int l = 0;l<2;l++){//side
            if(i < 10) iphi_ec << "00" << i;
            else if (i < 100) iphi_ec << "0" << i;
            else iphi_ec << i;
            if(j < 9) ieta_ec << "0" << j+1;
            else ieta_ec << j+1;
            EEAllPhiNoQetadrawnameLow << "ptEEAllPhiNoQetaResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            EEAllPhiNoQetadrawnameMiddle << "ptEEAllPhiNoQetaResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            EEAllPhiNoQetadrawnameHigh << "ptEEAllPhiNoQetaResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            SAAllPhiNoQetadrawnameLow << "ptSAAllPhiNoQetaResolutionLow" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            SAAllPhiNoQetadrawnameMiddle << "ptSAAllPhiNoQetaResolutionMiddle" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            SAAllPhiNoQetadrawnameHigh << "ptSAAllPhiNoQetaResolutionHigh" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            ptEEAllPhiNoQetaResolutionLow[i][j][k][l] = new TH1F(EEAllPhiNoQetadrawnameLow.str().c_str(),"",400,-2,2);
            ptEEAllPhiNoQetaResolutionMiddle[i][j][k][l] = new TH1F(EEAllPhiNoQetadrawnameMiddle.str().c_str(),"",400,-2,2);
            ptEEAllPhiNoQetaResolutionHigh[i][j][k][l] = new TH1F(EEAllPhiNoQetadrawnameHigh.str().c_str(),"",400,-2,2);
            ptSAAllPhiNoQetaResolutionLow[i][j][k][l] = new TH1F(SAAllPhiNoQetadrawnameLow.str().c_str(),"",400,-2,2);
            ptSAAllPhiNoQetaResolutionMiddle[i][j][k][l] = new TH1F(SAAllPhiNoQetadrawnameMiddle.str().c_str(),"",400,-2,2);
            ptSAAllPhiNoQetaResolutionHigh[i][j][k][l] = new TH1F(SAAllPhiNoQetadrawnameHigh.str().c_str(),"",400,-2,2);
            EEAllPhiNoQetadrawnameLow.str("");
            EEAllPhiNoQetadrawnameMiddle.str("");
            EEAllPhiNoQetadrawnameHigh.str("");
            SAAllPhiNoQetadrawnameLow.str("");
            SAAllPhiNoQetadrawnameMiddle.str("");
            SAAllPhiNoQetadrawnameHigh.str("");
            ieta_ec.str("");
            iphi_ec.str("");
          }
        }
      }
    }

    //barrel
    for(int i = 0;i<2;i++){//charge
      for(int j = 0;j<4;j++){//chamber
        for(int k = 0;k<30;k++){//etabin
          for(int l = 0;l<30;l++){//phibin
            if(k < 10) ieta_br << "0" << k;
            else ieta_br << k;
            if(l < 10) iphi_br << "0" << l;
            else iphi_br << l;
            BarreldrawnameLow << "ptBarrelResolutionLow" << i << j << ieta_br.str().c_str() << iphi_br.str().c_str();
            BarreldrawnameMiddle << "ptBarrelResolutionMiddle" << i << j << ieta_br.str().c_str() << iphi_br.str().c_str();
            BarreldrawnameHigh << "ptBarrelResolutionHigh" << i << j << ieta_br.str().c_str() << iphi_br.str().c_str();
            ptBarrelResolutionLow[i][j][k][l] = new TH1F(BarreldrawnameLow.str().c_str(),"",100,-1,1);
            ptBarrelResolutionMiddle[i][j][k][l] = new TH1F(BarreldrawnameMiddle.str().c_str(),"",100,-1,1);
            ptBarrelResolutionHigh[i][j][k][l] = new TH1F(BarreldrawnameHigh.str().c_str(),"",100,-1,1);
            BarreldrawnameLow.str("");
            BarreldrawnameMiddle.str("");
            BarreldrawnameHigh.str("");
            iphi_br.str("");
            ieta_br.str("");
          }
        }
      }
    }

    /*
    //EE
    for(int i = 0;i<8;i++){//etabin
    for(int j = 0;j<24;j++){//phibin
    for(int k = 0;k<2;k++){//qeta
    for(int l = 0;l<2;l++){//sl
    if(i < 9) ieta_ee << "0" << i+1;
    else ieta_ee << i+1;
    if(j < 10) iphi_ee << "0" << j;
    else iphi_ee << j;
    EEdrawnameLow << "ptEEResolutionLow" << iphi_ee.str().c_str() << ieta_ee.str().c_str() << k << l;
    EEdrawnameMiddle << "ptEEResolutionMiddle" << iphi_ee.str().c_str() << ieta_ee.str().c_str() << k << l;
    EEdrawnameHigh << "ptEEResolutionHigh" << iphi_ee.str().c_str() << ieta_ee.str().c_str() << k << l;
    ptEEResolutionLow[i][j][k][l] = new TH1F(EEdrawnameLow.str().c_str(),"",100,-1,1);
    ptEEResolutionMiddle[i][j][k][l] = new TH1F(EEdrawnameMiddle.str().c_str(),"",100,-1,1);
    ptEEResolutionHigh[i][j][k][l] = new TH1F(EEdrawnameHigh.str().c_str(),"",100,-1,1);
    EEdrawnameLow.str("");
    EEdrawnameMiddle.str("");
    EEdrawnameHigh.str("");
    iphi_ee.str("");
    ieta_ee.str("");
    }
    }
    }
    }
    */
  }
  
  TH1 *etaBin = new TH1I("etaBin",";etaBin;Events",30,0,30);
  TH1 *endcapOfflinePT = new TH1F("endcapOfflinePT",";endcap offline pt(GeV);Events",100,0,100);
  TH2* invPtVsAlpha02101 = new TH2F("invPtVsAlpha02101",";1/pt;alpha",250,0,0.25,800,0,0.4);
  TH2* invPtVsAlpha00101 = new TH2F("invPtVsAlpha00101",";1/pt;alpha",250,0,0.25,800,0,0.4);
  TH2* invPtVsAlpha00021 = new TH2F("invPtVsAlpha00021",";1/pt;alpha",250,0,0.25,800,0,0.4);
  TH2* invPtVsTgcAlpha01061 = new TH2F("invPtVsTgcAlpha01061",";1/pt;tgcalpha",250,0,0.25,800,0,0.4);
  TH1* phiBinBad = new TH1F("phiBinBad",";phibin;Events",192,0,192);
  TH2* EESAvsEEOffReso = new TH2F("EESAvsEEOffReso",";1-pt_{SA}/pt_{EE};1-pt{EE}/pt_{off}",400,-2,2,400,-2,2);
  TH2* invPtVsAlpha01021 = new TH2F("invPtVsAlpha01021",";1/pt;alpha",250,0,0.25,800,0,0.4);
  TH2* SAEEPtRatioVsOfflinePt = new TH2F("SAEEPtRatioVsOfflinePt",";offline pt(GeV);pt_{SA}/pt_{EE}",200,0,20,200,0,2);
  TH2* SAEEPtRatioVsEEResolution = new TH2F("SAEEPtRatioVsEEResolution",";EE pt resolution;pt_{SA}/pt_{EE}",400,-10,10,200,0,2);

  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;

  int aho1=0,aho2=0,aho3=0;
  for (int jentry=0; jentry<nentries;jentry++) {
    fChain->GetEntry(jentry);   
    etaBin->Fill(etaBinEC);
    if (jentry%10000000==0) cout << "entry=" << jentry << endl;
    //if (jentry>1000) break;
    ptAlphaCheck->Fill(L2_pt_alpha_rerun,L2_pt_alpha);
    ptTgcAlphaCheck->Fill(fabs(tgcPt_rerun),fabs(tgcPt));
    ptBetaCheck->Fill(L2_pt_beta_rerun,L2_pt_beta);
    ptBarrelCheck->Fill(L2_pt_br_radius,L2_pt_br_radius_rerun);
    ptECRadiusCheck->Fill(L2_pt_ec_radius_rerun,L2_pt_ec_radius);
    if (L2_saddress<0) endcapOfflinePT->Fill(offline_pt);
    int qeta = (L2_eta*L2_charge>0)? 1 : 0;
    int icharge = (L2_charge>0)? 1 : 0;
    float usedPhi = (fabs(tgcMid1_phi)>1e-5)? tgcMid1_phi : L1_phi;
    int phibinall = GetBinNumberAllPhi(usedPhi);
    int iside = (L2_etaMap>0)? 1 : 0;
    bool bad = isBadPhi(phibinall);
    float tgc_intercept = (fabs(tgcMid1_z)>1e-5 && fabs(tgcMid2_z)>1e-5)? 
      calcIntercept(tgcMid1_r,tgcMid1_z,tgcMid2_r,tgcMid2_z) : 0.;
    int tgc_charge = (tgc_intercept*tgcMid2_z<0)? -1 : 1;
    int tgc_qeta = (tgc_charge*tgcMid1_eta<0)? 0 : 1;
    float tgc_alpha = (fabs(tgcMid1_z)>1e-5 && fabs(tgcMid2_z)>1e-5)?
      calcAlpha(tgcMid1_r,tgcMid1_z,tgcMid2_r,tgcMid2_z) : 0;
    
    if (L2_saddress==-1){//endcap
      float L2_pt_rerun = 
        whichSAPT(L2_pt_alpha_rerun,L2_pt_beta_rerun,tgcPt_rerun,SP_inner_z,SP_middle_z,SP_outer_z,tgcMid1_z,L2_ec_beta);
      float L2_pt_new = 
        whichSAPT(L2_pt_alpha_new,L2_pt_beta_new,tgcPt_new,SP_inner_z,SP_middle_z,SP_outer_z,tgcMid1_z,L2_ec_beta);
      if (useEE){
        if (L2_pt_ee_allphi_noqeta_rerun>1e-5 && L2_pt_ee_allphi_noqeta_rerun<499 && !bad) 
          L2_pt_rerun=L2_pt_ee_allphi_noqeta_rerun;
        //if (L2_pt_ec_radius_noSL_rerun>1e-5 && L2_pt_ec_radius_noSL_rerun<499 && !bad) 
        //L2_pt_rerun=L2_pt_ec_radius_noSL_rerun;
      }
      if (fabs(L2_pt_rerun-L2_pt)>1) {
        /*cout << "*******************" << endl;
        cout << "offline pt=" << offline_pt << endl;
        cout << "SA pt default/rerun/new=" << L2_pt << "/" << L2_pt_rerun << "/" << L2_pt_new << endl;
        cout << "default alpha/beta/tgc=" << L2_pt_alpha << "/" << L2_pt_beta << "/" << tgcPt << endl;
        cout << "rerun alpha/beta/tgc=" << L2_pt_alpha_rerun << "/" << L2_pt_beta_rerun << "/" << tgcPt_rerun << endl;
        cout << "new alpha/beta/tgc=" << L2_pt_alpha_new << "/" << L2_pt_beta_new << "/" << tgcPt_new << endl;
        */
      }
      float reso8=1000;
      if (useDefaultLUT){//default LUT
        if(L2_pt_rerun>ZERO_LIMIT) {
          reso8 = 1 - offline_pt/L2_pt_rerun;
          if (seeResolution){
            //if (L2_pt_rerun<20) ptSAResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            //else if (L2_pt_rerun>=20 && L2_pt_rerun<40) ptSAResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            //else if (L2_pt_rerun>=40 && L2_pt_rerun<60) ptSAResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            if (offline_pt<20) ptSAResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            else if (offline_pt>=20 && offline_pt<40) ptSAResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            else if (offline_pt>=40 && offline_pt<60) ptSAResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso8);
          }
        }
      }
      else {//new LUT
        if(L2_pt_new>ZERO_LIMIT) {
          reso8 = 1 - offline_pt/L2_pt_new;
          if (seeResolution){
            //if (L2_pt_new<20) ptSAResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            //else if (L2_pt_new>=20 && L2_pt_new<40) ptSAResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            //else if (L2_pt_new>=40 && L2_pt_new<60) ptSAResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            if (offline_pt<20) ptSAResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            else if (offline_pt>=20 && offline_pt<40) ptSAResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso8);
            else if (offline_pt>=40 && offline_pt<60) ptSAResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso8);
          }
        }
      }

      if (L2_ec_alpha>1e-5){
        float reso1=1000;
        if (useDefaultLUT){//default LUT
          if (L2_pt_alpha_rerun>ZERO_LIMIT) {
            reso1 = 1 - offline_pt/L2_pt_alpha_rerun;
            if (seeResolution){
              //if (L2_pt_alpha_rerun<20) ptAlphaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              //else if (L2_pt_alpha_rerun>=20 && L2_pt_alpha_rerun<40) ptAlphaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              //else if (L2_pt_alpha_rerun>=40 && L2_pt_alpha_rerun<60) ptAlphaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              if (offline_pt<20) ptAlphaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              else if (offline_pt>=20 && offline_pt<40) ptAlphaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              else if (offline_pt>=40 && offline_pt<60) ptAlphaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso1);
            }
          }
        }
        else{//new LUT
          if (L2_pt_alpha_new>ZERO_LIMIT) {
            reso1 = 1 - offline_pt/L2_pt_alpha_new;
            if (seeResolution){
              //if (L2_pt_alpha_new<20) ptAlphaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              //else if (L2_pt_alpha_new>=20 && L2_pt_alpha_new<40) ptAlphaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              //else if (L2_pt_alpha_new>=40 && L2_pt_alpha_new<60) ptAlphaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              if (offline_pt<20) ptAlphaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              else if (offline_pt>=20 && offline_pt<40) ptAlphaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso1);
              else if (offline_pt>=40 && offline_pt<60) ptAlphaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso1);
            }
          }
        }
      }

      if (L2_ec_beta>1e-5){
        float reso2=1000;
        if (useDefaultLUT){//default LUT
          if(L2_pt_beta_rerun>ZERO_LIMIT) {
            reso2 = 1 - offline_pt/L2_pt_beta_rerun;
            if (seeResolution){
              //if (L2_pt_beta_rerun<20) ptBetaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              //else if (L2_pt_beta_rerun>=20 && L2_pt_beta_rerun<40) ptBetaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              //else if (L2_pt_beta_rerun>=40 && L2_pt_beta_rerun<60) ptBetaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              if (offline_pt<20) ptBetaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              else if (offline_pt>=20 && offline_pt<40) ptBetaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              else if (offline_pt>=40 && offline_pt<60) ptBetaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso2);
            }
          }
        }
        else{//new LUT
          if(L2_pt_beta_new>ZERO_LIMIT) {
            reso2 = 1 - offline_pt/L2_pt_beta_new;
            if (seeResolution){
              //if (L2_pt_beta_new<20) ptBetaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              //else if (L2_pt_beta_new>=20 && L2_pt_beta_new<40) ptBetaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              //else if (L2_pt_beta_new>=40 && L2_pt_beta_new<60) ptBetaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              if (offline_pt<20) ptBetaResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              else if (offline_pt>=20 && offline_pt<40) ptBetaResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso2);
              else if (offline_pt>=40 && offline_pt<60) ptBetaResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso2);
            }
          }
        }
      }

      /*if (L2_pt_ec_radius_rerun>1e-5 && L2_pt_ec_radius_rerun<499){
        float reso3 = 1 - offline_pt/L2_pt_ec_radius_rerun;
        if (L2_pt_ec_radius_rerun<20) ptEEResolutionLow[etaBinEC][phiBinEE][qeta][L2_SL]->Fill(reso3);
        else if (L2_pt_ec_radius_rerun>=20 && L2_pt_ec_radius_rerun<40) ptEEResolutionMiddle[etaBinEC][phiBinEE][qeta][L2_SL]->Fill(reso3);
        else if (L2_pt_ec_radius_rerun>=40 && L2_pt_ec_radius_rerun<60) ptEEResolutionHigh[etaBinEC][phiBinEE][qeta][L2_SL]->Fill(reso3);
      }*/

      /*if (L2_pt_ec_radius_noSL_rerun>1e-5 && L2_pt_ec_radius_noSL_rerun<499){
        float reso7 = 1 - offline_pt/L2_pt_ec_radius_noSL_rerun;
        if (L2_pt_ec_radius_noSL_rerun<20) 
          ptEENoSLResolutionLow[phiBinEC][etaBinEC][qeta]->Fill(reso7);
        else if (L2_pt_ec_radius_noSL_rerun>=20 && L2_pt_ec_radius_noSL_rerun<40) 
          ptEENoSLResolutionMiddle[phiBinEC][etaBinEC][qeta]->Fill(reso7);
        else if (L2_pt_ec_radius_noSL_rerun>=40 && L2_pt_ec_radius_noSL_rerun<60) 
          ptEENoSLResolutionHigh[phiBinEC][etaBinEC][qeta]->Fill(reso7);
      }*/

      if (tgc_alpha>1e-5){//tgcalpha is calculated
        float reso6 = 1000.;
        if (useDefaultLUT){//default LUT
          if(tgcPt_rerun>1e-5) {
            if (seeResolution){
              reso6 = 1 - offline_pt/tgcPt_rerun;
              //if (tgcPt_rerun<20) ptTgcResolutionLow[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              //else if (tgcPt_rerun>=20 && tgcPt_rerun<40) ptTgcResolutionMiddle[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              //else if (tgcPt_rerun>=40 && tgcPt_rerun<60) ptTgcResolutionHigh[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              if (offline_pt<20) ptTgcResolutionLow[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              else if (offline_pt>=20 && offline_pt<40) ptTgcResolutionMiddle[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              else if (offline_pt>=40 && offline_pt<60) ptTgcResolutionHigh[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
            }
          }
        }
        else{//new LUT
          if(tgcPt_new>1e-5) {
            if (seeResolution){
              reso6 = 1 - offline_pt/tgcPt_new;
              //if (tgcPt_new<20) ptTgcResolutionLow[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              //else if (tgcPt_new>=20 && tgcPt_new<40) ptTgcResolutionMiddle[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              //else if (tgcPt_new>=40 && tgcPt_new<60) ptTgcResolutionHigh[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              if (offline_pt<20) ptTgcResolutionLow[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              else if (offline_pt>=20 && offline_pt<40) ptTgcResolutionMiddle[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
              else if (offline_pt>=40 && offline_pt<60) ptTgcResolutionHigh[phiBinTGC][etaBinTGC][tgc_qeta]->Fill(reso6);
            }
          }
        }
      }

      if (L2_ec_radius>1e-5 && etaBinEC<8){
        if (seeResolution){
          float reso9=1000;
          if (L2_pt_ee_allphi_noqeta_rerun>ZERO_LIMIT && L2_pt_ee_allphi_noqeta_rerun<499) {
            reso9 = 1 - offline_pt/L2_pt_ee_allphi_noqeta_rerun;
            if (offline_pt<20){
              ptEEAllPhiNoQetaResolutionLow[phibinall][etaBinEC][icharge][iside]->Fill(reso9);
              ptSAAllPhiNoQetaResolutionLow[phibinall][etaBinEC][icharge][iside]->Fill(reso8);
            }
            else if (offline_pt>=20 && offline_pt<40) {
              ptEEAllPhiNoQetaResolutionMiddle[phibinall][etaBinEC][icharge][iside]->Fill(reso9);
              ptSAAllPhiNoQetaResolutionMiddle[phibinall][etaBinEC][icharge][iside]->Fill(reso8);
            }
            else if (offline_pt>=40 && offline_pt<60) {
              ptEEAllPhiNoQetaResolutionHigh[phibinall][etaBinEC][icharge][iside]->Fill(reso9);
              ptSAAllPhiNoQetaResolutionHigh[phibinall][etaBinEC][icharge][iside]->Fill(reso8);
            }
            if (etaBinEC==4 && phiBinEC==10 && qeta==1 && !bad){
              if (reso9<-0.5)
                phiBinBad->Fill(phibinall);
            }
            float reso=1-L2_pt_rerun/L2_pt_ee_allphi_noqeta_rerun;
            EESAvsEEOffReso->Fill(reso,reso9);

          }
        }
        if (fabs(L2_pt_ee_allphi_noqeta_rerun-L2_pt_rerun)>1e-5){
          float ratio = L2_pt_ee_allphi_noqeta_rerun/L2_pt_rerun;
          float resolution = 1 - offline_pt/L2_pt_ee_allphi_noqeta_rerun;
          //cout << "offpt/ratio/resolution=" << offline_pt << "/" << ratio << "/" << resolution << endl;
          SAEEPtRatioVsOfflinePt->Fill(offline_pt,ratio);
          SAEEPtRatioVsEEResolution->Fill(resolution,ratio);
        }
      }
    }

    else{//Barrel
      if (L2_pt_br_radius_rerun>1e-5){
        pair<int,int> barrelBin = GetBinNumberBarrel(L2_etaMap,L2_phiMap,L2_saddress);
        int etabinBarrel = barrelBin.first;
        int phibinBarrel = barrelBin.second;
        float reso5 = 1000.;
        if (useDefaultLUT){//default LUT
          if (seeResolution){
            if (L2_pt_br_radius_rerun>1e-5) reso5 = 1 - offline_pt/L2_pt_br_radius_rerun;
            if (L2_pt_br_radius_rerun<20) 
              ptBarrelResolutionLow[icharge][L2_saddress][etabinBarrel][phibinBarrel]->Fill(reso5);
            else if (L2_pt_br_radius_rerun>=20 && L2_pt_br_radius_rerun<40) 
              ptBarrelResolutionMiddle[icharge][L2_saddress][etabinBarrel][phibinBarrel]->Fill(reso5);
            else if (L2_pt_br_radius_rerun>=40 && L2_pt_br_radius_rerun<60) 
              ptBarrelResolutionHigh[icharge][L2_saddress][etabinBarrel][phibinBarrel]->Fill(reso5);
          }
        }
        else{//new LUT
          if (seeResolution){
            if (L2_pt_br_radius_new>1e-5) reso5 = 1 - offline_pt/L2_pt_br_radius_new;
            if (L2_pt_br_radius_new<20) 
              ptBarrelResolutionLow[icharge][L2_saddress][etabinBarrel][phibinBarrel]->Fill(reso5);
            else if (L2_pt_br_radius_new>=20 && L2_pt_br_radius_new<40) 
              ptBarrelResolutionMiddle[icharge][L2_saddress][etabinBarrel][phibinBarrel]->Fill(reso5);
            else if (L2_pt_br_radius_new>=40 && L2_pt_br_radius_new<60) 
              ptBarrelResolutionHigh[icharge][L2_saddress][etabinBarrel][phibinBarrel]->Fill(reso5);
          }
        }
      }
    }

  }

  file->Write();
}

