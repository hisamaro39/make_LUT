#define validationT_cxx
#include "macro/validationT.h"
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
double const PI_OVER_16 = PI/16.0;
double const PHI_RANGE = 12.0/PI_OVER_8;
double const ZERO_LIMIT = 1e-5;

float calcCosDphi(float xi,float yi,float xm,float ym){
  float bunsi=fabs(ym-yi);
  float bunbo=sqrt(pow(xm-xi,2)+pow(ym-yi,2));
  float ans=bunsi/bunbo;
  return ans;

}

float calcSlope(float eta){
  float bunbo = exp(eta)-exp(-eta);
  float ans = 2/bunbo;
  return ans;
}

float calcDis(float r1,float z1,float r2,float z2){
  float dr = sqrt((r1-r2)*(r1-r2)+(z1-z2)*(z1-z2));
  return dr;
}

float calcDisR(float segr,float segz,float spr,float spz){
  float a = spr/spz;
  float r = a * spz;
  float dr=segr-r;
  return dr;
}

float calcRfromTrackZ(float slope, float z){
  float r = slope*z;
  return r;
}

float calcZfromTrackR(float slope, float r){
  float z = r/slope;
  return z;
}

bool isBadPhi(int phibinall){
  int badNum[18] = {67,68,69,70,76,77,78,113,114,115,118,119,120,167,168,172,173,174};
  bool ans=false;
  for (int i=0;i<18;i++){
    if (phibinall==badNum[i])
      ans = true;
  }
  return ans;
}

bool isSmall(float eez){
  bool small = false;
  if (fabs(eez)>10000 && fabs(eez)<10600) small =true;
  else if(fabs(eez)>10600 && fabs(eez)<12000) small = false;
  return small;
}

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

pair<float,float> calcCenter(float x1,float y1,float x2,float y2,float x3,float y3){
  float xm1=(x1+x2)/2;
  float xm2=(x2+x3)/2;
  float ym1=(y1+y2)/2;
  float ym2=(y2+y3)/2;
  float a1=(x1-x2)/(y2-y1);
  float a2=(x2-x3)/(y3-y2);
  float x0=(a2*xm2-a1*xm1-ym2+ym1)/(a2-a1);//center of circle
  float y0=a1*(x0-xm1)+ym1;//center of circle
  return make_pair(x0,y0);
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

int GetBinNumberEE(float phi,int sl){
  int phiBin24=-1;
  int Octant = (int)(phi/PI_OVER_4);
  double PhiInOctant = fabs(phi - Octant * PI_OVER_4);
  if (PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  if ( sl==0 ){//Small
    int OctantSmall = Octant;
    double PhiInOctantSmall = PhiInOctant;
    if(phi<0) PhiInOctantSmall = fabs(phi - (OctantSmall-1)*PI_OVER_4);
    phiBin24 = PhiInOctantSmall * PHI_RANGE;
  }
  else {//Large
    //phi = phi + PI_OVER_8;
    int OctantLarge = (int)(phi / PI_OVER_4);
    double PhiInOctantLarge = fabs(phi - OctantLarge * PI_OVER_4);
    if (phi<0) PhiInOctantLarge = fabs(phi - (OctantLarge-1)*PI_OVER_4);
    phiBin24 = PhiInOctantLarge * PHI_RANGE;
  }
  return phiBin24;
}
//////////////////////////////////

int GetBinNumberEE2(float phi){
  int phiBin=-1;
  float Phi = phi + PI_OVER_16;
  int doubleOctant = (int)(Phi/PI_OVER_8);
  double PhiInDoubleOctant = fabs(Phi - doubleOctant * PI_OVER_8);
  phiBin = PhiInDoubleOctant * PHI_RANGE;
  return phiBin;
}

int GetBinNumberAllPhi(float phi){
  if (phi<0) phi = phi + 2*PI;
  float phibin = (int) (phi * 96/PI);

  return phibin;
}

float calcAlpha(float r1,float z1,float r2,float z2){
  float slope1 = r1/z1;
  float slope2 = (r2 - r1)/(z2 - z1);
  float alpha = 0;
  alpha = fabs(atan(slope1) - atan(slope2));
  return alpha;
}
////////////////////
float calcIntercept(float r1,float z1,float r2,float z2){
  float slope = (r1-r2)/(z1-z2);
  float intercept = r1 - slope*z1;
  return intercept;
}

////////////////////////
double calcSagitta(double InnerZ, double InnerR,
    double EEZ, double EER,
    double MiddleZ, double MiddleR){
  double a = (MiddleZ - InnerZ)/(MiddleR - InnerR);
  double sagitta = EEZ - EER*a - InnerZ + InnerR*a;
  return sagitta;
}

///////////////////////////////
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


void validationT::Loop()
{
  //TFile *fout = new TFile("inputLUT/data15_13TeV_run2_ee.root","recreate");
  //TFile *fout = new TFile("inputLUT/data15_13TeV_version9.root","recreate");
  //TFile *fout = new TFile("inputLUT/data15_13TeV_version9_separate.root","recreate");
  //TFile *fout = new TFile("inputLUT/data15_13TeV_version9_radius.root","recreate");
  //TFile *fout = new TFile("inputLUT/data15_13TeV_version9_radius_linear.root","recreate");
  //TFile *fout = new TFile("inputLUT/mc15_13TeV_version9_zmumu.root","recreate");
  //TFile *fout = new TFile("inputLUT/mc15_13TeV_version9_zmumu_segment.root","recreate");
  //TFile *fout = new TFile("inputLUT/data16_13TeV_version9_segment.root","recreate");
  //TFile *fout = new TFile("inputLUT/data15_13TeV_correct_sagitta.root","recreate");
  TFile *fout = new TFile("inputLUT/mc15_13TeV_correct_sagitta_offline.root","recreate");
  //TFile *fout = new TFile("inputLUT/aho.root","recreate");
  bool makeLUT=false;
  bool correctSagitta=true;

  ///////////////set region ///////////////////
  float phiminimum = -3.15;
  float phimaximum = 3.15;

  int phibinmax_ec = 12;
  int phibinmin_ec = 0;
  int phinumber = phibinmax_ec - phibinmin_ec + 1;

  int etabinmax_ec = 30;
  int etabinmin_ec = 0;
  int etanumber = etabinmax_ec - etabinmin_ec + 1;

  int phibinmax_br = 30;
  int phibinmin_br = 0;

  int etabinmax_br = 30;
  int etabinmin_br = 0;

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  //gStyle->SetOptStat(0);

  ////////////////////////////////////////////
  //endcap
  TH2F *hMDT_invpTvsAlpha[phinumber][etanumber][2];
  TH2F *hMDT_invpTvsTgcAlpha[phinumber][etanumber][2];
  TH2F *hMDT_invpTvsBeta[phinumber][etanumber][2];
  TH2F *hMDT_invpTvsEndcapRadius[24][8][2][2];//EE
  TH2F *hMDT_invpTvsEndcapRadiusNoSL[12][8][2];//EE no SL
  TH2F *hMDT_invpTvsEndcapRadiusAllPhiNoQeta[192][8][2][2];//EE all phi no qeta
  TH2F *hMDT_invpTvsEndcapRadiusLinear[192][8][2][2];//EE R vs pt
  TH2F *hMDT_invpTvsAlphaAnormal[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_invpTvsBetaAnormal[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_invpTvsEndcapRadiusAnormal[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_pTvsEndcapSagittaAnormal[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_invpTvsEndcapRadius680200[10];//Check EE anormal region
  TH1F *hMDT_EEDistanceR680200[10];//Check EE anormal region
  TH1F *hMDT_EEDistanceEta680200[10];//Check EE anormal region
  TH1F *hMDT_EEDistance[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EIDistance[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EMDistance[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_BIDistance[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EEDistanceR[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EIDistanceR[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EMDistanceR[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_BIDistanceR[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EEDistanceEta[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EIDistanceEta[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EMDistanceEta[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_BIDistanceEta[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EEDistancePhi[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EIDistancePhi[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_EMDistancePhi[192][8][2][2];//Check EE anormal region
  TH1F *hMDT_BIDistancePhi[192][8][2][2];//Check EE anormal region
  stringstream Alphadrawname;
  stringstream TgcAlphadrawname;
  stringstream Betadrawname;
  stringstream EndcapRadiusdrawname;
  stringstream EndcapRadiusNoSLdrawname;
  stringstream EndcapRadiusAllPhiNoQetadrawname;
  stringstream EndcapRadiusLineardrawname;
  stringstream EndcapRadiusAnormaldrawname;
  stringstream EndcapRadius680200drawname;
  stringstream EEDistanceR680200drawname;
  stringstream EEDistanceEta680200drawname;
  stringstream EndcapSagittaAnormaldrawname;
  stringstream AlphaAnormaldrawname;
  stringstream BetaAnormaldrawname;
  stringstream EEDistancedrawname;
  stringstream EIDistancedrawname;
  stringstream EMDistancedrawname;
  stringstream BIDistancedrawname;
  stringstream EEDistanceRdrawname;
  stringstream EIDistanceRdrawname;
  stringstream EMDistanceRdrawname;
  stringstream BIDistanceRdrawname;
  stringstream EEDistanceEtadrawname;
  stringstream EIDistanceEtadrawname;
  stringstream EMDistanceEtadrawname;
  stringstream BIDistanceEtadrawname;
  stringstream EEDistancePhidrawname;
  stringstream EIDistancePhidrawname;
  stringstream EMDistancePhidrawname;
  stringstream BIDistancePhidrawname;
  stringstream iphi_ec,ieta_ec;
  stringstream iphi_ee,ieta_ee;
  //barrel
  TH2F *hMDT_invpTvsinvRadius[2][4][etabinmax_br][phibinmax_br];
  stringstream BarrelRadiusdrawname;
  stringstream iphi_br,ieta_br;

  if (makeLUT){ 
    //endcap
    for(int i = 0;i<phibinmax_ec;i++){
      for(int j = 0;j<etabinmax_ec;j++){
        for(int k = 0;k<2;k++){//Qeta
          if(i < 10) iphi_ec << "0" << i;
          else iphi_ec << i;
          if(j < 9) ieta_ec << "0" << j+1;
          else ieta_ec << j+1;
          Alphadrawname << "Alpha" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          TgcAlphadrawname << "TgcAlpha" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          Betadrawname << "Beta" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          EndcapRadiusNoSLdrawname << "EndcapRadiusNoSL" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k ;
          hMDT_invpTvsAlpha[i][j][k] = new TH2F(Alphadrawname.str().c_str(),"",250,0,0.25,800,0,0.4);
          hMDT_invpTvsBeta[i][j][k] = new TH2F(Betadrawname.str().c_str(),"",250,0,0.25,800,0,0.4);
          hMDT_invpTvsTgcAlpha[i][j][k] = new TH2F(TgcAlphadrawname.str().c_str(),"",250,0,0.25,800,0,0.4);
          hMDT_invpTvsAlpha[i][j][k]->SetMarkerStyle(8);
          hMDT_invpTvsTgcAlpha[i][j][k]->SetMarkerStyle(8);
          hMDT_invpTvsBeta[i][j][k]->SetMarkerStyle(8);
          hMDT_invpTvsAlpha[i][j][k]->GetXaxis()->SetRangeUser(0,0.25);
          hMDT_invpTvsAlpha[i][j][k]->GetYaxis()->SetRangeUser(0,0.4);
          hMDT_invpTvsTgcAlpha[i][j][k]->GetXaxis()->SetRangeUser(0,0.25);
          hMDT_invpTvsTgcAlpha[i][j][k]->GetYaxis()->SetRangeUser(0,0.4);
          hMDT_invpTvsBeta[i][j][k]->GetXaxis()->SetRangeUser(0,0.25);
          hMDT_invpTvsBeta[i][j][k]->GetYaxis()->SetRangeUser(0,0.4);
          if (j<8){
            hMDT_invpTvsEndcapRadiusNoSL[i][j][k] = new TH2F(EndcapRadiusNoSLdrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
            hMDT_invpTvsEndcapRadiusNoSL[i][j][k]->SetMarkerStyle(8);
            hMDT_invpTvsEndcapRadiusNoSL[i][j][k]->GetXaxis()->SetRangeUser(0,0.25);
            hMDT_invpTvsEndcapRadiusNoSL[i][j][k]->GetYaxis()->SetRangeUser(0,0.00004);
          }
          Alphadrawname.str("");
          TgcAlphadrawname.str("");
          Betadrawname.str("");
          EndcapRadiusNoSLdrawname.str("");
          iphi_ec.str("");
          ieta_ec.str("");
        }
      }
    }


    //ee
    for(int i = 0;i<24;i++){//phi
      for(int j = 0;j<8;j++){//eta
        for(int k = 0;k<2;k++){//Qeta
          for(int l=0; l<2; l++){//samall large
            if(i < 10) iphi_ee << "0" << i;
            else iphi_ee << i;
            if(j < 9) ieta_ee << "0" << j+1;
            else ieta_ee << j+1;
            EndcapRadiusdrawname << "EndcapRadius" << iphi_ee.str().c_str() << ieta_ee.str().c_str() << k << l;
            hMDT_invpTvsEndcapRadius[i][j][k][l] = new TH2F(EndcapRadiusdrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
            hMDT_invpTvsEndcapRadius[i][j][k][l]->SetMarkerStyle(8);
            hMDT_invpTvsEndcapRadius[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
            hMDT_invpTvsEndcapRadius[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.00004);
            EndcapRadiusdrawname.str("");
            iphi_ee.str("");
            ieta_ee.str("");
          }
        }
      }
    }

    //barrel
    for(int k=0; k<2; k++){//charge
      for (int chamber=0; chamber<4; chamber++){//saddress
        for(int i=0; i<etabinmax_br; i++){
          for(int j=0; j<phibinmax_br; j++){
            if(i < 10) ieta_br << "0" << i;
            else ieta_br << i;
            if(j < 10) iphi_br << "0" << j;
            else iphi_br << j;
            BarrelRadiusdrawname << "Radius" << k << chamber << ieta_br.str().c_str() << iphi_br.str().c_str();
            hMDT_invpTvsinvRadius[k][chamber][i][j] = new TH2F(BarrelRadiusdrawname.str().c_str(),"",500,0,50000,120,0,60);
            hMDT_invpTvsinvRadius[k][chamber][i][j]->SetMarkerStyle(8);
            hMDT_invpTvsinvRadius[k][chamber][i][j]->GetXaxis()->SetRangeUser(0,50000);
            hMDT_invpTvsinvRadius[k][chamber][i][j]->GetYaxis()->SetRangeUser(0,60);
            BarrelRadiusdrawname.str("");
            iphi_br.str("");
            ieta_br.str("");
          }
        }
      }
    }

    //ee all phi charge eta
    for(int i = 0;i<192;i++){
      for(int j = 0;j<8;j++){
        for(int k = 0;k<2;k++){//charge
          for (int l=0;l<2;l++){//eta
            if(i < 10) iphi_ec << "0" << i;
            else iphi_ec << i;
            if(j < 9) ieta_ec << "0" << j+1;
            else ieta_ec << j+1;
            EndcapRadiusAllPhiNoQetadrawname << "EndcapRadiusAllPhiNoQeta" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            EndcapRadiusLineardrawname << "EndcapRadiusLinear" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            hMDT_invpTvsEndcapRadiusAllPhiNoQeta[i][j][k][l] = new TH2F(EndcapRadiusAllPhiNoQetadrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
            hMDT_invpTvsEndcapRadiusAllPhiNoQeta[i][j][k][l]->SetMarkerStyle(8);
            hMDT_invpTvsEndcapRadiusAllPhiNoQeta[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
            hMDT_invpTvsEndcapRadiusAllPhiNoQeta[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.00004);
            EndcapRadiusAllPhiNoQetadrawname.str("");
            hMDT_invpTvsEndcapRadiusLinear[i][j][k][l] = new TH2F(EndcapRadiusLineardrawname.str().c_str(),"",120,0,60,500,0,10000000);
            hMDT_invpTvsEndcapRadiusLinear[i][j][k][l]->SetMarkerStyle(8);
            hMDT_invpTvsEndcapRadiusLinear[i][j][k][l]->GetXaxis()->SetRangeUser(0,60);
            hMDT_invpTvsEndcapRadiusLinear[i][j][k][l]->GetYaxis()->SetRangeUser(0,1000000);
            EndcapRadiusLineardrawname.str("");
            iphi_ec.str("");
            ieta_ec.str("");
          }
        }
      }
    }

  }//make LUT 

  //Search anormal EE region
  for(int i = 0;i<192;i++){
    for(int j = 0;j<8;j++){
      if(i < 10) iphi_ec << "0" << i;
      else iphi_ec << i;
      if(j < 9) ieta_ec << "0" << j+1;
      else ieta_ec << j+1;
      for (int l=0;l<2;l++){//eta
        for(int k = 0;k<2;k++){//charge
        EndcapSagittaAnormaldrawname << "EndcapSagittaAnormal" << iphi_ec.str().c_str() << ieta_ec.str().c_str() <<  k << l;
        hMDT_pTvsEndcapSagittaAnormal[i][j][k][l] = new TH2F(EndcapSagittaAnormaldrawname.str().c_str(),"",200,10,50,400,-200,200);
        hMDT_pTvsEndcapSagittaAnormal[i][j][k][l]->SetMarkerStyle(8);
        hMDT_pTvsEndcapSagittaAnormal[i][j][k][l]->GetXaxis()->SetRangeUser(10,50);
        hMDT_pTvsEndcapSagittaAnormal[i][j][k][l]->GetYaxis()->SetRangeUser(-200,200);
        EndcapSagittaAnormaldrawname.str("");
        //if (correctSagitta) continue;
          EndcapRadiusAnormaldrawname << "EndcapRadiusAnormal" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          AlphaAnormaldrawname << "AlphaAnormal" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          BetaAnormaldrawname << "BetaAnormal" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EEDistancedrawname << "EEdistance" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EIDistancedrawname << "EIdistance" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EMDistancedrawname << "EMdistance" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          BIDistancedrawname << "BIdistance" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EEDistanceRdrawname << "EEdistanceR" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EIDistanceRdrawname << "EIdistanceR" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EMDistanceRdrawname << "EMdistanceR" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          BIDistanceRdrawname << "BIdistanceR" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EEDistanceEtadrawname << "EEdistanceEta" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EIDistanceEtadrawname << "EIdistanceEta" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EMDistanceEtadrawname << "EMdistanceEta" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          BIDistanceEtadrawname << "BIdistanceEta" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EEDistancePhidrawname << "EEdistancePhi" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EIDistancePhidrawname << "EIdistancePhi" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          EMDistancePhidrawname << "EMdistancePhi" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          BIDistancePhidrawname << "BIdistancePhi" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
          hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l] = new TH2F(EndcapRadiusAnormaldrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
          hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l]->SetMarkerStyle(8);
          hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
          hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.00004);
          hMDT_invpTvsAlphaAnormal[i][j][k][l] = new TH2F(AlphaAnormaldrawname.str().c_str(),"",250,0,0.25,800,0,0.4);
          hMDT_invpTvsAlphaAnormal[i][j][k][l]->SetMarkerStyle(8);
          hMDT_invpTvsAlphaAnormal[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
          hMDT_invpTvsAlphaAnormal[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.4);
          hMDT_invpTvsBetaAnormal[i][j][k][l] = new TH2F(BetaAnormaldrawname.str().c_str(),"",250,0,0.25,800,0,0.4);
          hMDT_invpTvsBetaAnormal[i][j][k][l]->SetMarkerStyle(8);
          hMDT_invpTvsBetaAnormal[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
          hMDT_invpTvsBetaAnormal[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.4);
          hMDT_EEDistance[i][j][k][l] = new TH1F(EEDistancedrawname.str().c_str(),";distance(mm);",100,0,200);
          hMDT_EIDistance[i][j][k][l] = new TH1F(EIDistancedrawname.str().c_str(),";distance(mm);",100,0,200);
          hMDT_EMDistance[i][j][k][l] = new TH1F(EMDistancedrawname.str().c_str(),";distance(mm);",100,0,200);
          hMDT_BIDistance[i][j][k][l] = new TH1F(BIDistancedrawname.str().c_str(),";distance(mm);",100,0,200);
          hMDT_EEDistanceR[i][j][k][l] = new TH1F(EEDistanceRdrawname.str().c_str(),";distance(mm);",200,-200,200);
          hMDT_EIDistanceR[i][j][k][l] = new TH1F(EIDistanceRdrawname.str().c_str(),";distance(mm);",200,-200,200);
          hMDT_EMDistanceR[i][j][k][l] = new TH1F(EMDistanceRdrawname.str().c_str(),";distance(mm);",200,-200,200);
          hMDT_BIDistanceR[i][j][k][l] = new TH1F(BIDistanceRdrawname.str().c_str(),";distance(mm);",200,-200,200);
          hMDT_EEDistanceEta[i][j][k][l] = new TH1F(EEDistanceEtadrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_EIDistanceEta[i][j][k][l] = new TH1F(EIDistanceEtadrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_EMDistanceEta[i][j][k][l] = new TH1F(EMDistanceEtadrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_BIDistanceEta[i][j][k][l] = new TH1F(BIDistanceEtadrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_EEDistancePhi[i][j][k][l] = new TH1F(EEDistancePhidrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_EIDistancePhi[i][j][k][l] = new TH1F(EIDistancePhidrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_EMDistancePhi[i][j][k][l] = new TH1F(EMDistancePhidrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          hMDT_BIDistancePhi[i][j][k][l] = new TH1F(BIDistancePhidrawname.str().c_str(),";distance(mm);",200,-0.1,0.1);
          EndcapRadiusAnormaldrawname.str("");
          AlphaAnormaldrawname.str("");
          BetaAnormaldrawname.str("");
          EEDistancedrawname.str("");
          EIDistancedrawname.str("");
          EMDistancedrawname.str("");
          BIDistancedrawname.str("");
          EEDistanceRdrawname.str("");
          EIDistanceRdrawname.str("");
          EMDistanceRdrawname.str("");
          BIDistanceRdrawname.str("");
          EEDistanceEtadrawname.str("");
          EIDistanceEtadrawname.str("");
          EMDistanceEtadrawname.str("");
          BIDistanceEtadrawname.str("");
          EEDistancePhidrawname.str("");
          EIDistancePhidrawname.str("");
          EMDistancePhidrawname.str("");
          BIDistancePhidrawname.str("");
        }
      }
      iphi_ec.str("");
      ieta_ec.str("");
    }
  }

  //Search EE region 680200 shifted
  for(int i = 0;i<11;i++){
    int shift = 10*i;
    EndcapRadius680200drawname << "EndcapRadius680200_" << shift ; 
    EEDistanceR680200drawname << "EEDistanceR680200_" << shift ; 
    EEDistanceEta680200drawname << "EEDistanceEta680200_" << shift ; 
    hMDT_invpTvsEndcapRadius680200[i] = new TH2F(EndcapRadius680200drawname.str().c_str(),";;",250,0,0.25,400,0,0.00004);
    hMDT_EEDistanceR680200[i] = new TH1F(EEDistanceR680200drawname.str().c_str(),";distance(mm);",200,-200,200);
    hMDT_EEDistanceEta680200[i] = new TH1F(EEDistanceEta680200drawname.str().c_str(),";#Delta#eta;",200,-0.1,0.1);
    hMDT_invpTvsEndcapRadius680200[i]->SetMarkerStyle(8);
    hMDT_invpTvsEndcapRadius680200[i]->GetXaxis()->SetRangeUser(0,0.25);
    hMDT_invpTvsEndcapRadius680200[i]->GetYaxis()->SetRangeUser(0,0.00004);
    EndcapRadius680200drawname.str("");
    EEDistanceR680200drawname.str("");
    EEDistanceEta680200drawname.str("");
  }

  TH2 *phiVsInvR = new TH2F("phiVsEndcapRadius",";phi;1/R",200,-3.2,3.2,400,0,0.00004);
  TH2 *PhiVsEESPR = new TH2F("PhiVsEESPR",";#phi;EE SP R(mm)",200,-3.2,3.2,200,5000,10000);
  TH2 *PhiVsEESPZ = new TH2F("PhiVsEESPZ",";#phi;EE SP Z(mm)",200,-3.2,3.2,200,10000,12000);
  TH2 *EtaVsEESPR = new TH2F("EtaVsEESPR",";#eta;EE SP R(mm)",500,-2.5,2.5,200,5000,10000);
  TH2 *EtaVsEESPZ = new TH2F("EtaVsEESPZ",";#eta;EE SP Z(mm)",500,-2.5,2.5,200,10000,12000);
  TH2 *PhiVsEMSPR = new TH2F("PhiVsEMSPR",";#phi;EM SP R(mm)",200,-3.2,3.2,500,0,12000);
  TH2 *PhiVsEMSPZ = new TH2F("PhiVsEMSPZ",";#phi;EM SP Z(mm)",200,-3.2,3.2,500,12000,15000);
  TH1 *Distance = new TH1F("Distance",";distance(mm);",500,0,500);
  TH2 *DistanceVsPt = new TH2F("DistanceVsPt",";offline pt(GeV);distance(mm)",100,0,50,500,0,500);
  TH2 *DistanceVsPhi = new TH2F("DistanceVsPhi",";offline #phi;distance(mm)",64,-3.2,3.2,500,0,500);
  TH2 *EESegmentXY = new TH2F("EESegmentXY",";offline X(mm);offline Y(mm)",1000,-10000,10000,1000,-10000,10000);
  TH2 *EISegmentXY = new TH2F("EISegmentXY",";offline X(mm);offline Y(mm)",1000,-8000,8000,1000,-8000,8000);
  TH2 *EMSegmentXY = new TH2F("EMSegmentXY",";offline X(mm);offline Y(mm)",1000,-12000,12000,1000,-12000,12000);
  TH2 *EESPXY = new TH2F("EESPXY",";SP X(mm);SP Y(mm)",1000,-10000,10000,1000,-10000,10000);
  TH2 *EISPXY = new TH2F("EISPXY",";SP X(mm);SP Y(mm)",1000,-8000,8000,1000,-8000,8000);
  TH2 *EMSPXY = new TH2F("EMSPXY",";SP X(mm);SP Y(mm)",1000,-12000,12000,1000,-12000,12000);

  ////////////////////////////////////////
  float m_offline_pt,m_L2_pt,m_tgcPt;
  //float dinvptinvradius=0.;
  float m_L2_eta,m_L2_phi;
  int   m_L2_saddress;
  float m_L2_etaMap,m_L2_phiMap;
  float alpha,beta,ec_radius,br_radius,tgcalpha;
  float m_tgcMid1_r,m_tgcMid1_z,m_tgcMid2_r,m_tgcMid2_z;
  int m_L2_charge;

  ////////////////////////////////////////////////
  float EtaMin[4] = {-1.145, -1.150, -1.050, -1.050};
  float PhiMin[4] = {-0.230, -0.233, -0.181, -0.181};
  float EtaMax[4] = {1.145, 1.150, 1.050, 1.050};
  float PhiMax[4] = {0.230, 0.233, 0.181, 0.181};
  float EtaStep[4], PhiStep[4];
  for (int i=0; i<4; i++){
    EtaStep[i] = (EtaMax[i]-EtaMin[i])/30;
    PhiStep[i] = (PhiMax[i]-PhiMin[i])/30;
  }

  /////////////////////////////////////
  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;
  int aho=0;
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000000==0) cout << "entry=" << jentry << endl;
    //cout << "entry=" << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry); 

    m_offline_pt = offline_pt;
    alpha = L2_ec_alpha;
    beta = L2_ec_beta;
    br_radius = L2_br_radius/10.;//cm
    ec_radius = L2_ec_radius;//mm
    m_L2_pt = L2_pt;
    m_L2_eta = L2_eta;
    m_L2_phi = L2_phi;
    m_L2_saddress = L2_saddress;
    m_L2_etaMap = L2_etaMap;
    m_L2_phiMap = L2_phiMap;
    m_L2_charge = L2_charge;
    m_tgcMid1_r = tgcMid1_r;
    m_tgcMid1_z = tgcMid1_z;
    m_tgcMid2_r = tgcMid2_r;
    m_tgcMid2_z = tgcMid2_z;
    //m_tgcPt = tgcPt;
    tgcalpha = (fabs(m_tgcMid1_z)>1e-5 && fabs(m_tgcMid2_z)>1e-5) ? calcAlpha(m_tgcMid1_r,m_tgcMid1_z,m_tgcMid2_r,m_tgcMid2_z) : 0;
    bool isEndcap = (m_L2_saddress<-0.5)? true : false;
    float deta = L1_eta-offline_eta;
    float dphi = acos(cos(L1_phi-offline_phi));
    float dr = sqrt(deta*deta+dphi*dphi);

    ////////////////////////////////////////////////////////////////////////
    if (isEndcap){ //endcap
      bool isTgcFailure = (fabs(tgcMid1_phi)<1e-5)? true : false;
      float usedEta = (isTgcFailure)? L1_eta : tgcMid1_eta;
      float usedPhi = (isTgcFailure)? L1_phi : tgcMid1_phi;
      pair<int, int> LUTbinNumbers = GetBinNumber( usedPhi, usedEta );
      pair<int, int> LUTbinNumbersTgc = GetBinNumber( tgcMid1_phi, tgcMid1_eta );
      int PhibinNumber = LUTbinNumbers.second;
      int EtabinNumber = LUTbinNumbers.first;
      int PhibinNumberTgc = LUTbinNumbersTgc.second;
      int EtabinNumberTgc = LUTbinNumbersTgc.first;
      int PhibinNumberEE2 = GetBinNumberEE2(usedPhi);
      int i = PhibinNumber;
      int j = EtabinNumber;
      int PhibinAll  = GetBinNumberAllPhi(usedPhi);
      bool badPhibin = isBadPhi(PhibinAll);
      float tgc_intercept = (fabs(tgcMid1_z)>1e-5 && fabs(tgcMid2_z)>1e-5)? 
        calcIntercept(tgcMid1_r,tgcMid1_z,tgcMid2_r,tgcMid2_z) : 0.;
      int tgc_charge = (tgc_intercept*tgcMid2_z<0)? -1 : 1;
      int tgc_qeta = (tgc_charge*tgcMid1_eta<0)? 0 : 1;
      int new_L2_SL = -1;
      float slope = calcSlope(offline_eta);
      float SP_ec_inner_z = sp_z->at(3);
      float SP_br_inner_z = sp_z->at(0);
      float SP_ec_ee_z = sp_z->at(6);
      float SP_ec_middle_z = sp_z->at(4);
      float SP_ec_inner_r = sp_r->at(3);
      float SP_ec_inner_r_shift = sp_r->at(3)-100;
      float SP_br_inner_r = sp_r->at(0);
      float SP_ec_ee_r = sp_r->at(6);
      float SP_ec_middle_r = sp_r->at(4);
      float SP_ec_ee_r_shift[10];
      for (int ishi=0;ishi<11;ishi++)
        SP_ec_ee_r_shift[ishi] = sp_r->at(6)+10*ishi; 
      if(SP_ec_ee_r>1e-5) EESPXY->Fill(SP_ec_ee_r*cos(L2_phi),SP_ec_ee_r*sin(L2_phi));
      if(SP_ec_inner_r>1e-5) EISPXY->Fill(SP_ec_inner_r*cos(L2_phi),SP_ec_inner_r*sin(L2_phi));
      if(SP_ec_middle_r>1e-5) EMSPXY->Fill(SP_ec_middle_r*cos(L2_phi),SP_ec_middle_r*sin(L2_phi));

      
      float ee_dis=10000,inner_dis=10000,middle_dis=10000,br_inner_dis=10000;
      float ee_dis_r=10000,inner_dis_r=10000,middle_dis_r=10000,br_inner_dis_r=10000;
      float ee_dis_r_shift[10],inner_dis_r_shift=10000,ee_dis_eta_shift[10];
      float ee_dis_eta=10000,inner_dis_eta=10000,middle_dis_eta=10000,br_inner_dis_eta=10000;
      float ee_dis_phi=10000,inner_dis_phi=10000,middle_dis_phi=10000,br_inner_dis_phi=10000;
      float ee_dis_phi_ip=10000,inner_dis_phi_ip=10000,middle_dis_phi_ip=10000,br_inner_dis_phi_ip=10000;

      float calc_ec_sagitta_offline=0.;
      //if (!correctSagitta){
        int nSegments = segment_r->size();
        float EIS_segr=0,EIL_segr=0,BIS_segr=0,BIL_segr=0;
        float EIS_segeta=0,EIL_segeta=0,BIS_segeta=0,BIL_segeta=0;
        float EIS_segz=0,EIL_segz=0,BIS_segz=0,BIL_segz=0;
        float EIS_segphi=0,EIL_segphi=0,BIS_segphi=0,BIL_segphi=0;
        float EES_segr=0,EEL_segr=0,EMS_segr=0,EML_segr=0;
        float EES_segeta=0,EEL_segeta=0,EMS_segeta=0,EML_segeta=0;
        float EES_segz=0,EEL_segz=0,EMS_segz=0,EML_segz=0;
        float EES_segphi=0,EEL_segphi=0,EMS_segphi=0,EML_segphi=0;
        float EIL_segx=0,BIS_segx=0,EMS_segx=0,EML_segx=0;
        float EIL_segy=0,BIS_segy=0,EMS_segy=0,EML_segy=0;
        for (int seg=0;seg<nSegments;seg++){
          float segr = segment_r->at(seg);
          float segz = segment_z->at(seg);
          int segch = segment_chamber->at(seg);
          float segeta = segment_eta->at(seg);
          float segphi = segment_phi->at(seg);
          if(segch==13 || segch==14) EESegmentXY->Fill(segr*cos(segphi),segr*sin(segphi));
          if(segch==7 || segch==8) EISegmentXY->Fill(segr*cos(segphi),segr*sin(segphi));
          if(segch==9 || segch==10) EMSegmentXY->Fill(segr*cos(segphi),segr*sin(segphi));
          switch(segch){
            case 0:
              BIS_segr=segr;
              BIS_segz=segz;
              BIS_segeta=segeta;
              BIS_segphi=segphi;
              BIS_segx=segr*cos(segphi);
              BIS_segy=segr*sin(segphi);
              break;
            case 1:
              BIL_segr=segr;
              BIL_segz=segz;
              BIL_segeta=segeta;
              BIL_segphi=segphi;
              break;
            case 7:
              EIS_segr=segr;
              EIS_segz=segz;
              EIS_segeta=segeta;
              EIS_segphi=segphi;
              break;
            case 8:
              EIL_segr=segr;
              EIL_segz=segz;
              EIL_segeta=segeta;
              EIL_segphi=segphi;
              EIL_segx=segr*cos(segphi);
              EIL_segy=segr*sin(segphi);
              break;
            case 13:
              EES_segr=segr;
              EES_segz=segz;
              EES_segeta=segeta;
              EES_segphi=segphi;
              break;
            case 14:
              EEL_segr=segr;
              EEL_segz=segz;
              EEL_segeta=segeta;
              EEL_segphi=segphi;
              break;
            case 9:
              EMS_segr=segr;
              EMS_segz=segz;
              EMS_segeta=segeta;
              EMS_segphi=segphi;
              EMS_segx=segr*cos(segphi);
              EMS_segy=segr*sin(segphi);
              break;
            case 10:
              EML_segr=segr;
              EML_segz=segz;
              EML_segeta=segeta;
              EML_segphi=segphi;
              EML_segx=segr*cos(segphi);
              EML_segy=segr*sin(segphi);
              break;
          }
        }
        for (int ishi=0;ishi<11;ishi++) {
          ee_dis_r_shift[ishi]=10000;
          ee_dis_eta_shift[ishi]=10000;
        }
        if (EES_segr>1e-5 && SP_ec_ee_r>1e-5){
          for (int shi=0;shi<11;shi++){
            float EES_segr_shift=EES_segr + 10*shi;
            float l_ee = sqrt(EES_segr_shift*EES_segr_shift+EES_segz*EES_segz); 
            float tan_ee = sqrt((l_ee-EES_segz)/(l_ee+EES_segz));
            float calc_eta_ee = -log(tan_ee);
            ee_dis_eta_shift[shi] = calc_eta_ee-ext_eta;
          }

          ee_dis = calcDis(EES_segr,EES_segz,SP_ec_ee_r,SP_ec_ee_z);
          ee_dis_r = calcDisR(EES_segr,EES_segz,SP_ec_ee_r,SP_ec_ee_z);
          for(int ishi=0;ishi<11;ishi++)
            ee_dis_r_shift[ishi] = calcDisR(EES_segr,EES_segz,SP_ec_ee_r_shift[ishi],SP_ec_ee_z);
          ee_dis_eta = EES_segeta-ext_eta;
          ee_dis_phi = acos(cos(EES_segphi-ext_phi));
          ee_dis_phi_ip = acos(cos(EES_segphi-offline_phi));
        }
        if (EEL_segr>1e-5 && SP_ec_ee_r>1e-5){
          for (int shi=0;shi<11;shi++){
            float EEL_segr_shift=EEL_segr + 10*shi;
            float l_ee = sqrt(EEL_segr_shift*EEL_segr_shift+EEL_segz*EEL_segz); 
            float tan_ee = sqrt((l_ee-EEL_segz)/(l_ee+EEL_segz));
            float calc_eta_ee = -log(tan_ee);
            ee_dis_eta_shift[shi] = calc_eta_ee-ext_eta;
          }
          ee_dis = calcDis(EEL_segr,EEL_segz,SP_ec_ee_r,SP_ec_ee_z);
          ee_dis_r = calcDisR(EEL_segr,EEL_segz,SP_ec_ee_r,SP_ec_ee_z);
          for(int ishi=0;ishi<11;ishi++)
            ee_dis_r_shift[ishi] = calcDisR(EEL_segr,EEL_segz,SP_ec_ee_r_shift[ishi],SP_ec_ee_z);
          ee_dis_eta = EEL_segeta-ext_eta;
          ee_dis_phi = acos(cos(EEL_segphi-ext_phi));
          ee_dis_phi_ip = acos(cos(EEL_segphi-offline_phi));
        }
        if (EIS_segr>1e-5 && SP_ec_inner_r>1e-5){
          inner_dis = calcDis(EIS_segr,EIS_segz,SP_ec_inner_r,SP_ec_inner_z);
          inner_dis_r = calcDisR(EIS_segr,EIS_segz,SP_ec_inner_r,SP_ec_inner_z);
          inner_dis_r_shift = calcDisR(EIS_segr,EIS_segz,SP_ec_inner_r_shift,SP_ec_inner_z);
          inner_dis_eta = EIS_segeta-ext_eta;
          inner_dis_phi = acos(cos(EIS_segphi-ext_phi));
          inner_dis_phi_ip = acos(cos(EIS_segphi-offline_phi));
        }
        if (EIL_segr>1e-5 && SP_ec_inner_r>1e-5){
          inner_dis = calcDis(EIL_segr,EIL_segz,SP_ec_inner_r,SP_ec_inner_z);
          inner_dis_r = calcDisR(EIL_segr,EIL_segz,SP_ec_inner_r,SP_ec_inner_z);
          inner_dis_r_shift = calcDisR(EIL_segr,EIL_segz,SP_ec_inner_r_shift,SP_ec_inner_z);
          inner_dis_eta = EIL_segeta-ext_eta;
          inner_dis_phi = acos(cos(EIL_segphi-ext_phi));
          inner_dis_phi_ip = acos(cos(EIL_segphi-offline_phi));
        }
        if (EMS_segr>1e-5 && SP_ec_middle_r>1e-5){
          middle_dis = calcDis(EMS_segr,EMS_segz,SP_ec_middle_r,SP_ec_middle_z);
          middle_dis_r = calcDisR(EMS_segr,EMS_segz,SP_ec_middle_r,SP_ec_middle_z);
          middle_dis_eta = EMS_segeta-ext_eta;
          middle_dis_phi = acos(cos(EMS_segphi-ext_phi));
          middle_dis_phi_ip = acos(cos(EMS_segphi-offline_phi));
        }
        if (EML_segr>1e-5 && SP_ec_middle_r>1e-5){
          middle_dis = calcDis(EML_segr,EML_segz,SP_ec_middle_r,SP_ec_middle_z);
          middle_dis_r = calcDisR(EML_segr,EML_segz,SP_ec_middle_r,SP_ec_middle_z);
          middle_dis_eta = EML_segeta-ext_eta;
          middle_dis_phi = acos(cos(EML_segphi-ext_phi));
          middle_dis_phi_ip = acos(cos(EML_segphi-offline_phi));
        }
        if (BIS_segr>1e-5 && SP_br_inner_r>1e-5){
          br_inner_dis = calcDis(BIS_segr,BIS_segz,SP_br_inner_r,SP_br_inner_z);
          br_inner_dis_r = calcDisR(BIS_segr,BIS_segz,SP_br_inner_r,SP_br_inner_z);
          br_inner_dis_eta = BIS_segeta-ext_eta;
          br_inner_dis_phi = acos(cos(BIS_segphi-ext_phi));
          br_inner_dis_phi_ip = acos(cos(BIS_segphi-offline_phi));
        }
        if (BIL_segr>1e-5 && SP_br_inner_r>1e-5){
          br_inner_dis = calcDis(BIL_segr,BIL_segz,SP_br_inner_r,SP_br_inner_z);
          br_inner_dis_r = calcDisR(BIL_segr,BIL_segz,SP_br_inner_r,SP_br_inner_z);
          br_inner_dis_eta = BIL_segeta-ext_eta;
          br_inner_dis_phi = acos(cos(BIL_segphi-ext_phi));
          br_inner_dis_phi_ip = acos(cos(BIL_segphi-offline_phi));
        }

        float dphi_EILEML=(fabs(EIL_segphi)>1e-5 && fabs(EML_segphi)>1e-5)? calcCosDphi(EIL_segx,EIL_segy,EML_segx,EML_segy):0;
        float dphi_EISEMS=(fabs(EIS_segphi)>1e-5 && fabs(EMS_segphi)>1e-5)? acos(cos(EIS_segphi-EMS_segphi)):0;
        float dphi_BISEMS=(fabs(BIS_segphi)>1e-5 && fabs(EMS_segphi)>1e-5)? calcCosDphi(BIS_segx,BIS_segy,EML_segx,EML_segy):0;

        float calc_ec_radius_offline=0.;
        float calc_off_ec_radius=0.;
        float calc_off_ec_radius_phi=0.;

        //calculate radius from offline segment
        if (BIS_segr>1e-5 && EES_segr>1e-5 && EMS_segr>1e-5){
          //EES_segr-=100;
          calc_ec_radius_offline = computeRadius3Points(BIS_segz,BIS_segr,
              EES_segz,EES_segr,
              EMS_segz,EMS_segr);
          calc_ec_sagitta_offline = calcSagitta(BIS_segz,BIS_segr,
              EES_segz,EES_segr,
              EMS_segz,EMS_segr);
        }
        if (EIL_segr>1e-5 && EEL_segr>1e-5 && EML_segr>1e-5){
          //EEL_segr-=100;
          calc_ec_radius_offline = computeRadius3Points(EIL_segz,EIL_segr,
              EEL_segz,EEL_segr,
              EML_segz,EML_segr);
          calc_ec_sagitta_offline = calcSagitta(EIL_segz,EIL_segr,
              EEL_segz,EEL_segr,
              EML_segz,EML_segr);
        }
        if (BIS_segr>1e-5 && EES_segr>1e-5 && EMS_segr>1e-5 ){
          calc_off_ec_radius = computeRadius3Points(BIS_segr,BIS_segz,EES_segr,EES_segz,EMS_segr,EMS_segz);
          calc_off_ec_radius_phi = calc_off_ec_radius/cos(dphi_BISEMS);
        }
        if (EIL_segr>1e-5 && EEL_segr>1e-5 && EML_segr>1e-5 ){
          calc_off_ec_radius = computeRadius3Points(EIL_segr,EIL_segz,EEL_segr,EEL_segz,EML_segr,EML_segz);
          calc_off_ec_radius_phi = calc_off_ec_radius/cos(dphi_EILEML);
        }

      //}//!correct sagitta

      float calc_ec_radius=0.;
      float calc_ec_radius_shift[10];
      float calc_ec_sagitta=0.;
      float distance=0.;
      pair<float,float> center;
      float vali_radius=0;
      float R0=0,Z0=0;
      bool small = isSmall(SP_ec_ee_z);
      int L2_SL = (isSmall)? 0 : 1;
      //calculate radius from SP
      if(small){
        if (SP_br_inner_r>1e-5 && SP_ec_ee_r>1e-5 && SP_ec_middle_r>1e-5){
          calc_ec_radius = computeRadius3Points(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          for(int ishi=0;ishi<11;ishi++){
            calc_ec_radius_shift[ishi] = computeRadius3Points(SP_br_inner_z,SP_br_inner_r,
                SP_ec_ee_z,SP_ec_ee_r_shift[ishi],
                SP_ec_middle_z,SP_ec_middle_r);
          }
          distance = calcDistance(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          center = calcCenter(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          calc_ec_sagitta = calcSagitta(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          Z0 = center.first;
          R0 = center.second;
          vali_radius = sqrt(pow((R0-SP_br_inner_r),2)+pow((Z0-SP_br_inner_z),2));
        }
      }
      else{
        if (SP_ec_inner_r>1e-5 && SP_ec_ee_r>1e-5 && SP_ec_middle_r>1e-5){
          calc_ec_radius = computeRadius3Points(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          for (int ishi=0;ishi<11;ishi++){
            calc_ec_radius_shift[ishi] = computeRadius3Points(SP_ec_inner_z,SP_ec_inner_r,
                SP_ec_ee_z,SP_ec_ee_r_shift[ishi],
                SP_ec_middle_z,SP_ec_middle_r);
          }
          distance = calcDistance(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          center = calcCenter(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          calc_ec_sagitta = calcSagitta(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          Z0 = center.first;
          R0 = center.second;
          vali_radius = sqrt(pow((R0-SP_ec_inner_r),2)+pow((Z0-SP_ec_inner_z),2));
        }
      }

      if (calc_ec_radius>1e-5) phiVsInvR->Fill(usedPhi,1./calc_ec_radius);
      int PhibinNumberEE = GetBinNumberEE(usedPhi, L2_SL);

      if (i<0 || j<0) continue;

      int Qeta = (m_L2_charge*m_L2_etaMap>0)? 1 : 0; 
      int ieta = (m_L2_etaMap>0) ? 1 : 0;
      int icharge = (m_L2_charge>0) ? 1 : 0;
      if (alpha>ZERO_LIMIT){//alpha was caluclated
        if (makeLUT) hMDT_invpTvsAlpha[i][j][Qeta]->Fill(fabs(1/m_offline_pt),alpha);
        /*if(!correctSagitta)*/ hMDT_invpTvsAlphaAnormal[PhibinAll][j][icharge][ieta]->Fill(fabs(1/m_offline_pt),alpha);
        float deta = offline_eta-L1_eta;
        float dphi = acos(cos(offline_phi-L1_phi));
        float dr = sqrt(deta*deta+dphi*dphi);
      }
      if (tgcalpha>ZERO_LIMIT){
        if(makeLUT) hMDT_invpTvsTgcAlpha[PhibinNumberTgc][EtabinNumberTgc][tgc_qeta]->Fill(fabs(1/m_offline_pt),tgcalpha);
      }
      if (beta>ZERO_LIMIT){//beta was caluclated
        /*if(!correctSagitta)*/ hMDT_invpTvsBetaAnormal[PhibinAll][j][icharge][ieta]->Fill(fabs(1/m_offline_pt),beta);
        if (makeLUT) hMDT_invpTvsBeta[i][j][Qeta]->Fill(fabs(1/m_offline_pt),beta);
      }
      if (calc_ec_radius>ZERO_LIMIT){
        Distance->Fill(distance);
        DistanceVsPt->Fill(offline_pt,distance);
        DistanceVsPhi->Fill(offline_phi,distance);

        if(j<8){
          hMDT_pTvsEndcapSagittaAnormal[PhibinAll][j][icharge][ieta]->Fill(m_offline_pt,calc_ec_sagitta_offline);
          //if(correctSagitta)continue;
            if(PhibinAll==68 && j==1 && icharge==0 && ieta==0){
              for(int ishi=0;ishi<11;ishi++){
                hMDT_invpTvsEndcapRadius680200[ishi]->Fill(fabs(1/m_offline_pt),1./calc_ec_radius_shift[ishi]);
                if(ee_dis_r_shift[ishi]<10000){
                  //hMDT_EEDistanceR680200[ishi]->Fill(ee_dis_r_shift[ishi]);
                }
                if(ee_dis_eta_shift[ishi]<10000){
                  //hMDT_EEDistanceEta680200[ishi]->Fill(ee_dis_eta_shift[ishi]);
                }
              }
            }
            hMDT_invpTvsEndcapRadiusAnormal[PhibinAll][j][icharge][ieta]->Fill(fabs(1/m_offline_pt),1./calc_ec_radius);
            //hMDT_invpTvsEndcapRadiusAnormal[PhibinAll][j][icharge][ieta]->Fill(fabs(1/m_offline_pt),1./calc_ec_radius_offline);
            if(inner_dis<10000) hMDT_EIDistance[PhibinAll][j][icharge][ieta]->Fill(inner_dis);
            if(ee_dis<10000) hMDT_EEDistance[PhibinAll][j][icharge][ieta]->Fill(ee_dis);
            if(br_inner_dis<10000) hMDT_BIDistance[PhibinAll][j][icharge][ieta]->Fill(br_inner_dis);
            if(middle_dis<10000) hMDT_EMDistance[PhibinAll][j][icharge][ieta]->Fill(middle_dis);
            if(inner_dis_r<10000) hMDT_EIDistanceR[PhibinAll][j][icharge][ieta]->Fill(inner_dis_r);
            if(ee_dis_r<10000) hMDT_EEDistanceR[PhibinAll][j][icharge][ieta]->Fill(ee_dis_r);
            if(br_inner_dis_r<10000) hMDT_BIDistanceR[PhibinAll][j][icharge][ieta]->Fill(br_inner_dis_r);
            if(middle_dis_r<10000) hMDT_EMDistanceR[PhibinAll][j][icharge][ieta]->Fill(middle_dis_r);
            if(inner_dis_eta<10000) hMDT_EIDistanceEta[PhibinAll][j][icharge][ieta]->Fill(inner_dis_eta);
            if(ee_dis_eta<10000) hMDT_EEDistanceEta[PhibinAll][j][icharge][ieta]->Fill(ee_dis_eta);
            if(br_inner_dis_eta<10000) hMDT_BIDistanceEta[PhibinAll][j][icharge][ieta]->Fill(br_inner_dis_eta);
            if(middle_dis_eta<10000) hMDT_EMDistanceEta[PhibinAll][j][icharge][ieta]->Fill(middle_dis_eta);
            if(inner_dis_phi<10000) hMDT_EIDistancePhi[PhibinAll][j][icharge][ieta]->Fill(inner_dis_phi);
            if(ee_dis_phi<10000) hMDT_EEDistancePhi[PhibinAll][j][icharge][ieta]->Fill(ee_dis_phi);
            if(br_inner_dis_phi<10000) hMDT_BIDistancePhi[PhibinAll][j][icharge][ieta]->Fill(br_inner_dis_phi);
            if(middle_dis_phi<10000) hMDT_EMDistancePhi[PhibinAll][j][icharge][ieta]->Fill(middle_dis_phi);
          }
          if (makeLUT){
            if (!badPhibin) {
              hMDT_invpTvsEndcapRadiusNoSL[i][j][Qeta]->Fill(fabs(1/m_offline_pt),1./calc_ec_radius);
              hMDT_invpTvsEndcapRadius[PhibinNumberEE][j][Qeta][L2_SL]->Fill(fabs(1/m_offline_pt),1./calc_ec_radius);
            }
            hMDT_invpTvsEndcapRadiusAllPhiNoQeta[PhibinAll][j][icharge][ieta]->Fill(fabs(1/m_offline_pt),1./calc_ec_radius);
            hMDT_invpTvsEndcapRadiusLinear[PhibinAll][j][icharge][ieta]->Fill(m_offline_pt,calc_ec_radius);
          }
        }
      }
    else{//barrel
      int etaBin_br = (int)((m_L2_etaMap - EtaMin[L2_saddress])/EtaStep[L2_saddress]);
      int phiBin_br = (int)((m_L2_phiMap - PhiMin[L2_saddress])/PhiStep[L2_saddress]);
      int icharge = (m_L2_charge>0)? 1 : 0;
      if(etaBin_br<=-1) etaBin_br = 0;
      if(etaBin_br>=30) etaBin_br = 29;
      if(phiBin_br<=-1) phiBin_br = 0;
      if(phiBin_br>=30) phiBin_br = 29;
      if (br_radius > 1e-5) {
        if (makeLUT) hMDT_invpTvsinvRadius[icharge][L2_saddress][etaBin_br][phiBin_br]->Fill(br_radius, fabs(offline_pt));
      }
    }
  }
  fout->Write();
}
