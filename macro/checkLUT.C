#define checkLUT_cxx
#include "macro/checkLUT.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TFile.h"
#include "TGraphErrors.h"
using namespace std;

#define PI 3.14159265258979

double const PI_OVER_4 = PI/4.0;
double const PI_OVER_8 = PI/8.0;
double const PHI_RANGE = 12.0/PI_OVER_8;
double const ZERO_LIMIT = 10e-5;

float calcAlpha(float r1,float z1,float r2,float z2){
  float slope1 = r1/z1;
  float slope2 = (r2 - r1)/(z2 - z1);
  float alpha = 0;
  alpha = fabs(atan(slope1) - atan(slope2));
  return alpha;
}

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

void checkLUT::Loop()
{
  if (fChain == 0) return;
  TFile *file = new TFile("outputCheckLUT/aho_check.root","recreate");
  //TFile *file = new TFile("outputCheckLUT/data15_13TeV_checkLUT.root","recreate");
  TH1 *SAPTEndcap = new TH1F("SAPTEndacp",";sa pt endcap(GeV);Events",60,0,60);
  TH1 *SAPTBarrel = new TH1F("SAPTBarrel",";sa pt endcap(GeV);Events",60,0,60);
  TH1 *ptResoEndcap = new TH1F("ptResoEndcap",";sa pt resolution;Events",50,-0.5,0.5);
  TH1 *ptResoBarrel = new TH1F("ptResoBarrel",";sa pt resolution;Events",50,-0.5,0.5);
  TH1 *ptAlphaReso = new TH1F("ptAlphaReso",";sa ptAlpha resolution;Events",220,-1.1,1.1);
  TH1 *ptBetaReso = new TH1F("ptBetaReso",";sa ptBeta resolution;Events",220,-1.1,1.1);
  TH1 *ptEcRadiusReso = new TH1F("ptEcRadiusReso",";sa ptEcRadius resolution;Events",22,-1.1,1.1);
  TH2 *ptEcRadiusResoPhi = new TH2F("ptEcRadiusResoPhi",";sa ptEcRadius resolution;Phi",22,-1.1,1.1,64,-3.2,3.2);
  TH1 *endcapRadius = new TH1F("endcapRadius",";Endcap Radius(mm);Events",100,0,1000000);
  TH2 *endcapRadiusPhi = new TH2F("endcapRadiusPhi",";Endcap Radius(mm);Phi",100,0,1000000,64,-3.2,3.2);
  TH2 *endcapRadiusPt = new TH2F("endcapRadiusPt",";Endcap Radius(mm);PT(GeV)",100,0,1000000,60,0,60);
  TH1 *yamaPhi = new TH1F("yamaPhi",";phi;Events",64,-3.2,3.2);
  TH2 *L2AlphaVsPhi = new TH2F("L2AlphaVsPhi",";offline #phi;L2 #alpha",64,-3.2,3.2,100,0,0.1);
  TH2 *OfflineAlphaVsPhi = new TH2F("OfflineAlphaVsPhi",";offline #phi;offline #alpha",64,-3.2,3.2,100,0,0.1);
  TH1 *drEE1 = new TH1F("drEE1",";dR;Events",100,-0.1,0.1);
  TH1 *drEE2 = new TH1F("drEE2",";dR;Events",100,-0.1,0.1);
  TH1 *drInner1 = new TH1F("drInner1",";dR;Events",100,-0.1,0.1);
  TH1 *drInner2 = new TH1F("drInner2",";dR;Events",100,-0.1,0.1);
  TH2 *fitEELUT050201 = new TH2F("fitEELUT050201",";1/offline pt;1/R",250,0,0.25,400,0,0.00004);
  TH2 *fitEELUT050201_ue = new TH2F("fitEELUT050201_ue",";1/offline pt;1/R",250,0,0.25,400,0,0.00004);
  TH1 *resoRinner050201_ue = new TH1F("resoRinner050201_ue",";resolution #eta;Events",100,-0.1,0.1);
  TH1 *resoRinner050201_sita = new TH1F("resoRinner050201_sita",";resolution #eta;Events",100,-0.1,0.1);
  TH1 *resoRee050201_ue = new TH1F("resoRee050201_ue",";resolution #eta;Events",100,-0.1,0.1);
  TH1 *resoRee050201_sita = new TH1F("resoRee050201_sita",";resolution #eta;Events",100,-0.1,0.1);
  TH1 *EESPZ_ue = new TH1F("EESPZ_ue",";SP z(mm);Events",100,10000,12000);
  TH1 *EESPZ_sita = new TH1F("EESPZ_sita",";SP z(mm);Events",100,10000,12000);
  TH1 *EESPZ_01031 = new TH1F("EESPZ_01031",";SP z(mm);Events",100,10000,12000);
  TH1 *resolutionR[8], *resolutionZ[8];
  string resolutionRName[8] = {"resorBI","resorBM","resorBO","resorEI","resorEM","resorEO","resorEE","resorEBI"};
  string resolutionZName[8] = {"resozBI","resozBM","resozBO","resozEI","resozEM","resozEO","resozEE","resozEBI"};
  for (int i=0; i<8; i++){
    resolutionR[i] = new TH1F(resolutionRName[i].c_str(),";resolution R;Events",100,-0.1,0.1);
    resolutionZ[i] = new TH1F(resolutionZName[i].c_str(),";resolution Z;Events",100,-0.1,0.1);
  }

  int nentries = fChain->GetEntries();
  cout << "number of events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000000==0) cout << jentry << endl;
    //cout << "jentry=" << jentry << endl;
    fChain->GetEntry(jentry); 
    if (L2_pt_alpha>10e-5){
      float ptResoAlpha =  1-offline_pt/L2_pt_alpha;
      /*if (L2_eta>-1.25 && L2_eta<-1.15)
        if (L2_phi>-0.07 && L2_phi<0)
          if (offline_pt>30 && offline_pt<32.5)*/ ptAlphaReso->Fill(ptResoAlpha);
    }
    if (L2_pt_beta>10e-5){
      float ptResoBeta =  1-offline_pt/L2_pt_beta;
      ptBetaReso->Fill(ptResoBeta);
    }
    if (L2_pt_ec_radius>10e-5 && L2_pt_ec_radius<500-10e-5){
      float ptResoEcRadius =  1-offline_pt/L2_pt_ec_radius;
      ptEcRadiusReso->Fill(ptResoEcRadius);
      ptEcRadiusResoPhi->Fill(ptResoEcRadius,L2_phi);
      //cout << "offpt/ecrapt/reso=" << offline_pt << "/" << L2_pt_ec_radius << "/" << ptResoEcRadius << endl;
    }
    if (L2_ec_radius>10e-5) {
      endcapRadius->Fill(L2_ec_radius);
      endcapRadiusPhi->Fill(L2_ec_radius,L2_phi);
      if (L2_pt_ec_radius>0) endcapRadiusPt->Fill(L2_ec_radius,L2_pt_ec_radius);
    }
    if(L2_saddress<-0.5){//endcap
      float offline_segment_middle_r = offline_segment_r->at(4);
      float offline_segment_middle_z = offline_segment_z->at(4);
      float offline_segment_outer_r = offline_segment_r->at(5);
      float offline_segment_outer_z = offline_segment_z->at(5);
      //cout << "offline segment z middle/outer=" << offline_segment_middle_z << "/" << offline_segment_outer_z << endl;
      //cout << "offline segment phi middle/outer=" << offline_segment_middle_phi << "/" << offline_segment_outer_phi << endl;
      pair<int,int> ansBin = GetBinNumber(L2_eta,L2_phi);
      int etaBin = ansBin.first;
      int phiBin = ansBin.second;
      int qeta = L2_eta*L2_charge<0 ? 0 : 1 ;
      SAPTEndcap->Fill(L2_pt);
      float resoEc = 1-offline_pt/L2_pt;
      ptResoEndcap->Fill(resoEc);
      float speez = fabs(L2_spz->at(6));
      if (phiBin==1 && etaBin+1==3 && qeta==1){
        EESPZ_01031->Fill(speez);
      }
      if (phiBin==5 && etaBin==1 && qeta==0){
        float offline_segment_ee_r = offline_segment_r->at(6);
        float offline_segment_inner_r = offline_segment_r->at(3);
        float l2_segment_ee_r = L2_spr->at(6);
        float l2_segment_inner_r = L2_spr->at(3);
        float reso_r_inner = (offline_segment_inner_r-l2_segment_inner_r)/offline_segment_inner_r;
        float reso_r_ee = (offline_segment_ee_r-l2_segment_ee_r)/offline_segment_ee_r;
        if (L2_spr->at(3)>1e-5 && L2_spr->at(0)<1e-5){
          fitEELUT050201->Fill(1/offline_pt,1./L2_ec_radius);
          if(1./L2_ec_radius > 7.5*1e-5/offline_pt-1e-6) {
            fitEELUT050201_ue->Fill(1./offline_pt,1./L2_ec_radius);
            EESPZ_ue->Fill(speez);
            resoRinner050201_ue->Fill(reso_r_inner);
            resoRee050201_ue->Fill(reso_r_ee);
          }
          else {
            resoRinner050201_sita->Fill(reso_r_inner);
            resoRee050201_sita->Fill(reso_r_ee);
            EESPZ_sita->Fill(speez);
          }
        }
      }
      if (fabs(offline_segment_middle_z)>1e-5 && fabs(offline_segment_outer_z)>1e-5){
        float offline_alpha = 
          calcAlpha(offline_segment_middle_r,offline_segment_middle_z,offline_segment_outer_r,offline_segment_outer_z);
        //cout << "offline alpha=" << offline_alpha << endl;
        //OfflineAlphaVsPhi->Fill(offline_alpha,offline_phi);

        //cout << "etabin/phibin/qeta=" << etaBin+1 << "/" << phiBin << "/" << qeta << endl;
        if (phiBin==6 && etaBin+1==3 && qeta==1){
          float offline_segment_ee_r = offline_segment_r->at(6);
          float offline_segment_inner_r = offline_segment_r->at(3);
          float l2_segment_ee_r = L2_spr->at(6);
          float l2_segment_inner_r = L2_spr->at(3);
          if (l2_segment_ee_r>1e-5 && offline_segment_ee_r>1e-5){
            float dree = (offline_segment_ee_r-l2_segment_ee_r)/offline_segment_ee_r; 
            //cout << "drEE=" << dree << endl;
            if (1./L2_ec_radius > 1e-4/offline_pt) drEE1->Fill(dree);
            else drEE2->Fill(dree);
          }
          if (l2_segment_inner_r>1e-5 && offline_segment_inner_r>1e-5){
            float drinner = (offline_segment_inner_r-l2_segment_inner_r)/offline_segment_inner_r; 
            //cout << "drInner=" << drinner << endl;
            if (1./L2_ec_radius > 1e-4/offline_pt) drInner1->Fill(drinner);
            else drInner2->Fill(drinner);
          }

        }
        
      }
    }else{//barrel
      float resoBr = 1-offline_pt/L2_pt;
      ptResoBarrel->Fill(resoBr);
      SAPTBarrel->Fill(L2_pt);
    }
    int saddress = L2_saddress;
    /*for (int sta=0; sta<7; sta++){
      float spz = L2_spz->at(sta);
      float spr = L2_spr->at(sta);
      //cout << "spr/spz/sta=" << spr << "/" << spz << "/" << sta << endl;
      for (int offsta=0; offsta<offline_segment_r->size(); offsta++){
        int offstaChamber = offline_segment_chamber->at(offsta);
        float offstaR = offline_segment_r->at(offsta);
        float offstaZ = offline_segment_z->at(offsta);
        //cout << "offstaR/offstaZ/chamber=" << offstaR << "/" << offstaZ << "/" << offstaChamber << endl;
        //float offstaEta = offline_segment_eta->at(offsta);
        //float offstaPhi = offline_segment_phi->at(offsta);
        if (offstaChamber<0) continue;
        if(fabs(spz)>10e-5){
          float resoz = (fabs(offstaZ)-fabs(spz))/fabs(offstaZ);
          if (offstaChamber==0){
            if(saddress<-0.5) resolutionZ[7]->Fill(resoz);
            else resolutionZ[0]->Fill(resoz);
          }
          else if (offstaChamber==sta) resolutionZ[offstaChamber]->Fill(resoz);
          if (offstaChamber==5)
            if (resoz>-0.03 && resoz<-0.01) yamaPhi->Fill(L2_phi);
        }
        if(fabs(spr)>10e-5){
          float resor = (offstaR-spr)/offstaR;
          if (offstaChamber==0){
            if(saddress<-0.5) resolutionR[7]->Fill(resor);
            else resolutionR[0]->Fill(resor);
          }
          else if (offstaChamber==sta) {
            resolutionR[offstaChamber]->Fill(resor);
          }
        }
      }
    }*/
  }
  file->Write();
}
