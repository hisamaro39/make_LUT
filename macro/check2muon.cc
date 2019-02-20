#define check2muon_cxx
#include "check2muon.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TEllipse.h>
#include <TLine.h>
#include <iostream>
#include <sstream>
#include <TLorentzVector.h>
using namespace std;

float calcDr(float eta1,float phi1,float eta2,float phi2){
  float deta=eta1-eta2;
  float dphi=acos(cos(phi1-phi2));
  float dr=sqrt(deta*deta+dphi*dphi);
  return dr;
}

void check2muon::Loop()
{
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  //gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.4);
  gStyle->SetPadBottomMargin(0.4);
  gStyle->SetOptStat(0);

  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLTmu4.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLTmu4msonly.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLT2mu4.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLT2mu4_default.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLT2mu4_add_road_separation.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLT2mu4_add_road_separation_newest_muComb.root","recreate");
  TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLTmu4msonly_default.root","recreate");
  //TFile *file = new TFile("outputEfficiency/2muon/efficiency_2muon_jpsi_mu4_mu20_HLTmu4msonly_add_road_separation.root","recreate");
  //TFile *file = new TFile("aho.root","recreate");
  int nentries = fChain->GetEntries();

  cout << "Number of events is " << nentries << endl;

  TH1 *OfflineDeltaR = new TH1F("OfflineDeltaR",";#DeltaR;Number of events",100,0,1);
  TH2 *DrSAEF = new TH2F("DrSAEF",";#DeltaR_{SA};#DeltaR_{EF}",100,0,0.5,100,0,0.5);
  TH1 *JpsiPt = new TH1F("JpsiPt",";J/#psi p_{T}(GeV);Number of events",100,0,100);
  TH2 *JpsiPtDeltaR = new TH2F("JpsiPtDeltaR",";J/#psi p_{T}(GeV);#DeltaR",100,0,100,100,0,0.5);
  stringstream name_OfflineEtaPhi1,name_OfflineEtaPhi2;
  stringstream name_ExtEtaPhi1,name_ExtEtaPhi2;
  stringstream name_RoIEtaPhi1,name_RoIEtaPhi2;
  stringstream name_SAEtaPhi1,name_SAEtaPhi2;
  stringstream name_SAEtaPhiMS1,name_SAEtaPhiMS2;
  stringstream name_EFEtaPhi1,name_EFEtaPhi2;
  stringstream name_MdtZR1,name_MdtZR2;
  stringstream canvas,display,canvas2,canvas3,display2,display3;
  int count=0;

  int num[3][7][40];//level: dr
  int num2[3][7][20];//level: jpsi pt
  for (int si=0;si<3;si++){
    for (int i=0;i<7;i++){
      for (int j=0;j<40;j++) num[si][i][j]=0;
      for (int k=0;k<40;k++) num2[si][i][k]=0;
    }
  }

  for (int jentry=0; jentry<nentries;jentry++) {
   // cout << "jentry=" << jentry << endl;
    //if (jentry==500) break;
    fChain->GetEntry(jentry); 
    float muon1_pt = tag_offline_pt;
    float muon2_pt = probe_offline_pt;
    float muon1_eta = tag_offline_eta;
    float muon2_eta = probe_offline_eta;
    float muon1_phi = tag_offline_phi;
    float muon2_phi = probe_offline_phi;
    float muon1_ext_eta = tag_ext_eta;
    float muon2_ext_eta = probe_ext_eta;
    float muon1_ext_phi = tag_ext_phi;
    float muon2_ext_phi = probe_ext_phi;
    float muon1_roi_eta = tag_roi_eta;
    float muon2_roi_eta = probe_roi_eta;
    float muon1_roi_phi = tag_roi_phi;
    float muon2_roi_phi = probe_roi_phi;
    float muon1_sa_pt = tag_sa_pt;
    float muon2_sa_pt = probe_sa_pt;
    float muon1_sa_eta = tag_sa_eta;
    float muon2_sa_eta = probe_sa_eta;
    float muon1_sa_etaMS = tag_sa_etaMS;
    float muon2_sa_etaMS = probe_sa_etaMS;
    float muon1_sa_phi = tag_sa_phi;
    float muon2_sa_phi = probe_sa_phi;
    float muon1_sa_phiMS = tag_sa_phiMS;
    float muon2_sa_phiMS = probe_sa_phiMS;
    float muon1_ef_pt = tag_ef_pt;
    float muon2_ef_pt = probe_ef_pt;
    float muon1_ef_eta = tag_ef_eta;
    float muon2_ef_eta = probe_ef_eta;
    float muon1_ef_phi = tag_ef_phi;
    float muon2_ef_phi = probe_ef_phi;
    bool muon1_passL1 = tag_passL1;
    bool muon2_passL1 = probe_passL1;
    bool muon1_passSA = tag_passSA;
    bool muon2_passSA = probe_passSA;
    bool muon1_passComb = tag_passComb;
    bool muon2_passComb = probe_passComb;
    bool muon1_passEF = tag_passEF;
    bool muon2_passEF = probe_passEF;
    bool muon1_Combexist = tag_Combexist;
    bool muon2_Combexist = probe_Combexist;
    bool muon1_EFexist = tag_EFexist;
    bool muon2_EFexist = probe_EFexist;
    bool dimuon_passL1 = (muon1_passL1 && muon2_passL1)? true : false;
    bool dimuon_passSA = (muon1_passSA && muon2_passSA && dimuon_passL1)? true : false;
    bool dimuon_passSAovrmv = (muon1_Combexist && muon2_Combexist && dimuon_passSA)? true : false;
    bool dimuon_passComb = (muon1_passComb && muon2_passComb && dimuon_passSAovrmv)? true : false;
    bool dimuon_passCombovrmv = (muon1_EFexist && muon2_EFexist && dimuon_passComb)? true : false;
    bool dimuon_passEF = (muon1_passEF && muon2_passEF && dimuon_passCombovrmv)? true : false;
    float dr=calcDr(muon1_eta,muon1_phi,muon2_eta,muon2_phi);
    float dr_sa=calcDr(muon1_sa_eta,muon1_sa_phi,muon2_sa_eta,muon2_sa_phi);
    float dr_ef=calcDr(muon1_ef_eta,muon1_ef_phi,muon2_ef_eta,muon2_ef_phi);
    TLorentzVector muon1,muon2;
    muon1.SetPtEtaPhiM(muon1_pt,muon1_eta,muon1_phi,0.1056);
    muon2.SetPtEtaPhiM(muon2_pt,muon2_eta,muon2_phi,0.1056);
    float jpsi_pt = (muon1+muon2).Pt();
    float jpsi_eta = (muon1+muon2).Eta();
    JpsiPt->Fill(jpsi_pt);
    JpsiPtDeltaR->Fill(jpsi_pt,dr);
    //cout << "****************" << endl;
    //cout << "jpsi pt=" << jpsi_pt << endl;
    //cout << "offline pt1/pt2/dr=" << muon1_pt << "/" << muon2_pt << "/" << dr << endl;
    //cout << "pass L1/SA/ovrmv/Comb/ovrmv/EF=" 
      //<< dimuon_passL1 << "/" << dimuon_passSA << "/" << dimuon_passSAovrmv << "/" << dimuon_passComb 
      //<< "/" << dimuon_passCombovrmv << "/" << dimuon_passEF << endl;
    vector<float> *muon1_tgc_hit_eta = tag_tgc_hit_eta;
    int muon1_nTgcHit = muon1_tgc_hit_eta->size();
    if (muon1_sa_pt>1e-5 && muon2_sa_pt>1e-5 && muon1_ef_pt>1e-5 && muon2_ef_pt>1e-5){
      if (dr_ef>0.1 && dr_sa<0.05){
        //if (count==15) break;
        count++;
        name_OfflineEtaPhi1 << "OfflineEtaPhi1_event" << count;
        name_OfflineEtaPhi2 << "OfflineEtaPhi2_event" << count;
        name_ExtEtaPhi1 << "ExtEtaPhi1_event" << count;
        name_ExtEtaPhi2 << "ExtEtaPhi2_event" << count;
        name_RoIEtaPhi1 << "RoIEtaPhi1_event" << count;
        name_RoIEtaPhi2 << "RoIEtaPhi2_event" << count;
        name_SAEtaPhi1 << "SAEtaPhi1_event" << count;
        name_SAEtaPhi2 << "SAEtaPhi2_event" << count;
        name_SAEtaPhiMS1 << "SAEtaPhiMS1_event" << count;
        name_SAEtaPhiMS2 << "SAEtaPhiMS2_event" << count;
        name_EFEtaPhi1 << "EFEtaPhi1_event" << count;
        name_EFEtaPhi2 << "EFEtaPhi2_event" << count;
        name_MdtZR1 << "MdtZR1_event" << count;
        name_MdtZR2 << "MdtZR2_event" << count;
        float etamin1=muon1_eta-0.3, etamax1=muon1_eta+0.3;
        float phimin1=muon1_phi-0.3, phimax1=muon1_phi+0.3;
        float etamin2=muon2_eta-0.3, etamax2=muon2_eta+0.3;
        float phimin2=muon2_phi-0.3, phimax2=muon2_phi+0.3;
        TH2 *OfflineEtaPhi1 = new TH2F(name_OfflineEtaPhi1.str().c_str(),";#eta;#phi",100,etamin1,etamax1,100,phimin1,phimax1);
        TH2 *OfflineEtaPhi2 = new TH2F(name_OfflineEtaPhi2.str().c_str(),";#eta;#phi",100,etamin2,etamax2,100,phimin2,phimax2);
        TH2 *ExtEtaPhi1 = new TH2F(name_ExtEtaPhi1.str().c_str(),";#eta;#phi",100,etamin1,etamax1,100,phimin1,phimax1);
        TH2 *ExtEtaPhi2 = new TH2F(name_ExtEtaPhi2.str().c_str(),";#eta;#phi",100,etamin2,etamax2,100,phimin2,phimax2);
        TH2 *RoIEtaPhi1 = new TH2F(name_RoIEtaPhi1.str().c_str(),";#eta;#phi",100,etamin1,etamax1,100,phimin1,phimax1);
        TH2 *RoIEtaPhi2 = new TH2F(name_RoIEtaPhi2.str().c_str(),";#eta;#phi",100,etamin2,etamax2,100,phimin2,phimax2);
        TH2 *SAEtaPhi1 = new TH2F(name_SAEtaPhi1.str().c_str(),";#eta;#phi",100,etamin1,etamax1,100,phimin1,phimax1);
        TH2 *SAEtaPhi2 = new TH2F(name_SAEtaPhi2.str().c_str(),";#eta;#phi",100,etamin2,etamax2,100,phimin2,phimax2);
        TH2 *SAEtaPhiMS1 = new TH2F(name_SAEtaPhiMS1.str().c_str(),";#eta;#phi",100,etamin1,etamax1,100,phimin1,phimax1);
        TH2 *SAEtaPhiMS2 = new TH2F(name_SAEtaPhiMS2.str().c_str(),";#eta;#phi",100,etamin2,etamax2,100,phimin2,phimax2);
        TH2 *EFEtaPhi1 = new TH2F(name_EFEtaPhi1.str().c_str(),";#eta;#phi",100,etamin1,etamax1,100,phimin1,phimax1);
        TH2 *EFEtaPhi2 = new TH2F(name_EFEtaPhi2.str().c_str(),";#eta;#phi",100,etamin2,etamax2,100,phimin2,phimax2);

        TLegend *leg = new TLegend(0.2,0.1,0.5,0.4);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        OfflineEtaPhi1->Fill(muon1_eta,muon1_phi);
        OfflineEtaPhi2->Fill(muon2_eta,muon2_phi);
        ExtEtaPhi1->Fill(muon1_ext_eta,muon1_ext_phi);
        ExtEtaPhi2->Fill(muon2_ext_eta,muon2_ext_phi);
        RoIEtaPhi1->Fill(muon1_roi_eta,muon1_roi_phi);
        RoIEtaPhi2->Fill(muon2_roi_eta,muon2_roi_phi);
        SAEtaPhi1->Fill(muon1_sa_eta,muon1_sa_phi);
        SAEtaPhi2->Fill(muon2_sa_eta,muon2_sa_phi);
        SAEtaPhiMS1->Fill(muon1_sa_etaMS,muon1_sa_phiMS);
        SAEtaPhiMS2->Fill(muon2_sa_etaMS,muon2_sa_phiMS);
        EFEtaPhi1->Fill(muon1_ef_eta,muon1_ef_phi);
        EFEtaPhi2->Fill(muon2_ef_eta,muon2_ef_phi);
        //leg->AddEntry(OfflineEtaPhi1,"muon1 offline","p");
        //leg->AddEntry(OfflineEtaPhi2,"muon2 offline","p");
        leg->AddEntry(ExtEtaPhi1,"muon1 offlineMS","p");
        leg->AddEntry(ExtEtaPhi2,"muon2 offlineMS","p");
        leg->AddEntry(RoIEtaPhi1,"muon1 RoI","p");
        leg->AddEntry(RoIEtaPhi2,"muon2 RoI","p");
        //leg->AddEntry(SAEtaPhi1,"muon1 SA","p");
        //leg->AddEntry(SAEtaPhi2,"muon2 SA","p");
        leg->AddEntry(SAEtaPhiMS1,"muon1 SAMS","p");
        leg->AddEntry(SAEtaPhiMS2,"muon2 SAMS","p");
        //leg->AddEntry(EFEtaPhi1,"muon1 EF","p");
        //leg->AddEntry(EFEtaPhi2,"muon2 EF","p");
        OfflineEtaPhi1->SetMarkerColor(2);
        ExtEtaPhi1->SetMarkerColor(2);
        RoIEtaPhi1->SetMarkerColor(2);
        SAEtaPhi1->SetMarkerColor(2);
        SAEtaPhiMS1->SetMarkerColor(2);
        EFEtaPhi1->SetMarkerColor(2);
        OfflineEtaPhi1->SetMarkerSize(2);
        ExtEtaPhi1->SetMarkerSize(2);
        RoIEtaPhi1->SetMarkerSize(2);
        SAEtaPhi1->SetMarkerSize(2);
        SAEtaPhiMS1->SetMarkerSize(2);
        EFEtaPhi1->SetMarkerSize(2);
        OfflineEtaPhi1->SetMarkerStyle(24);
        ExtEtaPhi1->SetMarkerStyle(28);
        RoIEtaPhi1->SetMarkerStyle(25);
        SAEtaPhi1->SetMarkerStyle(26);
        SAEtaPhiMS1->SetMarkerStyle(30);
        EFEtaPhi1->SetMarkerStyle(27);
        OfflineEtaPhi2->SetMarkerColor(4);
        ExtEtaPhi2->SetMarkerColor(4);
        RoIEtaPhi2->SetMarkerColor(4);
        SAEtaPhi2->SetMarkerColor(4);
        SAEtaPhiMS2->SetMarkerColor(4);
        EFEtaPhi2->SetMarkerColor(4);
        OfflineEtaPhi2->SetMarkerSize(2);
        ExtEtaPhi2->SetMarkerSize(2);
        RoIEtaPhi2->SetMarkerSize(2);
        SAEtaPhi2->SetMarkerSize(2);
        SAEtaPhiMS2->SetMarkerSize(2);
        EFEtaPhi2->SetMarkerSize(2);
        OfflineEtaPhi2->SetMarkerStyle(24);
        ExtEtaPhi2->SetMarkerStyle(28);
        RoIEtaPhi2->SetMarkerStyle(25);
        SAEtaPhi2->SetMarkerStyle(26);
        SAEtaPhiMS2->SetMarkerStyle(30);
        EFEtaPhi2->SetMarkerStyle(27);
        
        canvas << "canvas_" << count;
        //display << "outputEfficiency/2muon/EventDisplay/default/event" << count << ".pdf";
        display << "outputEfficiency/2muon/EventDisplay/add_road_separation/event" << count << ".pdf";
        //TCanvas *c1 = new TCanvas(canvas.str().c_str(),"",600,600);
        //OfflineEtaPhi1->Draw();
        //OfflineEtaPhi2->Draw("same");
        //ExtEtaPhi1->Draw();
        //ExtEtaPhi2->Draw("same");
        //RoIEtaPhi1->Draw("same");
        //RoIEtaPhi2->Draw("same");
        //SAEtaPhi1->Draw("same");
        //SAEtaPhi2->Draw("same");
        //SAEtaPhiMS1->Draw("same");
        //SAEtaPhiMS2->Draw("same");
        //EFEtaPhi1->Draw("same");
        //EFEtaPhi2->Draw("same");
        //leg->Draw();
        //c1->SaveAs(display.str().c_str());
        
        name_OfflineEtaPhi1.str("");
        name_OfflineEtaPhi2.str("");
        name_ExtEtaPhi1.str("");
        name_ExtEtaPhi2.str("");
        name_RoIEtaPhi1.str("");
        name_RoIEtaPhi2.str("");
        name_SAEtaPhi1.str("");
        name_SAEtaPhi2.str("");
        name_SAEtaPhiMS1.str("");
        name_SAEtaPhiMS2.str("");
        name_EFEtaPhi1.str("");
        name_EFEtaPhi2.str("");
        name_MdtZR1.str("");
        name_MdtZR2.str("");
        canvas.str("");
        display.str("");
        
        if (fabs(muon1_eta)>1.05) continue;
        cout << "***************" << endl;
        cout << "Close-by Event" << count << endl;
        cout << "OF muon1 pt/eta/phi=" << muon1_pt << "/" << muon1_eta << "/" << muon1_phi << endl;
        cout << "OF muon2 pt/eta/phi=" << muon2_pt << "/" << muon2_eta << "/" << muon2_phi << endl;
        cout << "SA muon1 pt/eta/phi=" << muon1_sa_pt << "/" << muon1_sa_eta << "/" << muon1_sa_phi << endl;
        cout << "SA muon2 pt/eta/phi=" << muon2_sa_pt << "/" << muon2_sa_eta << "/" << muon2_sa_phi << endl;
        cout << "EF muon1 pt/eta/phi=" << muon1_ef_pt << "/" << muon1_ef_eta << "/" << muon1_ef_phi << endl;
        cout << "EF muon2 pt/eta/phi=" << muon2_ef_pt << "/" << muon2_ef_eta << "/" << muon2_ef_phi << endl;
        cout << "Ext muon1 eta/phi=" << muon1_ext_eta << "/" << muon1_ext_phi << endl;
        cout << "Ext muon2 eta/phi=" << muon2_ext_eta << "/" << muon2_ext_phi << endl;
        cout << "dr off/sa/ef=" << dr << "/" << dr_sa << "/" << dr_ef << endl;
        float road_br_inner_aw1=tag_road_aw->at(0);
        float road_br_inner_bw1=tag_road_bw->at(0);
        float road_br_inner_aw2=probe_road_aw->at(0);
        float road_br_inner_bw2=probe_road_bw->at(0);
        float road_br_middle_aw1=tag_road_aw->at(1);
        float road_br_middle_bw1=tag_road_bw->at(1);
        float road_br_middle_aw2=probe_road_aw->at(1);
        float road_br_middle_bw2=probe_road_bw->at(1);
        float road_br_outer_aw1=tag_road_aw->at(2);
        float road_br_outer_bw1=tag_road_bw->at(2);
        float road_br_outer_aw2=probe_road_aw->at(2);
        float road_br_outer_bw2=probe_road_bw->at(2);
        cout << "barrel inner aw muon1/muon2=" << road_br_inner_aw1 << "/" << road_br_inner_aw2 << endl;
        cout << "barrel middle aw muon1/muon2=" << road_br_middle_aw1 << "/" << road_br_middle_aw2 << endl;
        cout << "barrel outer aw muon1/muon2=" << road_br_outer_aw1 << "/" << road_br_outer_aw2 << endl;
        float muon1_ec_inner_spr=tag_sp_r->at(3);
        float muon1_ec_inner_spz=tag_sp_z->at(3);
        float muon2_ec_inner_spr=probe_sp_r->at(3);
        float muon2_ec_inner_spz=probe_sp_z->at(3);
        float muon1_ec_middle_spr=tag_sp_r->at(4);
        float muon1_ec_middle_spz=tag_sp_z->at(4);
        float muon2_ec_middle_spr=probe_sp_r->at(4);
        float muon2_ec_middle_spz=probe_sp_z->at(4);
        float muon1_ec_inner_eta=-log(tan(atan(muon1_ec_inner_spr/fabs(muon1_ec_inner_spz))/2))*muon1_ec_inner_spz/fabs(muon1_ec_inner_spz);
        float muon2_ec_inner_eta=-log(tan(atan(muon2_ec_inner_spr/fabs(muon2_ec_inner_spz))/2))*muon2_ec_inner_spz/fabs(muon2_ec_inner_spz);
        float muon1_ec_middle_eta=-log(tan(atan(muon1_ec_middle_spr/fabs(muon1_ec_middle_spz))/2))*muon1_ec_middle_spz/fabs(muon1_ec_middle_spz);
        float muon2_ec_middle_eta=-log(tan(atan(muon2_ec_middle_spr/fabs(muon2_ec_middle_spz))/2))*muon2_ec_middle_spz/fabs(muon2_ec_middle_spz);
        //cout << "ec inner eta muon1/muon2=" << muon1_ec_inner_eta << "/" << muon2_ec_inner_eta << endl;
        //cout << "ec middle eta muon1/muon2=" << muon1_ec_middle_eta << "/" << muon2_ec_middle_eta << endl;
        float muon1_br_inner_spr=tag_sp_r->at(0);
        float muon1_br_inner_spz=tag_sp_z->at(0);
        float muon2_br_inner_spr=probe_sp_r->at(0);
        float muon2_br_inner_spz=probe_sp_z->at(0);
        float muon1_br_middle_spr=tag_sp_r->at(1);
        float muon1_br_middle_spz=tag_sp_z->at(1);
        float muon2_br_middle_spr=probe_sp_r->at(1);
        float muon2_br_middle_spz=probe_sp_z->at(1);
        float muon1_br_outer_spr=tag_sp_r->at(2);
        float muon1_br_outer_spz=tag_sp_z->at(2);
        float muon2_br_outer_spr=probe_sp_r->at(2);
        float muon2_br_outer_spz=probe_sp_z->at(2);
        float muon1_br_inner_eta=-log(tan(atan(muon1_br_inner_spr/fabs(muon1_br_inner_spz))/2))*muon1_br_inner_spz/fabs(muon1_br_inner_spz);
        float muon2_br_inner_eta=-log(tan(atan(muon2_br_inner_spr/fabs(muon2_br_inner_spz))/2))*muon2_br_inner_spz/fabs(muon2_br_inner_spz);
        float muon1_br_middle_eta=-log(tan(atan(muon1_br_middle_spr/fabs(muon1_br_middle_spz))/2))*muon1_br_middle_spz/fabs(muon1_br_middle_spz);
        float muon2_br_middle_eta=-log(tan(atan(muon2_br_middle_spr/fabs(muon2_br_middle_spz))/2))*muon2_br_middle_spz/fabs(muon2_br_middle_spz);
        float muon1_br_outer_eta=-log(tan(atan(muon1_br_outer_spr/fabs(muon1_br_outer_spz))/2))*muon1_br_outer_spz/fabs(muon1_br_outer_spz);
        float muon2_br_outer_eta=-log(tan(atan(muon2_br_outer_spr/fabs(muon2_br_outer_spz))/2))*muon2_br_outer_spz/fabs(muon2_br_outer_spz);
        cout << "br inner eta muon1/muon2=" << muon1_br_inner_eta << "/" << muon2_br_inner_eta << endl;
        cout << "br middle eta muon1/muon2=" << muon1_br_middle_eta << "/" << muon2_br_middle_eta << endl;
        cout << "br outer eta muon1/muon2=" << muon1_br_outer_eta << "/" << muon2_br_outer_eta << endl;
        
        float rmin1=9999999,rmax1=-9999999,zmin1=9999999,zmax1=-9999999;
        for (int mdt=0;mdt<tag_mdt_hit_r->size();mdt++){
          float mdtr=tag_mdt_hit_r->at(mdt);
          float mdtz=tag_mdt_hit_z->at(mdt);
          int mdtofflineid=tag_mdt_hit_offlineId->at(mdt);
          int mdtchamber=tag_mdt_hit_chamber->at(mdt);
          if (mdtchamber!=1) continue;
          cout << "muon1 mdt hits chamber/r/z/id=" 
            << mdtchamber << "/" << mdtr << "/" << mdtz << "/" << mdtofflineid << endl;
          if (mdtr<rmin1) rmin1=mdtr;
          if (mdtr>rmax1) rmax1=mdtr;
          if (mdtz<zmin1) zmin1=mdtz;
          if (mdtz>zmax1) zmax1=mdtz;
        }
        rmin1=rmin1-100;rmax1=rmax1+100;zmin1=zmin1-100;zmax1=zmax1+100;
        float road_ec_middle_aw1=tag_road_aw->at(4);
        float road_ec_middle_bw1=tag_road_bw->at(4);
        float road_rmin1=road_ec_middle_aw1*zmin1+road_ec_middle_bw1;
        float road_rmax1=road_ec_middle_aw1*zmax1+road_ec_middle_bw1;
        TLine *road1 = new TLine(zmin1,road_rmin1,zmax1,road_rmax1);
        //cout << "minr/maxr/minz/maxz=" << rmin1 << "/" << rmax1 << "/" << zmin1 << "/" << zmax1 << endl;
        canvas2 << "canvas2_" << count;
        //TCanvas *c2 = new TCanvas(canvas2.str().c_str(),"",800,600);
        TH2 *MdtZR1 = new TH2F(name_MdtZR1.str().c_str(),";Z(mm);R(mm)",100,zmin1,zmax1,100,rmin1,rmax1);
        MdtZR1->GetXaxis()->SetLabelSize(0.04);
        //MdtZR1->Draw();
        for (int mdt=0;mdt<tag_mdt_hit_r->size();mdt++){
          float mdtr=tag_mdt_hit_r->at(mdt);
          float mdtz=tag_mdt_hit_z->at(mdt);
          float mdtspace=tag_mdt_hit_space->at(mdt);
          int mdtchamber=tag_mdt_hit_chamber->at(mdt); 
          int mdtofflineid=tag_mdt_hit_offlineId->at(mdt);
          if (mdtchamber!=4) continue;
          TEllipse *circle = new TEllipse(mdtz,mdtr,15,15);
          circle->SetLineColor(2);
          if (mdtofflineid==1) circle->SetLineStyle(2);
          //circle->Draw();
          TEllipse *circle2 = new TEllipse(mdtz,mdtr,mdtspace,mdtspace);
          circle2->SetLineColor(4);
          if (mdtofflineid==1) circle2->SetLineStyle(2);
          //circle2->Draw();
        }
        //road1->Draw();
        display2 << "outputEfficiency/2muon/EventDisplayMDT/tag_event" << count << ".pdf";
        //c2->SaveAs(display2.str().c_str());
        canvas2.str("");
        display2.str("");

        float rmin2=9999999,rmax2=-9999999,zmin2=9999999,zmax2=-9999999;
        for (int mdt=0;mdt<probe_mdt_hit_r->size();mdt++){
          float mdtr=probe_mdt_hit_r->at(mdt);
          float mdtz=probe_mdt_hit_z->at(mdt);
          int mdtofflineid=probe_mdt_hit_offlineId->at(mdt);
          int mdtchamber=probe_mdt_hit_chamber->at(mdt);
          if (mdtchamber!=1) continue;
          cout << "muon2 mdt hits chamber/r/z/id=" 
            << mdtchamber << "/" << mdtr << "/" << mdtz << "/" << mdtofflineid << endl;
          if (mdtr<rmin2) rmin2=mdtr;
          if (mdtr>rmax2) rmax2=mdtr;
          if (mdtz<zmin2) zmin2=mdtz;
          if (mdtz>zmax2) zmax2=mdtz;
        }
        rmin2=rmin2-100;rmax2=rmax2+100;zmin2=zmin2-100;zmax2=zmax2+100;
        float road_ec_middle_aw2=probe_road_aw->at(4);
        float road_ec_middle_bw2=probe_road_bw->at(4);
        float road_rmin2=road_ec_middle_aw2*zmin2+road_ec_middle_bw2;
        float road_rmax2=road_ec_middle_aw2*zmax2+road_ec_middle_bw2;
        TLine *road2 = new TLine(zmin2,road_rmin2,zmax2,road_rmax2);
        //cout << "minr/maxr/minz/maxz=" << rmin2 << "/" << rmax2 << "/" << zmin2 << "/" << zmax2 << endl;
        canvas3 << "canvas3_" << count;
        //TCanvas *c3 = new TCanvas(canvas3.str().c_str(),"",800,600);
        TH2 *MdtZR2 = new TH2F(name_MdtZR2.str().c_str(),";Z(mm);R(mm)",100,zmin2,zmax2,100,rmin2,rmax2);
        MdtZR2->GetXaxis()->SetLabelSize(0.04);
        //MdtZR2->Draw();
        for (int mdt=0;mdt<probe_mdt_hit_r->size();mdt++){
          float mdtr=probe_mdt_hit_r->at(mdt);
          float mdtz=probe_mdt_hit_z->at(mdt);
          float mdtspace=probe_mdt_hit_space->at(mdt);
          int mdtchamber=probe_mdt_hit_chamber->at(mdt); 
          int mdtofflineid=probe_mdt_hit_offlineId->at(mdt);
          if (mdtchamber!=4) continue;
          TEllipse *circle = new TEllipse(mdtz,mdtr,15,15);
          circle->SetLineColor(2);
          if (mdtofflineid==1) circle->SetLineStyle(2);
          //circle->Draw();
          TEllipse *circle2 = new TEllipse(mdtz,mdtr,mdtspace,mdtspace);
          circle2->SetLineColor(4);
          if (mdtofflineid==1) circle2->SetLineStyle(2);
          //circle2->Draw();
        }
        //road2->Draw();
        display3 << "outputEfficiency/2muon/EventDisplayMDT/probe_event" << count << ".pdf";
        //c3->SaveAs(display3.str().c_str());
        canvas3.str("");
        display3.str("");
        
      }
      if (fabs(muon1_ef_pt-muon2_ef_pt)>1e-5) DrSAEF->Fill(dr_sa,dr_ef);
    }
    OfflineDeltaR->Fill(dr);
    if (jpsi_pt>15 && jpsi_pt<20){
      //cout << "jpsi pt=" << jpsi_pt << endl;
      //cout << "muon1 pt/eta/phi=" << muon1_pt << "/" << muon1_eta << "/" << muon1_phi << endl; 
      //cout << "muon2 pt/eta/phi=" << muon2_pt << "/" << muon2_eta << "/" << muon2_phi << endl; 
    }
    for (int side=0;side<3;side++){//0 All: 1 Barrel: 2 Endcap
      if(side==1 && fabs(jpsi_eta)>1.05) continue;
      if(side==2 && fabs(jpsi_eta)<1.05) continue;
      for (int idr=0;idr<40;idr++){
        if (dr>idr*0.01 && dr<(idr+1)*0.01){
          num[side][0][idr]++;
          if (dimuon_passL1) num[side][1][idr]++;
          if (dimuon_passSA) num[side][2][idr]++;
          if (dimuon_passSAovrmv) num[side][3][idr]++;
          if (dimuon_passComb) num[side][4][idr]++;
          if (dimuon_passCombovrmv) num[side][5][idr]++;
          if (dimuon_passEF) num[side][6][idr]++;
        }
      }
      for (int ipt=0;ipt<20;ipt++){
        if (jpsi_pt>ipt*5 && jpsi_pt<(ipt+1)*5){
          num2[side][0][ipt]++;
          if (dimuon_passL1) num2[side][1][ipt]++;
          if (dimuon_passSA) num2[side][2][ipt]++;
          if (dimuon_passSAovrmv) num2[side][3][ipt]++;
          if (dimuon_passComb) num2[side][4][ipt]++;
          if (dimuon_passCombovrmv) num2[side][5][ipt]++;
          if (dimuon_passEF) num2[side][6][ipt]++;
        }
      }
    }
  }

  float eff_dr[3][6][40],efferr_dr[3][6][40],xbin_dr[40],xbinerr_dr[40];
  float eff_pt[3][6][20],efferr_pt[3][6][20],xbin_pt[20],xbinerr_pt[20];
  TGraphErrors *EffDr[3][6],*EffJpsiPt[3][6];
  stringstream effName,effName2;
  string LEVEL[6] = {"L1","SA","SAovrmv","Comb","Combovrmv","EF"};
  string SIDE[3] = {"All","Barrel","Endcap"}; 
  for (int side=0;side<3;side++){
    for (int level=0;level<6;level++){
      for (int idr=0;idr<40;idr++){
        eff_dr[side][level][idr]=0;efferr_dr[side][level][idr]=0;
        xbin_dr[idr]=0.005+0.01*idr;
        xbinerr_dr[idr]=0.005;
        if (num[side][level][idr]>0){
          eff_dr[side][level][idr] = 1.*num[side][level+1][idr]/num[side][level][idr];
          efferr_dr[side][level][idr] = sqrt(eff_dr[side][level][idr]*(1-eff_dr[side][level][idr])/num[side][level][idr]);
        }
      }
      for (int ipt=0;ipt<20;ipt++){
        eff_pt[side][level][ipt]=0;efferr_pt[side][level][ipt]=0;
        xbin_pt[ipt]=2.5+5*ipt;
        xbinerr_pt[ipt]=2.5;
        //cout << "side/level/ipt/num=" << side << "/" << level << "/" << ipt << "/" << num2[side][level][ipt] << endl;
        if (num2[side][level][ipt]>0){
          eff_pt[side][level][ipt] = 1.*num2[side][level+1][ipt]/num2[side][level][ipt];
          efferr_pt[side][level][ipt] = sqrt(eff_pt[side][level][ipt]*(1-eff_pt[side][level][ipt])/num2[side][level][ipt]);
        }
      }
      effName << "EfficiencyDeltaR" << LEVEL[level] << SIDE[side] ;
      EffDr[side][level] = new TGraphErrors(40,xbin_dr,eff_dr[side][level],xbinerr_dr,efferr_dr[side][level]);
      EffDr[side][level]->SetName(effName.str().c_str());
      EffDr[side][level]->SetTitle("");
      EffDr[side][level]->GetXaxis()->SetLabelSize(0.07);
      EffDr[side][level]->GetYaxis()->SetLabelSize(0.07);
      file->Add(EffDr[side][level]);
      effName.str("");
      effName2 << "EfficiencyJpsiPt" << LEVEL[level] << SIDE[side];
      EffJpsiPt[side][level] = new TGraphErrors(20,xbin_pt,eff_pt[side][level],xbinerr_pt,efferr_pt[side][level]);
      EffJpsiPt[side][level]->SetName(effName2.str().c_str());
      EffJpsiPt[side][level]->SetTitle("");
      EffJpsiPt[side][level]->GetXaxis()->SetLabelSize(0.07);
      EffJpsiPt[side][level]->GetYaxis()->SetLabelSize(0.07);
      file->Add(EffJpsiPt[side][level]);
      effName2.str("");
    }
  }

  file->Write();
}
