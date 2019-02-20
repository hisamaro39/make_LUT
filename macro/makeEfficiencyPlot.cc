#include "macro/makeEfficiencyPlot.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <sstream>
#include <iostream>
#include "TFile.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TH2.h"
using namespace std;
  
#define PI  3.14159265258979

//// eta bins setup 
const double eta_bin[]  = {-2.4, -2.159, -1.7705, -1.4855, -1.29045, -1.1904, -1.05, -0.979, -0.8495, -0.7215, -0.564, -0.40, -0.228, -0.066, 0., 0.066, 0.228, 0.40, 0.564, 0.7215, 0.8495, 0.979, 1.05, 1.1904, 1.29045, 1.4855, 1.7705, 2.159, 2.4};
const int eta_nbin      = sizeof(eta_bin)/sizeof(double);

//// phi bins setup 
const double phi_bin[]  = 
{-16./16.*PI, -14./16.*PI, -12./16.*PI,  -10./16.*PI, -8./16.*PI,-6./16.*PI,-4./16.*PI, -2./16.*PI, 0./16.*PI, 2./16.*PI, 4./16.*PI, 6./16.*PI, 8./16.*PI, 10./16.*PI, 12./16.*PI, 14./16.*PI, 16./16.*PI};
const int phi_nbin      = sizeof(phi_bin)/sizeof(double);

//// ptvarcone30/pt bins setup
const double ptvarcone_bin[] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2,
  0.22,0.24,0.26,0.28,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1};
const int ptvarcone_nbin = sizeof(ptvarcone_bin)/sizeof(double);

const string side[2] = {"Barrel","Endcap"};
const string level[5] = {"L1","SA","Comb","EF","HLT"};

makeEffPlot::makeEffPlot(){}

makeEffPlot::~makeEffPlot(){}

void makeEffPlot::makePtEfficiencyPlot(const char* input_data,const char* input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period){

  int nChain = chain_list.size();

  TGraphErrors *ptEfficiency[5][2][26],*ptEfficiencyMC[5][2][26];
  TCanvas *c1[5][2][26];
  stringstream pteff_name,c1_name,pt_plot_name,pt_plot_title;

  TFile *file1 = TFile::Open(input_data);
  TFile *file2 = TFile::Open(input_mc);
  vector<int> divpt;
  vector<float> vec_ptcenter, vec_pterr;
  int nbin;
  for (int iside=0;iside<2;iside++){
    for (int ich=0;ich<nChain;ich++){
      divpt.clear();vec_ptcenter.clear();vec_pterr.clear();
      divpt = m_util.devidePt(chain_list[ich]);
      vec_ptcenter = m_util.getPtCenter(chain_list[ich]);
      vec_pterr = m_util.getPtError(chain_list[ich]);
      nbin = divpt.size()-1;
      for(int lv2=0;lv2<5;lv2++){ 
        pteff_name << level[lv2] << "PtEfficiency_" << chain_name_list[ich] << "_" << side[iside]; 
        ptEfficiency[lv2][iside][ich] = dynamic_cast<TGraphErrors*>(file1->Get(pteff_name.str().c_str()));
        ptEfficiencyMC[lv2][iside][ich] = dynamic_cast<TGraphErrors*>(file2->Get(pteff_name.str().c_str()));
        pteff_name.str("");
        ptEfficiency[lv2][iside][ich]->Draw("ap");

        c1_name << "canvas" << iside << ich << lv2 ;
        c1[lv2][iside][ich] = new TCanvas(c1_name.str().c_str(),"",800,1000);
        c1_name.str("");
        TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
        p1->Draw();
        p1->SetNumber(1);
        p1->SetBottomMargin(0.1);
        TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
        p2->SetBottomMargin(0.2);
        p2->Draw();
        p2->SetNumber(2);

        Double_t *pteffData = ptEfficiency[lv2][iside][ich]->GetY();
        Double_t *pteffMC = ptEfficiencyMC[lv2][iside][ich]->GetY();
        Double_t *ptefferrData = ptEfficiency[lv2][iside][ich]->GetEY();
        Double_t *ptefferrMC = ptEfficiencyMC[lv2][iside][ich]->GetEY();
        Double_t *ptbin = ptEfficiency[lv2][iside][ich]->GetX();
        Double_t *ptbinerr = ptEfficiency[lv2][iside][ich]->GetEX();
        
        float ptRatio[nbin],ptRatioErr[nbin], ptx[nbin], ptxerr[nbin];
        for (int i=0;i<nbin;i++){
          float pteff1 = pteffData[i];
          float pteff2 = pteffMC[i];
          float ptefferr1 = ptefferrData[i];
          float ptefferr2 = ptefferrMC[i];
          ptx[i]=ptbin[i];
          ptxerr[i]=ptbinerr[i];
          float ptratio = (pteffMC[i]>1e-5)? pteff1/pteff2 : 0;
          float ptratioerr = 
            (pteffMC[i]>1e-5)? sqrt(pteff2*pteff2*ptefferr1*ptefferr1 + pteff1*pteff1*ptefferr2*ptefferr2)/(pteff2*pteff2) : 0;
          ptRatio[i]=ptratio;
          ptRatioErr[i]=ptratioerr;
          
        }

        c1[lv2][iside][ich]->cd(1);
        pt_plot_title << chain_name_list[ich] << "  period " << period ;
        ptEfficiency[lv2][iside][ich]->SetTitle(pt_plot_title.str().c_str());
        ptEfficiency[lv2][iside][ich]->GetYaxis()->SetRangeUser(0,1.2);
        ptEfficiency[lv2][iside][ich]->GetYaxis()->SetTitle("Efficiency");
        ptEfficiency[lv2][iside][ich]->GetXaxis()->SetTitle("p_{T}(GeV)");
        ptEfficiency[lv2][iside][ich]->Draw("ap");
        ptEfficiencyMC[lv2][iside][ich]->SetLineColor(2);
        ptEfficiencyMC[lv2][iside][ich]->Draw("p");
        TLatex test;
        test.SetTextSize(0.03);
        test.DrawLatex(10,1.15,"ATLAS work in progress");
        if(iside==0) test.DrawLatex(10,1.08,"Z Tag&Probe |#eta|<1.05");
        else test.DrawLatex(10,1.08,"Z Tag&Probe |#eta|>1.05");
        test.SetTextSize(0.05);
        if(lv2==0) test.DrawLatex(70,1.1,"L1/probe");
        else if(lv2==1) test.DrawLatex(70,1.1,"SA/L1");
        else if(lv2==2) test.DrawLatex(70,1.1,"Comb/SA");
        else if(lv2==3) test.DrawLatex(70,1.1,"EF/Comb");
        else if(lv2==4) test.DrawLatex(70,1.1,"HLT/probe");
        test.DrawLatex(70,0.4,"Data16 25ns");
        test.SetTextColor(2);
        test.DrawLatex(70,0.3,"MC Zmumu");

        TGraphErrors *g_ptratio = new TGraphErrors(nbin,ptx,ptRatio,ptxerr,ptRatioErr);
        c1[lv2][iside][ich]->cd(2);
        g_ptratio->SetTitle("");
        g_ptratio->GetYaxis()->SetTitle("Ratio Data/MC");
        g_ptratio->GetYaxis()->SetTitleOffset(0.25);
        g_ptratio->GetYaxis()->SetTitleSize(0.13);
        g_ptratio->GetYaxis()->SetLabelSize(0.07);
        g_ptratio->GetYaxis()->SetRangeUser(0.7,1.3);
        g_ptratio->GetXaxis()->SetRangeUser(0,100);
        g_ptratio->GetXaxis()->SetLabelSize(0);
        g_ptratio->Draw("ap");
        TLine *line = new TLine(0,1,100,1);
        line->SetLineColor(4);
        line->SetLineStyle(2);
        line->Draw();
        pt_plot_name << "outputScaleFactor/result_pt_efficiency/" << level[lv2] << "PtEfficiency_" << chain_name_list[ich] << "_" << side[iside] << "_" << period << ".pdf"; 
        c1[lv2][iside][ich]->SaveAs(pt_plot_name.str().c_str());
        pt_plot_name.str("");
        pt_plot_title.str("");
      }
    }
  }
}
/////////////////////////////////////
void makeEffPlot::makeEtaEfficiencyPlot(const char* input_data, const char* input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period){

  int nChain = chain_list.size();

  TGraphErrors *etaEfficiency[5][26],*etaEfficiencyMC[5][26];
  TCanvas *c1[5][26];
  stringstream etaeff_name,c1_name,eta_plot_name,eta_plot_title;

  TFile *file1 = TFile::Open(input_data);
  TFile *file2 = TFile::Open(input_mc);

  for (int ich=0;ich<nChain;ich++){
    for(int lv2=0;lv2<5;lv2++){ 
      etaeff_name << level[lv2] << "etaEfficiency_" << chain_name_list[ich]; 
      etaEfficiency[lv2][ich] = dynamic_cast<TGraphErrors*>(file1->Get(etaeff_name.str().c_str()));
      etaEfficiencyMC[lv2][ich] = dynamic_cast<TGraphErrors*>(file2->Get(etaeff_name.str().c_str()));
      etaeff_name.str("");
      etaEfficiency[lv2][ich]->Draw("ap");

      c1_name << "canvas" << ich << lv2;
      c1[lv2][ich] = new TCanvas(c1_name.str().c_str(),"",800,1000);
      c1_name.str("");
      TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
      p1->Draw();
      p1->SetNumber(1);
      p1->SetBottomMargin(0.1);
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetBottomMargin(0.2);
      p2->Draw();
      p2->SetNumber(2);

      Double_t *etaeffData = etaEfficiency[lv2][ich]->GetY();
      Double_t *etaeffMC = etaEfficiencyMC[lv2][ich]->GetY();
      Double_t *etaefferrData = etaEfficiency[lv2][ich]->GetEY();
      Double_t *etaefferrMC = etaEfficiencyMC[lv2][ich]->GetEY();
      Double_t *etabin = etaEfficiency[lv2][ich]->GetX();
      Double_t *etabinerr = etaEfficiency[lv2][ich]->GetEX();

      float etaRatio[26],etaRatioErr[26], etax[26], etaxerr[26];
      for (int i=0;i<26;i++){
        float etaeff1 = etaeffData[i];
        float etaeff2 = etaeffMC[i];
        float etaefferr1 = etaefferrData[i];
        float etaefferr2 = etaefferrMC[i];
        etax[i]=etabin[i];
        etaxerr[i]=etabinerr[i];
        float etaratio = (etaeffMC[i]>1e-5)? etaeff1/etaeff2 : 0;
        float etaratioerr = 
          (etaeffMC[i]>1e-5)? sqrt(etaeff2*etaeff2*etaefferr1*etaefferr1 + etaeff1*etaeff1*etaefferr2*etaefferr2)/(etaeff2*etaeff2) : 0;
        etaRatio[i]=etaratio;
        etaRatioErr[i]=etaratioerr;

      }

      c1[lv2][ich]->cd(1);
      eta_plot_title << chain_name_list[ich]<< "  period " << period;
      etaEfficiency[lv2][ich]->SetTitle(eta_plot_title.str().c_str());
      etaEfficiency[lv2][ich]->GetYaxis()->SetRangeUser(0,1.2);
      etaEfficiency[lv2][ich]->GetYaxis()->SetTitle("Efficiency");
      etaEfficiency[lv2][ich]->GetXaxis()->SetTitle("#eta");
      etaEfficiency[lv2][ich]->Draw("ap");
      etaEfficiencyMC[lv2][ich]->SetLineColor(2);
      etaEfficiencyMC[lv2][ich]->Draw("p");
      TLatex test;
      test.SetTextSize(0.03);
      test.DrawLatex(-2.8,1.15,"ATLAS work in progress");
      test.DrawLatex(-2.8,1.08,"Z Tag&Probe");
      test.DrawLatex(-1.1,1.08,Form("Barrel p_{T}>%dGeV",m_util.returnPlateau(chain_list[ich],0)));
      test.DrawLatex(-1.1,1.0,Form("Endcap p_{T}>%dGeV",m_util.returnPlateau(chain_list[ich],1)));
      test.SetTextSize(0.05);
      if(lv2==0) test.DrawLatex(1,1.1,"L1/probe");
      else if(lv2==1) test.DrawLatex(1,1.1,"SA/L1");
      else if(lv2==2) test.DrawLatex(1,1.1,"Comb/SA");
      else if(lv2==3) test.DrawLatex(1,1.1,"EF/Comb");
      else if(lv2==4) test.DrawLatex(1,1.1,"HLT/probe");
      test.DrawLatex(1,0.4,"Data16 25ns");
      test.SetTextColor(2);
      test.DrawLatex(1,0.3,"MC Zmumu");

      TGraphErrors *g_etaratio = new TGraphErrors(26,etax,etaRatio,etaxerr,etaRatioErr);
      c1[lv2][ich]->cd(2);
      g_etaratio->SetTitle("");
      g_etaratio->GetYaxis()->SetTitle("Ratio Data/MC");
      g_etaratio->GetYaxis()->SetTitleOffset(0.25);
      g_etaratio->GetYaxis()->SetTitleSize(0.13);
      g_etaratio->GetYaxis()->SetLabelSize(0.07);
      g_etaratio->GetYaxis()->SetRangeUser(0.7,1.3);
      g_etaratio->GetXaxis()->SetRangeUser(-2.6,2.6);
      g_etaratio->GetXaxis()->SetLabelSize(0);
      g_etaratio->Draw("ap");
      TLine *line = new TLine(-2.6,1,2.6,1);
      line->SetLineColor(4);
      line->SetLineStyle(2);
      line->Draw();

      eta_plot_name << "outputScaleFactor/result_eta_efficiency/" << level[lv2] << "etaEfficiency_" << chain_name_list[ich] << "_" << period << ".pdf"; 
      c1[lv2][ich]->SaveAs(eta_plot_name.str().c_str());
      eta_plot_name.str("");
      eta_plot_title.str("");

    }
  }

}
//////////////////////////////////////////////////////
void makeEffPlot::makePhiEfficiencyPlot(const char* input_data, const char* input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period){

  int nChain = chain_list.size();

  TGraphErrors *phiEfficiency[5][2][26],*phiEfficiencyMC[5][2][26];
  TCanvas *c1[5][2][26];
  stringstream phieff_name,c1_name,phi_plot_name,phi_plot_title;

  TFile *file1 = TFile::Open(input_data);
  TFile *file2 = TFile::Open(input_mc);

  for (int iside=1;iside<2;iside++){
    for (int ich=0;ich<nChain;ich++){
      for(int lv2=4;lv2<5;lv2++){ 
        phieff_name << level[lv2] << "phiEfficiency_" << chain_name_list[ich] << "_" << side[iside]; 
        phiEfficiency[lv2][iside][ich] = dynamic_cast<TGraphErrors*>(file1->Get(phieff_name.str().c_str()));
        phiEfficiencyMC[lv2][iside][ich] = dynamic_cast<TGraphErrors*>(file2->Get(phieff_name.str().c_str()));
        phieff_name.str("");
        phiEfficiency[lv2][iside][ich]->Draw("ap");

        c1_name << "canvas" << iside << ich << lv2 ;
        c1[lv2][iside][ich] = new TCanvas(c1_name.str().c_str(),"",800,1000);
        c1_name.str("");
        TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
        p1->Draw();
        p1->SetNumber(1);
        p1->SetBottomMargin(0.1);
        TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
        p2->SetBottomMargin(0.2);
        p2->Draw();
        p2->SetNumber(2);

        Double_t *phieffData = phiEfficiency[lv2][iside][ich]->GetY();
        Double_t *phieffMC = phiEfficiencyMC[lv2][iside][ich]->GetY();
        Double_t *phiefferrData = phiEfficiency[lv2][iside][ich]->GetEY();
        Double_t *phiefferrMC = phiEfficiencyMC[lv2][iside][ich]->GetEY();
        Double_t *phibin = phiEfficiency[lv2][iside][ich]->GetX();
        Double_t *phibinerr = phiEfficiency[lv2][iside][ich]->GetEX();
        
        float phiRatio[32],phiRatioErr[32], phix[32], phixerr[32];
        for (int i=0;i<32;i++){
          float phieff1 = phieffData[i];
          float phieff2 = phieffMC[i];
          float phiefferr1 = phiefferrData[i];
          float phiefferr2 = phiefferrMC[i];
          phix[i]=phibin[i];
          phixerr[i]=phibinerr[i];
          float phiratio = (phieffMC[i]>1e-5)? phieff1/phieff2 : 0;
          float phiratioerr = 
            (phieffMC[i]>1e-5)? sqrt(phieff2*phieff2*phiefferr1*phiefferr1 + phieff1*phieff1*phiefferr2*phiefferr2)/(phieff2*phieff2) : 0;
          phiRatio[i]=phiratio;
          phiRatioErr[i]=phiratioerr;
          
        }

        c1[lv2][iside][ich]->cd(1);
        phi_plot_title << chain_name_list[ich] << "  period "  << period;
        phiEfficiency[lv2][iside][ich]->SetTitle(phi_plot_title.str().c_str());
        phiEfficiency[lv2][iside][ich]->GetYaxis()->SetRangeUser(0,1.2);
        phiEfficiency[lv2][iside][ich]->GetYaxis()->SetTitle("Efficiency");
        phiEfficiency[lv2][iside][ich]->GetXaxis()->SetTitle("#phi");
        phiEfficiency[lv2][iside][ich]->Draw("ap");
        phiEfficiencyMC[lv2][iside][ich]->SetLineColor(2);
        phiEfficiencyMC[lv2][iside][ich]->Draw("p");
        TLatex test;
        test.SetTextSize(0.03);
        test.DrawLatex(-3.4,0.3,"ATLAS work in progress");
        if(iside==0) test.DrawLatex(-3.4,0.23,"Z Tag&Probe |#eta|<1.05");
        else test.DrawLatex(-3.4,0.23,"Z Tag&Probe |#eta|>1.05");
        test.DrawLatex(-3.4,0.16,Form("p_{T}>%dGeV",m_util.returnPlateau(chain_list[ich],iside)));
        test.SetTextSize(0.05);
        if(lv2==0) test.DrawLatex(1.3,1.1,"L1/probe");
        else if(lv2==1) test.DrawLatex(1.2,1.1,"SA/L1");
        else if(lv2==2) test.DrawLatex(1.2,1.1,"Comb/SA");
        else if(lv2==3) test.DrawLatex(1.2,1.1,"EF/Comb");
        else if(lv2==4) test.DrawLatex(1.2,1.1,"HLT/probe");
        test.DrawLatex(1.2,0.4,"Data16 25ns");
        test.SetTextColor(2);
        test.DrawLatex(1.2,0.3,"MC Zmumu");

        TGraphErrors *g_phiratio = new TGraphErrors(32,phix,phiRatio,phixerr,phiRatioErr);
        c1[lv2][iside][ich]->cd(2);
        g_phiratio->SetTitle("");
        g_phiratio->GetYaxis()->SetTitle("Ratio Data/MC");
        g_phiratio->GetYaxis()->SetTitleOffset(0.25);
        g_phiratio->GetYaxis()->SetTitleSize(0.13);
        g_phiratio->GetYaxis()->SetLabelSize(0.07);
        g_phiratio->GetYaxis()->SetRangeUser(0.7,1.3);
        g_phiratio->GetXaxis()->SetRangeUser(-3.2,3.2);
        g_phiratio->GetXaxis()->SetLabelSize(0);
        g_phiratio->Draw("ap");
        TLine *line = new TLine(-3.2,1,3.2,1);
        line->SetLineColor(4);
        line->SetLineStyle(2);
        line->Draw();
        
        phi_plot_name << "outputScaleFactor/result_phi_efficiency/" << level[lv2] << "phiEfficiency_" << chain_name_list[ich] << "_" << side[iside]  << "_" << period << ".pdf"; 
        c1[lv2][iside][ich]->SaveAs(phi_plot_name.str().c_str());
        phi_plot_name.str("");
        phi_plot_title.str("");
        
      }
    }
  }
}
//////////////////////////////////////////////////////
void makeEffPlot::makeEtaPhiEfficiencyPlot(const char* input_data, const char* input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period){
  gStyle->SetOptStat(0);

  int nChain = chain_list.size();

  TFile *file1 = TFile::Open(input_data);
  TFile *file2 = TFile::Open(input_mc);

  TH2 *etaphiEfficiency[5][26],*etaphiEfficiencyMC[5][26],*ScaleFactor[5][26];
  TCanvas *c1[5][26],*c2[5][26],*c3[5][26];
  stringstream etaphieff_name,c1_name,c2_name,c3_name,etaphi_plot_name,title,sf_plot_name;

  for (int ich=0;ich<nChain;ich++){
    for(int lv2=0;lv2<5;lv2++){ 
      etaphieff_name << level[lv2] << "etaphiEfficiency_" << chain_name_list[ich]; 
      etaphiEfficiency[lv2][ich] = dynamic_cast<TH2F*>(file1->Get(etaphieff_name.str().c_str()));
      etaphiEfficiencyMC[lv2][ich] = dynamic_cast<TH2F*>(file2->Get(etaphieff_name.str().c_str()));
      etaphieff_name.str("");

      c1_name << "canvas1" << ich << "_" << lv2 ;
      c1[lv2][ich] = new TCanvas(c1_name.str().c_str(),"",800,600);
      c1_name.str("");
      if(lv2==0) title << "Data   " << chain_name_list[ich] << "   L1/probe   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==1) title << "Data   " << chain_name_list[ich] << "   SA/L1   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==2) title << "Data   " << chain_name_list[ich] << "   Comb/SA   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==3) title << "Data   " << chain_name_list[ich] << "   EF/Comb   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==4) title << "Data   " << chain_name_list[ich] << "   HLT/probe   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;

      etaphiEfficiency[lv2][ich]->SetTitle(title.str().c_str());
      etaphiEfficiency[lv2][ich]->GetZaxis()->SetRangeUser(0,1);
      etaphiEfficiency[lv2][ich]->Draw("colz");
      etaphi_plot_name << "outputScaleFactor/result_etaphi_efficiency/" << level[lv2] << "etaphiEfficiency_" << chain_name_list[ich] << "_" << period << ".pdf"; 
      c1[lv2][ich]->SaveAs(etaphi_plot_name.str().c_str());
      etaphi_plot_name.str("");
      title.str("");

      c2_name << "canvas2" << ich << "_" << lv2 ;
      c2[lv2][ich] = new TCanvas(c2_name.str().c_str(),"",800,600);
      c2_name.str("");
      if(lv2==0) title << "MC   " << chain_name_list[ich] << "   L1/probe   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==1) title << "MC   " << chain_name_list[ich] << "   SA/L1   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==2) title << "MC   " << chain_name_list[ich] << "   Comb/SA   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==3) title << "MC   " << chain_name_list[ich] << "   EF/Comb   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      else if(lv2==4) title << "MC   " << chain_name_list[ich] << "   HLT/probe   " << "Barrel:p_{T}>" << m_util.returnPlateau(chain_list[ich],0) << "GeV Endcap:p_{T}>" << m_util.returnPlateau(chain_list[ich],1) << "  period " << period;
      
      etaphiEfficiencyMC[lv2][ich]->SetTitle(title.str().c_str());
      etaphiEfficiencyMC[lv2][ich]->GetZaxis()->SetRangeUser(0,1);
      etaphiEfficiencyMC[lv2][ich]->Draw("colz");
      etaphi_plot_name << "outputScaleFactor/result_etaphi_efficiency/" << level[lv2] << "etaphiEfficiencyMC_" << chain_name_list[ich] << "_" << period << ".pdf"; 
      c2[lv2][ich]->SaveAs(etaphi_plot_name.str().c_str());
      etaphi_plot_name.str("");
      title.str("");
      sf_plot_name << level[lv2] << "ScaleFactor_" << chain_name_list[ich] ;
      ScaleFactor[lv2][ich] = new TH2F(sf_plot_name.str().c_str(),";#eta;#phi",eta_nbin-1,eta_bin,phi_nbin-1,phi_bin);
      sf_plot_name.str("");
      for (int ibin=0;ibin<eta_nbin-1;ibin++){
        for (int jbin=0;jbin<phi_nbin-1;jbin++){
          int etaphibin = etaphiEfficiency[lv2][ich]->GetBin(ibin+1,jbin+1);
          float eff_data = etaphiEfficiency[lv2][ich]->GetBinContent(etaphibin); 
          float eff_mc = etaphiEfficiencyMC[lv2][ich]->GetBinContent(etaphibin); 
          double ratio = eff_data/eff_mc;
          ScaleFactor[lv2][ich]->SetBinContent(etaphibin,ratio);
        }
      }
      title << "SF Data/MC   " << chain_name_list[ich] << "  period " << period << "   " << level[lv2] << "  period " << period;
      ScaleFactor[lv2][ich]->SetTitle(title.str().c_str());
      title.str("");
      ScaleFactor[lv2][ich]->GetZaxis()->SetRangeUser(0,1.2);
      c3_name << "canvas3" << lv2 << "_" << ich; 
      c3[lv2][ich] = new TCanvas(c3_name.str().c_str(),"",800,600);
      ScaleFactor[lv2][ich]->Draw("colz");
      sf_plot_name << "outputScaleFactor/result_SF/" << level[lv2] << "ScaleFactor_" << chain_name_list[ich] << "_" << period <<  ".pdf" ;
      c3[lv2][ich]->SaveAs(sf_plot_name.str().c_str());
      sf_plot_name.str("");
      c3_name.str("");
    }
  }

}

////////////////////////////////////////////////////////////
void makeEffPlot::makeComPtEfficiencyPlot(const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, vector<string> com_per_list){

  int nChain = chain_list.size();
  int nPeriod = com_per_list.size();
  cout << "size of period is " << nPeriod << endl;

  TGraphErrors *ptEfficiency[nPeriod], *ptEfficiencyMC, *g_ptratio[nPeriod];
  TCanvas *c1;
  stringstream pteff_name,c1_name,pt_plot_name,pt_plot_title,input_name;
  TFile *fileMC = TFile::Open("outputScaleFactor/output/mc15_13TeV_Zmumu.root") ;
  TFile *file;
  vector<int> divpt;
  vector<float> vec_ptcenter, vec_pterr;
  int nbin;
  for (int iside=0;iside<2;iside++){
    for (int ich=0;ich<nChain;ich++){
      divpt.clear();vec_ptcenter.clear();vec_pterr.clear();
      divpt = m_util.devidePt(chain_list[ich]);
      vec_ptcenter = m_util.getPtCenter(chain_list[ich]);
      vec_pterr = m_util.getPtError(chain_list[ich]);
      nbin = divpt.size()-1;
      for(int lv2=0;lv2<5;lv2++){ 
        for (int per=0;per<nPeriod;per++){
          input_name << "outputScaleFactor/output/data16_13TeV_period_" << com_per_list[per] << ".root";
          file = TFile::Open(input_name.str().c_str());
          input_name.str("");
          pteff_name << level[lv2] << "PtEfficiency_" << chain_name_list[ich] << "_" << side[iside]; 
          ptEfficiency[per] = dynamic_cast<TGraphErrors*>(file->Get(pteff_name.str().c_str()));
          ptEfficiencyMC = dynamic_cast<TGraphErrors*>(fileMC->Get(pteff_name.str().c_str()));
          pteff_name.str("");

          Double_t *pteffData = ptEfficiency[per]->GetY();
          Double_t *pteffMC = ptEfficiencyMC->GetY();
          Double_t *ptefferrData = ptEfficiency[per]->GetEY();
          Double_t *ptefferrMC = ptEfficiencyMC->GetEY();
          Double_t *ptbin = ptEfficiency[per]->GetX();
          Double_t *ptbinerr = ptEfficiency[per]->GetEX();

          float ptRatio[nbin],ptRatioErr[nbin], ptx[nbin], ptxerr[nbin];
          for (int i=0;i<nbin;i++){
            float pteff1 = pteffData[i];
            float pteff2 = pteffMC[i];
            float ptefferr1 = ptefferrData[i];
            float ptefferr2 = ptefferrMC[i];
            ptx[i]=ptbin[i];
            ptxerr[i]=ptbinerr[i];
            float ptratio = (pteffMC[i]>1e-5)? pteff1/pteff2 : 0;
            float ptratioerr = 
              (pteffMC[i]>1e-5)? sqrt(pteff2*pteff2*ptefferr1*ptefferr1 + pteff1*pteff1*ptefferr2*ptefferr2)/(pteff2*pteff2) : 0;
            ptRatio[i]=ptratio;
            ptRatioErr[i]=ptratioerr;
          }
          g_ptratio[per] = new TGraphErrors(nbin,ptx,ptRatio,ptxerr,ptRatioErr);
        }

        c1_name << "canvas" << iside << "_" << ich << "_" << lv2 ;
        c1 = new TCanvas(c1_name.str().c_str(),"",800,1000);
        c1_name.str("");
        TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
        p1->Draw();
        p1->SetNumber(1);
        p1->SetBottomMargin(0.1);
        TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
        p2->SetBottomMargin(0.2);
        p2->Draw();
        p2->SetNumber(2);
        c1->cd(1);
        pt_plot_title << chain_name_list[ich] ;
        ptEfficiencyMC->SetTitle(pt_plot_title.str().c_str());
        ptEfficiencyMC->GetYaxis()->SetRangeUser(0,1.2);
        ptEfficiencyMC->GetYaxis()->SetTitle("Efficiency");
        ptEfficiencyMC->GetXaxis()->SetTitle("p_{T}(GeV)");
        ptEfficiencyMC->Draw("ap");
        ptEfficiency[0]->SetLineColor(2);
        ptEfficiency[0]->Draw("p");
        ptEfficiency[1]->SetLineColor(4);
        ptEfficiency[1]->Draw("p");
        TLatex test;
        test.SetTextSize(0.03);
        test.DrawLatex(10,1.15,"ATLAS work in progress");
        if(iside==0) test.DrawLatex(10,1.08,"Z Tag&Probe |#eta|<1.05");
        else test.DrawLatex(10,1.08,"Z Tag&Probe |#eta|>1.05");
        test.SetTextSize(0.05);
        if(lv2==0) test.DrawLatex(70,1.1,"L1/probe");
        else if(lv2==1) test.DrawLatex(70,1.1,"SA/L1");
        else if(lv2==2) test.DrawLatex(70,1.1,"Comb/SA");
        else if(lv2==3) test.DrawLatex(70,1.1,"EF/Comb");
        else if(lv2==4) test.DrawLatex(70,1.1,"HLT/probe");
        test.DrawLatex(70,0.4,"MC Zmumu");
        test.SetTextColor(2);
        test.DrawLatex(70,0.3,"Data16 periodA");
        test.SetTextColor(4);
        test.DrawLatex(70,0.2,"Data16 periodB");

        c1->cd(2);
        g_ptratio[0]->SetTitle("");
        g_ptratio[0]->GetYaxis()->SetTitle("Ratio Data/MC");
        g_ptratio[0]->GetYaxis()->SetTitleOffset(0.25);
        g_ptratio[0]->GetYaxis()->SetTitleSize(0.13);
        g_ptratio[0]->GetYaxis()->SetLabelSize(0.07);
        g_ptratio[0]->GetYaxis()->SetRangeUser(0.7,1.3);
        g_ptratio[0]->GetXaxis()->SetRangeUser(0,100);
        g_ptratio[0]->GetXaxis()->SetLabelSize(0);
        g_ptratio[0]->SetLineColor(2);
        g_ptratio[0]->Draw("ap");
        g_ptratio[1]->SetLineColor(4);
        g_ptratio[1]->Draw("p");
        TLine *line = new TLine(0,1,100,1);
        line->SetLineColor(1);
        line->SetLineStyle(2);
        line->Draw();
        pt_plot_name << "outputScaleFactor/result_compare_pt_efficiency/" << level[lv2] << "PtEfficiency_" << chain_name_list[ich] << "_" << side[iside] << ".pdf"; 
        c1->SaveAs(pt_plot_name.str().c_str());
        pt_plot_name.str("");
        pt_plot_title.str("");
        
      }
    }

  }
}

///////////////////////////////////////////////////////
void makeEffPlot::makeComEtaEfficiencyPlot(const char* input_mc, vector<int> &chain_list, vector<string> &chain_name_list, vector<string> com_per_list){

  int nChain = chain_list.size();
  int nPeriod = com_per_list.size();
  cout << "size of period is " << nPeriod << endl;
  TGraphErrors *etaEfficiency[nPeriod], *etaEfficiencyMC, *g_etaratio[nPeriod];
  TCanvas *c1;
  stringstream etaeff_name,c1_name,eta_plot_name,eta_plot_title,input_name;
  TFile *fileMC = TFile::Open("outputScaleFactor/output/mc15_13TeV_Zmumu.root") ;
  TFile *file;

  int nbin = eta_nbin-1;

  for (int ich=0;ich<nChain;ich++){
    for(int lv2=0;lv2<5;lv2++){ 
      for (int per=0;per<nPeriod;per++){
        input_name << "outputScaleFactor/output/data16_13TeV_period_" << com_per_list[per] << ".root";
        file = TFile::Open(input_name.str().c_str());
        input_name.str("");
        etaeff_name << level[lv2] << "etaEfficiency_" << chain_name_list[ich] ; 
        etaEfficiency[per] = dynamic_cast<TGraphErrors*>(file->Get(etaeff_name.str().c_str()));
        etaEfficiencyMC = dynamic_cast<TGraphErrors*>(fileMC->Get(etaeff_name.str().c_str()));
        etaeff_name.str("");

        Double_t *etaeffData = etaEfficiency[per]->GetY();
        Double_t *etaeffMC = etaEfficiencyMC->GetY();
        Double_t *etaefferrData = etaEfficiency[per]->GetEY();
        Double_t *etaefferrMC = etaEfficiencyMC->GetEY();
        Double_t *etabin = etaEfficiency[per]->GetX();
        Double_t *etabinerr = etaEfficiency[per]->GetEX();

        float etaRatio[nbin],etaRatioErr[nbin], etax[nbin], etaxerr[nbin];
        for (int i=0;i<nbin;i++){
          float etaeff1 = etaeffData[i];
          float etaeff2 = etaeffMC[i];
          float etaefferr1 = etaefferrData[i];
          float etaefferr2 = etaefferrMC[i];
          etax[i]=etabin[i];
          etaxerr[i]=etabinerr[i];
          float etaratio = (etaeffMC[i]>1e-5)? etaeff1/etaeff2 : 0;
          float etaratioerr = 
            (etaeffMC[i]>1e-5)? sqrt(etaeff2*etaeff2*etaefferr1*etaefferr1 + etaeff1*etaeff1*etaefferr2*etaefferr2)/(etaeff2*etaeff2) : 0;
          etaRatio[i]=etaratio;
          etaRatioErr[i]=etaratioerr;
        }
        g_etaratio[per] = new TGraphErrors(nbin,etax,etaRatio,etaxerr,etaRatioErr);
      }
      c1_name << "canvas" << "_" << ich << "_" << lv2 ;
      c1 = new TCanvas(c1_name.str().c_str(),"",800,1000);
      c1_name.str("");
      TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
      p1->Draw();
      p1->SetNumber(1);
      p1->SetBottomMargin(0.1);
      TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
      p2->SetBottomMargin(0.2);
      p2->Draw();
      p2->SetNumber(2);
      c1->cd(1);
      eta_plot_title << chain_name_list[ich] ;
      etaEfficiencyMC->SetTitle(eta_plot_title.str().c_str());
      etaEfficiencyMC->GetYaxis()->SetRangeUser(0,1.2);
      etaEfficiencyMC->GetYaxis()->SetTitle("Efficiency");
      etaEfficiencyMC->GetXaxis()->SetTitle("#eta");
      etaEfficiencyMC->Draw("ap");
      etaEfficiency[0]->SetLineColor(2);
      etaEfficiency[0]->Draw("p");
      etaEfficiency[1]->SetLineColor(4);
      etaEfficiency[1]->Draw("p");
      TLatex test;
      test.SetTextSize(0.03);
      test.DrawLatex(-2.8,1.15,"ATLAS work in progress");
      test.DrawLatex(-2.8,1.08,"Z Tag&Probe");
      test.DrawLatex(-1.1,1.08,Form("Barrel p_{T}>%dGeV",m_util.returnPlateau(chain_list[ich],0)));
      test.DrawLatex(-1.1,1.0,Form("Endcap p_{T}>%dGeV",m_util.returnPlateau(chain_list[ich],1)));
      test.SetTextSize(0.05);
      if(lv2==0) test.DrawLatex(1,1.1,"L1/probe");
      else if(lv2==1) test.DrawLatex(1,1.1,"SA/L1");
      else if(lv2==2) test.DrawLatex(1,1.1,"Comb/SA");
      else if(lv2==3) test.DrawLatex(1,1.1,"EF/Comb");
      else if(lv2==4) test.DrawLatex(1,1.1,"HLT/probe");
      test.DrawLatex(1,0.4,"MC Zmumu");
      test.SetTextColor(2);
      test.DrawLatex(1,0.3,"Data16 periodA");
      test.SetTextColor(4);
      test.DrawLatex(1,0.2,"Data16 periodB");

      c1->cd(2);
      g_etaratio[0]->SetTitle("");
      g_etaratio[0]->GetYaxis()->SetTitle("Ratio Data/MC");
      g_etaratio[0]->GetYaxis()->SetTitleOffset(0.25);
      g_etaratio[0]->GetYaxis()->SetTitleSize(0.13);
      g_etaratio[0]->GetYaxis()->SetLabelSize(0.07);
      g_etaratio[0]->GetYaxis()->SetRangeUser(0.7,1.3);
      g_etaratio[0]->GetXaxis()->SetRangeUser(-2.6,2.6);
      g_etaratio[0]->GetXaxis()->SetLabelSize(0);
      g_etaratio[0]->SetLineColor(2);
      g_etaratio[0]->Draw("ap");
      g_etaratio[1]->SetLineColor(4);
      g_etaratio[1]->Draw("p");
      TLine *line = new TLine(-2.6,1,2.6,1);
      line->SetLineColor(1);
      line->SetLineStyle(2);
      line->Draw();
      eta_plot_name << "outputScaleFactor/result_compare_eta_efficiency/" << level[lv2] << "etaEfficiency_" << chain_name_list[ich] << ".pdf"; 
      c1->SaveAs(eta_plot_name.str().c_str());
      eta_plot_name.str("");
      eta_plot_title.str("");
    }
  }
}
///////////////////////////////////////////////////////
void makeEffPlot::makeComPhiEfficiencyPlot(const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, vector<string> com_per_list){

  int nChain = chain_list.size();
  int nPeriod = com_per_list.size();
  cout << "size of period is " << nPeriod << endl;

  TGraphErrors *phiEfficiency[nPeriod], *phiEfficiencyMC, *g_phiratio[nPeriod];
  TCanvas *c1;
  stringstream phieff_name,c1_name,phi_plot_name,phi_plot_title,input_name;
  TFile *fileMC = TFile::Open("outputScaleFactor/output/mc15_13TeV_Zmumu.root") ;
  TFile *file;
  int nbin = phi_nbin-1;
  for (int iside=0;iside<2;iside++){
    for (int ich=0;ich<nChain;ich++){
      for(int lv2=0;lv2<5;lv2++){ 
        for (int per=0;per<nPeriod;per++){
          input_name << "outputScaleFactor/output/data16_13TeV_period_" << com_per_list[per] << ".root";
          file = TFile::Open(input_name.str().c_str());
          input_name.str("");
          phieff_name << level[lv2] << "phiEfficiency_" << chain_name_list[ich] << "_" << side[iside]; 
          phiEfficiency[per] = dynamic_cast<TGraphErrors*>(file->Get(phieff_name.str().c_str()));
          phiEfficiencyMC = dynamic_cast<TGraphErrors*>(fileMC->Get(phieff_name.str().c_str()));
          phieff_name.str("");

          Double_t *phieffData = phiEfficiency[per]->GetY();
          Double_t *phieffMC = phiEfficiencyMC->GetY();
          Double_t *phiefferrData = phiEfficiency[per]->GetEY();
          Double_t *phiefferrMC = phiEfficiencyMC->GetEY();
          Double_t *phibin = phiEfficiency[per]->GetX();
          Double_t *phibinerr = phiEfficiency[per]->GetEX();

          float phiRatio[nbin],phiRatioErr[nbin], phix[nbin], phixerr[nbin];
          for (int i=0;i<nbin;i++){
            float phieff1 = phieffData[i];
            float phieff2 = phieffMC[i];
            float phiefferr1 = phiefferrData[i];
            float phiefferr2 = phiefferrMC[i];
            phix[i]=phibin[i];
            phixerr[i]=phibinerr[i];
            float phiratio = (phieffMC[i]>1e-5)? phieff1/phieff2 : 0;
            float phiratioerr = 
              (phieffMC[i]>1e-5)? sqrt(phieff2*phieff2*phiefferr1*phiefferr1 + phieff1*phieff1*phiefferr2*phiefferr2)/(phieff2*phieff2) : 0;
            phiRatio[i]=phiratio;
            phiRatioErr[i]=phiratioerr;
          }
          g_phiratio[per] = new TGraphErrors(nbin,phix,phiRatio,phixerr,phiRatioErr);
        }

        c1_name << "canvas" << iside << "_" << ich << "_" << lv2 ;
        c1 = new TCanvas(c1_name.str().c_str(),"",800,1000);
        c1_name.str("");
        TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
        p1->Draw();
        p1->SetNumber(1);
        p1->SetBottomMargin(0.1);
        TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
        p2->SetBottomMargin(0.2);
        p2->Draw();
        p2->SetNumber(2);
        c1->cd(1);
        phi_plot_title << chain_name_list[ich] ;
        phiEfficiencyMC->SetTitle(phi_plot_title.str().c_str());
        phiEfficiencyMC->GetYaxis()->SetRangeUser(0,1.2);
        phiEfficiencyMC->GetYaxis()->SetTitle("Efficiency");
        phiEfficiencyMC->GetXaxis()->SetTitle("#phi");
        phiEfficiencyMC->Draw("ap");
        phiEfficiency[0]->SetLineColor(2);
        phiEfficiency[0]->Draw("p");
        phiEfficiency[1]->SetLineColor(4);
        phiEfficiency[1]->Draw("p");
        TLatex test;
        test.SetTextSize(0.03);
        test.DrawLatex(-3.4,0.3,"ATLAS work in progress");
        if(iside==0) test.DrawLatex(-3.4,0.23,"Z Tag&Probe |#eta|<1.05");
        else test.DrawLatex(-3.4,0.23,"Z Tag&Probe |#eta|>1.05");
        test.DrawLatex(-3.4,0.16,Form("p_{T}>%dGeV",m_util.returnPlateau(chain_list[ich],iside)));
        test.SetTextSize(0.05);
        if(lv2==0) test.DrawLatex(1.3,1.1,"L1/probe");
        else if(lv2==1) test.DrawLatex(1.2,1.1,"SA/L1");
        else if(lv2==2) test.DrawLatex(1.2,1.1,"Comb/SA");
        else if(lv2==3) test.DrawLatex(1.2,1.1,"EF/Comb");
        else if(lv2==4) test.DrawLatex(1.2,1.1,"HLT/probe");
        test.DrawLatex(1.2,0.4,"MC Zmumu");
        test.SetTextColor(2);
        test.DrawLatex(1.2,0.3,"Data16 periodA");
        test.SetTextColor(4);
        test.DrawLatex(1.2,0.2,"Data16 periodB");

        c1->cd(2);
        g_phiratio[0]->SetTitle("");
        g_phiratio[0]->GetYaxis()->SetTitle("Ratio Data/MC");
        g_phiratio[0]->GetYaxis()->SetTitleOffset(0.25);
        g_phiratio[0]->GetYaxis()->SetTitleSize(0.13);
        g_phiratio[0]->GetYaxis()->SetLabelSize(0.07);
        g_phiratio[0]->GetYaxis()->SetRangeUser(0.7,1.3);
        g_phiratio[0]->GetXaxis()->SetRangeUser(-3.2,3.2);
        g_phiratio[0]->GetXaxis()->SetLabelSize(0);
        g_phiratio[0]->SetLineColor(2);
        g_phiratio[0]->Draw("ap");
        g_phiratio[1]->SetLineColor(4);
        g_phiratio[1]->Draw("p");
        TLine *line = new TLine(-3.2,1,3.2,1);
        line->SetLineColor(1);
        line->SetLineStyle(2);
        line->Draw();
        phi_plot_name << "outputScaleFactor/result_compare_phi_efficiency/" << level[lv2] << "phiEfficiency_" << chain_name_list[ich] << "_" << side[iside] << ".pdf"; 
        c1->SaveAs(phi_plot_name.str().c_str());
        phi_plot_name.str("");
        phi_plot_title.str("");

      }
    }
  }
}

////////////////////////////////////////////////////////////
void makeEffPlot::makeComPtVarConeEfficiencyPlot(const char *input){

  TGraphErrors *ptEfficiency[2], *g_ptratio;
  TCanvas *c1;
  stringstream pteff_name1,pteff_name2,c1_name,pt_plot_name,pt_plot_title,input_name;
  TFile *file = TFile::Open(input) ;

  for (int iside=1;iside<2;iside++){
    pteff_name1 << "HLTptvarconeEfficiency_HLT_mu26_ivarmedium" << "_" << side[iside]; 
    pteff_name2 << "HLTptvarconeEfficiency_HLT_mu26" << "_" << side[iside]; 
    ptEfficiency[0] = dynamic_cast<TGraphErrors*>(file->Get(pteff_name1.str().c_str()));
    ptEfficiency[1] = dynamic_cast<TGraphErrors*>(file->Get(pteff_name2.str().c_str()));
    pteff_name1.str("");
    pteff_name2.str("");

    Double_t *pteff = ptEfficiency[0]->GetY();
    Double_t *pteffiso = ptEfficiency[1]->GetY();
    Double_t *ptefferr = ptEfficiency[0]->GetEY();
    Double_t *ptefferriso = ptEfficiency[1]->GetEY();
    Double_t *ptbin = ptEfficiency[0]->GetX();
    Double_t *ptbinerr = ptEfficiency[0]->GetEX();

    float ptRatio[ptvarcone_nbin],ptRatioErr[ptvarcone_nbin], ptx[ptvarcone_nbin], ptxerr[ptvarcone_nbin];
    for (int i=0;i<ptvarcone_nbin;i++){
      float pteff1 = pteff[i];
      float pteff2 = pteffiso[i];
      float ptefferr1 = ptefferr[i];
      float ptefferr2 = ptefferriso[i];
      ptx[i]=ptbin[i];
      ptxerr[i]=ptbinerr[i];
      float ptratio = (pteff[i]>1e-5)? pteff1/pteff2 : 0;
      float ptratioerr = 
        (pteff[i]>1e-5)? sqrt(pteff2*pteff2*ptefferr1*ptefferr1 + pteff1*pteff1*ptefferr2*ptefferr2)/(pteff2*pteff2) : 0;
      ptRatio[i]=ptratio;
      ptRatioErr[i]=ptratioerr;
    }
    g_ptratio = new TGraphErrors(ptvarcone_nbin,ptx,ptRatio,ptxerr,ptRatioErr);

    c1_name << "canvas" << iside ;
    c1 = new TCanvas(c1_name.str().c_str(),"",800,1000);
    c1_name.str("");
    TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
    p1->Draw();
    p1->SetNumber(1);
    p1->SetBottomMargin(0.05);
    TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
    p2->SetBottomMargin(0.2);
    p2->Draw();
    p2->SetNumber(2);
    c1->cd(1);
    ptEfficiency[0]->GetYaxis()->SetRangeUser(0,1.1);
    ptEfficiency[0]->GetXaxis()->SetRangeUser(0,0.4);
    //ptEfficiency[0]->GetYaxis()->SetTitle("HLT Efficiency");
    ptEfficiency[0]->GetXaxis()->SetLabelSize(0);
    //ptEfficiency[0]->GetXaxis()->SetTitle("ptvarcone30/pt");
    //ptEfficiency[0]->GetXaxis()->SetTitleSize(0.05);
    ptEfficiency[0]->SetLineColor(2);
    ptEfficiency[0]->SetTitle("");
    ptEfficiency[0]->Draw("ap");
    ptEfficiency[1]->Draw("p");
    /*TLatex test;
    test.SetTextSize(0.03);
    test.DrawLatex(10,1.15,"ATLAS work in progress");
    if(iside==0) test.DrawLatex(10,1.08,"Z Tag&Probe |#eta|<1.05");
    else test.DrawLatex(10,1.08,"Z Tag&Probe |#eta|>1.05");
    test.SetTextSize(0.05);
    if(lv2==0) test.DrawLatex(70,1.1,"L1/probe");
    else if(lv2==1) test.DrawLatex(70,1.1,"SA/L1");
    else if(lv2==2) test.DrawLatex(70,1.1,"Comb/SA");
    else if(lv2==3) test.DrawLatex(70,1.1,"EF/Comb");
    else if(lv2==4) test.DrawLatex(70,1.1,"HLT/probe");
    test.DrawLatex(70,0.4,"MC Zmumu");
    test.SetTextColor(2);
    test.DrawLatex(70,0.3,"Data16 periodA");
    test.SetTextColor(4);
    test.DrawLatex(70,0.2,"Data16 periodB");
    */

    c1->cd(2);
    g_ptratio->SetTitle("");
    //g_ptratio->GetYaxis()->SetTitle("mu26_ivarmeidum / mu26");
    //g_ptratio->GetXaxis()->SetTitle("ptvarcone30 / pt");
    g_ptratio->GetYaxis()->SetTitleOffset(0.25);
    g_ptratio->GetYaxis()->SetTitleSize(0.08);
    g_ptratio->GetYaxis()->SetLabelSize(0.1);
    g_ptratio->GetYaxis()->SetRangeUser(0.4,1.1);
    g_ptratio->GetXaxis()->SetRangeUser(0,0.4);
    g_ptratio->GetXaxis()->SetLabelSize(0.15);
    g_ptratio->Draw("ap");
    TLine *line = new TLine(0,1,0.4,1);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();
    pt_plot_name << "outputScaleFactor/result_compare_ptvarcone_efficiency/" << "PtVarCone30Efficiency_" << side[iside] << "_periodB.pdf"; 
    c1->SaveAs(pt_plot_name.str().c_str());
    pt_plot_name.str("");
    pt_plot_title.str("");
  }
        
}

