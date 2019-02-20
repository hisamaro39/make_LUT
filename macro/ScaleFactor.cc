#define ScaleFactor_cxx
#include "macro/ScaleFactor.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <iostream>
#include <sstream>
#include <string>
#include "TTree.h"
#include "TF1.h"
using namespace std;

#define PI 3.14159265258979
double const ZERO_LIMIT = 1e-5;

void ScaleFactor:: ScaleFactor::getThreshold(const char *list,vector<int> &chain_list,vector<string> &chain_name_list){
  TChain *fChain = new TChain("validationT");
  ifstream finlist(list);
  string file_rec;
  while(finlist>>file_rec) 
    fChain->Add(file_rec.c_str());
  TTree *tree = static_cast<TTree*>(fChain);
  Init(tree);
  fChain->GetEntry(0); 

  chain_list.clear();chain_name_list.clear();
  for (int i=0;i<probe_trigger_threshold->size();i++){ 
    chain_list.push_back(probe_trigger_threshold->at(i));
    chain_name_list.push_back(probe_trigger_chain->at(i));
  }
  
}

void ScaleFactor::Loop(vector<string> list,string output)
{
  TChain *fChain = new TChain("validationT");
  for (int in=0;in<list.size();in++) {
    cout << "input is " << list[in].c_str() << endl;
    fChain->Add(list[in].c_str());
  }
  TTree *tree = static_cast<TTree*>(fChain);
  Init(tree);

  TFile *fout = new TFile(output.c_str(),"recreate");
  fChain->GetEntry(0); 
  int nChain = probe_trigger_chain->size();
  string chain[nChain];
  vector<int> threshold;
  threshold.clear();
  for (int ith=0;ith<nChain;ith++) {
    threshold.push_back(probe_trigger_threshold->at(ith));
    chain[ith] = probe_trigger_chain->at(ith);
  }

 //// eta bins setup 
  const double eta_bin[]  = {-2.4, -2.159, -1.7705, -1.4855, -1.29045, -1.1904, -1.05, -0.979, -0.8495, -0.7215, -0.564, -0.40, -0.228, -0.066, 0., 0.066, 0.228, 0.40, 0.564, 0.7215, 0.8495, 0.979, 1.05, 1.1904, 1.29045, 1.4855, 1.7705, 2.159, 2.4};
  const int eta_nbin      = sizeof(eta_bin)/sizeof(double);
 //// phi bins setup 
  const double phi_bin[]  = {-16./16.*PI, -14./16.*PI, -12./16.*PI,  -10./16.*PI, -8./16.*PI,-6./16.*PI,-4./16.*PI, -2./16.*PI, 0./16.*PI, 2./16.*PI, 4./16.*PI, 6./16.*PI, 8./16.*PI, 10./16.*PI, 12./16.*PI, 14./16.*PI, 16./16.*PI};
  const int phi_nbin      = sizeof(phi_bin)/sizeof(double);

  const int divpt[]  = {0,20,30,40,50,60,80,100};
  const int pt_nbin      = sizeof(divpt)/sizeof(int);

  //for phi efficiency
  int nbinphi=32, divphi[nbinphi+1];
  float phicenter[nbinphi],phierr[nbinphi];
  for (int iphi=0;iphi<nbinphi+1;iphi++){
    divphi[iphi] = iphi;
    if(iphi==nbinphi) continue;
    phicenter[iphi] = -3.1 + 0.2*iphi;
    phierr[iphi] = 0.1;
  }

  //pt plot
  TH1 *probe_pt[2],*pass_pt[4][2][nChain];
  stringstream probe_name,pt_name[4][2][nChain];
  string level[5]={"L1","SA","Comb","EF","HLT"};
  string side[2]={"Barrel","Endcap"};
  for (int iside=0;iside<2;iside++){
    probe_name << "ProbePt_" << side[iside] ;
    probe_pt[iside] = new TH1F(probe_name.str().c_str(),";p_{T}(GeV);Number of events",100,0,100);
    probe_name.str("");
    for (int ich=0;ich<nChain;ich++){
      for (int lv1=0;lv1<4;lv1++) {
        pt_name[lv1][iside][ich] << level[lv1] << "PassPt_" << probe_trigger_chain->at(ich) << "_" << side[iside];
        pass_pt[lv1][iside][ich] = new TH1F(pt_name[lv1][iside][ich].str().c_str(),";p_{T}(GeV);Number of events",100,0,100);
        pt_name[lv1][iside][ich].str("");
      }
    }
  }

  //eta plot
  TH1 *probe_eta[nChain];
  TH1 *pass_eta[4][nChain];
  stringstream eta_name[4][nChain];
  for (int ch=0;ch<nChain;ch++){
    probe_name << "ProbeEta_" << probe_trigger_chain->at(ch);
    probe_eta[ch] = new TH1F(probe_name.str().c_str(),";#eta;Number of events",eta_nbin-1,eta_bin);
    probe_name.str("");
    for(int lev=0;lev<4;lev++){
      eta_name[lev][ch] << level[lev] << "PassEta_" << probe_trigger_chain->at(ch) ;
      pass_eta[lev][ch] = new TH1F(eta_name[lev][ch].str().c_str(),";#eta;Number of events",eta_nbin-1,eta_bin);
      eta_name[lev][ch].str("");
    }
  }

  //phi plot
  TH1 *probe_phi[2][nChain];
  TH1 *pass_phi[4][2][nChain];
  stringstream phi_name[4][2][nChain];
  for (int iside=0;iside<2;iside++){
    for (int ch=0;ch<nChain;ch++){
      probe_name << "Probephi_" << probe_trigger_chain->at(ch) << "_" << side[iside];
      probe_phi[iside][ch] = new TH1F(probe_name.str().c_str(),";#phi;Number of events",phi_nbin-1,phi_bin);
      probe_name.str("");
      for(int lev=0;lev<4;lev++){
        phi_name[lev][iside][ch] << level[lev] << "PassPhi_" << probe_trigger_chain->at(ch) << "_" << side[iside];
        pass_phi[lev][iside][ch] = new TH1F(phi_name[lev][iside][ch].str().c_str(),";#phi;Number of events",phi_nbin-1,phi_bin);
        phi_name[lev][iside][ch].str("");
      }
    }
  }

  //eta-phi plot
  TH2 *probe_eta_phi[nChain];
  TH2 *pass_eta_phi[4][nChain];
  stringstream eta_phi_name[4][nChain];
  for (int ch=0;ch<nChain;ch++){
    probe_name << "ProbeEtaPhi_" << probe_trigger_chain->at(ch);
    probe_eta_phi[ch] = new TH2F(probe_name.str().c_str(),";#eta;#phi",eta_nbin-1,eta_bin,phi_nbin-1,phi_bin);
    probe_name.str("");
    for(int lev=0;lev<4;lev++){
      eta_phi_name[lev][ch] << level[lev] << "PassEtaPhi_" << probe_trigger_chain->at(ch) ;
      pass_eta_phi[lev][ch] = new TH2F(eta_phi_name[lev][ch].str().c_str(),";#eta;#phi",eta_nbin-1,eta_bin,phi_nbin-1,phi_bin);
      eta_phi_name[lev][ch].str("");
    }
  }

  //Pt residual
  TH1 *PtResidualChain[3][2][pt_nbin];
  string Chain[3]={"SA","Comb","EF"};
  string Region[2]={"Barrel","Endcap"};
  stringstream name;
  for(int ip=0;ip<pt_nbin-1;ip++){
    for(int ii=0;ii<3;ii++){
      for (int re=0;re<2;re++){
        name << "PtResidual" << Chain[ii] << Region[re] <<  "Pt" << divpt[ip] << "_" << divpt[ip+1];
        if(ii==2) PtResidualChain[ii][re][ip] = new TH1F(name.str().c_str(),";p_{T}residual;Number of events",200,-0.2,0.2);
        else PtResidualChain[ii][re][ip] = new TH1F(name.str().c_str(),";p_{T}residual;Number of events",500,-1,1);
        name.str("");
      }
    }
  }

  //TH2 *PtVsOfflineL1Dr = new TH2F("PtVsOfflineL1Dr",";p_{T}(GeV);#DeltaR_{offline,L1}",100,0,50,100,0,0.5);
  
  //make pt,eta,phi,eta-phi plot
  int nentries = fChain->GetEntries();
  cout << "Number of events is " << nentries << endl;
  
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%100000==0) cout << "entry=" << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry); 
    float comb_pt=probe_comb_pt/1000;
    float pt_residual_sa = (probe_sa_pt>1e-5)? 1-probe_offline_pt/probe_sa_pt : 1000;
    float pt_residual_comb = (comb_pt>1e-5)? 1-probe_offline_pt/comb_pt : 1000;
    float pt_residual_ef = (probe_ef_pt>1e-5)? 1-probe_offline_pt/probe_ef_pt : 1000;
    //cout << "pt sa/comb/ef/offline=" << probe_sa_pt << "/" << probe_comb_pt << "/" << probe_ef_pt << "/" << probe_offline_pt << endl;
    //cout << "pt residual sa/comb/ef=" << pt_residual_sa << "/" << pt_residual_comb << "/" << pt_residual_ef << endl;
    if(pt_residual_sa<1000) {
      for (int jp=0;jp<pt_nbin-1;jp++){
        if(probe_offline_pt<divpt[jp] || probe_offline_pt>divpt[jp+1]) continue;
        if(fabs(probe_offline_eta)<1.05) PtResidualChain[0][0][jp]->Fill(pt_residual_sa);
        else  PtResidualChain[0][1][jp]->Fill(pt_residual_sa);
      }
    }
    if(pt_residual_comb<1000) {
      for (int jp=0;jp<pt_nbin-1;jp++){
        if(probe_offline_pt<divpt[jp] || probe_offline_pt>divpt[jp+1]) continue;
        if(fabs(probe_offline_eta)<1.05) PtResidualChain[1][0][jp]->Fill(pt_residual_comb);
        else  PtResidualChain[1][1][jp]->Fill(pt_residual_comb);
      }
    }
    if(pt_residual_ef<1000) {
      for (int jp=0;jp<pt_nbin-1;jp++){
        if(probe_offline_pt<divpt[jp] || probe_offline_pt>divpt[jp+1]) continue;
        if(fabs(probe_offline_eta)<1.05 )PtResidualChain[2][0][jp]->Fill(pt_residual_ef);
        else  PtResidualChain[2][1][jp]->Fill(pt_residual_ef);
      }
    }
/*
    if (*tag_trigger_chain!="HLT_mu24_ivarmedium") continue;
    //float dr = m_util.calcDr(probe_offline_eta,probe_roi_eta,probe_offline_phi,probe_roi_phi);
    //PtVsOfflineL1Dr->Fill(probe_offline_pt,dr);
    //if (*tag_trigger_chain!="HLT_mu20_iloose_L1MU15") continue;
    int br_or_ec = (fabs(probe_offline_eta)<1.05)? 0 : 1;
    //make pt plot
    probe_pt[br_or_ec]->Fill(probe_offline_pt);
    for(int ich=0;ich<nChain;ich++){
      if(probe_passL1->at(ich)) pass_pt[0][br_or_ec][ich]->Fill(probe_offline_pt);
      if(probe_passSA->at(ich)) pass_pt[1][br_or_ec][ich]->Fill(probe_offline_pt);
      if(probe_passComb->at(ich)) pass_pt[2][br_or_ec][ich]->Fill(probe_offline_pt);
      if(probe_passEF->at(ich)) pass_pt[3][br_or_ec][ich]->Fill(probe_offline_pt);
      //make eta plot
      //if(probe_offline_pt < probe_trigger_threshold->at(ich)) continue;
      if(probe_offline_pt < m_util.returnPlateau(threshold[ich],br_or_ec) ) continue;//temporal
      probe_eta[ich]->Fill(probe_offline_eta);
      if(probe_passL1->at(ich)) pass_eta[0][ich]->Fill(probe_offline_eta);
      if(probe_passSA->at(ich)) pass_eta[1][ich]->Fill(probe_offline_eta);
      if(probe_passComb->at(ich)) pass_eta[2][ich]->Fill(probe_offline_eta);
      if(probe_passEF->at(ich)) pass_eta[3][ich]->Fill(probe_offline_eta);
      //make phi plot
      probe_phi[br_or_ec][ich]->Fill(probe_offline_phi);
      if(probe_passL1->at(ich)) pass_phi[0][br_or_ec][ich]->Fill(probe_offline_phi);
      if(probe_passSA->at(ich)) pass_phi[1][br_or_ec][ich]->Fill(probe_offline_phi);
      if(probe_passComb->at(ich)) pass_phi[2][br_or_ec][ich]->Fill(probe_offline_phi);
      if(probe_passEF->at(ich)) pass_phi[3][br_or_ec][ich]->Fill(probe_offline_phi);
      
      //make eta-phi plot
      probe_eta_phi[ich]->Fill(probe_offline_eta,probe_offline_phi);
      if(probe_passL1->at(ich)) pass_eta_phi[0][ich]->Fill(probe_offline_eta,probe_offline_phi);
      if(probe_passSA->at(ich)) pass_eta_phi[1][ich]->Fill(probe_offline_eta,probe_offline_phi);
      if(probe_passComb->at(ich)) pass_eta_phi[2][ich]->Fill(probe_offline_eta,probe_offline_phi);
      if(probe_passEF->at(ich)) pass_eta_phi[3][ich]->Fill(probe_offline_eta,probe_offline_phi);
    }
    */
  }

  //calculate pt efficiency
  /*int numEvPt[5];
  TGraphErrors *ptEfficiency[5][2][nChain];
  stringstream pteff_name;
  int nbinpt;
  vector<int> divpt;
  vector<float> vec_ptcenter,vec_pterr;
  for (int ich=0;ich<nChain;ich++){
    divpt.clear();vec_ptcenter.clear();vec_pterr.clear();
    divpt = m_util.devidePt(threshold[ich]);
    vec_ptcenter = m_util.getPtCenter(threshold[ich]);
    vec_pterr = m_util.getPtError(threshold[ich]);
    nbinpt = divpt.size()-1;
    for (int iside=0;iside<2;iside++){
      float pteff[5][nbinpt],ptefferr[5][nbinpt];//lv/side/chain/ptbin
      float ptcenter[nbinpt],pterr[nbinpt];
      for (int ibin=0;ibin<nbinpt;ibin++){
        ptcenter[ibin] = vec_ptcenter[ibin];
        pterr[ibin] = vec_pterr[ibin];
        numEvPt[0] = probe_pt[iside]->Integral(divpt[ibin],divpt[ibin+1]);
        numEvPt[1] = pass_pt[0][iside][ich]->Integral(divpt[ibin],divpt[ibin+1]);
        numEvPt[2] = pass_pt[1][iside][ich]->Integral(divpt[ibin],divpt[ibin+1]);
        numEvPt[3] = pass_pt[2][iside][ich]->Integral(divpt[ibin],divpt[ibin+1]);
        numEvPt[4] = pass_pt[3][iside][ich]->Integral(divpt[ibin],divpt[ibin+1]);
        for (int lv=0;lv<4;lv++){
          pteff[lv][ibin] = (numEvPt[lv]>0)? 1.*numEvPt[lv+1]/numEvPt[lv] : 0;
          ptefferr[lv][ibin] = (numEvPt[lv]>0)? sqrt(pteff[lv][ibin]*(1-pteff[lv][ibin]) / numEvPt[lv]) : 0;
        }
        //HLT/probe
        pteff[4][ibin] = (numEvPt[0]>0)? 1.*numEvPt[4]/numEvPt[0] : 0;
        ptefferr[4][ibin] = (numEvPt[0]>0)? sqrt(pteff[4][ibin]*(1-pteff[4][ibin]) / numEvPt[0]) : 0;
      }
      for(int lv2=0;lv2<5;lv2++){ 
        pteff_name << level[lv2] << "PtEfficiency_" << probe_trigger_chain->at(ich) << "_" << side[iside]; 
        ptEfficiency[lv2][iside][ich] = new TGraphErrors(nbinpt,ptcenter,pteff[lv2],pterr,ptefferr[lv2]);
        ptEfficiency[lv2][iside][ich]->SetName(pteff_name.str().c_str());
        fout->Add(ptEfficiency[lv2][iside][ich]);
        pteff_name.str("");
      }
    }
  }

  //calculate eta efficiency
  int numEvEta[5];
  float etaeff[5][nChain][eta_nbin-1],etaefferr[5][nChain][eta_nbin-1];//lv/chain/etabin
  float eta_center[eta_nbin-1],eta_err[eta_nbin-1];
  TGraphErrors *etaEfficiency[5][nChain];
  stringstream etaeff_name;
  for (int ich=0;ich<nChain;ich++){
    for (int ibin=0;ibin<eta_nbin-1;ibin++){
      numEvEta[0] = probe_eta[ich]->GetBinContent(ibin+1);
      numEvEta[1] = pass_eta[0][ich]->GetBinContent(ibin+1);
      numEvEta[2] = pass_eta[1][ich]->GetBinContent(ibin+1);
      numEvEta[3] = pass_eta[2][ich]->GetBinContent(ibin+1);
      numEvEta[4] = pass_eta[3][ich]->GetBinContent(ibin+1);
      eta_center[ibin] = 0.5*(eta_bin[ibin]+eta_bin[ibin+1]);
      eta_err[ibin] = eta_bin[ibin+1]-eta_center[ibin];
      for (int lv=0;lv<4;lv++){
        etaeff[lv][ich][ibin] = (numEvEta[lv]>0)? 1.*numEvEta[lv+1]/numEvEta[lv] : 0;
        etaefferr[lv][ich][ibin] = (numEvEta[lv]>0)? sqrt(etaeff[lv][ich][ibin]*(1-etaeff[lv][ich][ibin]) / numEvEta[lv]) : 0;
      }
      //HLT/probe
      etaeff[4][ich][ibin] = (numEvEta[0]>0)? 1.*numEvEta[4]/numEvEta[0] : 0;
      etaefferr[4][ich][ibin] = (numEvEta[0]>0)? sqrt(etaeff[4][ich][ibin]*(1-etaeff[4][ich][ibin]) / numEvEta[0]) : 0;
    }
    for(int lv2=0;lv2<5;lv2++){ 
      etaeff_name << level[lv2] << "etaEfficiency_" << probe_trigger_chain->at(ich) ; 
      etaEfficiency[lv2][ich] = new TGraphErrors(eta_nbin-1,eta_center,etaeff[lv2][ich],eta_err,etaefferr[lv2][ich]);
      etaEfficiency[lv2][ich]->SetName(etaeff_name.str().c_str());
      fout->Add(etaEfficiency[lv2][ich]);
      etaeff_name.str("");
    }
  }

  
  //calculate phi efficiency
  int numEvPhi[5];
  float phieff[5][2][nChain][phi_nbin-1],phiefferr[5][2][nChain][phi_nbin-1];//lv/side/chain/phibin
  TGraphErrors *phiEfficiency[5][2][nChain];
  stringstream phieff_name;
  float phi_center[phi_nbin-1],phi_err[phi_nbin-1];
  for (int ich=0;ich<nChain;ich++){
    for (int iside=0;iside<2;iside++){
      for (int ibin=0;ibin<phi_nbin-1;ibin++){
        numEvPhi[0] = probe_phi[iside][ich]->GetBinContent(ibin+1);
        numEvPhi[1] = pass_phi[0][iside][ich]->GetBinContent(ibin+1);
        numEvPhi[2] = pass_phi[1][iside][ich]->GetBinContent(ibin+1);
        numEvPhi[3] = pass_phi[2][iside][ich]->GetBinContent(ibin+1);
        numEvPhi[4] = pass_phi[3][iside][ich]->GetBinContent(ibin+1);
        phi_center[ibin] = 0.5*(phi_bin[ibin]+phi_bin[ibin+1]);
        phi_err[ibin] = phi_bin[ibin+1]-phi_center[ibin];
        for (int lv=0;lv<4;lv++){
          phieff[lv][iside][ich][ibin] = (numEvPhi[lv]>0)? 1.*numEvPhi[lv+1]/numEvPhi[lv] : 0;
          phiefferr[lv][iside][ich][ibin] = (numEvPhi[lv]>0)? sqrt(phieff[lv][iside][ich][ibin]*(1-phieff[lv][iside][ich][ibin]) / numEvPhi[lv]) : 0;
        }
        //HLT/probe
        phieff[4][iside][ich][ibin] = (numEvPhi[0]>0)? 1.*numEvPhi[4]/numEvPhi[0] : 0;
        phiefferr[4][iside][ich][ibin] = (numEvPhi[0]>0)? sqrt(phieff[4][iside][ich][ibin]*(1-phieff[4][iside][ich][ibin]) / numEvPhi[0]) : 0;
      }
      for(int lv2=0;lv2<5;lv2++){ 
        phieff_name << level[lv2] << "phiEfficiency_" << probe_trigger_chain->at(ich) << "_" << side[iside]; 
        phiEfficiency[lv2][iside][ich] = new TGraphErrors(phi_nbin,phi_center,phieff[lv2][iside][ich],phi_err,phiefferr[lv2][iside][ich]);
        phiEfficiency[lv2][iside][ich]->SetName(phieff_name.str().c_str());
        fout->Add(phiEfficiency[lv2][iside][ich]);
        phieff_name.str("");
      }
    }
  }

  //calculate eta-phi efficiency
  int numEvEtaPhi[5];
  float etaphieff[5][nChain][eta_nbin-1][phi_nbin-1],etaphiefferr[5][nChain][eta_nbin-1][phi_nbin-1];//lv/chain/etabin/phibin
  TH2 *etaphiEfficiency[5][nChain];
  stringstream etaphieff_name;
  for (int ich=0;ich<nChain;ich++){
    for (int ibin=0;ibin<eta_nbin-1;ibin++){
      for (int jbin=0;jbin<phi_nbin-1;jbin++){
        int etaphibin = probe_eta_phi[ich]->GetBin(ibin+1,jbin+1);
        numEvEtaPhi[0] = probe_eta_phi[ich]->GetBinContent(etaphibin);
        numEvEtaPhi[1] = pass_eta_phi[0][ich]->GetBinContent(etaphibin);
        numEvEtaPhi[2] = pass_eta_phi[1][ich]->GetBinContent(etaphibin);
        numEvEtaPhi[3] = pass_eta_phi[2][ich]->GetBinContent(etaphibin);
        numEvEtaPhi[4] = pass_eta_phi[3][ich]->GetBinContent(etaphibin);
        for (int lv=0;lv<4;lv++){
          etaphieff[lv][ich][ibin][jbin] = (numEvEtaPhi[lv]>0)? 1.*numEvEtaPhi[lv+1]/numEvEtaPhi[lv] : 0;
          etaphiefferr[lv][ich][ibin][jbin] = (numEvEtaPhi[lv]>0)? sqrt(etaphieff[lv][ich][ibin][jbin]*(1-etaphieff[lv][ich][ibin][jbin]) / numEvEtaPhi[lv]) : 0;
        }
        //HLT/probe
        etaphieff[4][ich][ibin][jbin] = (numEvEtaPhi[0]>0)? 1.*numEvEtaPhi[4]/numEvEtaPhi[0] : 0;
        etaphiefferr[4][ich][ibin][jbin] = (numEvEtaPhi[0]>0)? sqrt(etaphieff[4][ich][ibin][jbin]*(1-etaphieff[4][ich][ibin][jbin]) / numEvEtaPhi[0]) : 0;
        
      }
    }
    for(int lv2=0;lv2<5;lv2++){ 
      etaphieff_name << level[lv2] << "etaphiEfficiency_" << probe_trigger_chain->at(ich) ; 
      etaphiEfficiency[lv2][ich] = new TH2F(etaphieff_name.str().c_str(),";#eta;#phi",eta_nbin-1,eta_bin,phi_nbin-1,phi_bin);
      etaphieff_name.str("");
      for (int etabin=0;etabin<eta_nbin-1;etabin++){
        for (int phibin=0;phibin<phi_nbin-1;phibin++){
          int etaphibin = etaphiEfficiency[lv2][ich]->GetBin(etabin+1,phibin+1);
          etaphiEfficiency[lv2][ich]->SetBinContent(etaphibin,etaphieff[lv2][ich][etabin][phibin]);
        }
      }
    }
  }
 */ 
  fout->Write();
  
}
//////////////////////////////////////////////////////////

void ScaleFactor::makeResolutionPlot(vector<string> list,string output){

  TChain *fChain = new TChain("validationT");
  for (int in=0;in<list.size();in++) fChain->Add(list[in].c_str());
  TTree *tree = static_cast<TTree*>(fChain);
  Init(tree);

  TFile *fout = new TFile(output.c_str(),"recreate");
  TH1 *PtResidual[5][25][8];//saddress/eta/pt
  string address[5]={"Endcap","Small","SmallSpecial","Large","LargeSpecial"};
  stringstream name;
  for(int i=0;i<5;i++){
    for(int j=0;j<25;j++){
      for(int k=0;k<8;k++){
        name << "PtResidual" << address[i] << "Eta" << j << "Pt" << k ;
        PtResidual[i][j][k] = new TH1F(name.str().c_str(),";p_{T} residual;Number of events",200,-1,1);
        name.str("");
      }
    }
  }

  int nentries = fChain->GetEntries();
  int binwake[9] = {0,20,30,40,50,60,70,80,100};
  float binwakeeta[26];
  for (int e=0;e<26;e++) binwakeeta[e] = 0.1*e;
  cout << "Number of events is " << nentries << endl;

  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000000==0) cout << "entry=" << jentry << endl;
    //if(jentry==1000000) break;
    fChain->GetEntry(jentry); 

    if (*tag_trigger_chain!="HLT_mu24_ivarmedium") continue;
    float comb_pt = probe_comb_pt/1000.;
    float sa_re = (probe_sa_pt>0)? 1 - probe_offline_pt/probe_sa_pt : 9999;
    float comb_re = (comb_pt>0)? 1 - probe_offline_pt/comb_pt : 9999;
    float ef_re = (probe_ef_pt>0)? 1 - probe_offline_pt/probe_ef_pt : 9999;
    int addnum = probe_sa_saddress+1;
    for (int j=0;j<25;j++){//eta
      for (int k=0;k<8;k++){//pt
        if (fabs(probe_offline_eta)<binwakeeta[j] || fabs(probe_offline_eta)>binwakeeta[j+1]) continue;
        if(probe_offline_pt<binwake[k] || probe_offline_pt>binwake[k+1]) continue;
        if(sa_re<9999) PtResidual[addnum][j][k]->Fill(sa_re);
      }
    }
  }
  
  fout->Write();
}

void ScaleFactor::fitResolution(string input, string output){  

  TFile *file = TFile::Open(input.c_str());
  TFile *fout = new TFile(output.c_str(),"recreate");

  TH1 *PtResidual[3][3][8];
  stringstream name;
  string level[3] = {"SA","Comb","EF"};
  string side[3] = {"Barrel","Endcap","All"};
  for(int i=0;i<3;i++){ 
    for(int j=0;j<3;j++){ 
      for(int k=0;k<8;k++){ 
        name << level[i] << "PtResidual" << side[j] << "_" << k; 
        PtResidual[i][j][k] = dynamic_cast<TH1*>(file->Get(name.str().c_str()));
        name.str("");
      }
    }
  }
  
  TF1 *func1 = new TF1("func1","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2)");
  TF1 *func2 = new TF1("func2","[0]*exp(-0.5*((x-[1])/[2])^2)");
  TF1 *func3 = new TF1("func3","[0]*exp(-0.5*((x-[1])/[2])^2)");
  TF1 *func4 = new TF1("func4","[2]*[0]/((x-[1])*(x-[1])+[0]/4)");
  float mean[3][3][8],mean_err[3][3][8],sigma[3][3][8],sigma_err[3][3][8];
  float ptbin[8] = {10,25,35,45,55,65,75,90};
  float ptbinerr[8] = {10,5,5,5,5,5,5,10};
  stringstream canvas_name,plot_name;
  TGraphErrors *graph[3][3];
  
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<8;k++){
        float tmp_mean = PtResidual[i][j][k]->GetMean();
        float tmp_sigma = PtResidual[i][j][k]->GetRMS();
        /*func4->SetParameter(0,tmp_sigma*tmp_sigma);
        func4->SetParameter(1,tmp_mean);
        func4->SetParameter(2,PtResidual[i][j][k]->GetMaximum());
        PtResidual[i][j][k]->Fit(func4,"","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        */
        func1->SetLineColor(2);
        func1->SetParameter(1,tmp_mean);
        func1->SetParameter(2,tmp_sigma);
        func1->SetParLimits(0,0,PtResidual[i][j][k]->GetMaximum());
        //func1->SetParLimits(2,0,1000);
        func1->SetParLimits(3,0,PtResidual[i][j][k]->GetMaximum());
        func1->SetParameter(4,tmp_mean);
        func1->SetParLimits(5,0,tmp_sigma);
        PtResidual[i][j][k]->Fit(func1,"","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        /*tmp_mean = func1->GetParameter(1);
        tmp_sigma = func1->GetParameter(2);
        func1->SetParameter(1,tmp_mean);
        func1->SetParameter(2,tmp_sigma);
        //func1->SetParLimits(2,0,1000);
        func1->SetParameter(4,func1->GetParameter(4));
        func1->SetParLimits(5,0,0.5*tmp_sigma);
        PtResidual[i][j][k]->Fit(func1,"","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        */
        /*tmp_mean = func1->GetParameter(1);
        tmp_sigma = func1->GetParameter(2);
        func1->SetParameter(1,tmp_mean);
        func1->SetParameter(2,tmp_sigma);
        //func1->SetParLimits(2,0,1000);
        func1->SetParameter(4,func1->GetParameter(4));
        func1->SetParLimits(5,0,0.5*tmp_sigma);
        PtResidual[i][j][k]->Fit(func1,"+","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        */
        sigma[i][j][k] = 0.5*(func1->GetParameter(2) + func1->GetParameter(5));
        //sigma[i][j][k] = func1->GetParameter(5);
        sigma_err[i][j][k] = sqrt( pow(func1->GetParError(2),2) + pow(func1->GetParError(5),2) );
        //sigma_err[i][j][k] = func1->GetParError(5);
        func2->FixParameter(0,func1->GetParameter(0));
        func2->FixParameter(1,func1->GetParameter(1));
        func2->FixParameter(2,func1->GetParameter(2));
        func3->FixParameter(0,func1->GetParameter(3));
        func3->FixParameter(1,func1->GetParameter(4));
        func3->FixParameter(2,func1->GetParameter(5));
        func2->SetLineColor(4);
        func3->SetLineColor(6);
        PtResidual[i][j][k]->Fit(func2,"+","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        PtResidual[i][j][k]->Fit(func3,"+","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        fout->Add(PtResidual[i][j][k]);
      }
      plot_name << level[i] << "PtResolution" << side[j] ;
      graph[i][j] = new TGraphErrors(8,ptbin,sigma[i][j],ptbinerr,sigma_err[i][j]);
      graph[i][j]->SetName(plot_name.str().c_str());
      fout->Add(graph[i][j]);
      plot_name.str("");
    }
  }
  TCanvas *c2 = new TCanvas("c2","",800,600);
  graph[0][2]->SetLineColor(1);//SA
  graph[0][2]->SetTitle("");
  graph[0][2]->GetYaxis()->SetRangeUser(0,0.1);//SA
  graph[0][2]->Draw("ap");
  graph[1][2]->SetLineColor(2);//Comb
  graph[1][2]->Draw("p");
  graph[2][2]->SetLineColor(4);//EF
  graph[2][2]->Draw("p");
  c2->SaveAs("plot2.pdf");
  
  fout->Write();

  
}
