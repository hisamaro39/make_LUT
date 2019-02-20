#define newValidationT_cxx
#include "macro/newValidationT.h"
#include "macro/readLUT.h"
#include <TH2.h>
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

double const ZERO_LIMIT = 1e-5;

//pt residual
const int divpt[] = {0,5,10,15,20,25,30,40,50,60,100};
const int ndivpt = sizeof(divpt)/sizeof(int);
const int nChamber=5;
const float diveta[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
const int ndiveta = sizeof(diveta)/sizeof(float);

///////////////////////////////////////////////////////////////////////////
void newValidationT::makePtResidual(string list,string output)
{
  TFile *fout = new TFile(output.c_str(),"recreate");
  
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetOptStat(0);
  
  TChain *chain = new TChain("validationT");
  ifstream finlist(list.c_str());
  string file_rec;
  //cout << "input sample" << endl;
  while(finlist>>file_rec) {
    //cout << file_rec.c_str() << endl;
    chain->Add(file_rec.c_str());
  }
  TTree *tree = static_cast<TTree*>(chain);
  Init(tree);

  TH1 *hMDT_PtResidual[nChamber][ndiveta-1][ndivpt-1];//saddress/eta/pt
  stringstream PtResidualdrawname;

  for(int sad=0;sad<nChamber;sad++){
    for(int ieta=0;ieta<ndiveta-1;ieta++){
      for(int ipt=0;ipt<ndivpt-1;ipt++){
        PtResidualdrawname << "PtResidualSaddress" << sad << "Eta" << ieta << "Pt" << ipt;
        hMDT_PtResidual[sad][ieta][ipt] = new TH1F(PtResidualdrawname.str().c_str(),";p_{T} residual;Number of events",500,-1,1);
        PtResidualdrawname.str("");
      }
    }
  }

  int nentries = chain->GetEntries();
  cout << "Number of events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000000==0) cout << "entry=" << jentry << endl;
    chain->GetEntry(jentry); 

    float sa_residual = (L2_pt>0)? 1-offline_pt/L2_pt : 10000;
    //if(L2_saddress==1 && fabs(offline_eta)>0.2 && fabs(offline_eta)<0.3){

    for (int ieta=0;ieta<ndiveta-1;ieta++){
      for (int ipt=0;ipt<ndivpt-1;ipt++){
        if(fabs(offline_eta)>diveta[ieta] &&  fabs(offline_eta)<diveta[ieta+1] && 
            offline_pt>divpt[ipt] && offline_pt<divpt[ipt+1]) {
          if(L2_saddress==-1) hMDT_PtResidual[0][ieta][ipt]->Fill(sa_residual);
          if(L2_saddress==0)  hMDT_PtResidual[1][ieta][ipt]->Fill(sa_residual);
          if(L2_saddress==1)  hMDT_PtResidual[2][ieta][ipt]->Fill(sa_residual);
          if(L2_saddress==2)  hMDT_PtResidual[3][ieta][ipt]->Fill(sa_residual);
          if(L2_saddress==3)  hMDT_PtResidual[4][ieta][ipt]->Fill(sa_residual);
        }
      }
    }
  }
  fout->Write();
}

/////////////////////////////////////////////////
void newValidationT::Loop(string list,string output)
{
  TChain *chain = new TChain("validationT");
  ifstream finlist(list.c_str());
  string file_rec;
  while(finlist>>file_rec) chain->Add(file_rec.c_str());
  TTree *tree = static_cast<TTree*>(chain);
  Init(tree);

  TFile *fout = new TFile(output.c_str(),"recreate");
  
  bool makeLUTEndcap=false;
  bool makeLUTBarrel=false;
  bool makeLUTLargeSpecial=true;
  bool makeLUTLargeSpecialLight=false;
  bool makeLUTLargeSpecialNoSpr=false;
  bool makeLUTLargeSpecialInverse=false;
  bool correctEE=false;
  bool barrelSagitta=false;
  bool LSSagitta=false;
  bool barrelSagittaComb=false;
  
  ///Initialize
  m_readLUT.scanLUT();

  ///////////////set region ///////////////////
  int phibinmax_ec = 12;
  int phibinmin_ec = 0;
  int phinumber = phibinmax_ec - phibinmin_ec + 1;
  int divetamax_ec = 30;
  int divetamin_ec = 0;
  int etanumber = divetamax_ec - divetamin_ec + 1;
  int phibinmax_br = 30;
  int divetamax_br = 30;

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetTitleOffset(1.5,"y");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetOptStat(0);

  ////////////////////////////////////////////
  //endcap
  TH2F *hMDT_invpTvsAlpha[phinumber][etanumber][2];
  TH2F *hMDT_invpTvsTgcAlpha[phinumber][etanumber][2];
  TH2F *hMDT_invpTvsBeta[phinumber][etanumber][2];
  TH2F *hMDT_invpTvsEndcapRadius[24][8][2][2];//EE
  TH2F *hMDT_invpTvsEndcapRadiusNoSL[12][8][2];//EE no SL
  TH2F *hMDT_invpTvsEndcapRadiusAllPhiNoQeta[192][8][2][2];//EE all phi no qeta
  TH2F *hMDT_invpTvsEndcapRadiusLinear[192][8][2][2];//EE R vs pt
  TH2F *hMDT_invpTvsEndcapRadiusAnormal[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_invpTvsEndcapRadiusShift[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_pTvsEndcapSagittaAnormal[192][8][2][2];//Check EE anormal region
  TH2F *hMDT_pTvsBarrelSagitta[2][4][divetamax_br][phibinmax_br];//Check Barrel sagitta
  TH2F *hMDT_pTvsLSSagitta[2][2][divetamax_br][phibinmax_br];//Check LS sagitta
  TH2F *hMDT_pTvsBarrelSagittaComb[4][divetamax_br][phibinmax_br];//Check Barrel sagitta
  stringstream Alphadrawname,TgcAlphadrawname,Betadrawname;
  stringstream EndcapRadiusdrawname,EndcapRadiusNoSLdrawname,EndcapRadiusAllPhiNoQetadrawname;
  stringstream EndcapRadiusLineardrawname,EndcapRadiusAnormaldrawname,EndcapRadiusShiftdrawname;
  stringstream EndcapSagittaAnormaldrawname,BarrelSagittadrawname,LSSagittadrawname;
  stringstream iphi_ec,ieta_ec,iphi_ee,ieta_ee;
  //barrel
  TH2F *hMDT_invpTvsinvRadius[2][4][divetamax_br][phibinmax_br];
  TH2F *hMDT_invpTvsinvRadiusLS[2][2][divetamax_br][phibinmax_br];//charge/spr/etabin/phibin
  TH2F *hMDT_invpTvsinvRadiusLSInverse[2][2][divetamax_br][phibinmax_br];//charge/spr/etabin/phibin
  TH2F *hMDT_invpTvsinvRadiusLSNoSpr[2][divetamax_br][phibinmax_br];//charge/spr/etabin/phibin
  TH2F *hMDT_invpTvsinvRadiusLSLight[2][2][divetamax_br];//charge/spr/etabin/

  stringstream BarrelRadiusdrawname;
  stringstream iphi_br,ieta_br;


  if (makeLUTEndcap){ 
    //endcap
    for(int i = 0;i<phibinmax_ec;i++){
      for(int j = 0;j<divetamax_ec;j++){
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
  
  if(makeLUTBarrel){
    //barrel
    for(int k=0; k<2; k++){//charge
      for (int chamber=0; chamber<4; chamber++){//saddress
        for(int i=0; i<divetamax_br; i++){
          for(int j=0; j<phibinmax_br; j++){
            if(i < 10) ieta_br << "0" << i;
            else ieta_br << i;
            if(j < 10) iphi_br << "0" << j;
            else iphi_br << j;
            BarrelRadiusdrawname << "Radius" << k << chamber << ieta_br.str().c_str() << iphi_br.str().c_str();
            hMDT_invpTvsinvRadius[k][chamber][i][j] = new TH2F(BarrelRadiusdrawname.str().c_str(),"",500,0,500000,120,0,60);
            hMDT_invpTvsinvRadius[k][chamber][i][j]->SetMarkerStyle(8);
            hMDT_invpTvsinvRadius[k][chamber][i][j]->GetXaxis()->SetRangeUser(0,500000);
            hMDT_invpTvsinvRadius[k][chamber][i][j]->GetYaxis()->SetRangeUser(0,60);
            BarrelRadiusdrawname.str("");
            iphi_br.str("");
            ieta_br.str("");
          }
        }
      }
    }
  }

  //Large Special
  if(makeLUTLargeSpecial){
    for(int k=0; k<2; k++){//qeta
      for (int phi=0; phi<2; phi++){//phi=-1.5
        for(int i=0; i<divetamax_br; i++){//default 30
          for(int j=0; j<30; j++){//defualt 30
            if(i < 10) ieta_br << "0" << i;
            else ieta_br << i;
            if(j < 10) iphi_br << "0" << j;
            else iphi_br << j;
            BarrelRadiusdrawname << "RadiusLS" << k << phi << ieta_br.str().c_str() << iphi_br.str().c_str();
            hMDT_invpTvsinvRadiusLS[k][phi][i][j] = new TH2F(BarrelRadiusdrawname.str().c_str(),"",1000,0,1000000,200,0,100);
            hMDT_invpTvsinvRadiusLS[k][phi][i][j]->SetMarkerStyle(8);
            hMDT_invpTvsinvRadiusLS[k][phi][i][j]->GetXaxis()->SetRangeUser(0,1000000);
            hMDT_invpTvsinvRadiusLS[k][phi][i][j]->GetYaxis()->SetRangeUser(0,100);
            BarrelRadiusdrawname.str("");
            iphi_br.str("");
            ieta_br.str("");
          }
        }
      }
    }
  }
  
  //Large Special Inverse
  if(makeLUTLargeSpecialInverse){
    for(int k=0; k<2; k++){//qeta
      for (int spr=0; spr<2; spr++){//spr=-1.5
        for(int i=0; i<divetamax_br; i++){//default 30
          for(int j=0; j<30; j++){//defualt 30
            if(i < 10) ieta_br << "0" << i;
            else ieta_br << i;
            if(j < 10) iphi_br << "0" << j;
            else iphi_br << j;
            BarrelRadiusdrawname << "RadiusLS" << k << spr << ieta_br.str().c_str() << iphi_br.str().c_str();
            hMDT_invpTvsinvRadiusLSInverse[k][spr][i][j] = new TH2F(BarrelRadiusdrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
            hMDT_invpTvsinvRadiusLSInverse[k][spr][i][j]->SetMarkerStyle(8);
            hMDT_invpTvsinvRadiusLSInverse[k][spr][i][j]->GetXaxis()->SetRangeUser(0,0.25);
            hMDT_invpTvsinvRadiusLSInverse[k][spr][i][j]->GetYaxis()->SetRangeUser(0,0.00004);
            BarrelRadiusdrawname.str("");
            iphi_br.str("");
            ieta_br.str("");
          }
        }
      }
    }
  }

  //Large Special no spr
  if(makeLUTLargeSpecialNoSpr){
    for(int k=0; k<2; k++){//qeta
      for(int i=0; i<divetamax_br; i++){//default 30
        for(int j=0; j<30; j++){//defualt 30
          if(i < 10) ieta_br << "0" << i;
          else ieta_br << i;
          if(j < 10) iphi_br << "0" << j;
          else iphi_br << j;
          BarrelRadiusdrawname << "RadiusLSNoSpr" << k << ieta_br.str().c_str() << iphi_br.str().c_str();
          hMDT_invpTvsinvRadiusLSNoSpr[k][i][j] = new TH2F(BarrelRadiusdrawname.str().c_str(),"",1000,0,1000000,200,0,100);
          hMDT_invpTvsinvRadiusLSNoSpr[k][i][j]->SetMarkerStyle(8);
          hMDT_invpTvsinvRadiusLSNoSpr[k][i][j]->GetXaxis()->SetRangeUser(0,1000000);
          hMDT_invpTvsinvRadiusLSNoSpr[k][i][j]->GetYaxis()->SetRangeUser(0,100);
          BarrelRadiusdrawname.str("");
          iphi_br.str("");
          ieta_br.str("");
        }
      }
    }
  }

  //Large SpecialLight
  if(makeLUTLargeSpecialLight){
    for(int k=0; k<2; k++){//qeta
      for (int spr=0; spr<2; spr++){//SPR
        for(int i=0; i<divetamax_br; i++){//default 30
          if(i < 10) ieta_br << "0" << i;
          else ieta_br << i;
          BarrelRadiusdrawname << "RadiusLSLight" << k << spr << ieta_br.str().c_str() ;
          hMDT_invpTvsinvRadiusLSLight[k][spr][i] = new TH2F(BarrelRadiusdrawname.str().c_str(),"",1000,0,1000000,200,0,100);
          hMDT_invpTvsinvRadiusLSLight[k][spr][i]->SetMarkerStyle(8);
          hMDT_invpTvsinvRadiusLSLight[k][spr][i]->GetXaxis()->SetRangeUser(0,1000000);
          hMDT_invpTvsinvRadiusLSLight[k][spr][i]->GetYaxis()->SetRangeUser(0,100);
          BarrelRadiusdrawname.str("");
          iphi_br.str("");
          ieta_br.str("");
        }
      }
    }
  }

  if (correctEE){//Search anormal EE region
    for(int i = 0;i<192;i++){
      for(int j = 0;j<8;j++){
        if(i < 10) iphi_ec << "0" << i;
        else iphi_ec << i;
        if(j < 9) ieta_ec << "0" << j+1;
        else ieta_ec << j+1;
        for (int l=0;l<2;l++){//eta
          for(int k = 0;k<2;k++){//charge
            EndcapSagittaAnormaldrawname << "EndcapSagittaAnormal" << iphi_ec.str().c_str() << ieta_ec.str().c_str() <<  k << l;
            hMDT_pTvsEndcapSagittaAnormal[i][j][k][l] = new TH2F(EndcapSagittaAnormaldrawname.str().c_str(),"",400,0,80,600,-800,800);
            hMDT_pTvsEndcapSagittaAnormal[i][j][k][l]->SetMarkerStyle(8);
            hMDT_pTvsEndcapSagittaAnormal[i][j][k][l]->GetXaxis()->SetRangeUser(0,80);
            hMDT_pTvsEndcapSagittaAnormal[i][j][k][l]->GetYaxis()->SetRangeUser(-800,800);
            EndcapSagittaAnormaldrawname.str("");
            EndcapRadiusShiftdrawname << "EndcapRadiusShift" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            hMDT_invpTvsEndcapRadiusShift[i][j][k][l] = new TH2F(EndcapRadiusShiftdrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
            hMDT_invpTvsEndcapRadiusShift[i][j][k][l]->SetMarkerStyle(8);
            hMDT_invpTvsEndcapRadiusShift[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
            hMDT_invpTvsEndcapRadiusShift[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.00004);
            EndcapRadiusShiftdrawname.str("");
            EndcapRadiusAnormaldrawname << "EndcapRadiusAnormal" << iphi_ec.str().c_str() << ieta_ec.str().c_str() << k << l;
            hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l] = new TH2F(EndcapRadiusAnormaldrawname.str().c_str(),"",250,0,0.25,400,0,0.00004);
            hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l]->SetMarkerStyle(8);
            hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l]->GetXaxis()->SetRangeUser(0,0.25);
            hMDT_invpTvsEndcapRadiusAnormal[i][j][k][l]->GetYaxis()->SetRangeUser(0,0.00004);
            EndcapRadiusAnormaldrawname.str("");
          }
        }
        iphi_ec.str("");
        ieta_ec.str("");
      }
    }
  }

  if(barrelSagitta){//search barrel sagitta
    for (int chamber=0; chamber<4; chamber++){//saddress
      for(int i=0; i<divetamax_br; i++){
        for(int j=0; j<phibinmax_br; j++){
          if(i < 10) ieta_br << "0" << i;
          else ieta_br << i;
          if(j < 10) iphi_br << "0" << j;
          else iphi_br << j;
          for(int k=0; k<2; k++){//charge
            BarrelSagittadrawname << "BarrelSagitta" << k << chamber << ieta_br.str().c_str() << iphi_br.str().c_str();
            hMDT_pTvsBarrelSagitta[k][chamber][i][j] = new TH2F(BarrelSagittadrawname.str().c_str(),"",400,0,80,600,-800,800);
            hMDT_pTvsBarrelSagitta[k][chamber][i][j]->SetMarkerStyle(8);
            hMDT_pTvsBarrelSagitta[k][chamber][i][j]->GetXaxis()->SetRangeUser(0,80);
            hMDT_pTvsBarrelSagitta[k][chamber][i][j]->GetYaxis()->SetRangeUser(-800,800);
            BarrelSagittadrawname.str("");
          }
          iphi_br.str("");
          ieta_br.str("");
        }
      }
    }
  }

  if(LSSagitta){//search LS sagitta
    for (int spr=0; spr<2; spr++){//spr
      for(int i=0; i<divetamax_br; i++){
        for(int j=0; j<phibinmax_br; j++){
          if(i < 10) ieta_br << "0" << i;
          else ieta_br << i;
          if(j < 10) iphi_br << "0" << j;
          else iphi_br << j;
          for(int k=0; k<2; k++){//charge
            LSSagittadrawname << "LSSagitta" << k << spr << ieta_br.str().c_str() << iphi_br.str().c_str();
            hMDT_pTvsLSSagitta[k][spr][i][j] = new TH2F(LSSagittadrawname.str().c_str(),"",400,0,80,600,-800,800);
            hMDT_pTvsLSSagitta[k][spr][i][j]->SetMarkerStyle(8);
            hMDT_pTvsLSSagitta[k][spr][i][j]->GetXaxis()->SetRangeUser(0,80);
            hMDT_pTvsLSSagitta[k][spr][i][j]->GetYaxis()->SetRangeUser(-800,800);
            LSSagittadrawname.str("");
          }
          iphi_br.str("");
          ieta_br.str("");
        }
      }
    }
  }

  if(barrelSagittaComb){//search barrel sagitta
    for (int chamber=0; chamber<4; chamber++){//saddress
      for(int i=0; i<divetamax_br; i++){
        for(int j=0; j<phibinmax_br; j++){
          if(i < 10) ieta_br << "0" << i;
          else ieta_br << i;
          if(j < 10) iphi_br << "0" << j;
          else iphi_br << j;
          BarrelSagittadrawname << "BarrelSagittaComb" <<chamber << ieta_br.str().c_str() << iphi_br.str().c_str();
          hMDT_pTvsBarrelSagittaComb[chamber][i][j] = new TH2F(BarrelSagittadrawname.str().c_str(),"",400,0,80,600,-800,800);
          hMDT_pTvsBarrelSagittaComb[chamber][i][j]->SetMarkerStyle(8);
          hMDT_pTvsBarrelSagittaComb[chamber][i][j]->GetXaxis()->SetRangeUser(0,80);
          hMDT_pTvsBarrelSagittaComb[chamber][i][j]->GetYaxis()->SetRangeUser(-800,800);
          BarrelSagittadrawname.str("");
          iphi_br.str("");
          ieta_br.str("");
        }
      }
    }
  }
  TH2 *test1 = new TH2F("test1","",500,0,50000,120,0,60);
  TH2 *test2 = new TH2F("test2","",500,0,50000,120,0,60);
  TH2 *test3 = new TH2F("test3","",500,0,50000,120,0,60);
  TH2 *test4 = new TH2F("test4","",500,0,50000,120,0,60);
  TH2 *CorrectSagittaParPlus = new TH2I("CorrectSagittaParPlus","",30,0,30,30,0,30);
  TH2 *CorrectSagittaParMinus = new TH2I("CorrectSagittaParMinus","",30,0,30,30,0,30);
  TH1 *phiMapNormal = new TH1F("phiMapNormal",";phiMap;",300,-0.3,0.3);
  TH1 *phiMapAnormal = new TH1F("phiMapAnormal",";phiMap;",300,-0.3,0.3);
  TH2 *SPBIZvsR = new TH2F("SPBIZvsR",";Z(mm);R(mm)",400,0,-8000,400,4000,8000);
  TH2 *SPBIZvsRNormal = new TH2F("SPBIZvsRNormal",";Z(mm);R(mm)",400,0,-8000,400,4000,8000);
  TH2 *SPBIZvsRAnormal = new TH2F("SPBIZvsRAnormal",";Z(mm);R(mm)",400,0,-8000,400,4000,8000);
  TH1 *SPEIZ = new TH1F("SPEIZ",";Z(mm);",300,6000,9000);
  TH1 *SPEMZ = new TH1F("SPEMZ",";Z(mm);",300,12000,15000);
  TH1 *SPEOZ = new TH1F("SPEOZ",";Z(mm);",300,20000,23000);
  TH2 *SPEIZR = new TH2F("SPEIZR",";Z(mm);R(mm)",300,6000,9000,600,2000,8000);
  TH2 *RoIPhiVsEISPZ = new TH2F("RoIPhiVsEISPZ",";RoI #phi;Z(mm)",640,-3.2,3.2,300,6000,9000);
  TH2 *RoIPhiVsEMSPZ = new TH2F("RoIPhiVsEMSPZ",";RoI #phi;Z(mm)",640,-3.2,3.2,300,12000,15000);
  TH2 *RoIPhiVsEOSPZ = new TH2F("RoIPhiVsEOSPZ",";RoI #phi;Z(mm)",640,-3.2,3.2,300,20000,23000);
  TH2 *RoIPhiVsBISPR = new TH2F("RoIPhiVsBISPR",";RoI #phi;R(mm)",640,-3.2,3.2,300,4000,6000);
  TH2 *RoIPhiVsBMSPR = new TH2F("RoIPhiVsBMSPR",";RoI #phi;R(mm)",640,-3.2,3.2,200,6000,9000);
  TH2 *RoIPhiVsBOSPR = new TH2F("RoIPhiVsBOSPR",";RoI #phi;R(mm)",640,-3.2,3.2,300,9000,12000);
  TH1 *LSPhiMap = new TH1F("LSPhi",";phiMap;",200,-1,1);
  TH2 *invpTvsinvRadiusLS010309 = new TH2F("invpTvsinvRadiusLS010309","",250,0,0.25,400,0,0.00004);
  TH2 *invpTvsinvRadiusLSSPRlow010309 = new TH2F("invpTvsinvRadiusLSSPRlow010309","",250,0,0.25,400,0,0.00004);
  TH2 *invpTvsinvRadiusLSSPRhigh010309 = new TH2F("invpTvsinvRadiusLSSPRhigh010309","",250,0,0.25,400,0,0.00004);
  TH1 *etaMap010309 = new TH1F("etaMap010309","",200,-1,1);
  TH1 *phiMap010309 = new TH1F("phiMap010309","",400,-0.2,0.2);
  TH1 *SPBIR010309 = new TH1F("SPBIR010309","",500,4000,9000);
  TH2 *SPBIRvsPhi010309 = new TH2F("SPBIRvsPhi010309","",64,-3.2,3.2,500,4000,9000);
  TH1 *SPBILSR = new TH1F("SPBILSR","",500,4000,9000);

  ////////////////////////////////////////////////
  float EtaMin[4] = {-1.145, -1.150, -1.050, -1.050};
  float PhiMin[4] = {-0.230, -0.233, -0.181, -0.181};
  float EtaMax[4] = {1.145, 1.150, 1.050, 1.050};
  float PhiMax[4] = {0.230, 0.233, 0.181, 0.181};
  float EtaStep[4], PhiStep[4];
  for (int i=0; i<4; i++){
    EtaStep[i] = (EtaMax[i]-EtaMin[i])/30;
    PhiStep[i] = (PhiMax[i]-PhiMin[i])/30;
    //PhiStep[i] = (PhiMax[i]-PhiMin[i])/60;
  }

  /////////////////////////////////////
  int nentries = chain->GetEntries();
  cout << "Number of events is " << nentries << endl;
  for (int jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000000==0) cout << "entry=" << jentry << endl;
    //if(jentry==100000000) break;
    chain->GetEntry(jentry); 
    
    ////////////////////////////////////////

    float tgcalpha = (fabs(tgcMid1_z)>1e-5 && fabs(tgcMid2_z)>1e-5) ? m_util.calcAlpha(tgcMid1_r,tgcMid1_z,tgcMid2_r,tgcMid2_z) : 0;
    bool isEndcap = (L2_saddress<-0.5)? true : false;
    float center_phi = m_util.getChamberCenterPhi(tgcMid1_phi);
    float cosphidif = m_util.cosAminusB(tgcMid1_phi,center_phi);

    //SP
    float SP_br_inner_z = (L2_saddress!=-1)? sp_z->at(0) : 0;
    float SP_br_inner_r = (L2_saddress!=-1)? sp_r->at(0) : 0;
    float SP_br_middle_z = sp_z->at(1);
    float SP_br_middle_r = sp_r->at(1);
    float SP_br_outer_z = sp_z->at(2);
    float SP_br_outer_r = sp_r->at(2);
    float SP_ec_inner_z = sp_z->at(3);
    float SP_ec_inner_r = sp_r->at(3);
    float SP_ec_middle_z = sp_z->at(4);
    float SP_ec_middle_r = sp_r->at(4);
    float SP_ec_outer_z = sp_z->at(5);
    float SP_ec_outer_r = sp_r->at(5);
    float SP_ec_ee_z = sp_z->at(6);
    float SP_ec_ee_r = sp_r->at(6);
    float SP_ec_ee_r_shift = sp_r->at(6);
    SPEIZ->Fill(fabs(SP_ec_inner_z));
    SPEMZ->Fill(fabs(SP_ec_middle_z));
    SPEOZ->Fill(fabs(SP_ec_outer_z));
    SPEIZR->Fill(fabs(SP_ec_inner_z),SP_ec_inner_r);
    RoIPhiVsEISPZ->Fill(L1_phi,SP_ec_inner_z);
    RoIPhiVsEMSPZ->Fill(L1_phi,SP_ec_middle_z);
    RoIPhiVsEOSPZ->Fill(L1_phi,SP_ec_outer_z);
    RoIPhiVsBISPR->Fill(L1_phi,SP_br_inner_r);
    RoIPhiVsBMSPR->Fill(L1_phi,SP_br_middle_r);
    RoIPhiVsBOSPR->Fill(L1_phi,SP_br_outer_r);
    //-1.1 to -0.95
    //if(L1_phi>-1.1 && L1_phi<-0.95) SPEIZR->Fill(SP_ec_inner_z,SP_ec_inner_r);

    ////////////////////////////////////////////////////////////////////////
    if (isEndcap){ //endcap
      bool isTgcFailure = (fabs(tgcMid1_phi)<1e-5)? true : false;
      float usedEta = (isTgcFailure)? L1_eta : tgcMid1_eta;
      float usedPhi = (isTgcFailure)? L1_phi : tgcMid1_phi;
      pair<int, int> LUTbinNumbers = m_util.GetBinNumber( usedPhi, usedEta );
      pair<int, int> LUTbinNumbersTgc = m_util.GetBinNumber( tgcMid1_phi, tgcMid1_eta );
      int PhibinNumber = LUTbinNumbers.second;
      int divetaNumber = LUTbinNumbers.first;
      int PhibinNumberTgc = LUTbinNumbersTgc.second;
      int divetaNumberTgc = LUTbinNumbersTgc.first;
      int i = PhibinNumber;
      int j = divetaNumber;
      int PhibinAll  = m_util.GetBinNumberAllPhi(usedPhi);
      bool badPhibin = m_util.isBadPhi(PhibinAll);
      float tgc_intercept = (fabs(tgcMid1_z)>1e-5 && fabs(tgcMid2_z)>1e-5)? 
        m_util.calcIntercept(tgcMid1_r,tgcMid1_z,tgcMid2_r,tgcMid2_z) : 0.;
      int tgc_charge = (tgc_intercept*tgcMid2_z<0)? -1 : 1;
      int tgc_qeta = (tgc_charge*tgcMid1_eta<0)? 0 : 1;

      float calc_ec_radius=0.;
      float calc_ec_radius_shift=0.;
      float calc_ec_sagitta=0.;
      pair<float,float> center;
      bool small = m_util.isSmall(SP_ec_ee_z);
      int L2_SL = (small)? 0 : 1;
      int ieta = (L2_etaMap>0) ? 1 : 0;
      int icharge = (L2_charge>0) ? 1 : 0;

      //calculate radius from SP
      if(small){
        if (SP_br_inner_r>1e-5 && SP_ec_ee_r>1e-5 && SP_ec_middle_r>1e-5){
          //calc_ec_radius = m_util.computeRadius3Points(SP_br_inner_z,SP_br_inner_r,
            //  SP_ec_ee_z,SP_ec_ee_r,
              //SP_ec_middle_z,SP_ec_middle_r);
          calc_ec_radius = m_util.computeRadius3Points(SP_br_inner_z,SP_br_inner_r*cosphidif,
              SP_ec_ee_z,SP_ec_ee_r*cosphidif,
              SP_ec_middle_z,SP_ec_middle_r*cosphidif);
          center = m_util.calcCenter(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          calc_ec_sagitta = m_util.calcSagitta(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          SP_ec_ee_r_shift = m_readLUT.shiftEER(SP_br_inner_z, SP_br_inner_r,
              SP_ec_ee_z, SP_ec_ee_r,
              SP_ec_middle_z, SP_ec_middle_r,
              j, PhibinAll, icharge, ieta);
          calc_ec_radius_shift = m_util.computeRadius3Points(SP_br_inner_z,SP_br_inner_r,
              SP_ec_ee_z,SP_ec_ee_r_shift,
              SP_ec_middle_z,SP_ec_middle_r);
        }
      }
      else{
        if (SP_ec_inner_r>1e-5 && SP_ec_ee_r>1e-5 && SP_ec_middle_r>1e-5){
          calc_ec_radius = m_util.computeRadius3Points(SP_ec_inner_z,SP_ec_inner_r*cosphidif,
              SP_ec_ee_z,SP_ec_ee_r*cosphidif,
              SP_ec_middle_z,SP_ec_middle_r*cosphidif);
          center = m_util.calcCenter(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          calc_ec_sagitta = m_util.calcSagitta(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r,
              SP_ec_middle_z,SP_ec_middle_r);
          SP_ec_ee_r_shift = m_readLUT.shiftEER(SP_ec_inner_z, SP_ec_inner_r,
              SP_ec_ee_z, SP_ec_ee_r,
              SP_ec_middle_z, SP_ec_middle_r,
              j, PhibinAll, icharge, ieta);
          calc_ec_radius_shift = m_util.computeRadius3Points(SP_ec_inner_z,SP_ec_inner_r,
              SP_ec_ee_z,SP_ec_ee_r_shift,
              SP_ec_middle_z,SP_ec_middle_r);
        }
      }

      int PhibinNumberEE = m_util.GetBinNumberEE(usedPhi, L2_SL);

      if (i<0 || j<0) continue;

      int Qeta = (L2_charge*L2_etaMap>0)? 1 : 0; 
      if (L2_ec_alpha>ZERO_LIMIT){//alpha was caluclated
        if (makeLUTEndcap) hMDT_invpTvsAlpha[i][j][Qeta]->Fill(fabs(1/offline_pt),L2_ec_alpha);
      }
      if (tgcalpha>ZERO_LIMIT){
        if(makeLUTEndcap) hMDT_invpTvsTgcAlpha[PhibinNumberTgc][divetaNumberTgc][tgc_qeta]->Fill(fabs(1/offline_pt),tgcalpha);
      }
      if (L2_ec_beta>ZERO_LIMIT){//beta was caluclated
        if (makeLUTEndcap) hMDT_invpTvsBeta[i][j][Qeta]->Fill(fabs(1/offline_pt),L2_ec_beta);
      }
      if (calc_ec_radius>ZERO_LIMIT){
        if(j<8){
          if(correctEE){
            hMDT_pTvsEndcapSagittaAnormal[PhibinAll][j][icharge][ieta]->Fill(offline_pt,calc_ec_sagitta);
            hMDT_invpTvsEndcapRadiusAnormal[PhibinAll][j][icharge][ieta]->Fill(fabs(1/offline_pt),1./calc_ec_radius);
            hMDT_invpTvsEndcapRadiusShift[PhibinAll][j][icharge][ieta]->Fill(fabs(1/offline_pt),1./calc_ec_radius_shift);
          }
          if (makeLUTEndcap){
            if (!badPhibin) {
              hMDT_invpTvsEndcapRadiusNoSL[i][j][Qeta]->Fill(fabs(1/offline_pt),1./calc_ec_radius);
              hMDT_invpTvsEndcapRadius[PhibinNumberEE][j][Qeta][L2_SL]->Fill(fabs(1/offline_pt),1./calc_ec_radius);
            }
            hMDT_invpTvsEndcapRadiusAllPhiNoQeta[PhibinAll][j][icharge][ieta]->Fill(fabs(1/offline_pt),1./calc_ec_radius);
            hMDT_invpTvsEndcapRadiusLinear[PhibinAll][j][icharge][ieta]->Fill(offline_pt,calc_ec_radius);
          }
        }
      }
    }
    else{//barrel
      if(L2_saddress==1) LSPhiMap->Fill(L2_phiMap);
      int diveta_br = (int)((L2_etaMap - EtaMin[L2_saddress])/EtaStep[L2_saddress]);
      int phiBin_br = (int)((L2_phiMap - PhiMin[L2_saddress])/PhiStep[L2_saddress]);
      int icharge=-1;
      if(barrelSagitta || barrelSagittaComb || LSSagitta) icharge = (offline_charge>0)? 1 : 0;
      else icharge = (L2_charge>0)? 1 : 0;
      int iphi = (L2_phi < -1.5)? 0 : 1; 
      int ispr = (SP_br_inner_r<6000)? 0 : 1;
      //if(phiBin_br>=30) cout << "phi/bin=" << L2_phiMap << "/" << phiBin_br  << endl;
      if(diveta_br<=-1) diveta_br = 0;
      if(diveta_br>=30) diveta_br = 29;
      if(phiBin_br<=-1) phiBin_br = 0;
      if(phiBin_br>=30) {
        phiBin_br = 29;
        //continue;
      }
      int iqeta = ( (icharge==0 && diveta_br<=14) || (icharge==1 && diveta_br>14) )? 1 : 0; 
      float calc_shift_br_radius=L2_br_radius;
      if (SP_br_inner_r>1e-5 && SP_br_middle_r>1e-5 && SP_br_outer_r>1e-5){
        float calc_br_sagitta = m_util.calcSagitta(SP_br_inner_z,SP_br_inner_r,
            SP_br_middle_z,SP_br_middle_r,
            SP_br_outer_z,SP_br_outer_r);
        float dZ = m_readLUT.GetDeltaZ(L2_saddress, L2_etaMap, L2_phiMap, L2_phi, SP_br_inner_r);
        //if(L2_saddress==1) cout << "dZ=" << dZ << endl; 
        float calc_shift_br_sagitta = calc_br_sagitta;
        if(L2_saddress==1){
          float new_dZ = m_readLUT.GetNewDeltaZ(diveta_br, phiBin_br, L2_charge);
          float shift_SP_br_middle_z = SP_br_middle_z + dZ; 
          calc_shift_br_sagitta = m_util.calcSagitta(SP_br_inner_z,SP_br_inner_r,
              shift_SP_br_middle_z,SP_br_middle_r,
              SP_br_outer_z,SP_br_outer_r);
          calc_shift_br_radius = m_util.computeRadius3Points(SP_br_inner_z,SP_br_inner_r,
            shift_SP_br_middle_z,SP_br_middle_r,
            SP_br_outer_z,SP_br_outer_r);
        }
        if(barrelSagitta) hMDT_pTvsBarrelSagitta[icharge][L2_saddress][diveta_br][phiBin_br]->Fill(offline_pt, calc_br_sagitta);
        if(LSSagitta) hMDT_pTvsLSSagitta[icharge][ispr][diveta_br][phiBin_br]->Fill(offline_pt, calc_br_sagitta);
        if(barrelSagittaComb) hMDT_pTvsBarrelSagittaComb[L2_saddress][diveta_br][phiBin_br]->Fill(offline_pt, calc_br_sagitta);
      }
      if (L2_br_radius > 1e-5) {
        if (SP_br_inner_r>1e-5 && SP_br_middle_r>1e-5 && SP_br_outer_r>1e-5) ;
        else continue;
        if (makeLUTBarrel) hMDT_invpTvsinvRadius[icharge][L2_saddress][diveta_br][phiBin_br]->Fill(calc_shift_br_radius, offline_pt);
        if (L2_saddress==1) {
          SPBILSR->Fill(SP_br_inner_r);
          float new_radius = m_readLUT.shiftRadiusLS(L2_etaMap, L2_phiMap, SP_br_inner_r, L2_br_radius, L2_charge);
          //cout << "radius default/new=" << L2_br_radius << "/" << new_radius << endl;
          if(makeLUTLargeSpecial) {
            //hMDT_invpTvsinvRadiusLS[iqeta][ispr][diveta_br][phiBin_br]->Fill(L2_br_radius, offline_pt);
            hMDT_invpTvsinvRadiusLS[iqeta][ispr][diveta_br][phiBin_br]->Fill(calc_shift_br_radius, offline_pt);
            //hMDT_invpTvsinvRadiusLS[iqeta][ispr][diveta_br][phiBin_br]->Fill(new_radius, offline_pt);
          }
          if(makeLUTLargeSpecialLight) hMDT_invpTvsinvRadiusLSLight[iqeta][ispr][diveta_br]->Fill(L2_br_radius, offline_pt);
          if(makeLUTLargeSpecialNoSpr) hMDT_invpTvsinvRadiusLSNoSpr[iqeta][diveta_br][phiBin_br]->Fill(L2_br_radius, offline_pt);
          //if(makeLUTLargeSpecialInverse) hMDT_invpTvsinvRadiusLSInverse[iqeta][ispr][diveta_br][phiBin_br]->Fill(1./offline_pt, 1./L2_br_radius);
          //if(makeLUTLargeSpecialInverse) hMDT_invpTvsinvRadiusLSInverse[iqeta][ispr][diveta_br][phiBin_br]->Fill(1./offline_pt, 1./new_radius);
          if(makeLUTLargeSpecialInverse) hMDT_invpTvsinvRadiusLSInverse[iqeta][ispr][diveta_br][phiBin_br]->Fill(1./offline_pt, 1./calc_shift_br_radius);
          if(iqeta==0 && ispr==1 && diveta_br==3 && phiBin_br==9){
            if(L2_phiMap<-0.08) continue;
            invpTvsinvRadiusLS010309->Fill(1./offline_pt, 1./L2_br_radius);
            etaMap010309->Fill(L2_etaMap);
            phiMap010309->Fill(L2_phiMap);
            SPBIR010309->Fill(SP_br_inner_r);
            SPBIRvsPhi010309->Fill(L2_phi,SP_br_inner_r);
            if(SP_br_inner_r<6160)invpTvsinvRadiusLSSPRlow010309->Fill(1./offline_pt, 1./L2_br_radius);
            else invpTvsinvRadiusLSSPRhigh010309->Fill(1./offline_pt, 1./L2_br_radius);
            //cout << "SPR BI/BM/BO=" << SP_br_inner_r << "/" << SP_br_middle_r << "/" << SP_br_outer_r << endl;
          }
        }
      }
    }
  }
  fout->Write();
}

//////////////////////////////////////////////////////////////
void newValidationT::fitResidual(string input, string output){  
  
  TFile *file = TFile::Open(input.c_str());
  TFile *fout = new TFile(output.c_str(),"recreate");
  
  TH1 *PtResidual[nChamber][ndiveta-1][ndivpt-1];
  float ptcenter[ndivpt-1],pterr[ndivpt-1];
  stringstream name;
  for(int i=0;i<nChamber;i++){ 
    for(int j=0;j<ndiveta-1;j++){ 
      for(int k=0;k<ndivpt-1;k++){ 
        name << "PtResidualSaddress" << i << "Eta" << j << "Pt" << k;
        PtResidual[i][j][k] = dynamic_cast<TH1*>(file->Get(name.str().c_str()));
        name.str("");
        ptcenter[k] = 0.5*(divpt[k]+divpt[k+1]);
        pterr[k] = divpt[k+1]-ptcenter[k];
      }
    }
  }
  
  TF1 *func1 = new TF1("func1","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2)");
  TF1 *func2 = new TF1("func2","[0]*exp(-0.5*((x-[1])/[2])^2)");
  TF1 *func3 = new TF1("func3","[0]*exp(-0.5*((x-[1])/[2])^2)");
  float mean[nChamber][ndiveta-1][ndivpt-1],mean_err[nChamber][ndiveta-1][ndivpt-1];
  float sigma[nChamber][ndiveta-1][ndivpt-1],sigma_err[nChamber][ndiveta-1][ndivpt-1];
  float sigma2[nChamber][ndiveta-1][ndivpt-1],sigma_err2[nChamber][ndiveta-1][ndivpt-1];
  stringstream canvas_name,plot_name;
  TGraphErrors *graph[nChamber][ndiveta],*graph2[nChamber][ndiveta],*graph3[nChamber][ndiveta];
  
  for (int i=0;i<nChamber;i++){
    for (int j=0;j<ndiveta-1;j++){
      for (int k=0;k<ndivpt-1;k++){
        float tmp_mean = PtResidual[i][j][k]->GetMean();
        float tmp_sigma = PtResidual[i][j][k]->GetRMS();
        func1->SetLineColor(2);
        func1->SetParameter(1,tmp_mean);
        func1->SetParameter(2,tmp_sigma);
        func1->SetParLimits(0,0,PtResidual[i][j][k]->GetMaximum());
        func1->SetParLimits(3,0,PtResidual[i][j][k]->GetMaximum());
        func1->SetParameter(4,tmp_mean);
        func1->SetParLimits(5,0,tmp_sigma);
        PtResidual[i][j][k]->Fit(func1,"","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        if(func1->GetParameter(2) > func1->GetParameter(5)){
          func2->FixParameter(0,func1->GetParameter(0));
          func2->FixParameter(1,func1->GetParameter(1));
          func2->FixParameter(2,func1->GetParameter(2));
          func3->FixParameter(0,func1->GetParameter(3));
          func3->FixParameter(1,func1->GetParameter(4));
          func3->FixParameter(2,func1->GetParameter(5));
          sigma[i][j][k] = func1->GetParameter(5);
          sigma_err[i][j][k] = func1->GetParError(5);
          mean[i][j][k] = func1->GetParameter(4);
          mean_err[i][j][k] = func1->GetParError(4);
          sigma2[i][j][k] = func1->GetParameter(2);
          sigma_err2[i][j][k] = func1->GetParError(2);
        } else{
          func3->FixParameter(0,func1->GetParameter(0));
          func3->FixParameter(1,func1->GetParameter(1));
          func3->FixParameter(2,func1->GetParameter(2));
          func2->FixParameter(0,func1->GetParameter(3));
          func2->FixParameter(1,func1->GetParameter(4));
          func2->FixParameter(2,func1->GetParameter(5));
          sigma[i][j][k] = func1->GetParameter(2);
          sigma_err[i][j][k] = func1->GetParError(2);
          mean[i][j][k] = func1->GetParameter(1);
          mean_err[i][j][k] = func1->GetParError(1);
          sigma2[i][j][k] = func1->GetParameter(5);
          sigma_err2[i][j][k] = func1->GetParError(5);
        }
        func2->SetLineColor(4);
        func3->SetLineColor(6);
        PtResidual[i][j][k]->Fit(func2,"+","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        PtResidual[i][j][k]->Fit(func3,"+","",tmp_mean-3*tmp_sigma,tmp_mean+3*tmp_sigma);
        fout->Add(PtResidual[i][j][k]);
      }
      plot_name << "PtResolutionSharpChamber" << i << "Eta" << j ;
      graph[i][j] = new TGraphErrors(ndivpt-1,ptcenter,sigma[i][j],pterr,sigma_err[i][j]);
      graph[i][j]->SetName(plot_name.str().c_str());
      fout->Add(graph[i][j]);
      plot_name.str("");
      plot_name << "PtResolutionWideChamber" << i << "Eta" << j ;
      graph2[i][j] = new TGraphErrors(ndivpt-1,ptcenter,sigma2[i][j],pterr,sigma_err2[i][j]);
      graph2[i][j]->SetName(plot_name.str().c_str());
      fout->Add(graph2[i][j]);
      plot_name.str("");
      plot_name << "PtMeanSharpChamber" << i << "Eta" << j ;
      graph3[i][j] = new TGraphErrors(ndivpt-1,ptcenter,mean[i][j],pterr,mean_err[i][j]);
      graph3[i][j]->SetName(plot_name.str().c_str());
      fout->Add(graph3[i][j]);
      plot_name.str("");
    }
  }

  fout->Write();
}
