#include "macro/fittingCorrectBarrelSagittaLS.h"
#include "TFile.h"

fitting::fitting(  TString mode ){
  gausFunc=new TF1("gausFunc", "gaus");
  testFunc2=new TF1("testFunc2", "[0]*x+[1]*x*x+[2]", 0, 0.25);//default from 0 to 0.16
  //invFunc=new TF1("invFunc", "[0]/(x*x*x)+[1]/(x*x)+[2]/x+[3]", 0, 80);
  //invFunc=new TF1("invFunc", "[0]/(x*x)+[1]/x+[2]", 0, 80);
  invFunc=new TF1("invFunc", "[0]/x+[1]", 0, 80);
  flatFunc=new TF1("flatFunc", "[0]", 0, 80);
  expFunc=new TF1("expFunc","[0]*exp([1]*x)+[2]");
  Mode=mode;
}

fitting::~fitting(){}

void fitting::DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhigh_fit, TString outputtxtname,int charge,int chamber,int ieta,int iphi){
  double sag,sagerr;
  ofstream fout;
  fout.open(outputtxtname, std::ios::app);

  int numP=0,minpt=0,maxpt=0;
  TGraphErrors* g=MakeGraph( h2d, nbins, binlow, numP, minpt, maxpt );
  cout << "numP/minpt/maxpt==" << numP << "/" << minpt << "/" << maxpt << endl;

  if (numP<10) doFit=false;
  int nfitloop = (maxpt-minpt)/5;
  cout << "nfitloop=" << nfitloop << endl;
  
  float min_chi2=9999,par1=0,par2=0,par3=0;
  if(doFit){
    if (charge==0) invFunc->SetParLimits(0,0,100000000);
    else  invFunc->SetParLimits(0,-100000000,0);
    g->Fit(invFunc, "", "", 0, 80);
    par1 = invFunc->GetParameter(0);
  }

  /*if(doFit){
    for(int ifit=0;ifit<nfitloop;ifit++){
      if (charge==0) invFunc->SetParLimits(0,0,100000000);
      else  invFunc->SetParLimits(0,-100000000,0);
      if (charge==0) invFunc->SetParLimits(1,0,100000000);
      else  invFunc->SetParLimits(1,-100000000,0);
      invFunc->SetParLimits(2,-800,800);
      invFunc->SetLineColor(2);
      xlow_fit = minpt;
      xhigh_fit = minpt+5*(ifit+1);
      cout << "xlowfit/xhighfit=" << xlow_fit << "/" << xhigh_fit << endl;
      //g->Fit(invFunc, "", "", xlow_fit, xhigh_fit);
      g->Fit(invFunc, "", "", 0, 80);
      float chi2 = invFunc->GetChisquare()/(xhigh_fit-xlow_fit);
      cout << "chi2=" << chi2 << endl;
      if(chi2<min_chi2){
        min_chi2=chi2;
        par1 = invFunc->GetParameter(0);
        par2 = invFunc->GetParameter(1);
        par3 = invFunc->GetParameter(2);
      }
    }
  }*/

  if (fabs(par3)>799) par3=0;
  
  TCanvas* c=new TCanvas("c1", "", 800, 600);
  c->cd();
  TH1* frame=gPad->DrawFrame(0, -800, 80, 800);
  frame->SetTitle(";offline_pT  (GeV);Sagitta (mm)");
  frame->Draw();
  g->GetXaxis()->SetRangeUser(0,80);
  g->GetYaxis()->SetRangeUser(-800,800);
  g->SetMarkerColor(1);
  g->SetMarkerStyle(1);
  g->SetLineWidth(1);
  g->Draw("P");
  
  //fout << charge << " " << chamber << " " << ieta << " " << iphi << " " << par1 << " " << par2 << " " << par3 << endl;
  fout << charge << " " << chamber << " " << ieta << " " << iphi << " " << par1 << endl;
  
  TString pname="outputCorrectBarrelSagitta/pictureLS/";
  pname+=g->GetName();
  pname+=".root";
  TFile *gfile = new TFile(pname,"recreate");
  gfile->Add(g);
  gfile->Write();
  gfile->Close();
  delete gfile;
}

TGraphErrors* fitting::MakeGraph(TH2* h2d, int nbins, double* binlow, int &numP, int &minpt, int &maxpt){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    double ptcenter = 1 + 2*ibin;
    double pterror = 1;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 80 || fabs(ptcenter) < 0 ) continue;
    if (sag<-800 || sag>800) continue;

    g->GetXaxis()->SetRangeUser(0,80);
    g->GetYaxis()->SetRangeUser(-800,800);
    g->SetPoint( i, ptcenter, sag );
    g->SetPointError( i, pterror, sagerr );
    if(i==0) minpt=ptcenter; 
    maxpt=ptcenter;
    i++;
  }
  numP=i;
  TString gname="g_";
  gname+=h2d->GetName();
  gname+=Mode;
  g->SetName(gname);
  g->SetTitle(gname);

  return g;
}

TH1* fitting::GetTH1( TH2* h2d, int binNo, int ilow, int ihigh ){
  TString hname=h2d->GetName();
  hname+="_";
  hname+=Mode;
  hname+=binNo;
  TH1* h1d=((TH2*)h2d)->ProjectionY(hname, ilow, ihigh, "");
  return h1d;
}

pair<double, double> fitting::GetFitValues( TH1* h1d ){
  if( h1d->GetEntries() < 30 ) return make_pair( -1e30, -1e30 );
  double mean=h1d->GetMean();
  double sigma=h1d->GetRMS();
  gausFunc->SetParLimits(1,mean-1.2*sigma, mean+1.2*sigma);
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );
  mean=gausFunc->GetParameter(1);
  sigma=gausFunc->GetParameter(2);
  gausFunc->SetParLimits(1,mean-1.2*sigma, mean+1.2*sigma);
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );

  if(fabs(gausFunc->GetParameter(1) - h1d->GetMean())>50) 
    return make_pair( -1e30, -1e30 );

  if(gausFunc->GetParameter(2)>100) return make_pair(-1e30, -1e30);

  return make_pair( gausFunc->GetParameter(1),
		    gausFunc->GetParameter(2) );
}

