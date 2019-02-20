#include "macro/fittingCorrectSagitta.h"
#include "TFile.h"

fitting::fitting(  TString mode ){
  gausFunc=new TF1("gausFunc", "gaus");
  testFunc2=new TF1("testFunc2", "[0]*x+[1]*x*x+[2]", 0, 0.25);//default from 0 to 0.16
  //invFunc=new TF1("invFunc", "[0]/x+[1]", 1, 80);
  invFunc=new TF1("invFunc", "[0]/(x*x)+[1]/x+[2]", 0, 80);
  Mode=mode;
}

fitting::~fitting(){}

void fitting::DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhigh_fit, TString outputtxtname,int phi,int eta,int icharge,int ieta){
  double sag,sagerr;
  ofstream fout;
  fout.open(outputtxtname, std::ios::app);

  TGraphErrors* g=MakeGraph( h2d, nbins, binlow );

  int numP = GetPointsNum(h2d,nbins,binlow);
  if (numP<6) doFit=false;
  
  if(doFit){
    if (icharge==0) invFunc->SetParLimits(0,0,100000000);
    else  invFunc->SetParLimits(0,-100000000,0);
    if (icharge==0) invFunc->SetParLimits(1,0,100000000);
    else  invFunc->SetParLimits(1,-100000000,0);
    invFunc->SetParLimits(2,-800,800);
    invFunc->SetLineColor(2);
    g->Fit(invFunc, "", "", xlow_fit, xhigh_fit);
    //g->Fit(invFunc, "", "", 7, 80);
  }

  float par1 = invFunc->GetParameter(0);
  float par2 = invFunc->GetParameter(1);
  float par3 = invFunc->GetParameter(2);
  
  /*float x[6] = {0,0.05,0.1,0.15,0.2,0.25};
  float y[6] = {0,0,0,0,0,0};
  //cout << "par1/par2=" << par1 << "/" << par2 << endl;
  cout << "par1/par2/par3=" << par1 << "/" << par2 << "/" << par3 << endl;
  for (int i=0; i<6; i++) y[i] = par1*x[i] + par2*x[i]*x[i] + par3;*/
  
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
  //TGraph *fitpol2 = new TGraph(6,x,y);
  //fitpol2->SetName("fitpol2");
  //fitpol2->SetLineColor(4);
  //fitpol2->Draw("c");
  
  //TLegend* leg=new TLegend(0.5, 0.12, 0.9, 0.4);
  //leg->SetBorderSize(0);
  //leg->SetFillColor(0);
  //leg->AddEntry( g, "data", "lep");
  //fout << eta << " " << phi << " " << testFunc->GetParameter(0) << " " << testFunc->GetParameter(1) << " " << icharge << " " << ieta << endl;
  fout << eta << " " << phi << " " << icharge << " " << ieta << " " << par1 << " " << par2 << " " << par3 << endl;
  //leg->Draw();
  
  TString pname="outputCorrectSagitta/picture3/";
  pname+=g->GetName();
  pname+=".root";
  TFile *gfile = new TFile(pname,"recreate");
  gfile->Add(g);
  //gfile->Add(fitpol2);
  gfile->Write();
  gfile->Close();
  delete gfile;
}

TGraphErrors* fitting::MakeGraph(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );
    //if( h1d->GetEntries() < 50 ) continue;

    double ptcenter = 1 + 2*ibin;
    double pterror = 1;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 80 || fabs(ptcenter) < 0 ) continue;
    if (sag<-800 || sag>800) continue;
    //if (sagerr>0.000005) continue;

    g->GetXaxis()->SetRangeUser(0,80);
    g->GetYaxis()->SetRangeUser(-800,800);
    g->SetPoint( i, ptcenter, sag );
    g->SetPointError( i, pterror, sagerr );
    i++;
  }
  TString gname="g_";
  gname+=h2d->GetName();
  gname+=Mode;
  g->SetName(gname);
  g->SetTitle(gname);

  return g;
}

int fitting::GetPointsNum(TH2* h2d, int nbins, double* binlow){
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
    //if (sagerr>0.000005) continue;

    i++;
  }
  return i;
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

  return make_pair( gausFunc->GetParameter(1),
		    gausFunc->GetParameter(2) );
}

