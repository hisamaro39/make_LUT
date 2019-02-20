#include "macro/fittingEENoSL.h"
#include "TFile.h"

fitting::fitting(  TString mode ){
  gausFunc=new TF1("gausFunc", "gaus");
  invXFunc=new TF1("invXFunc", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0, 0.16);
  invX2Func=new TF1("invX2Func", "[0]+[1]*x+[2]*x*x", 0, 0.16);
  testFunc=new TF1("testFunc", "[0]*x+[1]*x*x", 0, 0.25);//default from 0 to 0.16
  testFunc2=new TF1("testFunc2", "[0]*x+[1]*x*x+[2]", 0, 0.25);//default from 0 to 0.16
  Mode=mode;
}

fitting::~fitting(){}

void fitting::DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhigh_fit, TString outputtxtname,int phi,int eta,int qeta){
  double sag,sagerr;
  ofstream fout;
  fout.open(outputtxtname, std::ios::app);

  TGraphErrors* g=MakeGraph( h2d, nbins, binlow );

  int numP = GetPointsNum(h2d,nbins,binlow);
  if (numP<5) doFit=false;
  float highLimit = GetHighLimit(h2d, nbins, binlow);
  xhigh_fit = highLimit;
  //float first_slope = GetSlope(h2d, nbins, binlow);
  
  if(doFit){
    //testFunc->SetParLimits(1,0,10);
    //testFunc->SetLineColor(4);
    //testFunc->SetLineWidth(1);
    testFunc2->SetParLimits(1,0,10);
    testFunc2->SetLineColor(4);
    testFunc2->SetLineWidth(1);
    g->Fit(testFunc2, "", "", xlow_fit, xhigh_fit);
    //g->Fit(testFunc, "", "", xlow_fit, xhigh_fit);
    //float par1 = testFunc->GetParameter(0);
    //float par2 = testFunc->GetParameter(1);
    if (testFunc2->GetParameter(0)<1e-5 && testFunc2->GetParameter(1)<1e-5 )
      g->Fit(testFunc2, "", "", xlow_fit,0.25);
  }
  float x[6] = {0,0.05,0.1,0.15,0.2,0.25};
  float y[6] = {0,0,0,0,0,0};
  //float par1 = testFunc->GetParameter(0);
  //float par2 = testFunc->GetParameter(1);
  float par1 = testFunc2->GetParameter(0);
  float par2 = testFunc2->GetParameter(1);
  float par3 = testFunc2->GetParameter(2);
  //cout << "par1/par2=" << par1 << "/" << par2 << endl;
  //for (int i=0; i<6; i++) y[i] = par1*x[i] + par2*x[i]*x[i];
  cout << "par1/par2/par3=" << par1 << "/" << par2 << "/" << par3 << endl;
  for (int i=0; i<6; i++) y[i] = par1*x[i] + par2*x[i]*x[i] + par3;
  
  TCanvas* c=new TCanvas("c1", "", 800, 600);
  c->cd();
  //TH1* frame=gPad->DrawFrame(xlow_fit, -1, xhigh_fit, 1);
  TH1* frame=gPad->DrawFrame(0, 0, 0.25, 0.00004);
  frame->SetTitle(";offline_pT  (GeV);Sagitta (mm)");
  frame->Draw();
  g->GetXaxis()->SetRangeUser(0,0.25);
  g->GetYaxis()->SetRangeUser(0,0.00004);
  g->SetMarkerColor(1);
  g->SetMarkerStyle(1);
  g->SetLineWidth(1);
  g->Draw("P");
  TGraph *fitpol2 = new TGraph(6,x,y);
  fitpol2->SetName("fitpol2");
  fitpol2->SetLineColor(4);
  fitpol2->Draw("c");
  
  TLegend* leg=new TLegend(0.5, 0.12, 0.9, 0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry( g, "data", "lep");
  //fout << eta << " " << phi << " " << testFunc->GetParameter(0) << " " << testFunc->GetParameter(1) << " " << qeta << endl;
  fout << eta << " " << phi << " " << par1 << " " << par2 << " " << par3 << " " << qeta << endl;
  leg->Draw();
  TString pname="pictureEENoSL/version7/";
  pname+=g->GetName();
  pname+=".root";
  TFile *gfile = new TFile(pname,"recreate");
  gfile->Add(g);
  gfile->Add(fitpol2);
  gfile->Write();
  gfile->Close();
  delete gfile;

  //g->SaveAs(pname);
}

TGraphErrors* fitting::MakeGraph(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );
    //if( h1d->GetEntries() < 1000 ){
    //  continue;
    //}

    double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;
    //cout << "sag/sagerr=" << sag << "/" << sagerr << endl;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.0001) continue;
    if (sagerr>0.000005) continue;

    g->GetXaxis()->SetRangeUser(0,0.25);
    g->GetYaxis()->SetRangeUser(0,0.00004);
    g->SetPoint( i, ptcenter, sag );
    g->SetPointError( i, pterror, sagerr );
    i++;
  }
  //cout << "i=" << i << endl;
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
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;
    //cout << "sag/sagerr=" << sag << "/" << sagerr << endl;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.0001) continue;
    if (sagerr>0.000005) continue;

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
//  h1d->Rebin(2);
  return h1d;
}

pair<double, double> fitting::GetFitValues( TH1* h1d ){
  if( h1d->GetEntries() < 500 ) return make_pair( -1e30, -1e30 );
  double mean=h1d->GetMean();
  double mean_first=h1d->GetMean();
  double sigma=h1d->GetRMS();
  double sigma_first=h1d->GetRMS();
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );
  mean=gausFunc->GetParameter(1);
  sigma=gausFunc->GetParameter(2);
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );
  mean=gausFunc->GetParameter(1);
  sigma=gausFunc->GetParameter(2);
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );
  mean=gausFunc->GetParameter(1);
  sigma=gausFunc->GetParameter(2);
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );

  if (fabs(mean_first-gausFunc->GetParameter(1))>1.2*sigma_first)
    return make_pair( -1, -1);

  return make_pair( gausFunc->GetParameter(1),
		    gausFunc->GetParError(1) );
}

float fitting::GetHighLimit(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  float x[nbins],y[nbins];
  float maxnum = 0;
  int maxbin=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    x[ibin]=0,y[ibin]=0;
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.4) continue;
    if (sagerr>0.02) continue;
    cout << "ibin/sag=" << ibin << "/" << sag << endl;
    if (sag>1e-5 && maxnum<1e-5) maxnum=sag;
    else if (sag<maxnum && sag>1e-5){
      maxnum=sag;
      maxbin=ibin;
      break;
    }
    maxnum=sag;
  }
  float high;
  high = 0.01*maxbin;
  
  cout << "high limit = " << high << endl;

  return high;
}

float fitting::GetSlope(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  float x[nbins+1],y[nbins+1];
  x[0]=0.;y[0]=0.;
  for( int ibin=1; ibin<nbins+1; ibin++ ){
    x[ibin]=0,y[ibin]=0;
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.4) continue;
    if (sagerr>0.02) continue;
    x[ibin]=ptcenter;
    y[ibin]=sag;

  }
  float this_slope = 0;
  for (int k=0;k<nbins;k++){
    cout << "k=" << k << endl;
    this_slope = (y[k+1]-y[k])/(x[k+1]-x[k]);
    cout << "slope = " << this_slope << endl;
    if (this_slope>1e-5) break;
  }

  return this_slope;
}

