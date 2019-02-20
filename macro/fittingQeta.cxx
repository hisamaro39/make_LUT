#include "macro/fittingQeta.h"
#include "TFile.h"

fitting::fitting(  TString mode ){
  gausFunc=new TF1("gausFunc", "gaus");
  testFunc=new TF1("testFunc", "[0]*x+[1]*x*x", 0, 0.25);//default from 0 to 0.16
  Mode=mode;
}

fitting::~fitting(){}

void fitting::DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhigh_fit, TString outputtxtname,int phi,int eta,int qeta){
  double sag,sagerr;
  ofstream fout;
  fout.open(outputtxtname, std::ios::app);

  TGraphErrors* g=MakeGraph( h2d, nbins, binlow );
  
  int numP = GetPointsNum(h2d, nbins, binlow);
  if (numP<5) doFit=false;
  float highLimit = GetHighLimit(h2d, nbins, binlow);//decide high limit
  xhigh_fit = highLimit;
  float first_slope = GetSlope(h2d, nbins, binlow);

  float par1=0,par2=0;
  if (doFit){
    testFunc->SetLineColor(4);
    testFunc->SetLineWidth(1);
    testFunc->SetParameter(0,first_slope);
    testFunc->SetParLimits(1,0,100);
    /*float chi2[10];
    float chi2_min=1000;
    int minnum=-1;
    for (int k=0;k<10;k++){
      xhigh_fit=0.16+0.01*k;
      g->Fit(testFunc,"","",xlow_fit,xhigh_fit);
      chi2[k]=testFunc->GetChisquare()/(16+k);
      cout << "xlow_fit/xhigh_fit/chi2=" << xlow_fit << "/" << xhigh_fit << "/" << chi2[k] << endl;
      if (chi2[k]<chi2_min){
        chi2_min=chi2[k];
        minnum=k;
      }
    }
    xhigh_fit=0.16+0.01*minnum;
    cout << "final xhigh_fit=" << xhigh_fit << endl;*/
    g->Fit(testFunc,"","",xlow_fit,xhigh_fit);
    //testFunc->SetParLimits(0,0,5);
    //g->Fit(testFunc, "", "", xlow_fit,0.14);
    //cout << "chi2=" << testFunc->GetChisquare()/numP << endl;
    /*int num1=15,num2=20,num3=25;
    float chi2_1=g->Fit(testFunc, "", "", xlow_fit,0.15)/15;
    float chi2_2=g->Fit(testFunc, "", "", xlow_fit,0.2)/20;
    float chi2_3=g->Fit(testFunc, "", "", xlow_fit,0.25)/25;
    if (chi2_1<chi2_2){
      if (chi2_3<chi2_1) xhigh_fit=0.25;
      else xhigh_fit=0.15;
    }
    else{
      if (chi2_3<chi2_2) xhigh_fit=0.25;
      else xhigh_fit = 0.2;
    }
    g->Fit(testFunc, "", "", xlow_fit,xhigh_fit);
    */
    
    /*int num=15;
    xhigh_fit=0.15;
    float chi2_temp=1000;
    g->Fit(testFunc, "", "", xlow_fit,0.15);
    float chi2;
    while(1){
      chi2 = testFunc->GetChisquare()/num;
      cout << "num/highlimit=" << num << "/" << xhigh_fit << endl;
      cout << "chi2/chi2_temp=" << chi2 << "/" << chi2_temp << endl;
      if(num==25)break;
      if (chi2<chi2_temp){
        chi2_temp=chi2;
        num++;
        xhigh_fit+=0.01;
        g->Fit(testFunc, "", "", xlow_fit,xhigh_fit);
        par1 = testFunc->GetParameter(0);
        par2 = testFunc->GetParameter(1);
        cout << "nonum chi2=" << testFunc->GetChisquare() << endl;
      }
      else break;
      
    }
    if (chi2>1){
      while (1){
        xhigh_fit = xhigh_fit - 0.01;
        cout << "xhigh_fit=" << xhigh_fit << endl;
        g->Fit(testFunc, "", "", xlow_fit,xhigh_fit);
        chi2 = testFunc->GetChisquare()/numP;
        cout << "chi2=" << testFunc->GetChisquare()/numP << endl;
        if (chi2<1)break;
      }
    }*/
  }
  
  par1 = testFunc->GetParameter(0);
  par2 = testFunc->GetParameter(1);

  float x[6] = {0,0.05,0.1,0.15,0.2,0.25};
  float y[6] = {0,0,0,0,0,0};
  //cout << "par1/par2=" << par1 << "/" << par2 << endl;
  for (int i=0; i<6; i++) {
    y[i] = par1*x[i] + par2*x[i]*x[i];
    //cout << "x=" << x[i] << "::" << "y=" << y[i] << endl;
  }
  
  TCanvas* c=new TCanvas("c1", "", 800, 600);
  c->cd();
  //TH1* frame=gPad->DrawFrame(xlow_fit, -1, xhigh_fit, 1);
  TH1* frame=gPad->DrawFrame(0, 0, 0.25, 0.4);
  frame->SetTitle(";offline_pT  (GeV);Sagitta (mm)");
  frame->Draw();
  g->GetXaxis()->SetRangeUser(0,0.25);
  g->GetYaxis()->SetRangeUser(0,0.4);
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
  testFunc->Draw("same");
  fout << eta << " " << phi << " " << par1 << " " << par2 << " " << qeta << endl;
  //cout << "eta/phi/qeta/par1/par2=" 
    //<< eta << "/" << phi << "/" << qeta << "/" << par1 << "/" << par2 << "/" << endl;
  leg->Draw();
  TString pname="pictureAlphaBeta/version9_separate/";
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
  float x[nbins],y[nbins];
  for( int ibin=0; ibin<nbins; ibin++ ){
    //cout << "ibin=" << ibin << endl;
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );
    //if( h1d->GetEntries() < 1000 ){
    //  continue;
    //}

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    //cout << "binlow[" << ibin << "]=" << binlow[ibin] << endl;
    //cout << "GetBinLowEdge(binlow[" << ibin << "])=" << h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) << endl;
    //cout << "ptcenter/pterror=" << ptcenter << "/" << pterror << endl;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;
    //cout << "sag/sagerr=" << sag << "/" << sagerr << endl;
    //cout << "sag_over_sagerr=" << sag/sagerr << endl;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.4) continue;
    if (sagerr>0.02) continue;

    g->GetXaxis()->SetRangeUser(0,0.25);
    g->GetYaxis()->SetRangeUser(0,0.4);
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
  if( h1d->GetEntries() < 1000 ) return make_pair( -1e30, -1e30 );
  double mean=h1d->GetMean();
  double mean_first=h1d->GetMean();
//  double mean=0;
  double sigma=h1d->GetRMS();
  double sigma_first=h1d->GetRMS();
  //if (sigma > 0.04) return make_pair( -1e30, -1e30 );
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

  //fit probably failed
  if (fabs(mean_first-gausFunc->GetParameter(1))>1.2*sigma_first)
    return make_pair( -1, -1);

  return make_pair( gausFunc->GetParameter(1),
		    gausFunc->GetParError(1) );
}

int fitting::GetPointsNum(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    //double radiuscenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double radiuserror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) - radiuscenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.4) continue;
    if (sagerr>0.02) continue;

    i++;
  }

  return i;
}

float fitting::GetHighLimit(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  float x[nbins],y[nbins];
  float maxnum = 0;
  int maxbin=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    x[ibin]=0,y[ibin]=0;
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.4) continue;
    if (sagerr>0.02) continue;
    if (ibin>9){
      cout << "ibin/sag=" << ibin << "/" << sag << endl;
      if (sag>1e-5 && maxnum<1e-5) {
        cout << "sag>0 && maxnum==0" << endl;
        maxnum=sag;
      }
      else if (sag>maxnum) {
        maxnum=sag;
        maxbin=ibin;
        cout << "sag>maxnum maxnum/maxbin=" << maxnum << "/" << maxbin << endl;
      }
      else if (sag<maxnum && sag>1e-5){
        cout << "sag<maxnum" << endl;
        break;
      }
    }
  }
  float high;
  high = 0.01*(maxbin+1);
  
  cout << "high limit = " << high << endl;

  return high;
}

float fitting::GetHighLimit2(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  float x[nbins],y[nbins];
  for( int ibin=0; ibin<nbins; ibin++ ){
    x[ibin]=0,y[ibin]=0;
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.4) continue;
    if (sagerr>0.02) continue;
    x[ibin]=ptcenter;y[ibin]=sag;
  }

  float max_slope = (y[2]-y[0])/(x[2]-x[0]);
  cout << "first slope=" << max_slope << endl;
  int max_bin=0;
  for (int j=0;j<7;j++){
    float this_slope=0;
    if(y[3*j+5]>1e-5 && y[3*j+2]>1e-5) 
      this_slope = (y[3*j+5]-y[3*j+3])/(x[3*j+5]-x[3*j+3]);
    cout << "No." << j+1 << " slope!" << endl;
    cout << "x1/y1=" << x[3*j+3] << "/" << y[3*j+3] << endl;
    cout << "x2/y2=" << x[3*j+5] << "/" << y[3*j+5] << endl;
    cout << "this slope=" << this_slope << endl;
    if (this_slope>max_slope && fabs(this_slope)>1e-5) {
      max_slope=this_slope;
      max_bin = 3*(j+1);
    }
    //else break;
  }

  float high;
  high = 0.01*max_bin;
  
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
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]-1) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
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

