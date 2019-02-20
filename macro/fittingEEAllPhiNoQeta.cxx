#include "macro/fittingEEAllPhiNoQeta.h"
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

void fitting::DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhigh_fit, TString outputtxtname,int phi,int eta,int icharge,int ieta){
  double sag,sagerr;
  ofstream fout;
  fout.open(outputtxtname, std::ios::app);

  TGraphErrors* g=MakeGraph( h2d, nbins, binlow );

  int numP = GetPointsNum(h2d,nbins,binlow);
  if (numP<5) doFit=false;
  
  if(doFit){
    testFunc2->SetParLimits(1,0,1000);
    testFunc2->FixParameter(2,0);//offset=0
    testFunc2->SetLineColor(4);
    xlow_fit=GetFirstPoint(h2d,nbins,binlow);
    cout << "xlow_fit=" << xlow_fit << endl;
    g->Fit(testFunc2, "", "", xlow_fit, xhigh_fit);
    float chi2 = testFunc2->GetChisquare()/numP;
    cout << "chi2=" << chi2 << endl;
    if (chi2>1){
      while (1){
        xhigh_fit = xhigh_fit - 0.01;
        cout << "xhigh_fit=" << xhigh_fit << endl;
        testFunc2->FixParameter(2,0);//offset=0
        g->Fit(testFunc2, "", "", xlow_fit,xhigh_fit);
        chi2 = testFunc2->GetChisquare()/numP;
        cout << "chi2=" << testFunc2->GetChisquare()/numP << endl;
        if (chi2<1)break;
      }
    }
  }

  //float par1 = testFunc->GetParameter(0);
  //float par2 = testFunc->GetParameter(1);
  float par1 = testFunc2->GetParameter(0);
  float par2 = testFunc2->GetParameter(1);
  float par3 = testFunc2->GetParameter(2);
  
  int badNum[19] = {68,69,70,71,75,76,77,78,113,114,115,118,119,120,167,168,169,173,174};//v2
  for (int phibin=0;phibin<19;phibin++){
    if (phi==badNum[phibin]){
      par1=0;
      par2=0;
      par3=0;
    }
  }
  
  float x[6] = {0,0.05,0.1,0.15,0.2,0.25};
  float y[6] = {0,0,0,0,0,0};
  //cout << "par1/par2=" << par1 << "/" << par2 << endl;
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
  //fout << eta << " " << phi << " " << testFunc->GetParameter(0) << " " << testFunc->GetParameter(1) << " " << icharge << " " << ieta << endl;
  fout << eta << " " << phi << " " << par1 << " " << par2 << " " << par3 << " " << icharge << ieta << endl;
  leg->Draw();
  TString pname="pictureEEAllPhiNoQeta/version9_final/";
  pname+=g->GetName();
  pname+=".root";
  TFile *gfile = new TFile(pname,"recreate");
  gfile->Add(g);
  gfile->Add(fitpol2);
  gfile->Write();
  gfile->Close();
  delete gfile;

}

TGraphErrors* fitting::MakeGraph(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );
    //if( h1d->GetEntries() < 1000 ){
    //  continue;
    //}

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
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
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
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

float fitting::GetFirstPoint(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    //double ptcenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]) )/2.;
    //double pterror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]) - ptcenter );
    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;
    //cout << "sag/sagerr=" << sag << "/" << sagerr << endl;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.0001) continue;
    if (sagerr>0.000005) continue;

    i=ibin;
    break;
  }

  float lowlimit=i*0.01;
  return lowlimit;
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
  if( h1d->GetEntries() < 100 ) return make_pair( -1e30, -1e30 );
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

