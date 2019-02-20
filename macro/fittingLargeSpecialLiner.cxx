#include "macro/fittingLargeSpecialLiner.h"
#include "TFile.h"

fitting::fitting(  TString mode ){
  gausFunc=new TF1("gausFunc", "gaus");
  testFunc=new TF1("testFunc", "[0]*x*x+[1]*x+[2]", 0, 1000000);
  testFunc2=new TF1("testFunc2", "[0]*sqrt(x)+[1]", 0, 1000000);
  linerFunc=new TF1("linerFunc", "[0]*x+[1]", 0, 1000000);
  Mode=mode;
}

fitting::~fitting(){}

void fitting::DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhigh_fit, TString outputtxtname,int qeta,int phi,int etabin,int phibin){
  ofstream fout;
  fout.open(outputtxtname, std::ios::app);

  TGraphErrors* g=MakeGraph( h2d, nbins, binlow );

  int numP = GetPointsNum(h2d, nbins, binlow);
  cout << "numP=" << numP << endl;
  if (numP<2) doFit=false;
  xlow_fit = GetRangeFirst(h2d, nbins, binlow);
  cout << "xlow_fit/xhigh_fit=" << xlow_fit << "/" << xhigh_fit << endl;

  if(doFit){
    linerFunc->SetParLimits(0,0,1000);
    linerFunc->SetLineColor(4);
    linerFunc->SetLineWidth(1);
    g->Fit(linerFunc, "", "", xlow_fit, xhigh_fit);
  }

  float x[3] = {0,500000,1000000};
  float y[3] = {0,0,0};
  float par1 = linerFunc->GetParameter(0);
  float par2 = linerFunc->GetParameter(1);
  for (int i=0; i<3; i++) y[i] = par1*x[i] + par2;
  cout << "par1/par2=" << par1 << "/" << par2 <<  endl;
  
  TCanvas* c=new TCanvas("c1", "", 800, 600);
  c->cd();
  //TH1* frame=gPad->DrawFrame(xlow_fit, -1, xhigh_fit, 1);
  TH1* frame=gPad->DrawFrame(0, 0, 1000000, 100);
  frame->SetTitle(";Radius(cm);offline pt(GeV)");
  frame->Draw();
  g->GetXaxis()->SetRangeUser(0,1000000);
  g->GetYaxis()->SetRangeUser(0,100);
  g->SetMarkerColor(1);
  g->SetMarkerStyle(1);
  g->SetLineWidth(1);
  g->Draw("P");
  TGraph *fitpol2 = new TGraph(3,x,y);
  fitpol2->SetName("fitpol2");
  fitpol2->SetLineColor(4);
  //fitpol2->Draw("l");
  
  TLegend* leg=new TLegend(0.5, 0.12, 0.9, 0.4);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry( g, "data", "lep");
  fout << qeta << " " << phi << " " << etabin << " " << phibin << " " << par1 << " " << par2 << endl;
  //if(qeta==0) fout << par1 << " " << par2 << " " << par3 << " ";
  //else fout << par4 << " " << par5 << " ";
  //if (etabin==29 && phibin==29) fout << endl;
  leg->Draw();
  TString pname="pictureBarrelRadius/data16_LS_liner/";
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
  float max_ptcenter=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    double radiuscenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]+1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]+1) )/2.;
    double radiuserror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]+1) - radiuscenter );
    cout << "binlow ibin/ibin+1=" << binlow[ibin] << "/" << binlow[ibin+1] << endl;
    cout << "radiuscenter/radiuserror=" << radiuscenter << "/" << radiuserror << endl;
    pair<double, double> par=GetFitValues( h1d );
    double ptcenter=par.first;
    double pterr=par.second;
    cout << "BinLowEdge bin/value=" << ibin << "/" << h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]+1) << endl;
    cout << "bin/center/err=" << i << "/" << ptcenter << "/" << pterr << endl;

    if( fabs(radiuscenter) > 1000000 || fabs(radiuscenter) < 0 ) continue;
    if (ptcenter<4 && pterr<0) continue;
    if (pterr>10) continue;
    //if (ptcenter<max_ptcenter) break;

    g->GetXaxis()->SetRangeUser(0,1000000);
    g->GetYaxis()->SetRangeUser(0,100);
    g->SetPoint( i, radiuscenter, ptcenter );
    g->SetPointError( i, radiuserror,pterr );
    i++;
    max_ptcenter=ptcenter;
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
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    double radiuscenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]+1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]+1) )/2.;
    double radiuserror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]+1) - radiuscenter );
    pair<double, double> par=GetFitValues( h1d );
    double ptcenter=par.first;
    double pterr=par.second;

    if( fabs(radiuscenter) > 1000000 || fabs(radiuscenter) < 0 ) continue;
    if (ptcenter<4 || pterr<0) continue;
    if (pterr>10) continue;

    i++;
  }

  return i;
}

float fitting::GetRangeFirst(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  float first=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    cout << "ibin=" << ibin << endl;
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin], binlow[ibin+1] );

    double radiuscenter=( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]+1) + h2d->ProjectionX()->GetBinLowEdge(binlow[ibin]+1) )/2.;
    double radiuserror =( h2d->ProjectionX()->GetBinLowEdge(binlow[ibin+1]+1) - radiuscenter );
    pair<double, double> par=GetFitValues( h1d );
    double ptcenter=par.first;
    double pterr=par.second;

    cout << "ptcenter/pterr=" << ptcenter << "/" << pterr << endl;
    cout << "radiuscenter/radiuserr=" << radiuscenter << "/" << radiuserror << endl;
    if( fabs(radiuscenter) > 1000000 || fabs(radiuscenter) < 0 ) continue;
    if (ptcenter<4 || pterr<0) continue;
    if (pterr>10) continue;
    first = ibin*5000 + 5000;
    break;
  }
  return first;
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
  
  if (fabs(mean_first-gausFunc->GetParameter(1))>2*sigma_first)
    return make_pair( -1e30, -1e30);

  return make_pair( gausFunc->GetParameter(1),
		    gausFunc->GetParError(1) );
}

