#include "macro/fittingEEAllPhiNoQetaShift.h"
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
  cout << "aho" << endl;
  
  int numP = GetPointsNum(h2d,nbins,binlow);
  if (numP<3) doFit=false;
  
  if(doFit){
    testFunc2->SetParLimits(0,0,1000);
    testFunc2->SetParLimits(1,0,1000);
    //testFunc2->FixParameter(2,0);//offset=0
    testFunc2->SetLineColor(4);
    g->Fit(testFunc2, "", "", xlow_fit, xhigh_fit);
  }

  cout << "aho2" << endl;

  float par1 = testFunc2->GetParameter(0);
  float par2 = testFunc2->GetParameter(1);
  float par3 = testFunc2->GetParameter(2);
  
  cout << "aho3" << endl;
  //TH1* frame = gPad->DrawFrame(0, 0, 0.25, 0.00004);
  cout << "aho4" << endl;
  //frame->SetTitle(";1/offline_pT(GeV);1/Radius(mm)");
  //frame->Draw();
  g->GetXaxis()->SetRangeUser(0,0.25);
  g->GetYaxis()->SetRangeUser(0,0.00004);
  g->SetMarkerColor(1);
  g->SetMarkerStyle(1);
  g->SetLineWidth(1);
  g->Draw("AP");
  

  fout << eta << " " << phi << " " << par1 << " " << par2 << " " << par3 << " " << icharge << ieta << endl;
  TString pname="pictureEEAllPhiNoQetaShift/data15/";
  pname+=g->GetName();
  pname+=".root";
  TFile *gfile = new TFile(pname,"recreate");
  gfile->Add(g);
  gfile->Write();
  gfile->Close();
  delete gfile;

}

TGraphErrors* fitting::MakeGraph(TH2* h2d, int nbins, double* binlow){
  TGraphErrors* g=new TGraphErrors;
  
  int i=0;
  for( int ibin=0; ibin<nbins; ibin++ ){
    TH1* h1d=GetTH1( h2d, ibin, binlow[ibin]+1, binlow[ibin+1] );

    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;
    cout << "sag/sagerr=" << sag << "/" << sagerr << endl;

    if( fabs(ptcenter) > 0.26 || fabs(ptcenter) < 0.001 ) continue;
    if (sag<0 || sag>0.0001) continue;
    if (sagerr>0.000005) continue;

    g->GetXaxis()->SetRangeUser(0,0.25);
    g->GetYaxis()->SetRangeUser(0,0.00004);
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

    double ptcenter = 0.005 + 0.01*ibin;
    double pterror = 0.005;
    pair<double, double> par=GetFitValues( h1d );
    double sag=par.first;
    double sagerr=par.second;

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
  return h1d;
}

pair<double, double> fitting::GetFitValues( TH1* h1d ){
  if( h1d->GetEntries() < 50 ) return make_pair( -1e30, -1e30 );
  double mean=h1d->GetMean();
  double sigma=h1d->GetRMS();
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );
  mean=gausFunc->GetParameter(1);
  sigma=gausFunc->GetParameter(2);
  h1d->Fit( gausFunc, "", "", mean-1.2*sigma, mean+1.2*sigma );
  
  if(fabs(gausFunc->GetParameter(1) - h1d->GetMean())>0.000005) 
    return make_pair( -1e30, -1e30 );

  if (gausFunc->GetParameter(2)>0.000005) 
    return make_pair(-1e30, -1e30);

  return make_pair( gausFunc->GetParameter(1),
		    gausFunc->GetParameter(2) );
}

