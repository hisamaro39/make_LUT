#include "fittingLargeSpecialLiner.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TGraphErrors.h"
#include <fstream>
void LoopLargeSpecialLiner(const int nbins, double* binlow, TString path, float xlow_fit, float xhigh_fit, TString inputfilename, TString outputtxtname, TString mode, TString histoname,int charge,int phi,int etabin ,int phibin){
  fitting* Fitting=new fitting( mode );    

  TFile* fin=new TFile(inputfilename);
  //////////Fill the name of histograms you like to analyze//////////
  const int nhist=1;
  TString hname[nhist]={histoname};
  ///////////////////////////////////////////////////////////////////

  for( int i=0; i<nhist; i++ ){
    gDirectory->Cd(path);

    TH2* h=(TH2F*)fin->Get(hname[i]);
    gDirectory->Append(h);

    Fitting->DrawGraph( h, true, nbins, binlow, xlow_fit, xhigh_fit, outputtxtname,charge,phi,etabin,phibin);

    TGraphErrors* g=Fitting->MakeGraph(h, nbins, binlow);
    gDirectory->Append(g);
  }
  fin->Close();
  delete fin;
}
