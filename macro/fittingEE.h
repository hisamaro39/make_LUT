#ifndef fittingEE_h
#define fittingEE_h

#include <map>
#include <fstream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TLegend.h"

class fitting{
 public:
  fitting( TString mode );
  ~fitting();
  
  void DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhish_fit, TString outputtxtname ,int phi ,int eta ,int qeta,int sl);
  TGraphErrors* MakeGraph( TH2* h2d, int nbins, double* binlow );
  TH1* GetTH1( TH2* h2d, int binNo, int ilow, int ihigh );
  pair<double, double> GetFitValues( TH1* h1d );

 private:
  TF1* gausFunc;
  TF1* invXFunc;
  TF1* invX2Func;
  TF1* testFunc;
  TString Mode;
};

#endif
