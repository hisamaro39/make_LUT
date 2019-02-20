#ifndef fittingLargeSpecialInv_h
#define fittingLargeSpecialInv_h

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
  
  void DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhish_fit, TString outputtxtname ,int qeta ,int phi ,int etabin,int phibin);
  TGraphErrors* MakeGraph( TH2* h2d, int nbins, double* binlow );
  int GetPointsNum( TH2* h2d, int nbins, double* binlow );
  float GetRangeFirst(TH2* h2d, int nbins, double* binlow);
  float GetRangeFinal(TH2* h2d, int nbins, double* binlow);
  TH1* GetTH1( TH2* h2d, int binNo, int ilow, int ihigh );
  pair<double, double> GetFitValues( TH1* h1d );

 private:
  TF1* gausFunc;
  TF1* invXFunc;
  TF1* invX2Func;
  TF1* testFunc;
  TF1* testFunc2;
  TString Mode;
};

#endif
