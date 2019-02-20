#ifndef fittingCorrectBarrelSagittaLS_h
#define fittingCorrectBarrelSagittaLS_h

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
  
  void DrawGraph( TH2* h2d, bool doFit, int nbins, double* binlow, float xlow_fit, float xhish_fit, TString outputtxtname ,int charge ,int chamber ,int ieta, int iphi);
  TGraphErrors* MakeGraph( TH2* h2d, int nbins, double* binlow, int &numP, int &minpt, int &maxpt );
  TH1* GetTH1( TH2* h2d, int binNo, int ilow, int ihigh );
  pair<double, double> GetFitValues( TH1* h1d );

 private:
  TF1* gausFunc;
  TF1* testFunc2;
  TF1* invFunc;
  TF1* flatFunc;
  TF1* expFunc;
  TString Mode;
};

#endif
