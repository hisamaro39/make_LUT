#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
using namespace std;

#define PI 3.14159265258979
const float PI_OVER_4=PI/4.;
const float PI_OVER_8=PI/4.;
const float PHI_RANGE = 12.0/PI_OVER_8;

pair<int,int> GetBinNumber(float m_tgcMid1_phi, float m_tgcMid1_eta){
  int Octant = (int)(m_tgcMid1_phi/PI_OVER_4);
  double PhiInOctant = fabs(m_tgcMid1_phi - Octant*PI_OVER_4);
  if(PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  int phiBin = static_cast<int>(PhiInOctant*PHI_RANGE);
  int etaBin = static_cast<int>((fabs(m_tgcMid1_eta)-1.)/0.05);
  if (etaBin == -1) etaBin =  0;
  if (etaBin == 30) etaBin = 29;
  if (etaBin < -0.5 || etaBin > 29.5 || phiBin < -0.5 || phiBin > 11.5) return make_pair(-1,-1);
  return make_pair(etaBin,phiBin);
}

void getBinDefault(){
  float eta,phi;
  cout << "enter eta" << endl;
  cin >> eta;
  cout << "enter phi" << endl;
  cin >> phi;

  pair<int, int> ans = GetBinNumber(phi,eta);
  cout << "etabin=" << ans.first+1 << endl;
  cout << "phibin=" << ans.second << endl;

}

