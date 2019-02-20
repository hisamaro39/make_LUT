#include "util.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
const double PI = 3.14159265258979;

double const PI_OVER_4 = PI/4.0;
double const PI_OVER_8 = PI/8.0;
double const PI_OVER_16 = PI/16.0;
double const PHI_RANGE = 12.0/PI_OVER_8;
double const ZERO_LIMIT = 1e-5;

util::util(){}

util::~util(){}

float util::calcCosDphi(float xi,float yi,float xm,float ym){
  float bunsi=fabs(ym-yi);
  float bunbo=sqrt(pow(xm-xi,2)+pow(ym-yi,2));
  float ans=bunsi/bunbo;
  return ans;
}
///////////////////
float util::calcSlope(float eta){
  float bunbo = exp(eta)-exp(-eta);
  float ans = 2/bunbo;
  return ans;
}
////////////////////
float util::calcDis(float r1,float z1,float r2,float z2){
  float dr = sqrt((r1-r2)*(r1-r2)+(z1-z2)*(z1-z2));
  return dr;
}
////////////////////
float util::calcDisR(float segr,float segz,float spr,float spz){
  float a = spr/spz;
  float r = a * spz;
  float dr=segr-r;
  return dr;
}
///////////////////
float util::calcRfromTrackZ(float slope, float z){
  float r = slope*z;
  return r;
}
////////////////////////
float util::calcZfromTrackR(float slope, float r){
  float z = r/slope;
  return z;
}
///////////////////////
bool util::isBadPhi(int phibinall){
  int badNum[18] = {67,68,69,70,76,77,78,113,114,115,118,119,120,167,168,172,173,174};
  bool ans=false;
  for (int i=0;i<18;i++){
    if (phibinall==badNum[i])
      ans = true;
  }
  return ans;
}
///////////////
bool util::isSmall(float eez){
  bool small = false;
  if (fabs(eez)>10000 && fabs(eez)<10600) small =true;
  else if(fabs(eez)>10600 && fabs(eez)<12000) small = false;
  return small;
}
/////////////////////////////////
pair<int,int> util::GetBinNumber(float m_tgcMid1_phi, float m_tgcMid1_eta){
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
////////////////////////////////
pair<float,float> util::calcCenter(float x1,float y1,float x2,float y2,float x3,float y3){
  float xm1=(x1+x2)/2;
  float xm2=(x2+x3)/2;
  float ym1=(y1+y2)/2;
  float ym2=(y2+y3)/2;
  float a1=(x1-x2)/(y2-y1);
  float a2=(x2-x3)/(y3-y2);
  float x0=(a2*xm2-a1*xm1-ym2+ym1)/(a2-a1);//center of circle
  float y0=a1*(x0-xm1)+ym1;//center of circle
  return make_pair(x0,y0);
}
/////////////////////////////////
float util::calcDistance(float x1,float y1,float x2,float y2,float x3,float y3){
  float xm1=(x1+x2)/2;
  float xm2=(x2+x3)/2;
  float ym1=(y1+y2)/2;
  float ym2=(y2+y3)/2;
  float a1=(x1-x2)/(y2-y1);
  float a2=(x2-x3)/(y3-y2);
  float x0=(a2*xm2-a1*xm1-ym2+ym1)/(a2-a1);//center of circle
  float y0=a1*(x0-xm1)+ym1;//center of circle
  float a = (x0-x1)/(y1-y0);//slope of sessen
  float b = y1+x1*(x1-x0)/(y1-y0);//intercept of sessen
  float d=fabs(b)/sqrt(a*a+1);
  return d;
}
//////////////////////////////////
int util::GetBinNumberEE(float phi,int sl){
  int phiBin24=-1;
  int Octant = (int)(phi/PI_OVER_4);
  double PhiInOctant = fabs(phi - Octant * PI_OVER_4);
  if (PhiInOctant > PI_OVER_8) PhiInOctant = PI_OVER_4 - PhiInOctant;
  if ( sl==0 ){//Small
    int OctantSmall = Octant;
    double PhiInOctantSmall = PhiInOctant;
    if(phi<0) PhiInOctantSmall = fabs(phi - (OctantSmall-1)*PI_OVER_4);
    phiBin24 = PhiInOctantSmall * PHI_RANGE;
  }
  else {//Large
    //phi = phi + PI_OVER_8;
    int OctantLarge = (int)(phi / PI_OVER_4);
    double PhiInOctantLarge = fabs(phi - OctantLarge * PI_OVER_4);
    if (phi<0) PhiInOctantLarge = fabs(phi - (OctantLarge-1)*PI_OVER_4);
    phiBin24 = PhiInOctantLarge * PHI_RANGE;
  }
  return phiBin24;
}
//////////////////////////////////
int util::GetBinNumberEE2(float phi){
  int phiBin=-1;
  float Phi = phi + PI_OVER_16;
  int doubleOctant = (int)(Phi/PI_OVER_8);
  double PhiInDoubleOctant = fabs(Phi - doubleOctant * PI_OVER_8);
  phiBin = PhiInDoubleOctant * PHI_RANGE;
  return phiBin;
}
///////////////////////////////
int util::GetBinNumberAllPhi(float phi){
  if (phi<0) phi = phi + 2*PI;
  float phibin = (int) (phi * 96/PI);

  return phibin;
}
/////////////////////////////////
float util::calcAlpha(float r1,float z1,float r2,float z2){
  float slope1 = r1/z1;
  float slope2 = (r2 - r1)/(z2 - z1);
  float alpha = 0;
  alpha = fabs(atan(slope1) - atan(slope2));
  return alpha;
}
//////////////////////////////////
float util::calcIntercept(float r1,float z1,float r2,float z2){
  float slope = (r1-r2)/(z1-z2);
  float intercept = r1 - slope*z1;
  return intercept;
}
/////////////////////////////////
double util::calcSagitta(double InnerZ, double InnerR,
    double EEZ, double EER,
    double MiddleZ, double MiddleR){
  double a = (MiddleZ - InnerZ)/(MiddleR - InnerR);
  double sagitta = EEZ - EER*a - InnerZ + InnerR*a;
  return sagitta;
}
////////////////////////////////
double util::computeRadius3Points(double InnerZ, double InnerR, 
    double EEZ, double EER,
    double MiddleZ, double MiddleR)
{
  double radius_EE;

  double a3;

  double m = 0.;
  double cost = 0.;
  double x0 = 0., y0 = 0., x2 = 0., y2 = 0., x3 = 0., y3 = 0.;
  double tm = 0.;

  a3 = ( MiddleZ - InnerZ ) / ( MiddleR - InnerR );

  m = a3;
  cost = cos(atan(m));
  x2 = EER - InnerR;
  y2 = EEZ - InnerZ;
  x3 = MiddleR - InnerR;
  y3 = MiddleZ - InnerZ;

  tm = x2;
  x2 = ( x2   + y2*m)*cost;
  y2 = (-tm*m + y2  )*cost;

  tm = x3;
  x3 = ( x3   + y3*m)*cost;
  y3 = (-tm*m + y3  )*cost;

  x0 = x3/2.;
  y0 = (y2*y2 + x2*x2 -x2*x3)/(2*y2);

  radius_EE = sqrt(x0*x0 + y0*y0);
  return radius_EE;
}
/////////////////////////////////////
double util::calcShiftEER(double InnerZ, double InnerR, double EEZ, double MiddleZ, double MiddleR, double sagitta){
  double a = (MiddleZ - InnerZ)/(MiddleR - InnerR);
  double EER = (EEZ-InnerZ-sagitta)/a + InnerR;
  return EER;
}
////////////////////////
int util::returnPlateau(int threshold,int barrel){//0:barrel 1:endcap
  int plateau=-1;
  if(barrel==0){
    if(threshold==4) plateau=20;
    else if(threshold==6) plateau=30;
    else if(threshold==10) plateau=30;
    else if(threshold==11) plateau=30;
    else if(threshold==14) plateau=30;
    else if(threshold==18) plateau=30;
    else if(threshold==20) plateau=30;
    else if(threshold==22) plateau=30;
    else if(threshold==24) plateau=30;
    else if(threshold==26) plateau=30;
    else if(threshold==28) plateau=35;
    else if(threshold==40) plateau=45;
    else if(threshold==50) plateau=55;
  }
  else{
    if(threshold==4) plateau=20;
    else if(threshold==6) plateau=20;
    else if(threshold==10) plateau=20;
    else if(threshold==11) plateau=20;
    else if(threshold==14) plateau=20;
    else if(threshold==18) plateau=25;
    else if(threshold==20) plateau=30;
    else if(threshold==22) plateau=30;
    else if(threshold==24) plateau=30;
    else if(threshold==26) plateau=30;
    else if(threshold==28) plateau=35;
    else if(threshold==40) plateau=45;
    else if(threshold==50) plateau=55;
  }

  return plateau;
}
//////////////////////////
const int divpt4[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,55,60,80,100};
const int divpt6[] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,35,40,45,50,55,60,80,100};
const int divpt10[] = {0,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30,35,40,45,50,55,60,80,100};
const int divpt11[] = {0,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30,35,40,45,50,55,60,80,100};
const int divpt14[] = {0,5,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,35,40,45,50,55,60,80,100};
const int divpt18[] = {0,5,10,13,14,15,16,17,18,19,20,21,22,24,26,28,30,35,40,45,50,55,60,80,100};
const int divpt20[] = {0,5,10,15,16,17,18,19,20,21,22,23,24,25,26,28,30,35,40,45,50,55,60,80,100};
const int divpt22[] = {0,5,10,15,17,18,19,20,21,22,23,24,25,26,27,28,30,35,40,45,50,55,60,80,100};
const int divpt24[] = {0,5,10,15,19,20,21,22,23,24,25,26,27,28,29,30,35,40,45,50,55,60,80,100};
const int divpt26[] = {0,5,10,15,20,22,23,24,25,26,27,28,29,30,31,32,35,40,45,50,55,60,80,100};
const int divpt28[] = {0,5,10,15,20,23,24,25,26,27,28,29,30,31,32,33,35,40,45,50,55,60,80,100};
const int divpt40[] = {0,5,10,15,20,25,30,35,36,37,38,39,40,41,42,43,44,45,50,55,60,80,100};
const int divpt50[] = {0,5,10,15,20,25,30,35,40,45,46,47,48,49,50,51,52,53,54,55,60,80,100};
//////////////////////////
vector<int> util::devidePt(int threshold){
  vector<int> divpt;
  divpt.clear();
  if(threshold==4) 
    for (unsigned int i=0;i<sizeof(divpt4)/sizeof(int);i++) divpt.push_back(divpt4[i]);
  if(threshold==6) 
    for (unsigned int i=0;i<sizeof(divpt6)/sizeof(int);i++) divpt.push_back(divpt6[i]);
  if(threshold==10) 
    for (unsigned int i=0;i<sizeof(divpt10)/sizeof(int);i++) divpt.push_back(divpt10[i]);
  if(threshold==11) 
    for (unsigned int i=0;i<sizeof(divpt11)/sizeof(int);i++) divpt.push_back(divpt11[i]);
  if(threshold==14) 
    for (unsigned int i=0;i<sizeof(divpt14)/sizeof(int);i++) divpt.push_back(divpt14[i]);
  if(threshold==18) 
    for (unsigned int i=0;i<sizeof(divpt18)/sizeof(int);i++) divpt.push_back(divpt18[i]);
  if(threshold==20) 
    for (unsigned int i=0;i<sizeof(divpt20)/sizeof(int);i++) divpt.push_back(divpt20[i]);
  if(threshold==22) 
    for (unsigned int i=0;i<sizeof(divpt22)/sizeof(int);i++) divpt.push_back(divpt22[i]);
  if(threshold==24) 
    for (unsigned int i=0;i<sizeof(divpt24)/sizeof(int);i++) divpt.push_back(divpt24[i]);
  if(threshold==26) 
    for (unsigned int i=0;i<sizeof(divpt26)/sizeof(int);i++) divpt.push_back(divpt26[i]);
  if(threshold==28) 
    for (unsigned int i=0;i<sizeof(divpt28)/sizeof(int);i++) divpt.push_back(divpt28[i]);
  if(threshold==40) 
    for (unsigned int i=0;i<sizeof(divpt40)/sizeof(int);i++) divpt.push_back(divpt40[i]);
  if(threshold==50) 
    for (unsigned int i=0;i<sizeof(divpt50)/sizeof(int);i++) divpt.push_back(divpt50[i]);

  return divpt;
}
//////////////////
vector<float> util::getPtCenter(int threshold){
  vector<int> divpt;
  divpt.clear();
  int size;
  if(threshold==4){ 
    size = sizeof(divpt4)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt4[i]);
  }
  if(threshold==6){ 
    size = sizeof(divpt6)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt6[i]);
  }
  if(threshold==10){ 
    size = sizeof(divpt10)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt10[i]);
  }
  if(threshold==11){ 
    size = sizeof(divpt11)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt11[i]);
  }
  if(threshold==14){ 
    size = sizeof(divpt14)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt14[i]);
  }
  if(threshold==18){ 
    size = sizeof(divpt18)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt18[i]);
  }
  if(threshold==20){ 
    size = sizeof(divpt20)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt20[i]);
  }
  if(threshold==22){ 
    size = sizeof(divpt22)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt22[i]);
  }
  if(threshold==24){ 
    size = sizeof(divpt24)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt24[i]);
  }
  if(threshold==26){ 
    size = sizeof(divpt26)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt26[i]);
  }
  if(threshold==28){ 
    size = sizeof(divpt28)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt28[i]);
  }
  if(threshold==40){ 
    size = sizeof(divpt40)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt40[i]);
  }
  if(threshold==50){ 
    size = sizeof(divpt50)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt50[i]);
  }
  vector<float> ptcenter;
  ptcenter.clear();
  float center;
  for (int j=0;j<size-1;j++){
    center = 0.5*(divpt[j]+divpt[j+1]);
    ptcenter.push_back(center);
  }
  return ptcenter;
}
///////////////////////////
vector<float> util::getPtError(int threshold){
  vector<int> divpt;
  divpt.clear();
  int size;
  if(threshold==4){ 
    size = sizeof(divpt4)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt4[i]);
  }
  if(threshold==6){ 
    size = sizeof(divpt6)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt6[i]);
  }
  if(threshold==10){ 
    size = sizeof(divpt10)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt10[i]);
  }
  if(threshold==11){ 
    size = sizeof(divpt11)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt11[i]);
  }
  if(threshold==14){ 
    size = sizeof(divpt14)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt14[i]);
  }
  if(threshold==18){ 
    size = sizeof(divpt18)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt18[i]);
  }
  if(threshold==20){ 
    size = sizeof(divpt20)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt20[i]);
  }
  if(threshold==22){ 
    size = sizeof(divpt22)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt22[i]);
  }
  if(threshold==24){ 
    size = sizeof(divpt24)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt24[i]);
  }
  if(threshold==26){ 
    size = sizeof(divpt26)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt26[i]);
  }
  if(threshold==28){ 
    size = sizeof(divpt28)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt28[i]);
  }
  if(threshold==40){ 
    size = sizeof(divpt40)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt40[i]);
  }
  if(threshold==50){ 
    size = sizeof(divpt50)/sizeof(int);
    for (int i=0;i<size;i++) divpt.push_back(divpt50[i]);
  }
  vector<float> pterr;
  pterr.clear();
  float center,err;
  for (int j=0;j<size-1;j++){
    center = 0.5*(divpt[j]+divpt[j+1]);
    err = divpt[j+1]-center;
    pterr.push_back(err);
  }
  return pterr;
}
////////////////////////////////////
vector<int> util::getRunNumber(string period){
  vector<int> runNumber;
  runNumber.clear();
  if(period=="All" || period=="A" || period=="A1"){
    runNumber.push_back(296939);
    runNumber.push_back(296942);
  };
  if(period=="All" || period=="A" || period == "A2"){
    runNumber.push_back(297041);
    runNumber.push_back(297170);
    runNumber.push_back(297447);
  };
  if(period=="All" || period=="A" || period=="A3"){
    runNumber.push_back(297730);
  };
  if(period=="All" || period=="A" || period=="A4"){
    runNumber.push_back(298591);
  };
  if(period=="All" || period=="A" || period=="A5"){
    runNumber.push_back(298595);
    runNumber.push_back(298609);
    runNumber.push_back(298633);
    runNumber.push_back(298687);
    runNumber.push_back(298690);
    runNumber.push_back(298771);
    runNumber.push_back(298773);
    runNumber.push_back(298862);
    runNumber.push_back(298967);
  };
  if(period=="All" || period=="A" || period=="A6"){
    runNumber.push_back(299055);
  };
  if(period=="All" || period=="A" || period=="A7"){
    runNumber.push_back(299144);
    runNumber.push_back(299147);
    runNumber.push_back(299184);
  };
  if(period=="All" || period=="A" || period=="A8"){
    runNumber.push_back(299241);
    runNumber.push_back(299243);
    runNumber.push_back(299278);
    runNumber.push_back(299288);
    runNumber.push_back(299340);
  };
  if(period=="All" || period=="A" || period=="A9"){
    runNumber.push_back(299390);
    runNumber.push_back(300287);
    runNumber.push_back(299315);
  };
  if(period=="All" || period=="A" || period=="A10"){
    runNumber.push_back(299584);
    runNumber.push_back(300279);
  };
  if(period=="All" || period=="B" || period=="B1"){
    runNumber.push_back(300345);
    runNumber.push_back(300415);
  };
  if(period=="All" || period=="B" || period=="B2"){
    runNumber.push_back(300418);
    runNumber.push_back(300487);
    runNumber.push_back(300540);
    runNumber.push_back(300571);
    runNumber.push_back(300600);
    runNumber.push_back(300655);
    runNumber.push_back(300687);
  };
  if(period=="All" || period=="B" || period=="B3"){
    runNumber.push_back(300784);
    runNumber.push_back(300800);
    runNumber.push_back(300863);
    runNumber.push_back(300908);
  };
  if(period=="All" || period=="C" || period=="C1"){
    runNumber.push_back(301912);
    runNumber.push_back(301915);
    runNumber.push_back(301918);
  };
  if(period=="All" || period=="C" || period=="C2"){
    runNumber.push_back(301932);
    runNumber.push_back(301973);
  };
  if(period=="All" || period=="C" || period=="C3"){
    runNumber.push_back(302053);
    runNumber.push_back(302137);
    runNumber.push_back(302265);
    runNumber.push_back(302269);
    runNumber.push_back(302300);
    runNumber.push_back(302347);
    runNumber.push_back(302380);
    runNumber.push_back(302391);
  }
  if(period=="All" || period=="C" || period=="C4"){
    runNumber.push_back(302393);
  }
  if(period=="All" || period=="D" || period=="D1"){
    runNumber.push_back(302737);
  }
  if(period=="All" || period=="D" || period=="D1"){
    runNumber.push_back(302737);
  }
  return runNumber;
}
/////////////////////////////////
vector<string> util::getDataList(vector<int> runNumber, const char* input){
  vector<string> list;
  list.clear();
  ifstream finlist(input);
  stringstream sample;
  string file_rec;
  while(finlist>>file_rec){
    //cout << file_rec << endl;
    for (int run=0;run<runNumber.size();run++){
      sample.str("");
      sample << "user.mtanaka.data16_13TeV.00" << runNumber[run];
      if(strstr(file_rec.c_str(),sample.str().c_str()) == NULL) continue;
      list.push_back(file_rec);
    }
  }
  return list;
}
///////////////////////////////////
vector<string> util::getMCList( const char* input){
  vector<string> list;
  list.clear();
  ifstream finlist(input);
  stringstream sample;
  string file_rec;
  while(finlist>>file_rec)
    list.push_back(file_rec);
  return list;
}
////////////////////////////////
float util::calcDr(float eta1,float eta2,float phi1,float phi2){
  float deta = eta1-eta2;
  float dphi = acos(cos(phi1-phi2));
  float dr = sqrt(deta*deta+dphi*dphi);
  return dr;
}
///////////////////////////////
float util::getChamberCenterPhi(float phi){
  if(-PI/16 < phi && phi <= PI/16) return 0; 
  else if(PI/16<phi && phi<=PI/16*3) return PI/16*2;
  else if(PI/16*3<phi && phi<=PI/16*5) return PI/16*4;
  else if(PI/16*5<phi && phi<=PI/16*7) return PI/16*6;
  else if(PI/16*7<phi && phi<=PI/16*9) return PI/16*8;
  else if(PI/16*9<phi && phi<=PI/16*11) return PI/16*10;
  else if(PI/16*11<phi && phi<=PI/16*13) return PI/16*12;
  else if(PI/16*13<phi && phi<=PI/16*15) return PI/16*14;
  else if(-PI/16>phi && phi>=-PI/16*3) return -PI/16*2;
  else if(-PI/16*3>phi && phi>=-PI/16*5) return -PI/16*4;
  else if(-PI/16*5>phi && phi>=-PI/16*7) return -PI/16*6;
  else if(-PI/16*7>phi && phi>=-PI/16*9) return -PI/16*8;
  else if(-PI/16*9>phi && phi>=-PI/16*11) return -PI/16*10;
  else if(-PI/16*11>phi && phi>=-PI/16*13) return -PI/16*12;
  else if(-PI/16*13>phi && phi>=-PI/16*15) return -PI/16*14;
  else if(PI/16*15<phi && phi<=PI/16*16) return PI;
  else if(-PI/16*15>phi && phi>=-PI/16*16) return -PI;
  else return 9999;
}
/////////////////////////////
float util::cosAminusB(float phi1, float phi2){
  float ans = cos(phi1)*cos(phi2)+sin(phi1)*sin(phi2);
  return ans;
}

