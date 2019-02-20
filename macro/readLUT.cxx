#include "readLUT.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TMath.h"
using namespace std;

readLUT::readLUT(){}

readLUT::~readLUT(){}

void readLUT::scanLUT(){

  //////////LUT for correcting EE potision//////////////////////////////////////////
  string correctEELUT = "outputCorrectSagitta/parameter/fitCorrectSagitta3.txt";
  cout << "Reading EE alignment LUT..." << endl;
  cout << correctEELUT << endl;
  for (int i=0;i<8;i++)//eta
    for (int j=0;j<192;j++)//phi
      for (int k=0;k<2;k++)//charge
        for (int l=0;l<2;l++)//side
        SagittaPar[i][j][k][l]=0.; 

  ifstream ifs(correctEELUT.c_str());
  
  string line;
  while (!ifs.eof()) {
    getline(ifs, line);
    if (line.empty()) continue;
    int iEta, iPhi, iCharge, iSide;
    double sagittapar1, sagittapar2, sagittapar3;
    sscanf(line.c_str(), "%d %d %d %d %lf %lf %lf", &iEta, &iPhi, &iCharge, &iSide, &sagittapar1, &sagittapar2, &sagittapar3);
    SagittaPar[iEta-1][iPhi][iCharge][iSide] = sagittapar3;
  }

  //////////////Alignment Barrel LUT/////////////////////////////////////////////
  string barrelAlignmentLUT = "LUT_alignment/dZ_barrel.lut";
  cout << "Reading Barrel alignment LUT..." << endl;

  int saddress, innerR;
  int N0, N1, N2;
  double A0, A1, A2;
  std::ifstream file;

  for(int i_saddress=0; i_saddress<4; i_saddress++) {
    for(int i_innerR=0; i_innerR<2; i_innerR++) {
      NbinEta[i_saddress][i_innerR]=0;
      EtaMin[i_saddress][i_innerR]=0;
      EtaMax[i_saddress][i_innerR]=0;
      EtaStep[i_saddress][i_innerR]=0;
      NbinPhi[i_saddress][i_innerR]=0;
      PhiMin[i_saddress][i_innerR]=0;
      PhiMax[i_saddress][i_innerR]=0;
      PhiStep[i_saddress][i_innerR]=0;
      for(int i_eta=0; i_eta<15; i_eta++) 
	for(int i_phi=0; i_phi<30; i_phi++) 
	  for(int i_etaQ=0; i_etaQ<2; i_etaQ++) 
	    dZ[i_saddress][i_innerR][i_eta][i_phi][i_etaQ] = 0;
    }
  }
  
  file.open(barrelAlignmentLUT.c_str());
  if (!file) {
    cout << "Failed to open barrel alignment LUT file" << endl;
  }

  for(int i_lut=0; i_lut<2; i_lut++) {

    file >> saddress >> innerR;
    file >> EtaMin[saddress][innerR] >> EtaMax[saddress][innerR]
	 >> PhiMin[saddress][innerR] >> PhiMax[saddress][innerR]
	 >> NbinEta[saddress][innerR] >> NbinPhi[saddress][innerR];

    EtaStep[saddress][innerR] = (EtaMax[saddress][innerR] - EtaMin[saddress][innerR]) / (float)NbinEta[saddress][innerR];
    PhiStep[saddress][innerR] = (PhiMax[saddress][innerR] - PhiMin[saddress][innerR]) / (float)NbinPhi[saddress][innerR];

    for (int i_eta=0; i_eta<15; i_eta++) {
      for (int i_phi=0; i_phi<30; i_phi++) {
	for (int i_etaQ=0; i_etaQ<2; i_etaQ++) {

	  file >> N0 >> N1 >> N2 >> A0 >> A1 >> A2;

	  dZ[saddress][innerR][i_eta][i_phi][i_etaQ] = A0;

	} // etaQ loop
      } // phi loop
    } // eta loop
  } // nlut loop

  //////////////Alignment Barrel LUT data16/////////////////////////////////////////////
  string barrelAlignmentLUTnew = "LUT_alignment/dZ_barrel_data16_v2.lut";
  cout << "Reading Barrel alignment LUT made by data16..." << endl;

  int M0, M1, M2;
  double B0, B1, B2;
  std::ifstream file2;

  for(int i_charge=0;i_charge<2;i_charge++)
    for(int i_eta=0; i_eta<30; i_eta++) 
      for(int i_phi=0; i_phi<30; i_phi++) 
        dZnew[i_charge][i_eta][i_phi] = 0;
  
  file2.open(barrelAlignmentLUTnew.c_str());
  if (!file2) {
    cout << "Failed to open barrel alignment LUT file made by data16" << endl;
  }

  for(int i_charge=0;i_charge<2;i_charge++){
    for (int i_eta=0; i_eta<30; i_eta++) {
      for (int i_phi=0; i_phi<30; i_phi++) {
	  file2 >> M0 >> M1 >> M2 >> B0 >> B1 >> B2;
	  dZnew[i_charge][i_eta][i_phi] = B2;
          //cout << "charge/eta/phi/par=" << i_charge << "/" << i_eta << "/" << i_phi << "/" << B2 << endl;
      }
    }
  }

  cout << "Reading new Large Special LUT pol2 & sqrt..." << endl;
  string newLargeSpecialLUT = "/data/maxi174/zp/mtanaka/ReprocessLUT/LUT/pt_barrel_LS_v4.lut";
  cout << newLargeSpecialLUT << endl;
  double D0,D1,D2,D3,D4,D5;
  int qeta,spr;
  double etamin,etamax,phimin,phimax,temp1,temp2;
  
  for (int i_charge=0; i_charge<2; i_charge++) 
    for(int i_spr=0; i_spr<2; i_spr++) 
      for(int i_eta=0; i_eta<30; i_eta++) 
	for(int i_phi=0; i_phi<30; i_phi++) 
	  for(int i_parms=0; i_parms<5; i_parms++) 
	    BarrelLSnewPar[i_charge][i_spr][i_eta][i_phi][i_parms] = 0.;

  std::ifstream fileLSnew;
  fileLSnew.open(newLargeSpecialLUT.c_str());
  
  for(int nlut=0;nlut<4;++nlut) {
    fileLSnew >> qeta >> spr;
    fileLSnew >> etamin >> etamax >> phimin >> phimax
	 >> temp1 >> temp2;
    //cout << "qeta/spr/etamin/etamax/phimin/phimax/temp1/temp2=" 
      //<< qeta << "/" << spr << "/" << etamin << "/" << etamax << "/" << phimin << "/" << phimax << "/" << temp1 << "/" << temp2 << endl;
    
    if(qeta==0){
      for(int i=0;i<30;i++) {
        for(int j=0;j<30;j++) {
          fileLSnew >> D2 >> D1 >> D0;
          //cout << "par1/par2/par3=" << D2 << "/" << D1 << "/" << D0 << endl;
          BarrelLSnewPar[qeta][spr][i][j][0] = D2;
          BarrelLSnewPar[qeta][spr][i][j][1] = D1; 
          BarrelLSnewPar[qeta][spr][i][j][2] = D0; 
        }
      }
    }
    else{
      for(int i=0;i<30;i++) {
        for(int j=0;j<30;j++) {
          fileLSnew >> D4 >> D3 ;
          //cout << "par1/par2=" << D4 << "/" << D3 << "/" << endl;
          BarrelLSnewPar[qeta][spr][i][j][3] = D4;
          BarrelLSnewPar[qeta][spr][i][j][4] = D3;
        }
      }
    }

  }	
  fileLSnew.close();


}

/////////////////////////////////////////////
float readLUT::shiftEER(float InnerZ, float InnerR, float EEZ, float EER, float MiddleZ, float MiddleR, int iEta, int iPhi, int iCharge, int iSide){
  float par = SagittaPar[iEta][iPhi][iCharge][iSide];
  float sagitta = m_util.calcSagitta(InnerZ, InnerR, EEZ, EER, MiddleZ, MiddleR);
  float new_sagitta = sagitta - par;
  float new_EER = m_util.calcShiftEER(InnerZ, InnerR, EEZ, MiddleZ, MiddleR, new_sagitta);
  return new_EER;
}

//////////////////////////////////////////////////
double readLUT::GetDeltaZ(int saddress, double etaMap, double phiMap, double MFphi, float  sp1R) 
{
  if (saddress == 1) {
    int innerR = (sp1R > 6000)? 1: 0;
    std::pair<int, int> bins = GetBinNumber(saddress, innerR, etaMap, phiMap);
    int iEta = bins.first;
    int iPhi = bins.second;
    int iChamber=( fabs(MFphi)>90*TMath::DegToRad() )  ? 1:0;
    int iEta_inv=29-iEta;
    int iPhi_inv=29-iPhi;
    int iEta_bin=iEta;
    int iPhi_bin=iPhi;
    if( iEta<15 ){
      if( iChamber == 0 ){}
      else iPhi_bin=iPhi_inv;
    }
    else{
      if( iChamber == 0 ) iEta_bin=iEta_inv;
      else{
	iEta_bin=iEta_inv;
	iPhi_bin=iPhi_inv;
      }
    }
    int sign_etam=(iEta>14)?-1:1;
    int sign_etap=-sign_etam;
    return 10*(dZ[saddress][innerR][iEta_bin][iPhi_bin][0]*sign_etam
	    + dZ[saddress][innerR][iEta_bin][iPhi_bin][1]*sign_etap) / 2.;

  } 
  else return 0;
}

//////////////////////////////////////////////////
double readLUT::GetNewDeltaZ(int etabin, int phibin, int charge){ 
  //value of Large Special
  /*float EtaMin = -1.150;
  float PhiMin = -0.233;
  float EtaMax = 1.150;
  float PhiMax = 0.233;
  float EtaStep = (EtaMax-EtaMin)/30; 
  float PhiStep = (PhiMax-PhiMin)/30;

  int diveta_br = (int)((etaMap - EtaMin)/EtaStep);
  int phiBin_br = (int)((phiMap - PhiMin)/PhiStep);
  if(diveta_br<=-1) diveta_br = 0;
  if(diveta_br>=30) diveta_br = 29;
  if(phiBin_br<=-1) phiBin_br = 0;
  if(phiBin_br>=30) phiBin_br = 29;
*/
  int icharge = (charge>0)? 1 : 0;
  return dZnew[icharge][etabin][phibin];

}

//////////////////////////////////////////////////////
std::pair<int, int> readLUT :: GetBinNumber(int saddress, int innerR, double etaMap, double phiMap) 
{
  if(saddress > 3 || saddress < 0 || innerR < 0 || innerR > 1)  return std::make_pair(-5,-5);
  int etaBin = (int)((etaMap - EtaMin[saddress][innerR])/EtaStep[saddress][innerR]);
  int phiBin = (int)((phiMap - PhiMin[saddress][innerR])/PhiStep[saddress][innerR]);
  if(etaBin <= -1) etaBin = 0;
  if(etaBin >= NbinEta[saddress][innerR]) etaBin = NbinEta[saddress][innerR] - 1;
  if(phiBin <= -1) phiBin = 0;
  if(phiBin >= NbinPhi[saddress][innerR]) phiBin = NbinPhi[saddress][innerR] - 1;
  return std::make_pair(etaBin,phiBin);
}

////////////////////////////////////////////
float readLUT :: shiftRadiusLS(float etaMap, float phiMap, float spr, float radius, int charge){
  float EtaMinLS = -1.150;
  float PhiMinLS = -0.233;
  float EtaMaxLS = 1.150;
  float PhiMaxLS = 0.233;
  float EtaStepLS = (EtaMaxLS-EtaMinLS)/30; 
  float PhiStepLS = (PhiMaxLS-PhiMinLS)/30;
  int diveta_br = (int)((etaMap - EtaMinLS)/EtaStepLS);
  int phiBin_br = (int)((phiMap - PhiMinLS)/PhiStepLS);
  if(diveta_br<=-1) diveta_br = 0;
  if(diveta_br>=30) diveta_br = 29;
  if(phiBin_br<=-1) phiBin_br = 0;
  if(phiBin_br>=30) phiBin_br = 29;

  int qeta = (charge*etaMap>0)? 1 : 0;
  int ispr = (spr<6000)? 0 : 1;
  //cout << "etaMap/phiMap=" << etaMap << "/" << phiMap << endl;
  //cout << "qeta/spr/etabin/phibin=" << qeta << "/" << ispr << "/" << diveta_br << "/" << phiBin_br << endl;

  float A0 = BarrelLSnewPar[qeta][ispr][diveta_br][phiBin_br][0];
  float A1 = BarrelLSnewPar[qeta][ispr][diveta_br][phiBin_br][1];
  float A2 = BarrelLSnewPar[qeta][ispr][diveta_br][phiBin_br][2];
  float A3 = BarrelLSnewPar[qeta][ispr][diveta_br][phiBin_br][3];
  float A4 = BarrelLSnewPar[qeta][ispr][diveta_br][phiBin_br][4];
  //cout << "A0/A1/A2/A3/A4=" << A0 << "/" << A1 << "/" << A2 << "/" << A3 << "/" << A4<< endl;
  //cout << "radius=" << radius << endl;

  float new_radius=0.;
  if ( qeta==0 )  new_radius = radius + (A0/A1)*radius*radius;
  else new_radius = radius;

  return new_radius;

}
