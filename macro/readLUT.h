#ifndef readLUT_h
#define readLUT_h
#include "macro/util.h"

class readLUT {
  public:
    
    readLUT();
    ~readLUT();

    util m_util;

    void scanLUT();
    double SagittaPar[8][192][2][2];
    float shiftEER(float InnerZ, float InnerR, float EEZ, float EER, float MiddleZ, float MiddleR, int iEta, int iPhi, int iCharge, int iSide);
    double GetDeltaZ(int saddress, double etaMap, double phiMap, double MFphi, float  sp1R);
    double GetNewDeltaZ(int etabin, int phibin, int charge); 
    std::pair<int, int> GetBinNumber(int saddress, int innerR, double etaMap, double phiMap); 
    float shiftRadiusLS(float etaMap, float phiMap, float spr, float radius, int charge);

    double dZ[4][2][15][30][2]; // [s_address][innerR][eta][phi][etaQ]
    double dZnew[2][30][30]; // [charge][eta][phi]
    int NbinEta[4][2]; // [s_address][innerR]
    float EtaMin[4][2];
    float EtaMax[4][2];
    float EtaStep[4][2];
    int NbinPhi[4][2];
    float PhiMin[4][2];
    float PhiMax[4][2];
    float PhiStep[4][2];
    double BarrelLSnewPar[4][2][30][30][2];

};
#endif //
