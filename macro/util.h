#ifndef util_h
#define util_h

class util {
  public:
    util();
    ~util();
    float calcCosDphi(float xi,float yi,float xm,float ym);
    float calcSlope(float eta);
    float calcDis(float r1,float z1,float r2,float z2);
    float calcDisR(float segr,float segz,float spr,float spz);
    float calcRfromTrackZ(float slope, float z);
    float calcZfromTrackR(float slope, float r);
    bool isBadPhi(int phibinall);
    bool isSmall(float eez);
    pair<int,int> GetBinNumber(float m_tgcMid1_phi, float m_tgcMid1_eta);
    pair<float,float> calcCenter(float x1,float y1,float x2,float y2,float x3,float y3);
    float calcDistance(float x1,float y1,float x2,float y2,float x3,float y3);
    int GetBinNumberEE(float phi,int sl);
    int GetBinNumberEE2(float phi);
    int GetBinNumberAllPhi(float phi);
    float calcAlpha(float r1,float z1,float r2,float z2);
    float calcIntercept(float r1,float z1,float r2,float z2);
    double calcSagitta(double InnerZ, double InnerR,
        double EEZ, double EER,
        double MiddleZ, double MiddleR);
    double computeRadius3Points(double InnerZ, double InnerR, 
        double EEZ, double EER,
        double MiddleZ, double MiddleR);
    double calcShiftEER(double InnerZ, double InnerR, double EEZ, double MiddleZ, double MiddleR, double sagitta);
    int returnPlateau(int threshold, int barrel);
    vector<int> devidePt(int threshold);
    vector<float> getPtCenter(int threshold);
    vector<float> getPtError(int threshold);
    vector<int> getRunNumber(string period);
    vector<string> getDataList(vector<int> runNumber, const char* input);
    vector<string> getMCList( const char* input);
    float calcDr(float eta1,float eta2,float phi1,float phi2);
    float getChamberCenterPhi(float phi);
    float cosAminusB(float phi1, float phi2);

};

#endif //
