#ifndef makeEfficiencyPlot_h
#define makeEfficiencyPlot_h
#include "macro/util.h"

class makeEffPlot {

  public:
    util m_util;

    makeEffPlot();
    ~makeEffPlot();
    void makePtEfficiencyPlot(const char *input_data, const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period); 
    void makeEtaEfficiencyPlot(const char *input_data, const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period); 
    void makePhiEfficiencyPlot(const char *input_data, const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period); 
    void makeEtaPhiEfficiencyPlot(const char *input_data, const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, string period); 
    void makeComPtEfficiencyPlot(const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, vector<string> com_per_list); 
    void makeComEtaEfficiencyPlot(const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, vector<string> com_per_list); 
    void makeComPhiEfficiencyPlot(const char *input_mc, vector<int> &chain_list, vector<string> &chain_name_list, vector<string> com_per_list); 
    void makeComPtVarConeEfficiencyPlot(const char *input_mc); 

};

#endif //
