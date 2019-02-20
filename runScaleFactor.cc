void runScaleFactor(){
  gROOT->LoadMacro("macro/util.cxx++");
  gROOT->LoadMacro("macro/makeEfficiencyPlot.cc++");
  gROOT->LoadMacro("macro/ScaleFactor.cc++");
  ScaleFactor m_ScaleFactor;
  makeEffPlot m_makeEffPlot;
  util m_util;
  const char *input_mc = "datalist/inputScaleFactorMCZmumu_v8.list" ;
  const char *output_mc = "outputScaleFactor/output/mc15_13TeV_Zmumu.root" ;
  vector<string> input_list_mc;
  input_list_mc.clear();
  input_list_mc = m_util.getMCList(input_mc);
  //m_ScaleFactor.Loop(input_list_mc,output_mc); 
  const char *output_mc_resolution = "outputScaleFactor/output/mc15_13TeV_Zmumu_resolution.root" ;
  const char *output_data_resolution = "outputScaleFactor/output/data16_13TeV_Zmumu_resolution.root" ;
  const char *output_fit = "outputScaleFactor/output/mc15_13TeV_Zmumu_fit_resolution.root" ;
  //m_ScaleFactor.makeResolutionPlot(input_list_mc,output_mc_resolution); 
  //m_ScaleFactor.fitResolution(output_mc_resolution,output_fit); 

  const char *input_data = "datalist/inputScaleFactor_v5.list" ;
  //const char *input_data = "datalist/inputScaleFactor_jpsi_test.list" ;
  vector<int> chain_list;
  vector<string> chain_name_list;
  m_ScaleFactor.getThreshold(input_data,chain_list,chain_name_list);
  
  //string period_list[] = {"A","A1","A2","A3","A4","A5","A6","A7","A8","A9","A10","B","B1","B2","B3","C","C1","C2","C3"};
  string period_list[] = {"All"};
  int per_size = sizeof(period_list)/sizeof(string);
  for (int per=0;per<per_size;per++){
    string period = period_list[per];
    cout << "Period" << period_list[per] << endl;
    vector<int> runNumber;
    runNumber.clear();
    runNumber = m_util.getRunNumber(period);

    stringstream output_data;
    output_data << "outputScaleFactor/output_mishitsu/data16_13TeV_period_" << period << ".root";
    //output_data << "outputScaleFactor/output/data15_13TeV_jpsi_test.root";
    vector<string> input_list_data;
    input_list_data.clear();
    input_list_data = m_util.getDataList(runNumber,input_data);
    if(input_list_data.size()==0) {
      cout << "no data sets found in period " << period_list[per] << endl;
      continue;
    }
    m_ScaleFactor.Loop(input_list_data,output_data.str().c_str()); 
    //m_makeEffPlot.makePtEfficiencyPlot(output_data.str().c_str(),output_mc,chain_list,chain_name_list,period);
    //m_makeEffPlot.makeEtaEfficiencyPlot(output_data.str().c_str(),output_mc,chain_list,chain_name_list,period);
    //m_makeEffPlot.makePhiEfficiencyPlot(output_data.str().c_str(),output_mc,chain_list,chain_name_list,period);
    //m_makeEffPlot.makeEtaPhiEfficiencyPlot(output_data.str().c_str(),output_mc,chain_list,chain_name_list,period);
    //m_ScaleFactor.makeResolutionPlot(input_list_data,output_data_resolution); 
    output_data.str("");
  }
  string com_period_list[] = {"A","B"};
  int com_per_size = sizeof(com_period_list)/sizeof(string);
  vector<string> vec_per_list;
  vec_per_list.clear();
  for (int j=0;j<com_per_size;j++) vec_per_list.push_back(com_period_list[j]);
  //m_makeEffPlot.makeComPtEfficiencyPlot(output_mc,chain_list,chain_name_list,vec_per_list);
  //m_makeEffPlot.makeComEtaEfficiencyPlot(output_mc,chain_list,chain_name_list,vec_per_list);
  //m_makeEffPlot.makeComPhiEfficiencyPlot(output_mc,chain_list,chain_name_list,vec_per_list);

}
