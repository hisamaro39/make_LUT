void runScaleFactorTrack(){
  gROOT->LoadMacro("macro/util.cxx++");
  gROOT->LoadMacro("macro/makeEfficiencyPlot.cc++");
  gROOT->LoadMacro("macro/ScaleFactorTrack.cc++");
  ScaleFactorTrack m_ScaleFactor;
  makeEffPlot m_makeEffPlot;
  util m_util;

  const char *input_data = "datalist/inputScaleFactorTrackPeriodB.list" ;
  const char *output_data = "outputScaleFactor/output_track/data16_13TeV_periodB.root";

  //m_ScaleFactor.Loop(input_data,output_data); 
  m_makeEffPlot.makeComPtVarConeEfficiencyPlot(output_data);
}
