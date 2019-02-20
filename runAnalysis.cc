void runAnalysis(){
  gROOT->LoadMacro("macro/analysis.cc++");
  analysis t;
  t.Loop(); 
}
