void runCheckEff(){
  gROOT->LoadMacro("macro/CheckEff.cc++");
  CheckEff t;
  t.Loop(); 
}
