void runCheck2muon(){
  gROOT->LoadMacro("macro/check2muon.cc++");
  check2muon t;
  t.Loop(); 
}
