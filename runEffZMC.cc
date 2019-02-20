void runEffZMC(){
  gROOT->LoadMacro("macro/effZMC.cc++");
  effZMC t;
  t.Loop(); 
}
