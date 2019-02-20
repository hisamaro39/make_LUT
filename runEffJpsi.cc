void runEffJpsi(){
  gROOT->LoadMacro("macro/effJpsi.cc++");
  effJpsi t;
  t.Loop(); 
}
