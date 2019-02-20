void runCheckLUT(){
  gROOT->LoadMacro("macro/checkLUT.C+");
  checkLUT t;
  t.Loop(); 
}
