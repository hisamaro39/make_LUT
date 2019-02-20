void runResolution(){
  gROOT->LoadMacro("macro/resolution.C+");
  resolution t;
  t.Loop(); 
}
