void runRerun(){
  gROOT->LoadMacro("macro/Rerun.C+");
  Rerun t;
  t.Loop(); 
}
