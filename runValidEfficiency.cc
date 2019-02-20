void runValidEfficiency(){
  gROOT->LoadMacro("macro/validEfficiency.C+");
  validEfficiency t;
  t.Loop(); 
}
