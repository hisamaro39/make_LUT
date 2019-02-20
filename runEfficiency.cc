void runEfficiency(){
  gROOT->LoadMacro("macro/efficiency.C+");
  efficiency t;
  t.Loop(); 
}
