void runEfficiencyMC(){
  gROOT->LoadMacro("macro/efficiency_mc.C+");
  efficiency_mc t;
  t.Loop(); 
}
