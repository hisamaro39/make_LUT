void runEffZ(){
  gROOT->LoadMacro("macro/effZ.cc++");
  effZ t;
  t.Loop(); 
}
