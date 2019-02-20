void runValidation(){
  gROOT->LoadMacro("macro/validationT.cc++");
  validationT t;
  t.Loop(); 
}
