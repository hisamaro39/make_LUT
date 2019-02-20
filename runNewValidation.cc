void runNewValidation(){
  gROOT->LoadMacro("macro/util.cxx++");
  gROOT->LoadMacro("macro/readLUT.cxx++");
  gROOT->LoadMacro("macro/newValidationT.cc++");
  newValidationT t;
  
  const char *input_data = "datalist/inputValidationData16.list";
  //const char *input_data = "datalist/input_mc_cheny.list";///mc test cheny
  //const char *output_data = "inputLUT/data16_13TeV_LS_sagitta.root";
  const char *output_data = "inputLUT/data16_13TeV_makeLUT_large_special_liner.root";
  //const char *output_data = "inputLUT/data16_13TeV_makeLUT_barrel_shift_v3.root";
  //const char *output_data = "inputLUT/data16_13TeV_makeLUT_large_special_inverse_invdefaultSP.root";
  //const char *output_data = "inputLUT/data16_13TeV_barrel_sagitta_large_special.root";
  //const char *output_data = "inputLUT/aho.root";
  t.Loop(input_data,output_data); 
  //t.makePtResidual(input_data,output_data_pt_residual); 
  
  const char *output_data_pt_residual = "inputLUT/data16_13TeV_pt_residual.root";
  const char *output_data_pt_residual_fitting = "inputLUT/data16_13TeV_pt_residual_fitting.root";
  //t.fitResidual(output_data_pt_residual,output_data_pt_residual_fitting); 
}
