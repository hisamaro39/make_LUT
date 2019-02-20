stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runMakeLUTEENoSL(TString inputfilename="inputLUT/data15_13TeV_version7.root"){
  gROOT->LoadMacro("macro/fittingEENoSL.cxx+");
  gROOT->LoadMacro("macro/LoopEENoSL.C+");
  int region,phimax,phimin,fitting;
  TString valiable;

  //for calculate parameter

  const int bins = 25;
  double bin[bins+1] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250};
  //const int bins = 12;
  //double bin[bins+1] = {0,20,40,60,80,100,120,140,160,180,200,220,240};
  //cout << "choose fitting mode" << endl;
  //cout << "1.ALL regions 2.each region" << endl;
  //cin >> fitting;

  //if(fitting == 1){
  //cout << "select valiables" << endl;
  //cin >> valiable;
  valiable = "EndcapRadiusNoSL";
  phimax = 12;
  phimin = 0;
  int etamax = 30;
  stringstream iphi,ieta;
  outputtxtname << "LUTEndcapRadiusNoSL/version7/fit" << valiable << ".txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for(int k = 0 ; k < 2;k++){//Qeta
    for(int j = 1 ; j <=8;j++){//eta from 1 to 8
      for(int i = 0 ; i < 12;i++){//phi from 0 to 11
        if (i < 10) iphi << "0" << i;
        else iphi << i;
        if (j < 10) ieta << "0" << j;
        else ieta << j;
        mode << valiable << iphi.str().c_str() << ieta.str().c_str() << k ;
        cout << "mode=" << mode.str() << endl;
        outputfilename << "outputMakeLUTEENoSL/version7/fit" << mode.str().c_str() << ".root";
        TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
        TString path=fout->GetPath();  
        LoopEENoSL(bins, bin, path, 0, 0.25, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),i,j,k);
        fout->Write();
        fout->Close();
        delete fout;
        mode.str("");
        outputfilename.str("");
        iphi.str("");
        ieta.str("");
      }
    }
  }
  outputtxtname.str("");
}
