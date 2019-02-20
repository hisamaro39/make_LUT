stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runMakeLUTLargeSpecial(TString inputfilename="inputLUT/data16_13TeV_makeLUT_large_special_shift.root"){
  gROOT->LoadMacro("macro/fittingLargeSpecial.cxx+");
  gROOT->LoadMacro("macro/LoopLargeSpecial.C+");
  int region,phimax,phimin,fitting;
  TString valiable;

  double bin[] = {0,20,40,60,80,100,120,140,160,180,200,300,400,500,600,700,800,900,1000}; 
  int bins = sizeof(bin)/sizeof(double)-1;

  //cout << "select valiables" << endl;
  //cin >> valiable;
  valiable = "RadiusLS";
  phimax = 30;
  phimin = 0;
  int etamax = 30;
  stringstream iphi,ieta;
  outputtxtname << "LUTBarrelRadius/data16_LS/fit" << valiable << "_liner.txt";
  //outputtxtname << "LUTBarrelRadius/data16_LS/fit" << valiable << "_test.txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for(int qeta = 1 ; qeta < 2;qeta++){//qeta
    for (int phi=0; phi<1; phi++){//phi
      for(int i = 18 ; i < 19;i++){//eta from 0 to 29
        for(int j = 16 ; j < 17;j++){//phi from 0 to 29
          if (i < 10) ieta << "0" << i;
          else ieta << i;
          if (j < 10) iphi << "0" << j;
          else iphi << j;
          mode << valiable << qeta << phi << ieta.str().c_str() << iphi.str().c_str();
          outputfilename << "outputMakeLUTBarrelRadius/data16_LS/fit" << mode.str().c_str() << "_liner.root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopLargeSpecial(bins, bin, path, 0, 1000000, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),qeta,phi,i,j);
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
  }
  outputtxtname.str("");
}
