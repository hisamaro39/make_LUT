stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runMakeLUTLargeSpecialInv(TString inputfilename="inputLUT/data16_13TeV_makeLUT_large_special_inverse_invdefaultSP.root"){
  gROOT->LoadMacro("macro/fittingLargeSpecialInv.cxx+");
  gROOT->LoadMacro("macro/LoopLargeSpecialInv.C+");
  int region,phimax,phimin,fitting;
  TString valiable;

  const int bins = 25;
  double bin[bins+1] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250};

  //cout << "select valiables" << endl;
  //cin >> valiable;
  valiable = "RadiusLS";
  phimax = 30;
  phimin = 0;
  int etamax = 30;
  stringstream iphi,ieta;
  outputtxtname << "LUTBarrelRadius/data16_LS/fit" << valiable << "_inv_v2.txt";
  //outputtxtname << "LUTBarrelRadius/data16_LS/fit" << valiable << "_test.txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for(int qeta = 0 ; qeta < 2;qeta++){//qeta
    for (int phi=0; phi<2; phi++){//phi
      for(int i = 0 ; i < 30;i++){//eta from 0 to 29
        for(int j = 0 ; j < 30;j++){//phi from 0 to 29
          if (i < 10) ieta << "0" << i;
          else ieta << i;
          if (j < 10) iphi << "0" << j;
          else iphi << j;
          mode << valiable << qeta << phi << ieta.str().c_str() << iphi.str().c_str();
          outputfilename << "outputMakeLUTBarrelRadius/data16_LS/fit" << mode.str().c_str() << "_inv.root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopLargeSpecialInv(bins, bin, path, 0, 0.25, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),qeta,phi,i,j);
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
