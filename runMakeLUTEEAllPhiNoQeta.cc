stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runMakeLUTEEAllPhiNoQeta(TString inputfilename="inputLUT/data15_13TeV_version9_separate.root"){
  gROOT->LoadMacro("macro/fittingEEAllPhiNoQeta.cxx+");
  gROOT->LoadMacro("macro/LoopEEAllPhiNoQeta.C+");
  int region,fitting;
  TString valiable;

  const int bins = 25;
  double bin[bins+1] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250};

  valiable = "EndcapRadiusAllPhiNoQeta";
  stringstream iphi,ieta;
  outputtxtname << "LUTEndcapRadiusAllPhiNoQeta/version9_final/fit" << valiable << ".txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for (int l=0; l<2; l++){//side
    for(int k = 0 ; k < 2;k++){//charge
      for(int j = 1 ; j <=8;j++){//eta from 1 to 8
        for(int i = 0 ; i < 192;i++){//phi from 0 to 191
          if (i < 10) iphi << "0" << i;
          else iphi << i;
          if (j < 10) ieta << "0" << j;
          else ieta << j;
          mode << valiable << iphi.str().c_str() << ieta.str().c_str() << k << l;
          cout << "mode=" << mode.str() << endl;
          outputfilename << "outputMakeLUTEEAllPhiNoQeta/version9_final/fit" << mode.str().c_str() << ".root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopEEAllPhiNoQeta(bins, bin, path, 0, 0.25, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),i,j,k,l);
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
