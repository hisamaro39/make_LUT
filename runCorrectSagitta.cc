stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runCorrectSagitta(TString inputfilename="inputLUT/data15_13TeV_shift_EER.root"){
  gROOT->LoadMacro("macro/fittingCorrectSagitta.cxx+");
  gROOT->LoadMacro("macro/LoopCorrectSagitta.C+");
  
  const int bins = 40;
  double bin[bins+1];
  for(int ibin=0;ibin<bins+1;ibin++) {
    bin[ibin] = 10*ibin;
    //cout << bin[ibin] << endl;
  }

  stringstream iphi,ieta;
  outputtxtname << "outputCorrectSagitta/parameter/fitCorrectSagitta3.txt";
  //outputtxtname << "aho.txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for (int l=0; l<2; l++){//side 0-2
    for(int k = 0 ; k < 2;k++){//charge 0-2
      for(int j = 2 ; j <=7;j++){//eta from 2-7
        for(int i = 0 ; i < 192;i++){//phi from 0-192
          if (i < 10) iphi << "0" << i;
          else iphi << i;
          if (j < 10) ieta << "0" << j;
          else ieta << j;
          mode << "EndcapSagittaAnormal" << iphi.str().c_str() << ieta.str().c_str() << k << l;
          cout << "mode=" << mode.str() << endl;
          outputfilename << "outputCorrectSagitta/fitResult3/fit" << mode.str().c_str() << ".root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopCorrectSagitta(bins, bin, path, 1, 80, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),i,j,k,l);
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
