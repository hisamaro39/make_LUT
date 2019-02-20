stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runCorrectBarrelSagitta(TString inputfilename="inputLUT/data16_13TeV_barrel_sagitta.root"){
  gROOT->LoadMacro("macro/fittingCorrectBarrelSagitta.cxx+");
  gROOT->LoadMacro("macro/LoopCorrectBarrelSagitta.C+");
  
  const int bins = 40;
  double bin[bins+1];
  for(int ibin=0;ibin<bins+1;ibin++) {
    bin[ibin] = 10*ibin;
    //cout << bin[ibin] << endl;
  }

  stringstream iphi,ieta;
  outputtxtname << "outputCorrectBarrelSagitta/parameter/fitCorrectBarrelSagitta.txt";
  //outputtxtname << "aho.txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for(int charge = 0 ; charge < 2;charge++){//charge
    for (int chamber=1; chamber<2; chamber++){//chamber
      for(int i = 0 ; i < 30;i++){//eta from 0 to 29
        for(int j = 0 ; j < 30;j++){//phi from 0 to 29
          if (i < 10) ieta << "0" << i;
          else ieta << i;
          if (j < 10) iphi << "0" << j;
          else iphi << j;
          mode << "BarrelSagitta" << charge << chamber << ieta.str().c_str() << iphi.str().c_str();
          cout << "mode=" << mode.str() << endl;
          outputfilename << "outputCorrectBarrelSagitta/fitResult/fit" << mode.str().c_str() << ".root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopCorrectBarrelSagitta(bins, bin, path, 0, 40, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),charge,chamber,i,j);
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
