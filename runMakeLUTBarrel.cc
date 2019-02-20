stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runMakeLUTBarrel(TString inputfilename="inputLUT/data16_13TeV_makeLUT_barrel_mm.root"){
  gROOT->LoadMacro("macro/fittingBarrel.cxx+");
  gROOT->LoadMacro("macro/LoopBarrel.C+");
  int region,phimax,phimin,fitting;
  TString valiable;

  const int bins = 25;
  double bin[bins+1]; 
  for (int i=0; i<bins+1; i++) bin[i]=20*i;

  //cout << "select valiables" << endl;
  //cin >> valiable;
  valiable = "Radius";
  phimax = 30;
  phimin = 0;
  int etamax = 30;
  stringstream iphi,ieta;
  outputtxtname << "LUTBarrelRadius/data16/fit3" << valiable << "_mm.txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for(int charge = 0 ; charge < 2;charge++){//charge
    for (int chamber=0; chamber<4; chamber++){//chamber
      for(int i = 0 ; i < 30;i++){//eta from 0 to 29
        for(int j = 0 ; j < 30;j++){//phi from 0 to 29
          if (i < 10) ieta << "0" << i;
          else ieta << i;
          if (j < 10) iphi << "0" << j;
          else iphi << j;
          mode << valiable << charge << chamber << ieta.str().c_str() << iphi.str().c_str();
          outputfilename << "outputMakeLUTBarrelRadius/data16/fit" << mode.str().c_str() << "_mm.root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopBarrel(bins, bin, path, 0, 500000, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),charge,chamber,i,j);
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
