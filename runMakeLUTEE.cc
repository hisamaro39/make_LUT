stringstream  mode,name;
stringstream outputfilename;
stringstream outputtxtname;

void runMakeLUTEE(TString inputfilename="inputLUT/data15_13TeV_final_again.root"){
  gROOT->LoadMacro("macro/fittingEE.cxx+");
  gROOT->LoadMacro("macro/LoopEE.C+");
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
  valiable = "EndcapRadius";
  phimax = 24;
  phimin = 0;
  int etamax = 30;
  stringstream iphi,ieta;
  outputtxtname << "LUTEndcapRadius/final_again/fit" << valiable << ".txt";
  ofstream fout_txt(outputtxtname.str().c_str());

  for(int k = 0 ; k < 2;k++){//Qeta
    for (int l=0; l<2; l++){//small large
      for(int j = 1 ; j <=8;j++){//eta from 1 to 8
        for(int i = 0 ; i < 24;i++){//phi from 0 to 23
          if (i < 10) iphi << "0" << i;
          else iphi << i;
          if (j < 10) ieta << "0" << j;
          else ieta << j;
          mode << valiable << iphi.str().c_str() << ieta.str().c_str() << k << l;
          cout << "mode=" << mode.str() << endl;
          outputfilename << "outputMakeLUTEE/final_again/fit" << mode.str().c_str() << ".root";
          TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
          TString path=fout->GetPath();  
          LoopEE(bins, bin, path, 0, 0.25, inputfilename, outputtxtname.str().c_str(),"",mode.str().c_str(),i,j,k,l);
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
  //}
  /*  else if(fitting == 2){
      cout << "input valiable" << endl;
      cin >> valiable;
      name << valiable;
      outputfilename << "output/fit" << name.str().c_str() << ".root";
      outputtxtname << "output/fit" << name.str().c_str() << ".txt";
      TFile* fout=new TFile(outputfilename.str().c_str(), "recreate");
      ofstream fout_txt(outputtxtname.str().c_str());
      TString path=fout->GetPath();
      Loop(bins, bin, path, 0, 1.0, inputfilename, outputtxtname.str().c_str(),"",name.str().c_str(),0,0,0);
      fout->Write();
      }
      */
}
