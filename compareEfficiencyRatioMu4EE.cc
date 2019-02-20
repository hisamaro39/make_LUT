void compareEfficiencyRatioMu4EE(){

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetPadLeftMargin(0.12);
  //gStyle->SetPadBottomMargin(0.14);
  gStyle->SetOptStat(0);
  
  TFile *file1 = TFile::Open("input/efficiency_jpsi_version9_ee_distance_all.root");
  TLatex test;

  TCanvas *c1 = new TCanvas("c1","",800,1000);
  TPad* p1 = new TPad("main","main",0.0,0.2,1.0,1.0,10,0,0);
  p1->Draw();
  p1->SetNumber(1);
  p1->SetBottomMargin(0.15);
  TPad* p2 = new TPad("ratio","ratio",0.0,0.0,1.0,0.2,10,0,0);
  //p2->SetBottomMargin(0.2);
  p2->Draw();
  p2->SetNumber(2);

  c1->cd(1);
  EfficiencyPtSAMu4EndcapDefault->SetTitle("");
  EfficiencyPtSAMu4EndcapDefault->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPtSAMu4EndcapDefault->GetXaxis()->SetTitle("p_{T,off}(GeV)");
  EfficiencyPtSAMu4EndcapDefault->GetXaxis()->SetTitleSize(0.05);
  EfficiencyPtSAMu4EndcapDefault->GetYaxis()->SetTitleSize(0.05);
  EfficiencyPtSAMu4EndcapDefault->GetYaxis()->SetTitleOffset(1.1);
  EfficiencyPtSAMu4EndcapDefault->GetXaxis()->SetLabelSize(0.05);
  EfficiencyPtSAMu4EndcapDefault->GetYaxis()->SetLabelSize(0.05);
  EfficiencyPtSAMu4EndcapDefault->GetXaxis()->SetRangeUser(0,15);
  EfficiencyPtSAMu4EndcapDefault->SetLineColor(2);
  EfficiencyPtSAMu4EndcapDefault->Draw("ap");
  EfficiencyPtSAMu4EndcapEE->SetLineColor(4);
  EfficiencyPtSAMu4EndcapEE->Draw("p");
  test.SetTextSize(0.03);
  test.DrawLatex(8,0.6,"ATLAS work in progress");
  test.DrawLatex(8,0.5,"#sqrt{s}=13TeV, #intLdt=3.26fb^{-1}");
  test.DrawLatex(8,0.4,"J/#psi Tag&Probe |#eta|>1.0");
  test.SetTextSize(0.04);
  test.SetTextColor(2);
  test.DrawLatex(8,0.85,"#alpha,#beta");
  test.SetTextColor(4);
  test.DrawLatex(8,0.77,"#alpha,#beta,EE Radiud");

  Double_t *effDefault = EfficiencyPtSAMu4EndcapDefault->GetY();
  Double_t *effEE = EfficiencyPtSAMu4EndcapEE->GetY();
  Double_t *efferrDefault = EfficiencyPtSAMu4EndcapDefault->GetEY();
  Double_t *efferrEE = EfficiencyPtSAMu4EndcapEE->GetEY();
  Double_t *ptbin = EfficiencyPtSAMu4EndcapDefault->GetX();
  Double_t *ptbinerr = EfficiencyPtSAMu4EndcapDefault->GetEX();
  float Ratio[60],RatioErr[60], x[60], xerr[60];
  for (int i=0;i<60;i++){
    float eff1 = effDefault[i];
    float eff2 = effEE[i];
    float efferr1 = efferrDefault[i];
    float efferr2 = efferrEE[i];
    float ratio=0, ratioerr=0;
    x[i]=ptbin[i];
    xerr[i]=ptbinerr[i];
    cout << "pt=" << x[i] << " +- " << xerr[i] << endl; 
    cout << "effDefault=" << eff1 << " +- " << efferr1 << endl;
    cout << "effEE=" << eff2 << " +- " << efferr2 << endl;
    if (effDefault[i]>1e-5){
      ratio = eff2/eff1;
      ratioerr = sqrt(eff1*eff1*efferr2*efferr2 + eff2*eff2*efferr1*efferr1)/(eff1*eff1);
    }
    cout << "ratio=" << ratio << " +- " << ratioerr << endl;
    Ratio[i]=ratio;
    RatioErr[i]=ratioerr;
  }
  
  TGraphErrors *g_ratio = new TGraphErrors(60,x,Ratio,xerr,RatioErr);
  c1->cd(2);
  g_ratio->SetTitle("");
  g_ratio->GetYaxis()->SetTitle("Ratio EE/Default");
  g_ratio->GetYaxis()->SetTitleSize(0.1);
  g_ratio->GetYaxis()->SetTitleOffset(0.5);
  g_ratio->GetYaxis()->SetLabelSize(0.07);
  g_ratio->GetYaxis()->SetRangeUser(0.9,1.1);
  g_ratio->GetXaxis()->SetRangeUser(0,15);
  g_ratio->GetXaxis()->SetLabelSize(0);
  g_ratio->Draw("ap");
  TLine *line = new TLine(0,1,15,1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw();

  c1->SaveAs("output/compareEfficiencyPtSAMu4EE.eps");

}
