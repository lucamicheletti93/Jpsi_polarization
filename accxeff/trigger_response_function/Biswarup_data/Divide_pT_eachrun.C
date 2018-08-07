void Divide_pT_eachrun(){
TFile *fLHC15o = new TFile("PbPb_PS_246991.root");
TList *chist0 = (TList*)fLHC15o->Get("chist0");
TH1D *h1_c = (TH1D*)chist0->FindObject("hAllpt_1M");

TH1D *h1_d = (TH1D*)chist0->FindObject("hLpt_1M");
h1_d->Draw();

TH1F *Divide_Pt = new TH1F("Divide_Pt", "Divide_Pt", 500, 0., 50.);
Divide_Pt->Sumw2();  
Divide_Pt->Divide(h1_d,h1_c,1,1);
  
 TCanvas *c1 = new TCanvas("c1","divide_pT",10,10,600,1000);
  c1->SetLogy();
  Divide_Pt->Draw();

   TFile *fout = new TFile("Data_divide_Pt.root","recreate");
  Divide_Pt->Write();
   fout->Close();

}
