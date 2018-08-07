Divide_pT(){
TFile *fLHC15o = new TFile("SumPbPb2015_MC_AOD_Allpt_Lpt.root");

TH1D *h1_c = (TH1D*)fLHC15o->Get("hAllpt_25y3_1M");
TH1D *h1_d = (TH1D*)fLHC15o->Get("hLpt_25y3_1M");

//TH1D *h1_c = (TH1D*)fLHC15o->Get("hAllpt_3y35_1M");
//TH1D *h1_d = (TH1D*)fLHC15o->Get("hLpt_3y35_1M");

//TH1D *h1_c = (TH1D*)fLHC15o->Get("hAllpt_35y4_1M");
//TH1D *h1_d = (TH1D*)fLHC15o->Get("hLpt_35y4_1M");

  TH1F *Divide_Pt = new TH1F("Divide_Pt", "Divide_Pt", 500, 0., 50.);
  Divide_Pt->Divide(h1_d,h1_c,1,1);

  for (Int_t ii=0;ii<h1_c->GetNbinsX();ii++){
    Double_t Nineff=h1_c->GetBinContent(ii+1)-h1_d->GetBinContent(ii+1);
//    printf("Nineff=%d\n",Nineff);
    Double_t ErrNineff = TMath::Sqrt(Nineff);
    if(ErrNineff==0){
    Divide_Pt->SetBinError(ii+1,0);}
    else{
    Divide_Pt->SetBinError(ii+1,ErrNineff/h1_c->GetBinContent(ii+1));
    }
  }

 TCanvas *c1 = new TCanvas("c1","divide_pT",10,10,600,1000);
  c1->SetLogy();
  Divide_Pt->Draw();

   TFile *fout = new TFile("MC_25y3_divide_Pt_new.root","recreate");
//   TFile *fout = new TFile("MC_3y35_divide_Pt_new.root","recreate");
//   TFile *fout = new TFile("MC_35y4_divide_Pt_new.root","recreate");
  Divide_Pt->Write();
   fout->Close();


}
