Divide_pT(){
TFile *fLHC15o = new TFile("Sum_LHC15o_PS.root");

TH1D *h0Apt = (TH1D*)fLHC15o->Get("hAllpt_1M");
TH1D *h0Lpt = (TH1D*)fLHC15o->Get("hLpt_1M");

TH1D *h1Apt = (TH1D*)fLHC15o->Get("hAllpt_25y3_1M");
TH1D *h1Lpt = (TH1D*)fLHC15o->Get("hLpt_25y3_1M");

TH1D *h2Apt = (TH1D*)fLHC15o->Get("hAllpt_3y35_1M");
TH1D *h2Lpt = (TH1D*)fLHC15o->Get("hLpt_3y35_1M");

TH1D *h3Apt = (TH1D*)fLHC15o->Get("hAllpt_35y4_1M");
TH1D *h3Lpt = (TH1D*)fLHC15o->Get("hLpt_35y4_1M");


  TH1F *hData25y4 = new TH1F("hData25y4", "hData25y4", 500, 0., 50.);
  hData25y4->Divide(h0Lpt,h0Apt,1,1);

  TH1F *hData25y3 = new TH1F("hData25y3", "hData25y3", 500, 0., 50.);
  hData25y3->Divide(h1Lpt,h1Apt,1,1);

  TH1F *hData3y35 = new TH1F("hData3y35", "hData3y35", 500, 0., 50.);
  hData3y35->Divide(h2Lpt,h2Apt,1,1);

  TH1F *hData35y4 = new TH1F("hData35y4", "hData35y4", 500, 0., 50.);
  hData35y4->Divide(h3Lpt,h3Apt,1,1);

  for (Int_t ii=0;ii<h0Apt->GetNbinsX();ii++){
    Double_t Nineff=h0Apt->GetBinContent(ii+1)-h0Lpt->GetBinContent(ii+1);
    Double_t Nineff1=h1Apt->GetBinContent(ii+1)-h1Lpt->GetBinContent(ii+1);
    Double_t Nineff2=h2Apt->GetBinContent(ii+1)-h2Lpt->GetBinContent(ii+1);
    Double_t Nineff3=h3Apt->GetBinContent(ii+1)-h3Lpt->GetBinContent(ii+1);
//    printf("Nineff=%d\n",Nineff);
    Double_t ErrNineff = TMath::Sqrt(Nineff);
    Double_t ErrNineff1 = TMath::Sqrt(Nineff1);
    Double_t ErrNineff2 = TMath::Sqrt(Nineff2);
    Double_t ErrNineff3 = TMath::Sqrt(Nineff3);
    if(ErrNineff==0 && ErrNineff1==0 && ErrNineff2==0 && ErrNineff3==0){
    hData25y4->SetBinError(ii+1,0);
    hData25y3->SetBinError(ii+1,0);
    hData3y35->SetBinError(ii+1,0);
    hData35y4->SetBinError(ii+1,0);
       }
    else{
    hData25y4->SetBinError(ii+1,ErrNineff/h0Apt->GetBinContent(ii+1));
    hData25y3->SetBinError(ii+1,ErrNineff1/h1Apt->GetBinContent(ii+1));
    hData3y35->SetBinError(ii+1,ErrNineff2/h2Apt->GetBinContent(ii+1));
    hData35y4->SetBinError(ii+1,ErrNineff3/h3Apt->GetBinContent(ii+1));
    }
  }

 TCanvas *c1 = new TCanvas("c1","divide_pT",10,10,600,1000);
  c1->SetLogy();
  hData25y4->Draw();

   TFile *fout = new TFile("Data_divide_Pt.root","recreate");
   hData25y4->Write();
   hData25y3->Write();
   hData3y35->Write();
   hData35y4->Write();
   fout->Close();


}
