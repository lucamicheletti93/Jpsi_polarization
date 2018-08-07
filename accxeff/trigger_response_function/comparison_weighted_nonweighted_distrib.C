#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

#include <TROOT.h>
#include <TMinuit.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TPaveText.h>
#include <TGaxis.h>

#include "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.h"
#endif

void comparison_weighted_nonweighted_distrib(int pt_min, int pt_max){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  double PI = TMath::Pi();

  ostringstream convert_pt_min;
  convert_pt_min << pt_min;
  string str_pt_min =  convert_pt_min.str();

  ostringstream convert_pt_max;
  convert_pt_max << pt_max;
  string str_pt_max =  convert_pt_max.str();

  string dataset = str_pt_min + "pt" + str_pt_max;
  string nameOption = "_test";

  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + nameOption + ".root";
  string fileInputName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/HistosFromOfficialTree_Jpsi_PbPb_Nopol_TH3rec.root";
  string fileInputWeightedName = "Histos_Jpsi_Rec_Weigthed.root";
  string fileAccxEffName = "accxeff_" + dataset + nameOption + ".root";

  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  TFile *fileInput = new TFile(fileInputName.c_str(),"READ");
  TFile *fileInputWeighted = new TFile(fileInputWeightedName.c_str(),"READ");

  //============================================================================
  printf("---> Reading the binning file ... \n");
  //============================================================================
  vector <double> CostValues;
  vector <double> PhiValues;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  //============================================================================
  printf("---> Defining the pT ranges ... \n");
  //============================================================================
  TH2D *histCostPhiGen = (TH2D*) fileInput -> Get(Form("hCostPhiHE_%ipt%i_2m_gen",pt_min,pt_min+1)); // Gen distrib taken from non-weighted distribution
  TH3D *histCostPhiMassRec = (TH3D*) fileInput -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec",pt_min,pt_min+1)); // Gen distrib taken from noin-weighted distribution
  TH3D *histCostPhiMassRecWeighted = (TH3D*) fileInputWeighted -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec_weighted",pt_min,pt_min+1)); // Gen distrib taken from weighted distribution

  cout << pt_min << " - " << pt_min+1 << endl;

  for(int i = pt_min+1;i < pt_max;i++){
    cout << i << " - " << i+1 << endl;
    TH2D *tmpTH2HistoGen = (TH2D*) fileInput -> Get(Form("hCostPhiHE_%ipt%i_2m_gen",i,i+1));
    histCostPhiGen -> Add(tmpTH2HistoGen);

    TH3D *tmpTH3Histo = (TH3D*) fileInput -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec",i,i+1));
    histCostPhiMassRec -> Add(tmpTH3Histo);

    TH3D *tmpTH3HistoWeighted = (TH3D*) fileInputWeighted -> Get(Form("hCostPhiMassHE_%ipt%i_2m_rec_weighted",i,i+1));
    histCostPhiMassRecWeighted -> Add(tmpTH3HistoWeighted);
  }
  TH2D *histCostPhiRec = (TH2D*) histCostPhiMassRec -> Project3D("yx");
  TH2D *histCostPhiRecWeighted = (TH2D*) histCostPhiMassRecWeighted -> Project3D("yx");

  //============================================================================
  printf("---> Creating the TH2D with the binning tuned on data ... \n");
  //============================================================================
  TH2D *histCostPhiGenRebin = new TH2D("histCostPhiGenRebin","histCostPhiGenRebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  TH2D *histCostPhiRecRebin = new TH2D("histCostPhiRecRebin","histCostPhiRecRebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  TH2D *histCostPhiRecWeightedRebin = new TH2D("histCostPhiRecWeightedRebin","histCostPhiRecWeightedRebin",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      histCostPhiGenRebin -> SetBinContent(i+1,j+1,histCostPhiGen -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
      histCostPhiRecRebin -> SetBinContent(i+1,j+1,histCostPhiRec -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
      histCostPhiRecWeightedRebin -> SetBinContent(i+1,j+1,histCostPhiRecWeighted -> Integral(CostBinsMin[i],CostBinsMax[i],PhiBinsMin[j],PhiBinsMax[j]));
    }
  }

  //============================================================================
  printf("---> Computing Acc X Eff ... \n");
  //============================================================================
  cout << histCostPhiGenRebin -> Integral() << endl;
  cout << histCostPhiRecRebin -> Integral() << endl;
  cout << histCostPhiRecWeightedRebin -> Integral() << endl;

  histCostPhiGenRebin -> Sumw2();
  histCostPhiRecRebin -> Sumw2();
  histCostPhiRecWeightedRebin -> Sumw2();

  TH2D *histAccxEff = new TH2D("histAccxEff","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histAccxEff -> Divide(histCostPhiRecRebin,histCostPhiGenRebin,1,1,"B");

  TH2D *histAccxEffWeighted = new TH2D("histAccxEffWeighted","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histAccxEffWeighted -> Divide(histCostPhiRecWeightedRebin,histCostPhiGenRebin,1,1,"B");

  TCanvas *cAccxEff = new TCanvas("cAccxEff","cAccxEff",20,20,600,600);
  histAccxEff -> Draw("COLZtext");

  TCanvas *cAccxEffWeighted = new TCanvas("cAccxEffWeighted","cAccxEffWeighted",20,20,600,600);
  histAccxEffWeighted -> Draw("COLZtext");

  cout << histAccxEff -> GetBinContent(5,5) << " +- " << histAccxEff -> GetBinError(5,5) << endl;
  cout << histAccxEffWeighted -> GetBinContent(5,5) << " +- " << histAccxEffWeighted -> GetBinError(5,5) << endl;

  //////////
  TH1D *hhh = (TH1D*) histAccxEff -> ProjectionX("hhh");
  hhh -> SetLineColor(kRed);
  TH1D *hhhWeighted = (TH1D*) histAccxEffWeighted -> ProjectionX("hhhWeighted");
  hhhWeighted -> SetLineColor(kGreen);

  TCanvas *cccc = new TCanvas("cccc","cccc",20,20,600,600);
  hhh -> Draw("E");
  hhhWeighted -> Draw("Esame");
  ///////////

  TH2D *hist_accxeff_HE_rebin = (TH2D*) histAccxEffWeighted -> Clone("hist_accxeff_HE_rebin");

  //============================================================================
  printf("---> Saving Acc X Eff ... \n");
  //============================================================================
  printf("Saving Acc x Eff in %s ... \n",fileAccxEffName.c_str());
  TFile *fileOutput = new TFile(fileAccxEffName.c_str(),"RECREATE");
  hist_accxeff_HE_rebin -> Write();
  fileOutput -> Close();

}
