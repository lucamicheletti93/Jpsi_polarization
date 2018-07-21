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

void mohamad_cross_check(){

  //============================================================================
  printf("---> Comparing my results with Mohamad's ... \n");
  //============================================================================
  printf("Range 2 < pT < 4 GeV/c \n");
  //============================================================================
  TFile *fileBinning2pt4 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_2pt4_test.root","READ");
  vector <double> CostValues2pt4, CostWidth2pt4, PhiValues2pt4, PhiWidth2pt4;
  vector <int> CostBinsMin2pt4, CostBinsMax2pt4, PhiBinsMin2pt4, PhiBinsMax2pt4;

  Binning *binning2pt4 = (Binning*) fileBinning2pt4 -> Get("Binning");
  CostValues2pt4 = binning2pt4 -> GetCostValues();
  CostWidth2pt4 = binning2pt4 -> GetCostWidth();
  CostBinsMin2pt4 = binning2pt4 -> GetCostBinsMin();
  CostBinsMax2pt4 = binning2pt4 -> GetCostBinsMax();
  const int NCostBins2pt4 = CostValues2pt4.size() - 1;
  PhiValues2pt4 = binning2pt4 -> GetPhiValues();
  PhiWidth2pt4 = binning2pt4 -> GetPhiWidth();
  PhiBinsMin2pt4 = binning2pt4 -> GetPhiBinsMin();
  PhiBinsMax2pt4 = binning2pt4 -> GetPhiBinsMax();
  const int NPhiBins2pt4 = PhiValues2pt4.size() - 1;

  TFile *fileNJpsi1D2pt4 = new TFile("/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_2pt4_test/2pt4.root","READ");
  TH1D *histNJpsi1DCost2pt4 = (TH1D*) fileNJpsi1D2pt4 -> Get("histNJpsiCost");

  TH1D *histNJpsi1DCost2pt4Scaled = new TH1D("histNJpsi1DCost2pt4Scaled","",NCostBins2pt4,&CostValues2pt4[0]);
  histNJpsi1DCost2pt4Scaled -> SetLineColor(kRed);
  histNJpsi1DCost2pt4Scaled -> SetMarkerColor(kRed);
  for(int i = 0;i < NCostBins2pt4;i++){
    histNJpsi1DCost2pt4Scaled -> SetBinContent(i+1,(histNJpsi1DCost2pt4 -> GetBinContent(i+1))/CostWidth2pt4[i]);
    histNJpsi1DCost2pt4Scaled -> SetBinError(i+1,(histNJpsi1DCost2pt4 -> GetBinError(i+1))/CostWidth2pt4[i]);
  }

  double NJpsiMoham2pt3[7] = {11639,11917,10339,9561,7315,5480,2826};
  double statJpsiMoham2pt3[7] = {348,376,349,361,376,245,195};
  double systJpsiMoham2pt3[7] = {150,156,196,245,212,233};
  double errJpsiMoham2pt3[7];

  double NJpsiMoham3pt4[7] = {6731,6002,5737,5389,4534,3260,2162};
  double statJpsiMoham3pt4[7] = {210,197,206,224,223,219,164};
  double systJpsiMoham3pt4[7] = {120,71,155,114,174,113,25};
  double errJpsiMoham3pt4[7];

  TH1D *histNJpsiMoham2pt3 = new TH1D("histNJpsiMoham2pt3","",7,0,0.7);
  TH1D *histNJpsiMoham3pt4 = new TH1D("histNJpsiMoham3pt4","",7,0,0.7);

  for(int i = 0;i < 7;i++){
    errJpsiMoham2pt3[i] = TMath::Sqrt(statJpsiMoham2pt3[i]*statJpsiMoham2pt3[i] + systJpsiMoham2pt3[i]*systJpsiMoham2pt3[i]);
    histNJpsiMoham2pt3 -> SetBinContent(i+1,NJpsiMoham2pt3[i]);
    histNJpsiMoham2pt3 -> SetBinError(i+1,errJpsiMoham2pt3[i]);

    errJpsiMoham3pt4[i] = TMath::Sqrt(statJpsiMoham3pt4[i]*statJpsiMoham3pt4[i] + systJpsiMoham3pt4[i]*systJpsiMoham3pt4[i]);
    histNJpsiMoham3pt4 -> SetBinContent(i+1,NJpsiMoham3pt4[i]);
    histNJpsiMoham3pt4 -> SetBinError(i+1,errJpsiMoham3pt4[i]);
  }

  TH1D *histNJpsiMoham2pt4 = new TH1D("histNJpsiMoham2pt4","",7,0,0.7);
  histNJpsiMoham2pt4 -> Add(histNJpsiMoham2pt3,histNJpsiMoham3pt4);
  histNJpsiMoham2pt4 -> Scale(5);

  TCanvas *cNJpsiMoham2pt4 = new TCanvas("cNJpsiMoham2pt4","cNJpsiMoham2pt4",20,20,600,600);
  histNJpsiMoham2pt4 -> Draw();
  histNJpsi1DCost2pt4Scaled -> Draw("same");

  //============================================================================
  printf("Range 4 < pT < 6 GeV/c \n");
  //============================================================================
  TFile *fileBinning4pt6 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_4pt6_test.root","READ");
  vector <double> CostValues4pt6, CostWidth4pt6, PhiValues4pt6, PhiWidth4pt6;
  vector <int> CostBinsMin4pt6, CostBinsMax4pt6, PhiBinsMin4pt6, PhiBinsMax4pt6;

  Binning *binning4pt6 = (Binning*) fileBinning4pt6 -> Get("Binning");
  CostValues4pt6 = binning4pt6 -> GetCostValues();
  CostWidth4pt6 = binning4pt6 -> GetCostWidth();
  CostBinsMin4pt6 = binning4pt6 -> GetCostBinsMin();
  CostBinsMax4pt6 = binning4pt6 -> GetCostBinsMax();
  const int NCostBins4pt6 = CostValues4pt6.size() - 1;
  PhiValues4pt6 = binning4pt6 -> GetPhiValues();
  PhiWidth4pt6 = binning4pt6 -> GetPhiWidth();
  PhiBinsMin4pt6 = binning4pt6 -> GetPhiBinsMin();
  PhiBinsMax4pt6 = binning4pt6 -> GetPhiBinsMax();
  const int NPhiBins4pt6 = PhiValues4pt6.size() - 1;

  TFile *fileNJpsi1D4pt6 = new TFile("/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_4pt6_test/4pt6.root","READ");
  TH1D *histNJpsi1DCost4pt6 = (TH1D*) fileNJpsi1D4pt6 -> Get("histNJpsiCost");

  TH1D *histNJpsi1DCost4pt6Scaled = new TH1D("histNJpsi1DCost4pt6Scaled","",NCostBins4pt6,&CostValues4pt6[0]);
  histNJpsi1DCost4pt6Scaled -> SetLineColor(kRed);
  histNJpsi1DCost4pt6Scaled -> SetMarkerColor(kRed);
  for(int i = 0;i < NCostBins4pt6;i++){
    histNJpsi1DCost4pt6Scaled -> SetBinContent(i+1,(histNJpsi1DCost4pt6 -> GetBinContent(i+1))/CostWidth4pt6[i]);
    histNJpsi1DCost4pt6Scaled -> SetBinError(i+1,(histNJpsi1DCost4pt6 -> GetBinError(i+1))/CostWidth4pt6[i]);
  }

  double NJpsiMoham4pt5[7] = {3050,3124,2705,2968,2554,2007,1431};
  double statJpsiMoham4pt5[7] = {120,120,126,141,143,143,102};
  double systJpsiMoham4pt5[7] = {33,62,54,54,88,50,17};
  double errJpsiMoham4pt5[7];

  double NJpsiMoham5pt6[7] = {1624,1721,1597,1488,1239,1258,1014};
  double statJpsiMoham5pt6[7] = {77,76,81,77,83,79,61};
  double systJpsiMoham5pt6[7] = {28,14,21,17,45,39,7};
  double errJpsiMoham5pt6[7];

  TH1D *histNJpsiMoham4pt5 = new TH1D("histNJpsiMoham4pt5","",7,0,0.7);
  TH1D *histNJpsiMoham5pt6 = new TH1D("histNJpsiMoham5pt6","",7,0,0.7);

  for(int i = 0;i < 7;i++){
    errJpsiMoham4pt5[i] = TMath::Sqrt(statJpsiMoham4pt5[i]*statJpsiMoham4pt5[i] + systJpsiMoham4pt5[i]*systJpsiMoham4pt5[i]);
    histNJpsiMoham4pt5 -> SetBinContent(i+1,NJpsiMoham4pt5[i]);
    histNJpsiMoham4pt5 -> SetBinError(i+1,errJpsiMoham4pt5[i]);

    errJpsiMoham5pt6[i] = TMath::Sqrt(statJpsiMoham5pt6[i]*statJpsiMoham5pt6[i] + systJpsiMoham5pt6[i]*systJpsiMoham5pt6[i]);
    histNJpsiMoham5pt6 -> SetBinContent(i+1,NJpsiMoham5pt6[i]);
    histNJpsiMoham5pt6 -> SetBinError(i+1,errJpsiMoham5pt6[i]);
  }

  TH1D *histNJpsiMoham4pt6 = new TH1D("histNJpsiMoham4pt6","",7,0,0.7);
  histNJpsiMoham4pt6 -> Add(histNJpsiMoham4pt5,histNJpsiMoham5pt6);
  histNJpsiMoham4pt6 -> Scale(5);

  TCanvas *cNJpsiMoham4pt6 = new TCanvas("cNJpsiMoham4pt6","cNJpsiMoham4pt6",20,20,600,600);
  histNJpsiMoham4pt6 -> Draw();
  histNJpsi1DCost4pt6Scaled -> Draw("same");

  //============================================================================
  printf("Range 6 < pT < 10 GeV/c \n");
  //============================================================================
  TFile *fileBinning6pt10 = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_6pt10_test.root","READ");
  vector <double> CostValues6pt10, CostWidth6pt10, PhiValues6pt10, PhiWidth6pt10;
  vector <int> CostBinsMin6pt10, CostBinsMax6pt10, PhiBinsMin6pt10, PhiBinsMax6pt10;

  Binning *binning6pt10 = (Binning*) fileBinning6pt10 -> Get("Binning");
  CostValues6pt10 = binning6pt10 -> GetCostValues();
  CostWidth6pt10 = binning6pt10 -> GetCostWidth();
  CostBinsMin6pt10 = binning6pt10 -> GetCostBinsMin();
  CostBinsMax6pt10 = binning6pt10 -> GetCostBinsMax();
  const int NCostBins6pt10 = CostValues6pt10.size() - 1;
  PhiValues6pt10 = binning6pt10 -> GetPhiValues();
  PhiWidth6pt10 = binning6pt10 -> GetPhiWidth();
  PhiBinsMin6pt10 = binning6pt10 -> GetPhiBinsMin();
  PhiBinsMax6pt10 = binning6pt10 -> GetPhiBinsMax();
  const int NPhiBins6pt10 = PhiValues6pt10.size() - 1;

  TFile *fileNJpsi1D6pt10 = new TFile("/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_6pt10_test/6pt10.root","READ");
  TH1D *histNJpsi1DCost6pt10 = (TH1D*) fileNJpsi1D6pt10 -> Get("histNJpsiCost");

  TH1D *histNJpsi1DCost6pt10Scaled = new TH1D("histNJpsi1DCost6pt10Scaled","",NCostBins6pt10,&CostValues6pt10[0]);
  histNJpsi1DCost6pt10Scaled -> SetLineColor(kRed);
  histNJpsi1DCost6pt10Scaled -> SetMarkerColor(kRed);
  for(int i = 0;i < NCostBins6pt10;i++){
    histNJpsi1DCost6pt10Scaled -> SetBinContent(i+1,(histNJpsi1DCost6pt10 -> GetBinContent(i+1))/CostWidth6pt10[i]);
    histNJpsi1DCost6pt10Scaled -> SetBinError(i+1,(histNJpsi1DCost6pt10 -> GetBinError(i+1))/CostWidth6pt10[i]);
  }
}
