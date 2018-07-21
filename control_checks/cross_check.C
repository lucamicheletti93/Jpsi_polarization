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

void cross_check(int pt_min, int pt_max){
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
  string fileNJpsiFixedSigmaName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "_fixed_sigma" + nameOption + ".root";
  string fileNJpsiFreeSigmaName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + nameOption + ".root";
  string fileNJpsi1DName = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_" + dataset + nameOption + "/" + dataset + ".root";

  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  TFile *fileNJpsiFixedSigma = new TFile(fileNJpsiFixedSigmaName.c_str(),"READ");
  TFile *fileNJpsiFreeSigma = new TFile(fileNJpsiFreeSigmaName.c_str(),"READ");
  TFile *fileNJpsi1D = new TFile(fileNJpsi1DName.c_str(),"READ");

  //============================================================================
  printf("---> Reading the binning file ... \n");
  //============================================================================
  vector <double> CostValues, CostWidth;
  vector <double> PhiValues, PhiWidth;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;

  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostWidth = binning -> GetCostWidth();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiWidth = binning -> GetPhiWidth();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  //============================================================================
  printf("---> Getting N J/psi for fixed and free J/psi sigma and drawing ... \n");
  //============================================================================
  TH2D *histNJpsiFixedSigma = (TH2D*) fileNJpsiFixedSigma -> Get("histNJpsi");
  histNJpsiFixedSigma -> SetTitle("N_{J/#psi} #sigma_{J/#psi} fixed");
  TH2D *histNJpsiFreeSigma = (TH2D*) fileNJpsiFreeSigma -> Get("histNJpsi");
  histNJpsiFreeSigma -> SetTitle("N_{J/#psi} #sigma_{J/#psi} free");

  TCanvas *cNJpsiComparison = new TCanvas("cNJpsiComparison","cNJpsiComparison",20,20,600,600);
  cNJpsiComparison -> Divide(2,2);

  cNJpsiComparison -> cd(1);
  histNJpsiFixedSigma -> Draw("COLZtext");

  cNJpsiComparison -> cd(2);
  histNJpsiFreeSigma -> Draw("COLZtext");

  cNJpsiComparison -> cd(3);
  TH2D *histRatioNJpsiFreeFixedSigma = new TH2D("histRatioNJpsiFreeFixedSigma","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  histRatioNJpsiFreeFixedSigma -> Divide(histNJpsiFreeSigma,histNJpsiFixedSigma,1,1);
  histRatioNJpsiFreeFixedSigma -> Draw("COLZtext");
  gStyle -> SetPaintTextFormat("2.2f");

  //============================================================================
  printf("---> Comparing the N J/psi in 1D and drawing ... \n");
  //============================================================================
  TH1D *histNJpsi1DCost = (TH1D*) fileNJpsi1D -> Get("histNJpsiCost");
  histNJpsi1DCost -> SetMarkerStyle(20);
  histNJpsi1DCost -> SetMarkerColor(kBlack);
  histNJpsi1DCost -> SetLineColor(kBlack);

  TH1D *histNJpsi2DFreeSigmaCost = (TH1D*) histNJpsiFreeSigma -> ProjectionX("histNJpsi2DFreeSigmaCost");
  histNJpsi2DFreeSigmaCost -> SetMarkerStyle(20);
  histNJpsi2DFreeSigmaCost -> SetMarkerColor(kRed);
  histNJpsi2DFreeSigmaCost -> SetLineColor(kRed);

  TH1D *histNJpsi2DFixedSigmaCost = (TH1D*) histNJpsiFixedSigma -> ProjectionX("histNJpsi2DFixedSigmaCost");
  histNJpsi2DFixedSigmaCost -> SetMarkerStyle(20);
  histNJpsi2DFixedSigmaCost -> SetMarkerColor(kGreen+2);
  histNJpsi2DFixedSigmaCost -> SetLineColor(kGreen+2);

  TCanvas *cNJpsi1DComparison = new TCanvas("cNJpsi1DComparison","cNJpsi1DComparison",20,20,600,600);
  TH2D *hNJpsi1DComparison = new TH2D("N_{J/#psi} along cos(#it{#theta}) comparison","",100,-1,1,100,0,1000);
  hNJpsi1DComparison -> Draw();
  histNJpsi2DFreeSigmaCost -> Draw("same");
  histNJpsi2DFixedSigmaCost -> Draw("same");
  histNJpsi1DCost -> Draw("same");

  TH1D *histNJpsi1DCostScaled = new TH1D("histNJpsi1DCostScaled","",NCostBins,&CostValues[0]);
  histNJpsi1DCostScaled -> SetLineColor(kRed);
  histNJpsi1DCostScaled -> SetMarkerColor(kRed);
  for(int i = 0;i < NCostBins;i++){
    histNJpsi1DCostScaled -> SetBinContent(i+1,(histNJpsi1DCost -> GetBinContent(i+1))/CostWidth[i]);
    histNJpsi1DCostScaled -> SetBinError(i+1,(histNJpsi1DCost -> GetBinError(i+1))/CostWidth[i]);
  }

}
