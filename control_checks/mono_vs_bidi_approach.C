#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TMinuit.h>
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
#include <TCollection.h>
#include <TKey.h>
#include <TGaxis.h>

#include "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.h"
#endif

void mono_vs_bidi_appoach(){

  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  string dataset = "6pt10";

  string fileInputName2D = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "free_sigma.root";
  TFile *fileInput2D = new TFile(fileInputName2D.c_str());
  TH2D *histNJpsi = (TH2D*) fileInput2D -> Get("histNJpsi");
  //TCanvas *ccc = new TCanvas("ccc","ccc",20,20,600,600);
  //histNJpsi -> Draw("COLZtext");

  //string fileInputName2DRebin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "free_sigma_rebin.root";
  //string fileInputName2DRebin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "fixed_sigma.root";
  string fileInputName2DRebin = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + "_test.root";
  TFile *fileInput2DRebin = new TFile(fileInputName2DRebin.c_str());
  TH2D *histNJpsiRebin = (TH2D*) fileInput2DRebin -> Get("histNJpsi");
  //TCanvas *cccc = new TCanvas("cccc","cccc",20,20,600,600);
  //histNJpsiRebin -> Draw("COLZtext");

  //string fileInputName1D = "../signal_extraction/slice_fit_" + dataset + "/" + dataset + ".root";
  //string fileInputName1DBinned = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_" + dataset + "/" + dataset + ".root";
  string fileInputName1DBinned = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_" + dataset + "_test/" + dataset + ".root";
  TFile *fileInput1DBinned = new TFile(fileInputName1DBinned.c_str());

  string fileInputName1DUnbinned = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/unbinned_1D_" + dataset + "/" + dataset + ".root";
  TFile *fileInput1DUnbinned = new TFile(fileInputName1DUnbinned.c_str());
  fileInput1DUnbinned -> ls();

  //============================================================================
  printf("---> Getting the number of J/psi ... \n");
  //============================================================================
  TH1D *projNJpsiCost = (TH1D*) histNJpsi -> ProjectionX("projNJpsiCost");
  projNJpsiCost -> SetMarkerColor(kRed);
  projNJpsiCost -> SetLineColor(kRed);

  TH1D *projNJpsiCostRebin = (TH1D*) histNJpsiRebin -> ProjectionX("projNJpsiCostRebin");
  projNJpsiCostRebin -> SetMarkerColor(kGreen+2);
  projNJpsiCostRebin -> SetLineColor(kGreen+2);

  TH1D *histNJpsiCostBinned = (TH1D*) fileInput1DBinned -> Get("histNJpsiCost");
  histNJpsiCostBinned -> SetMarkerColor(kBlue);
  histNJpsiCostBinned -> SetLineColor(kBlue);

  TH1D *histSigmaJpsiCostBinned = (TH1D*) fileInput1DBinned -> Get("histSigmaJpsiCost");
  histSigmaJpsiCostBinned -> SetMarkerColor(kBlue);
  histSigmaJpsiCostBinned -> SetLineColor(kBlue);

  TH1D *histNJpsiCostUnbinned = (TH1D*) fileInput1DUnbinned -> Get("histNJpsiCost");
  histNJpsiCostUnbinned -> SetMarkerColor(kMagenta);
  histNJpsiCostUnbinned -> SetLineColor(kMagenta);

  TH1D *histSigmaJpsiCostUnbinned = (TH1D*) fileInput1DUnbinned -> Get("histSigmaJpsiCost");
  histSigmaJpsiCostUnbinned -> SetMarkerColor(kMagenta);
  histSigmaJpsiCostUnbinned -> SetLineColor(kMagenta);

  TH2D *hNJpsiCost = new TH2D("hNJpsiCost"," ",100,-1.,1.,100,0,1000);
  hNJpsiCost -> GetYaxis() -> SetTitle("N_{J/#psi}");
  hNJpsiCost -> GetYaxis() -> SetTitleOffset(0.95);
  hNJpsiCost -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisNJpsiCost = new TGaxis(-1.,10,-1.,900,10,900,510,"");
  axisNJpsiCost -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisNJpsiCost -> SetLabelSize(15);

  TH2D *hSigmaJpsiCost = new TH2D("hSigmaJpsiCost"," ",100,-1.,1.,100,0,0.15);
  hSigmaJpsiCost -> GetYaxis() -> SetTitle("#sigma_{J/#psi}");
  hSigmaJpsiCost -> GetYaxis() -> SetTitleSize(25);
  hSigmaJpsiCost -> GetYaxis() -> SetTitleOffset(0.95);
  hSigmaJpsiCost -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisSigmaJpsiCost = new TGaxis(-1.,0.01,-1.,0.14,0.01,0.14,510,"");
  axisSigmaJpsiCost -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisSigmaJpsiCost -> SetLabelSize(15);

  TCanvas *cJpsiCost = new TCanvas("cJpsiCost","cJpsiCost",800,800);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
  pad1 -> SetBottomMargin(0);
  pad1 -> Draw();
  pad1 -> cd();
  projNJpsiCost -> SetTitle("N_{J/#psi} vs cos#it{#theta}");
  hNJpsiCost -> Draw();
  axisNJpsiCost -> Draw("same");
  projNJpsiCost -> Draw("same");
  projNJpsiCostRebin -> Draw("same");
  histNJpsiCostBinned -> Draw("same");
  histNJpsiCostUnbinned -> Draw("same");

  cJpsiCost -> cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
  pad2 -> SetTopMargin(0);
  pad2 -> SetBottomMargin(0.2);
  pad2 -> Draw();
  pad2 -> cd();
  hSigmaJpsiCost -> Draw();
  axisSigmaJpsiCost -> Draw("same");
  histSigmaJpsiCostBinned -> Draw("same");
  histSigmaJpsiCostUnbinned -> Draw("same");

  return;

  TCanvas *cNJpsiPhi = new TCanvas("cNJpsiPhi","cNJpsiPhi",20,20,600,600);

  TH1D *projNJpsiPhi = (TH1D*) histNJpsi -> ProjectionY("projNJpsiPhi");
  projNJpsiPhi -> SetMarkerColor(kRed);
  projNJpsiPhi -> SetLineColor(kRed);

  TH1D *projNJpsiPhiRebin = (TH1D*) histNJpsiRebin -> ProjectionY("projNJpsiPhiRebin");
  projNJpsiPhiRebin -> SetMarkerColor(kGreen+2);
  projNJpsiPhiRebin -> SetLineColor(kGreen+2);

  TH1D *histNJpsiPhi = (TH1D*) fileInput1DBinned -> Get("histNJpsiPhi");
  histNJpsiPhi -> SetMarkerColor(kBlue);
  histNJpsiPhi -> SetLineColor(kBlue);

  projNJpsiPhi -> SetTitle("N_{J/#psi} vs #it{#varphi}");
  projNJpsiPhi -> Draw();
  projNJpsiPhiRebin -> Draw("same");
  histNJpsiPhi -> Draw("same");

}
