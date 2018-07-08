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

  string fileInputName1D = "../signal_extraction/slice_fit_" + dataset + "/" + dataset + ".root";
  TFile *fileInput1D = new TFile(fileInputName1D.c_str());

  //============================================================================
  printf("---> Getting the number of J/psi ... \n");
  //============================================================================
  TH2D *histNJpsi = (TH2D*) fileInput2D -> Get("histNJpsi");

  TCanvas *cNJpsiCost = new TCanvas("cNJpsiCost","cNJpsiCost",20,20,600,600);

  TH1D *projNJpsiCost = (TH1D*) histNJpsi -> ProjectionX();
  projNJpsiCost -> SetMarkerColor(kRed);
  projNJpsiCost -> SetLineColor(kRed);

  TH1D *histNJpsiCost = (TH1D*) fileInput1D -> Get("histNJpsiCost");
  histNJpsiCost -> SetMarkerColor(kBlue);
  histNJpsiCost -> SetLineColor(kBlue);

  projNJpsiCost -> SetTitle("N_{J/#psi} vs cos#it{#theta}");
  projNJpsiCost -> Draw();
  histNJpsiCost -> Draw("same");

  TCanvas *cNJpsiPhi = new TCanvas("cNJpsiPhi","cNJpsiPhi",20,20,600,600);

  TH1D *projNJpsiPhi = (TH1D*) histNJpsi -> ProjectionY();
  projNJpsiPhi -> SetMarkerColor(kRed);
  projNJpsiPhi -> SetLineColor(kRed);

  TH1D *histNJpsiPhi = (TH1D*) fileInput1D -> Get("histNJpsiPhi");
  histNJpsiPhi -> SetMarkerColor(kBlue);
  histNJpsiPhi -> SetLineColor(kBlue);

  projNJpsiPhi -> SetTitle("N_{J/#psi} vs #it{#varphi}");
  projNJpsiPhi -> Draw();
  histNJpsiPhi -> Draw("same");

}
