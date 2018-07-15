#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

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

void mono_vs_bidi_appoach(int pt_min, int pt_max){
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

  string fileInput2DName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + nameOption + ".root";
  string fileInput1DBinnedName = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/binned_1D_" + dataset + nameOption + "/" + dataset + ".root";
  string fileInput1DUnbinnedName = "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/signal_extraction/1D_fit/unbinned_1D_" + dataset + "/" + dataset + ".root";

  nameOption = "fixed_sigma_test";
  string fileInput2DFixedSigmaName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/NEW_GIT_OUTPUT/N_Jpsi_" + dataset + nameOption + ".root";

  TFile *fileInput2D = new TFile(fileInput2DName.c_str(),"READ");
  TFile *fileInput1DBinned = new TFile(fileInput1DBinnedName.c_str(),"READ");
  TFile *fileInput1DUnbinned = new TFile(fileInput1DUnbinnedName.c_str(),"READ");
  TFile *fileInput2DFixedSigma = new TFile(fileInput2DFixedSigmaName.c_str(),"READ");

  //============================================================================
  printf("---> Getting the number of J/psi ... \n");
  //============================================================================
  TH1D *histNJpsi1DCostBinned = (TH1D*) fileInput1DBinned -> Get("histNJpsiCost");
  histNJpsi1DCostBinned -> SetMarkerColor(kBlack);
  histNJpsi1DCostBinned -> SetLineColor(kBlack);

  TH1D *histSigmaJpsi1DCostBinned = (TH1D*) fileInput1DBinned -> Get("histSigmaJpsiCost");
  histSigmaJpsi1DCostBinned -> SetMarkerColor(kBlack);
  histSigmaJpsi1DCostBinned -> SetLineColor(kBlack);

  TH2D *histNJpsiFreeSigma = (TH2D*) fileInput2D -> Get("histNJpsi");

  TH1D *histNJpsi2DFreeSigmaCost = (TH1D*) histNJpsiFreeSigma -> ProjectionX("histNJpsi2DFreeSigmaCost");
  histNJpsi2DFreeSigmaCost -> SetMarkerColor(kRed);
  histNJpsi2DFreeSigmaCost -> SetLineColor(kRed);

  TH1D *histNJpsi2DFreeSigmaPhi = (TH1D*) histNJpsiFreeSigma -> ProjectionY("histNJpsi2DFreeSigmaPhi");
  histNJpsi2DFreeSigmaPhi -> SetMarkerColor(kRed);
  histNJpsi2DFreeSigmaPhi -> SetLineColor(kRed);

  TH1D *histNJpsi1DCostUnbinned = (TH1D*) fileInput1DUnbinned -> Get("histNJpsiCost");
  histNJpsi1DCostUnbinned -> SetMarkerColor(kMagenta);
  histNJpsi1DCostUnbinned -> SetLineColor(kMagenta);

  TH1D *histSigmaJpsi1DCostUnbinned = (TH1D*) fileInput1DUnbinned -> Get("histSigmaJpsiCost");
  histSigmaJpsi1DCostUnbinned -> SetMarkerColor(kMagenta);
  histSigmaJpsi1DCostUnbinned -> SetLineColor(kMagenta);

  TH2D *histNJpsiFixedSigma = (TH2D*) fileInput2DFixedSigma -> Get("histNJpsi");

  TH1D *histNJpsi2DFixedSigmaCost = (TH1D*) histNJpsiFixedSigma -> ProjectionX("histNJpsi2DFixedSigmaCost");
  histNJpsi2DFixedSigmaCost -> SetMarkerColor(kGreen+2);
  histNJpsi2DFixedSigmaCost -> SetLineColor(kGreen+2);

  TH1D *histNJpsi2DFixedSigmaPhi = (TH1D*) histNJpsiFixedSigma -> ProjectionY("histNJpsi2DFixedSigmaPhi");
  histNJpsi2DFixedSigmaPhi -> SetMarkerColor(kGreen+2);
  histNJpsi2DFixedSigmaPhi -> SetLineColor(kGreen+2);

  //============================================================================
  printf("---> Computing the ratio of the between the number of J/psi ... \n");
  //============================================================================
  TH1D *histRatioNJpsi1DBinned2DFreeSigma = (TH1D*) histNJpsi2DFreeSigmaCost -> Clone("histRatioNJpsi1DBinned2DFreeSigma");
  histRatioNJpsi1DBinned2DFreeSigma -> Divide(histNJpsi1DCostBinned);

  TH1D *histRatioNJpsi1DBinned2DFixedSigma = (TH1D*) histNJpsi2DFixedSigmaCost -> Clone("histRatioNJpsi1DBinned2DFixedSigma");
  histRatioNJpsi1DBinned2DFixedSigma -> Divide(histNJpsi1DCostBinned);

  //============================================================================
  printf("---> Comparing the number of J/psi ... \n");
  //============================================================================
  TLine *lUnity = new TLine(-1,1,1,1);
  lUnity -> SetLineColor(kGray+2);
  lUnity -> SetLineStyle(2);

  TCanvas *cJpsiCost = new TCanvas("cJpsiCost","cJpsiCost",800,800);

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
  pad1 -> SetBottomMargin(0);
  pad1 -> Draw();
  pad1 -> cd();

  TH2D *hNJpsiCost = new TH2D("hNJpsiCost"," ",100,-1.,1.,100,0,1000);
  hNJpsiCost -> GetYaxis() -> SetTitle("N_{J/#psi}");
  hNJpsiCost -> GetYaxis() -> SetTitleOffset(0.95);
  hNJpsiCost -> GetYaxis() -> SetLabelSize(0);

  TLegend *legNJpsiCost = new TLegend(0.25,0.7,.55,.89,"","brNDC");
  legNJpsiCost -> SetBorderSize(0);
  legNJpsiCost -> SetFillColor(10);
  legNJpsiCost -> SetFillStyle(1);
  legNJpsiCost -> SetLineStyle(0);
  legNJpsiCost -> SetLineColor(0);
  legNJpsiCost -> SetTextFont(42);
  legNJpsiCost -> SetTextSize(0.035);
  legNJpsiCost -> AddEntry(histNJpsi2DFreeSigmaCost,"N_{J/#psi} free #sigma_{J/#psi} 2D","L");
  legNJpsiCost -> AddEntry(histNJpsi2DFixedSigmaCost,"N_{J/#psi} fixed #sigma_{J/#psi} 2D","L");
  legNJpsiCost -> AddEntry(histNJpsi1DCostBinned,"N_{J/#psi} 1D","L");

  TGaxis *axisNJpsiCost = new TGaxis(-1.,10,-1.,900,10,900,510,"");
  axisNJpsiCost -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisNJpsiCost -> SetLabelSize(15);

  hNJpsiCost -> Draw();
  axisNJpsiCost -> Draw("same");
  histNJpsi2DFreeSigmaCost -> Draw("same");
  histNJpsi1DCostBinned -> Draw("same");
  histNJpsi2DFixedSigmaCost -> Draw("same");
  legNJpsiCost -> Draw("same");

  cJpsiCost -> cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
  pad2 -> SetTopMargin(0);
  pad2 -> SetBottomMargin(0.2);
  pad2 -> Draw();
  pad2 -> cd();

  TH2D *hRatioNJpsiCost = new TH2D("hRatioNJpsiCost"," ",100,-1.,1.,100,0.4,1.6);
  hRatioNJpsiCost -> GetYaxis() -> SetTitle("Ratio");
  hRatioNJpsiCost -> GetYaxis() -> SetTitleSize(25);
  hRatioNJpsiCost -> GetYaxis() -> SetTitleOffset(0.95);
  hRatioNJpsiCost -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisSigmaJpsiCost = new TGaxis(-1.,0.45,-1.,1.55,0.45,1.55,510,"");
  axisSigmaJpsiCost -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisSigmaJpsiCost -> SetLabelSize(15);

  hRatioNJpsiCost -> Draw();
  axisSigmaJpsiCost -> Draw("same");
  lUnity -> Draw("same");
  histRatioNJpsi1DBinned2DFreeSigma -> Draw("same");
  histRatioNJpsi1DBinned2DFixedSigma -> Draw("same");

}
