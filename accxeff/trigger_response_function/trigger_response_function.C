#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string.h>

#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TStopwatch.h>
#endif

void trigger_response_function(){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);

  TFile *fileHistMCFullStat = new TFile("MC_histograms/Hist_MC_Full_Stat.root","READ");
  TH1D *histMCLowpt = (TH1D*) fileHistMCFullStat -> Get("histoLowPt");
  TH1D *histMCAllpt = (TH1D*) fileHistMCFullStat -> Get("histoAllPt");

  TH1D *histTriggerResponseFunctionMC = new TH1D("histTriggerResponseFunctionMC","",100,0,10);
  histTriggerResponseFunctionMC -> Divide(histMCLowpt,histMCAllpt,1,1,"B");
  histTriggerResponseFunctionMC -> SetLineColor(kBlue);
  histTriggerResponseFunctionMC -> SetMarkerColor(kBlue);

  TFile *fileHistDataFullStat = new TFile("data_histograms/Hist_Data_Full_Stat.root","READ");
  TH1D *histDataLowpt = (TH1D*) fileHistDataFullStat -> Get("histoLowPt");
  TH1D *histDataAllpt = (TH1D*) fileHistDataFullStat -> Get("histoAllPt");

  TH1D *histTriggerResponseFunctionData = new TH1D("histTriggerResponseFunctionData","",100,0,10);
  histTriggerResponseFunctionData -> Divide(histDataLowpt,histDataAllpt,1,1,"B");
  histTriggerResponseFunctionData -> SetLineColor(kRed);
  histTriggerResponseFunctionData -> SetMarkerColor(kRed);

  TCanvas *cTriggerResponseFunctionComparison = new TCanvas("cTriggerResponseFunctionComparison","cTriggerResponseFunctionComparison",20,20,600,600);
  histTriggerResponseFunctionMC -> Draw("E");
  histTriggerResponseFunctionData -> Draw("Esame");

  //============================================================================
  printf("---> Ratio of the trigger response functions ... \n");
  //============================================================================
  TH1D *histRatioTriggerResponseFunction = new TH1D("histRatioTriggerResponseFunction","",100,0,10);
  histRatioTriggerResponseFunction -> Divide(histTriggerResponseFunctionData,histTriggerResponseFunctionMC,1,1);
  histRatioTriggerResponseFunction -> SetLineColor(kBlack);
  histRatioTriggerResponseFunction -> SetMarkerColor(kBlack);

  TCanvas *cRatioTriggerResponseFunctionComparison = new TCanvas("cRatioTriggerResponseFunctionComparison","cRatioTriggerResponseFunctionComparison",20,20,600,600);
  histRatioTriggerResponseFunction -> Draw("E");

  //============================================================================
  printf("---> Saving histograms ... \n");
  //============================================================================
  TFile *fileTriggerResponseFunction = new TFile("Trigger_Response_Function.root","RECREATE");
  histTriggerResponseFunctionData -> Write();
  histTriggerResponseFunctionMC -> Write();
  histRatioTriggerResponseFunction -> Write();
  fileTriggerResponseFunction -> Close();

}
