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

void create_trigger_response_function_from_data(){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  string fileTreeMCFullStatName = "~/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root";
  TFile *fileTreeMCFullStat = new TFile(fileTreeMCFullStatName.c_str(),"READ");

  Int_t NMuons_rec;
  Double_t Pt_rec[3000];
  Int_t MatchTrig_rec[3000];

  TTree *treeMCFullStat = (TTree*) fileTreeMCFullStat -> Get("MCTree");
  treeMCFullStat -> SetBranchAddress("NMuons_rec",&NMuons_rec);
  treeMCFullStat -> SetBranchAddress("Pt_rec",Pt_rec);
  treeMCFullStat -> SetBranchAddress("MatchTrig_rec",MatchTrig_rec);

  TH1D *histoLowPt = new TH1D("histoLowPt","",100,0,10);
  histoLowPt -> Sumw2();
  TH1D *histoAllPt = new TH1D("histoAllPt","",100,0,10);
  histoAllPt -> Sumw2();

  double percentage = 0;
  int nEntries = treeMCFullStat -> GetEntries();

  for(int i = 0;i < nEntries;i++){
    treeMCFullStat -> GetEntry(i);
    percentage = ((double) i)/((double) nEntries)*100;
    printf("Progress %2.1f %% ... \r",percentage);
    //printf("%i - %i \n",MatchTrig_rec[0],MatchTrig_rec[1]);
    if(MatchTrig_rec[0] >= 1){histoAllPt -> Fill(Pt_rec[0]);}
    if(MatchTrig_rec[1] >= 1){histoAllPt -> Fill(Pt_rec[1]);}
    if(MatchTrig_rec[0] >= 2){histoLowPt -> Fill(Pt_rec[0]);}
    if(MatchTrig_rec[1] >= 2){histoLowPt -> Fill(Pt_rec[1]);}
  }
  printf("\n");

  //TH1D *histTriggerResponseFunction = new TH1D("histTriggerResponseFunction","",100,0,10);
  //histTriggerResponseFunction -> Divide(histoLowPt,histoAllPt,1,1,"B");
  //histTriggerResponseFunction -> Draw("E");

  TFile *fileHistMCFullStat = new TFile("MC_histograms/Hist_MC_Full_Stat.root","RECREATE");
  histoLowPt -> Write();
  histoAllPt -> Write();
  fileHistMCFullStat -> Close();
  fileTreeMCFullStat -> Close();
}
