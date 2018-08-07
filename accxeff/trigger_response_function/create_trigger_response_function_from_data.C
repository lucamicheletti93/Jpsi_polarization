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

//------------------------------------------------------------------------------
// To run this macro use loop_on_trees.sh
// File are stored in "/media/luca/488AE2208AE20A70/PbPb_2015_Trees/"
//------------------------------------------------------------------------------
void create_trigger_response_function_from_data(int runNumber = 246994){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  TFile *fileTreeData = new TFile(Form("/media/luca/488AE2208AE20A70/PbPb_2015_Trees/Tree_%i.root",runNumber),"READ");

  char TrigClass[500];
  Float_t PercV0M, PercCL0, PercCL1;
  Int_t NMuons, NTracklets, NContributors;
  Double_t Vertex[3];
  Double_t Pt[300],E[300], Px[300], Py[300], Pz[300], Y[300], Eta[300];
  Double_t TrackChi2[300], MatchTrigChi2[300], DCA[300], RAtAbsEnd[300];
  Int_t Charge[300], MatchTrig[300];
  Int_t NDimu;
  Double_t DimuPt[3000], DimuPx[3000], DimuPy[3000],DimuPz[3000], DimuY[3000];
  Double_t DimuMass[3000];
  Int_t DimuCharge[3000], DimuMatch[3000];
  Int_t DimuMu[3000][2];
  Double_t CostHE[3000], PhiHE[3000], CostCS[3000], PhiCS[3000];
  UInt_t inpmask;
  Bool_t IsPhysSelected;

  TTree *treeData = (TTree*) fileTreeData -> Get("PbPbTree");
  treeData -> SetBranchAddress("FiredTriggerClasses",TrigClass);
  treeData -> SetBranchAddress("NMuons",&NMuons);
  /*treeData -> SetBranchAddress("Vertex",Vertex);
  treeData -> SetBranchAddress("PercentV0M",&PercV0M);*/
  treeData -> SetBranchAddress("Pt",Pt);
  /*treeData -> SetBranchAddress("E",E);
  treeData -> SetBranchAddress("Px",Px);
  treeData -> SetBranchAddress("Py",Py);
  treeData -> SetBranchAddress("Pz",Pz);
  treeData -> SetBranchAddress("Y",Y);
  treeData -> SetBranchAddress("Eta",Eta);*/
  treeData -> SetBranchAddress("MatchTrig",MatchTrig);
  /*treeData -> SetBranchAddress("MatchTrigChi2",MatchTrigChi2);
  treeData -> SetBranchAddress("Charge",Charge);
  treeData -> SetBranchAddress("RAtAbsEnd",RAtAbsEnd);
  treeData -> SetBranchAddress("NDimu",&NDimu);
  treeData -> SetBranchAddress("DimuPt",DimuPt);
  treeData -> SetBranchAddress("DimuPx",DimuPx);
  treeData -> SetBranchAddress("DimuPy",DimuPy);
  treeData -> SetBranchAddress("DimuPz",DimuPz);
  treeData -> SetBranchAddress("DimuY",DimuY);
  treeData -> SetBranchAddress("DimuMass",DimuMass);
  treeData -> SetBranchAddress("DimuCharge",DimuCharge);
  treeData -> SetBranchAddress("DimuMatch",DimuMatch);
  treeData -> SetBranchAddress("DimuMu",DimuMu);
  treeData -> SetBranchAddress("CostHE",CostHE);
  treeData -> SetBranchAddress("PhiHE",PhiHE);
  treeData -> SetBranchAddress("CostCS",CostCS);
  treeData -> SetBranchAddress("PhiCS",PhiCS);*/
  treeData -> SetBranchAddress("IsPhysSelected",&IsPhysSelected);

  TH1D *histoLowPt = new TH1D("histoLowPt","",100,0,10);
  histoLowPt -> Sumw2();
  TH1D *histoAllPt = new TH1D("histoAllPt","",100,0,10);
  histoAllPt -> Sumw2();

  double percentage = 0;
  int nEntries = treeData -> GetEntries();

  for(int i = 0;i < nEntries;i++){
    treeData -> GetEntry(i);
    TString Trigger = TrigClass;
    Bool_t TriggerSelected = kFALSE;
    percentage = ((double) i)/((double) nEntries)*100;
    printf("Run %i - progress %2.1f %% ... \r",runNumber,percentage);
    if(IsPhysSelected){
      if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")){TriggerSelected = kTRUE;}
      if(TriggerSelected){
        for(int j = 0;j < NMuons;j++){
          if(MatchTrig[j] >= 1){histoAllPt -> Fill(Pt[j]);}
          if(MatchTrig[j] >= 2){histoLowPt -> Fill(Pt[j]);}
        }
      }
    }
  }
  printf("\n");

  TH1D *histTriggerResponseFunctionData = new TH1D("histTriggerResponseFunctionData","",100,0,10);
  histTriggerResponseFunctionData -> Divide(histoLowPt,histoAllPt,1,1,"B");
  histTriggerResponseFunctionData -> SetLineColor(kRed);
  histTriggerResponseFunctionData -> SetMarkerColor(kRed);
  histTriggerResponseFunctionData -> Draw("E");

  TFile *fileHistData = new TFile(Form("data_histograms/Hist_%i.root",runNumber),"RECREATE");
  histoLowPt -> Write();
  histoAllPt -> Write();
  fileHistData -> Close();
  fileTreeData -> Close();
}
