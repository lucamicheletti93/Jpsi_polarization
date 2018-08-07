#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

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

void create_weighted_rec_distribution(){

  double PI = TMath::Pi();
  //============================================================================
  printf("---> Setting output histograms ... \n");
  //============================================================================
  TH3D *hCostPhiMassHE_0pt1_2m_rec_weighted = new TH3D("hCostPhiMassHE_0pt1_2m_rec_weighted","hCostPhiMassHE_0pt1_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_1pt2_2m_rec_weighted = new TH3D("hCostPhiMassHE_1pt2_2m_rec_weighted","hCostPhiMassHE_1pt2_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_2pt3_2m_rec_weighted = new TH3D("hCostPhiMassHE_2pt3_2m_rec_weighted","hCostPhiMassHE_2pt3_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_3pt4_2m_rec_weighted = new TH3D("hCostPhiMassHE_3pt4_2m_rec_weighted","hCostPhiMassHE_3pt4_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_4pt5_2m_rec_weighted = new TH3D("hCostPhiMassHE_4pt5_2m_rec_weighted","hCostPhiMassHE_4pt5_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_5pt6_2m_rec_weighted = new TH3D("hCostPhiMassHE_5pt6_2m_rec_weighted","hCostPhiMassHE_5pt6_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_6pt7_2m_rec_weighted = new TH3D("hCostPhiMassHE_6pt7_2m_rec_weighted","hCostPhiMassHE_6pt7_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_7pt8_2m_rec_weighted = new TH3D("hCostPhiMassHE_7pt8_2m_rec_weighted","hCostPhiMassHE_7pt8_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_8pt9_2m_rec_weighted = new TH3D("hCostPhiMassHE_8pt9_2m_rec_weighted","hCostPhiMassHE_8pt9_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_9pt10_2m_rec_weighted = new TH3D("hCostPhiMassHE_9pt10_2m_rec_weighted","hCostPhiMassHE_9pt10_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_10pt11_2m_rec_weighted = new TH3D("hCostPhiMassHE_10pt11_2m_rec_weighted","hCostPhiMassHE_10pt11_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_11pt12_2m_rec_weighted = new TH3D("hCostPhiMassHE_11pt12_2m_rec_weighted","hCostPhiMassHE_11pt12_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_12ptinf_2m_rec_weighted = new TH3D("hCostPhiMassHE_12ptinf_2m_rec_weighted","hCostPhiMassHE_12ptinf_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);

  TH3D *hCostPhiMassCS_0pt1_2m_rec_weighted = new TH3D("hCostPhiMassCS_0pt1_2m_rec_weighted","hCostPhiMassCS_0pt1_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_1pt2_2m_rec_weighted = new TH3D("hCostPhiMassCS_1pt2_2m_rec_weighted","hCostPhiMassCS_1pt2_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_2pt3_2m_rec_weighted = new TH3D("hCostPhiMassCS_2pt3_2m_rec_weighted","hCostPhiMassCS_2pt3_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_3pt4_2m_rec_weighted = new TH3D("hCostPhiMassCS_3pt4_2m_rec_weighted","hCostPhiMassCS_3pt4_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_4pt5_2m_rec_weighted = new TH3D("hCostPhiMassCS_4pt5_2m_rec_weighted","hCostPhiMassCS_4pt5_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_5pt6_2m_rec_weighted = new TH3D("hCostPhiMassCS_5pt6_2m_rec_weighted","hCostPhiMassCS_5pt6_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_6pt7_2m_rec_weighted = new TH3D("hCostPhiMassCS_6pt7_2m_rec_weighted","hCostPhiMassCS_6pt7_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_7pt8_2m_rec_weighted = new TH3D("hCostPhiMassCS_7pt8_2m_rec_weighted","hCostPhiMassCS_7pt8_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_8pt9_2m_rec_weighted = new TH3D("hCostPhiMassCS_8pt9_2m_rec_weighted","hCostPhiMassCS_8pt9_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_9pt10_2m_rec_weighted = new TH3D("hCostPhiMassCS_9pt10_2m_rec_weighted","hCostPhiMassCS_9pt10_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_10pt11_2m_rec_weighted = new TH3D("hCostPhiMassCS_10pt11_2m_rec_weighted","hCostPhiMassCS_10pt11_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_11pt12_2m_rec_weighted = new TH3D("hCostPhiMassCS_11pt12_2m_rec_weighted","hCostPhiMassCS_11pt12_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_12ptinf_2m_rec_weighted = new TH3D("hCostPhiMassCS_12ptinf_2m_rec_weighted","hCostPhiMassCS_12ptinf_2m_rec_weighted",100,-1,1,50,0,PI,120,2,5);

  //============================================================================
  printf("---> Opening Trigger_Response_Function.root ... \n");
  //============================================================================
  TFile *fileTriggerResponseFunction = new TFile("Trigger_Response_Function.root","READ");
  TH1D *histRatioTriggerResponseFunction = (TH1D*) fileTriggerResponseFunction -> Get("histRatioTriggerResponseFunction");
  histRatioTriggerResponseFunction -> Draw("E");

  //============================================================================
  printf("---> Opening MC_official_tree_Jpsi_PbPb_Nopol.root ... \n");
  //============================================================================
  Int_t NDimu_rec;
  Double_t DimuPt_rec[3000], DimuY_rec[3000];
  Double_t DimuMass_rec[3000];
  Int_t DimuMatch_rec[3000];
  Double_t Px_rec[3000];
  Double_t Pt_rec[3000];
  Double_t DimuPx_rec[3000];
  Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];

  TFile *fileAccxEff = new TFile("~/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root","READ");
  TTree *treeAccxEff = (TTree*) fileAccxEff -> Get("MCTree");
  treeAccxEff -> SetBranchAddress("NDimu_rec",&NDimu_rec);
  treeAccxEff -> SetBranchAddress("DimuPt_rec",DimuPt_rec);
  treeAccxEff -> SetBranchAddress("DimuY_rec",DimuY_rec);
  treeAccxEff -> SetBranchAddress("DimuMass_rec",DimuMass_rec);
  treeAccxEff -> SetBranchAddress("DimuMatch_rec",DimuMatch_rec);
  treeAccxEff -> SetBranchAddress("CostHE_rec",CostHE_rec);
  treeAccxEff -> SetBranchAddress("PhiHE_rec",PhiHE_rec);
  treeAccxEff -> SetBranchAddress("CostCS_rec",CostCS_rec);
  treeAccxEff -> SetBranchAddress("PhiCS_rec",PhiCS_rec);

  treeAccxEff -> SetBranchAddress("Pt_rec",Pt_rec);
  treeAccxEff -> SetBranchAddress("Px_rec",Px_rec);
  treeAccxEff -> SetBranchAddress("DimuPx_rec",DimuPx_rec);

  //============================================================================
  printf("---> Weighting the AccxEff ... \n");
  //============================================================================
  double weightMu1 = 0, weightMu2 = 0;
  double weight = 0;

  Int_t nEntries = treeAccxEff -> GetEntries();
  for(int i = 0;i < nEntries;i++){
    printf("GENERATED %i -> %i : %3.2f%% \r",nEntries,i,(double) i/ (double) nEntries*100);
    treeAccxEff -> GetEntry(i);
    for(int k = 0;k < NDimu_rec;k++){
      if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
        if(DimuMatch_rec[k] == 2){
          if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
            weightMu1 = histRatioTriggerResponseFunction -> GetBinContent(histRatioTriggerResponseFunction -> FindBin(Pt_rec[0]));
            weightMu2 = histRatioTriggerResponseFunction -> GetBinContent(histRatioTriggerResponseFunction -> FindBin(Pt_rec[1]));
            weight = weightMu1*weightMu2;

            if(DimuPt_rec[k] > 0 && DimuPt_rec[k] <= 1){hCostPhiMassHE_0pt1_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_0pt1_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 1 && DimuPt_rec[k] <= 2){hCostPhiMassHE_1pt2_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_1pt2_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 3){hCostPhiMassHE_2pt3_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_2pt3_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 3 && DimuPt_rec[k] <= 4){hCostPhiMassHE_3pt4_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_3pt4_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 4 && DimuPt_rec[k] <= 5){hCostPhiMassHE_4pt5_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_4pt5_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 5 && DimuPt_rec[k] <= 6){hCostPhiMassHE_5pt6_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_5pt6_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 6 && DimuPt_rec[k] <= 7){hCostPhiMassHE_6pt7_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_6pt7_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 7 && DimuPt_rec[k] <= 8){hCostPhiMassHE_7pt8_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_7pt8_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 8 && DimuPt_rec[k] <= 9){hCostPhiMassHE_8pt9_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_8pt9_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 9 && DimuPt_rec[k] <= 10){hCostPhiMassHE_9pt10_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_9pt10_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 10 && DimuPt_rec[k] <= 11){hCostPhiMassHE_10pt11_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_10pt11_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 11 && DimuPt_rec[k] <= 12){hCostPhiMassHE_11pt12_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_11pt12_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
            if(DimuPt_rec[k] > 12){hCostPhiMassHE_12ptinf_2m_rec_weighted -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k],weight); hCostPhiMassCS_12ptinf_2m_rec_weighted -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k],weight);}
          }
        }
      }
    }
  }
  printf("\n");

  //==========================================================================
  printf("---> Saving histos ... \n");
  //==========================================================================
  TFile *fileJpsiRecWeighted = new TFile("Histos_Jpsi_Rec_Weigthed.root","RECREATE");
  fileJpsiRecWeighted -> cd();
  hCostPhiMassHE_0pt1_2m_rec_weighted -> Write();
  hCostPhiMassHE_1pt2_2m_rec_weighted -> Write();
  hCostPhiMassHE_2pt3_2m_rec_weighted -> Write();
  hCostPhiMassHE_3pt4_2m_rec_weighted -> Write();
  hCostPhiMassHE_4pt5_2m_rec_weighted -> Write();
  hCostPhiMassHE_5pt6_2m_rec_weighted -> Write();
  hCostPhiMassHE_6pt7_2m_rec_weighted -> Write();
  hCostPhiMassHE_7pt8_2m_rec_weighted -> Write();
  hCostPhiMassHE_8pt9_2m_rec_weighted -> Write();
  hCostPhiMassHE_9pt10_2m_rec_weighted -> Write();
  hCostPhiMassHE_10pt11_2m_rec_weighted -> Write();
  hCostPhiMassHE_11pt12_2m_rec_weighted -> Write();
  hCostPhiMassHE_12ptinf_2m_rec_weighted -> Write();

  hCostPhiMassCS_0pt1_2m_rec_weighted -> Write();
  hCostPhiMassCS_1pt2_2m_rec_weighted -> Write();
  hCostPhiMassCS_2pt3_2m_rec_weighted -> Write();
  hCostPhiMassCS_3pt4_2m_rec_weighted -> Write();
  hCostPhiMassCS_4pt5_2m_rec_weighted -> Write();
  hCostPhiMassCS_5pt6_2m_rec_weighted -> Write();
  hCostPhiMassCS_6pt7_2m_rec_weighted -> Write();
  hCostPhiMassCS_7pt8_2m_rec_weighted -> Write();
  hCostPhiMassCS_8pt9_2m_rec_weighted -> Write();
  hCostPhiMassCS_9pt10_2m_rec_weighted -> Write();
  hCostPhiMassCS_10pt11_2m_rec_weighted -> Write();
  hCostPhiMassCS_11pt12_2m_rec_weighted -> Write();
  hCostPhiMassCS_12ptinf_2m_rec_weighted -> Write();
  fileJpsiRecWeighted -> Close();


}
