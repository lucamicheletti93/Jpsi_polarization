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

void create_TH2_from_MCtree(){

  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  char *PATH_IN = "~/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT";
  char FILE_NAME_IN[400];
  sprintf(FILE_NAME_IN,"%s/MC_official_tree_Jpsi_PbPb_Nopol.root",PATH_IN);
  
  //============================================================================
  printf("---> Defining TH2 (gen) and TH3 (rec) histogram ... \n");
  //============================================================================
  double PI = TMath::Pi();

  TH2D *hCostPhiHE_0pt1_2m_gen = new TH2D("hCostPhiHE_0pt1_2m_gen","hCostPhiHE_0pt1_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_1pt2_2m_gen = new TH2D("hCostPhiHE_1pt2_2m_gen","hCostPhiHE_1pt2_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_2pt3_2m_gen = new TH2D("hCostPhiHE_2pt3_2m_gen","hCostPhiHE_2pt3_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_3pt4_2m_gen = new TH2D("hCostPhiHE_3pt4_2m_gen","hCostPhiHE_3pt4_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_4pt5_2m_gen = new TH2D("hCostPhiHE_4pt5_2m_gen","hCostPhiHE_4pt5_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_5pt6_2m_gen = new TH2D("hCostPhiHE_5pt6_2m_gen","hCostPhiHE_5pt6_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_6pt7_2m_gen = new TH2D("hCostPhiHE_6pt7_2m_gen","hCostPhiHE_6pt7_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_7pt8_2m_gen = new TH2D("hCostPhiHE_7pt8_2m_gen","hCostPhiHE_7pt8_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_8pt9_2m_gen = new TH2D("hCostPhiHE_8pt9_2m_gen","hCostPhiHE_8pt9_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_9pt10_2m_gen = new TH2D("hCostPhiHE_9pt10_2m_gen","hCostPhiHE_9pt10_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_10pt11_2m_gen = new TH2D("hCostPhiHE_10pt11_2m_gen","hCostPhiHE_10pt11_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_11pt12_2m_gen = new TH2D("hCostPhiHE_11pt12_2m_gen","hCostPhiHE_11pt12_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiHE_12ptinf_2m_gen = new TH2D("hCostPhiHE_12ptinf_2m_gen","hCostPhiHE_12ptinf_2m_gen",100,-1,1,50,0,PI);

  TH2D *hCostPhiCS_0pt1_2m_gen = new TH2D("hCostPhiCS_0pt1_2m_gen","hCostPhiCS_0pt1_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_1pt2_2m_gen = new TH2D("hCostPhiCS_1pt2_2m_gen","hCostPhiCS_1pt2_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_2pt3_2m_gen = new TH2D("hCostPhiCS_2pt3_2m_gen","hCostPhiCS_2pt3_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_3pt4_2m_gen = new TH2D("hCostPhiCS_3pt4_2m_gen","hCostPhiCS_3pt4_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_4pt5_2m_gen = new TH2D("hCostPhiCS_4pt5_2m_gen","hCostPhiCS_4pt5_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_5pt6_2m_gen = new TH2D("hCostPhiCS_5pt6_2m_gen","hCostPhiCS_5pt6_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_6pt7_2m_gen = new TH2D("hCostPhiCS_6pt7_2m_gen","hCostPhiCS_6pt7_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_7pt8_2m_gen = new TH2D("hCostPhiCS_7pt8_2m_gen","hCostPhiCS_7pt8_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_8pt9_2m_gen = new TH2D("hCostPhiCS_8pt9_2m_gen","hCostPhiCS_8pt9_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_9pt10_2m_gen = new TH2D("hCostPhiCS_9pt10_2m_gen","hCostPhiCS_9pt10_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_10pt11_2m_gen = new TH2D("hCostPhiCS_10pt11_2m_gen","hCostPhiCS_10pt11_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_11pt12_2m_gen = new TH2D("hCostPhiCS_11pt12_2m_gen","hCostPhiCS_11pt12_2m_gen",100,-1,1,50,0,PI);
  TH2D *hCostPhiCS_12ptinf_2m_gen = new TH2D("hCostPhiCS_12ptinf_2m_gen","hCostPhiCS_12ptinf_2m_gen",100,-1,1,50,0,PI);

  TH3D *hCostPhiMassHE_0pt1_2m_rec = new TH3D("hCostPhiMassHE_0pt1_2m_rec","hCostPhiMassHE_0pt1_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_1pt2_2m_rec = new TH3D("hCostPhiMassHE_1pt2_2m_rec","hCostPhiMassHE_1pt2_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_2pt3_2m_rec = new TH3D("hCostPhiMassHE_2pt3_2m_rec","hCostPhiMassHE_2pt3_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_3pt4_2m_rec = new TH3D("hCostPhiMassHE_3pt4_2m_rec","hCostPhiMassHE_3pt4_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_4pt5_2m_rec = new TH3D("hCostPhiMassHE_4pt5_2m_rec","hCostPhiMassHE_4pt5_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_5pt6_2m_rec = new TH3D("hCostPhiMassHE_5pt6_2m_rec","hCostPhiMassHE_5pt6_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_6pt7_2m_rec = new TH3D("hCostPhiMassHE_6pt7_2m_rec","hCostPhiMassHE_6pt7_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_7pt8_2m_rec = new TH3D("hCostPhiMassHE_7pt8_2m_rec","hCostPhiMassHE_7pt8_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_8pt9_2m_rec = new TH3D("hCostPhiMassHE_8pt9_2m_rec","hCostPhiMassHE_8pt9_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_9pt10_2m_rec = new TH3D("hCostPhiMassHE_9pt10_2m_rec","hCostPhiMassHE_9pt10_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_10pt11_2m_rec = new TH3D("hCostPhiMassHE_10pt11_2m_rec","hCostPhiMassHE_10pt11_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_11pt12_2m_rec = new TH3D("hCostPhiMassHE_11pt12_2m_rec","hCostPhiMassHE_11pt12_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassHE_12ptinf_2m_rec = new TH3D("hCostPhiMassHE_12ptinf_2m_rec","hCostPhiMassHE_12ptinf_2m_rec",100,-1,1,50,0,PI,120,2,5);

  TH3D *hCostPhiMassCS_0pt1_2m_rec = new TH3D("hCostPhiMassCS_0pt1_2m_rec","hCostPhiMassCS_0pt1_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_1pt2_2m_rec = new TH3D("hCostPhiMassCS_1pt2_2m_rec","hCostPhiMassCS_1pt2_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_2pt3_2m_rec = new TH3D("hCostPhiMassCS_2pt3_2m_rec","hCostPhiMassCS_2pt3_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_3pt4_2m_rec = new TH3D("hCostPhiMassCS_3pt4_2m_rec","hCostPhiMassCS_3pt4_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_4pt5_2m_rec = new TH3D("hCostPhiMassCS_4pt5_2m_rec","hCostPhiMassCS_4pt5_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_5pt6_2m_rec = new TH3D("hCostPhiMassCS_5pt6_2m_rec","hCostPhiMassCS_5pt6_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_6pt7_2m_rec = new TH3D("hCostPhiMassCS_6pt7_2m_rec","hCostPhiMassCS_6pt7_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_7pt8_2m_rec = new TH3D("hCostPhiMassCS_7pt8_2m_rec","hCostPhiMassCS_7pt8_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_8pt9_2m_rec = new TH3D("hCostPhiMassCS_8pt9_2m_rec","hCostPhiMassCS_8pt9_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_9pt10_2m_rec = new TH3D("hCostPhiMassCS_9pt10_2m_rec","hCostPhiMassCS_9pt10_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_10pt11_2m_rec = new TH3D("hCostPhiMassCS_10pt11_2m_rec","hCostPhiMassCS_10pt11_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_11pt12_2m_rec = new TH3D("hCostPhiMassCS_11pt12_2m_rec","hCostPhiMassCS_11pt12_2m_rec",100,-1,1,50,0,PI,120,2,5);
  TH3D *hCostPhiMassCS_12ptinf_2m_rec = new TH3D("hCostPhiMassCS_12ptinf_2m_rec","hCostPhiMassCS_12ptinf_2m_rec",100,-1,1,50,0,PI,120,2,5);

  //============================================================================
  printf("---> Setting the chain ... \n");
  //============================================================================
  TChain *chain = new TChain("MCTree");
  Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;
  if(gSystem -> GetPathInfo(FILE_NAME_IN,dummy1,dummy2,dummy3,dummy4) != 0){
    printf("NOT EXISTING FILE\n");
    return;
  }

  printf("Opening %s\n",FILE_NAME_IN);
  chain -> Add(FILE_NAME_IN);

  //============================================================================
  printf("---> Setting tree variable (gen) ... \n");
  //============================================================================
  Int_t NDimu_gen;
  Double_t DimuPt_gen[3000], DimuY_gen[3000];
  Double_t CostHE_gen[3000], PhiHE_gen[3000], CostCS_gen[3000], PhiCS_gen[3000];

  chain -> SetBranchAddress("NDimu_gen",&NDimu_gen);
  chain -> SetBranchAddress("DimuPt_gen",DimuPt_gen);
  chain -> SetBranchAddress("DimuY_gen",DimuY_gen);
  chain -> SetBranchAddress("CostHE_gen",CostHE_gen);
  chain -> SetBranchAddress("PhiHE_gen",PhiHE_gen);
  chain -> SetBranchAddress("CostCS_gen",CostCS_gen);
  chain -> SetBranchAddress("PhiCS_gen",PhiCS_gen);

  //============================================================================
  printf("---> Filling TH2D (gen) ... \n");
  //============================================================================
  Int_t NEntries = chain -> GetEntries();
  for(int i = 0;i < NEntries;i++){
    printf("GENERATED %i -> %i : %3.2f%\r",NEntries,i,(double) i/NEntries*100);
    chain -> GetEntry(i);
    for(int k = 0;k < NDimu_gen;k++){
      if(DimuY_gen[k] > -4. && DimuY_gen[k] < -2.5){
          if(DimuPt_gen[k] > 0 && DimuPt_gen[k] <= 1){hCostPhiHE_0pt1_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_0pt1_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 1 && DimuPt_gen[k] <= 2){hCostPhiHE_1pt2_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_1pt2_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 2 && DimuPt_gen[k] <= 3){hCostPhiHE_2pt3_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_2pt3_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 3 && DimuPt_gen[k] <= 4){hCostPhiHE_3pt4_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_3pt4_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 4 && DimuPt_gen[k] <= 5){hCostPhiHE_4pt5_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_4pt5_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 5 && DimuPt_gen[k] <= 6){hCostPhiHE_5pt6_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_5pt6_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 6 && DimuPt_gen[k] <= 7){hCostPhiHE_6pt7_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_6pt7_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 7 && DimuPt_gen[k] <= 8){hCostPhiHE_7pt8_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_7pt8_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 8 && DimuPt_gen[k] <= 9){hCostPhiHE_8pt9_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_8pt9_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 9 && DimuPt_gen[k] <= 10){hCostPhiHE_9pt10_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_9pt10_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 10 && DimuPt_gen[k] <= 11){hCostPhiHE_10pt11_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_10pt11_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 11 && DimuPt_gen[k] <= 12){hCostPhiHE_11pt12_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_11pt12_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
          if(DimuPt_gen[k] > 12){hCostPhiHE_12ptinf_2m_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k])); hCostPhiCS_12ptinf_2m_gen -> Fill(CostCS_gen[k],TMath::Abs(PhiCS_gen[k]));}
        }
      }
    }

    //==========================================================================
    printf("---> Setting tree variable (rec) ... \n");
    //==========================================================================
    Int_t NDimu_rec;
    Double_t DimuPt_rec[3000], DimuY_rec[3000];
    Double_t DimuMass_rec[3000];
    Int_t DimuMatch_rec[3000];
    Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];
    //Bool_t IsPhysSelected;

    //chain -> SetBranchAddress("FiredTriggerClasses",TrigClass); //???
    chain -> SetBranchAddress("NDimu_rec",&NDimu_rec);
    chain -> SetBranchAddress("DimuPt_rec",DimuPt_rec);
    chain -> SetBranchAddress("DimuY_rec",DimuY_rec);
    chain -> SetBranchAddress("DimuMass_rec",DimuMass_rec);
    chain -> SetBranchAddress("DimuMatch_rec",DimuMatch_rec);
    chain -> SetBranchAddress("CostHE_rec",CostHE_rec);
    chain -> SetBranchAddress("PhiHE_rec",PhiHE_rec);
    chain -> SetBranchAddress("CostCS_rec",CostCS_rec);
    chain -> SetBranchAddress("PhiCS_rec",PhiCS_rec);
    //chain -> SetBranchAddress("IsPhysSelected",&IsPhysSelected); //???

    //==========================================================================
    printf("---> Filling TH3D (rec) ... \n");
    //==========================================================================
    NEntries = chain -> GetEntries();
    for(int i = 0;i < NEntries;i++){
      printf("RECONSTRUCTED %i -> %i : %3.2f%\r",NEntries,i,(double) i/NEntries*100);
      chain -> GetEntry(i);
      for(int k = 0;k < NDimu_rec;k++){

        //if(IsPhysSelected){
          //TString Trigger = TrigClass;
          //Bool_t TriggerSelected = kFALSE;
          //if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
          if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
            //if(TriggerSelected){
              if(DimuMatch_rec[k] == 2){
                if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
                  if(DimuPt_rec[k] > 0 && DimuPt_rec[k] <= 1){hCostPhiMassHE_0pt1_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_0pt1_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 1 && DimuPt_rec[k] <= 2){hCostPhiMassHE_1pt2_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_1pt2_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 3){hCostPhiMassHE_2pt3_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_2pt3_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 3 && DimuPt_rec[k] <= 4){hCostPhiMassHE_3pt4_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_3pt4_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 4 && DimuPt_rec[k] <= 5){hCostPhiMassHE_4pt5_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_4pt5_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 5 && DimuPt_rec[k] <= 6){hCostPhiMassHE_5pt6_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_5pt6_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 6 && DimuPt_rec[k] <= 7){hCostPhiMassHE_6pt7_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_6pt7_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 7 && DimuPt_rec[k] <= 8){hCostPhiMassHE_7pt8_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_7pt8_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 8 && DimuPt_rec[k] <= 9){hCostPhiMassHE_8pt9_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_8pt9_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 9 && DimuPt_rec[k] <= 10){hCostPhiMassHE_9pt10_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_9pt10_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 10 && DimuPt_rec[k] <= 11){hCostPhiMassHE_10pt11_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_10pt11_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 11 && DimuPt_rec[k] <= 12){hCostPhiMassHE_11pt12_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_11pt12_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                  if(DimuPt_rec[k] > 12){hCostPhiMassHE_12ptinf_2m_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),DimuMass_rec[k]); hCostPhiMassCS_12ptinf_2m_rec -> Fill(CostCS_rec[k],TMath::Abs(PhiCS_rec[k]),DimuMass_rec[k]);}
                }
              }
            //}
          }
        //}
      }
    }

    //==========================================================================
    printf("---> Saving histos ... \n");
    //==========================================================================
    char FILE_NAME_OUT[400];
    sprintf(FILE_NAME_OUT,"~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/HistosFromOfficialTree_Jpsi_PbPb_Nopol_TH3rec.root");
    TFile *file_out = new TFile(FILE_NAME_OUT,"RECREATE");
    file_out -> cd();
    hCostPhiHE_0pt1_2m_gen -> Write();
    hCostPhiHE_1pt2_2m_gen -> Write();
    hCostPhiHE_2pt3_2m_gen -> Write();
    hCostPhiHE_3pt4_2m_gen -> Write();
    hCostPhiHE_4pt5_2m_gen -> Write();
    hCostPhiHE_5pt6_2m_gen -> Write();
    hCostPhiHE_6pt7_2m_gen -> Write();
    hCostPhiHE_7pt8_2m_gen -> Write();
    hCostPhiHE_8pt9_2m_gen -> Write();
    hCostPhiHE_9pt10_2m_gen -> Write();
    hCostPhiHE_10pt11_2m_gen -> Write();
    hCostPhiHE_11pt12_2m_gen -> Write();
    hCostPhiHE_12ptinf_2m_gen -> Write();

    hCostPhiCS_0pt1_2m_gen -> Write();
    hCostPhiCS_1pt2_2m_gen -> Write();
    hCostPhiCS_2pt3_2m_gen -> Write();
    hCostPhiCS_3pt4_2m_gen -> Write();
    hCostPhiCS_4pt5_2m_gen -> Write();
    hCostPhiCS_5pt6_2m_gen -> Write();
    hCostPhiCS_6pt7_2m_gen -> Write();
    hCostPhiCS_7pt8_2m_gen -> Write();
    hCostPhiCS_8pt9_2m_gen -> Write();
    hCostPhiCS_9pt10_2m_gen -> Write();
    hCostPhiCS_10pt11_2m_gen -> Write();
    hCostPhiCS_11pt12_2m_gen -> Write();
    hCostPhiCS_12ptinf_2m_gen -> Write();

    hCostPhiMassHE_0pt1_2m_rec -> Write();
    hCostPhiMassHE_1pt2_2m_rec -> Write();
    hCostPhiMassHE_2pt3_2m_rec -> Write();
    hCostPhiMassHE_3pt4_2m_rec -> Write();
    hCostPhiMassHE_4pt5_2m_rec -> Write();
    hCostPhiMassHE_5pt6_2m_rec -> Write();
    hCostPhiMassHE_6pt7_2m_rec -> Write();
    hCostPhiMassHE_7pt8_2m_rec -> Write();
    hCostPhiMassHE_8pt9_2m_rec -> Write();
    hCostPhiMassHE_9pt10_2m_rec -> Write();
    hCostPhiMassHE_10pt11_2m_rec -> Write();
    hCostPhiMassHE_11pt12_2m_rec -> Write();
    hCostPhiMassHE_12ptinf_2m_rec -> Write();

    hCostPhiMassCS_0pt1_2m_rec -> Write();
    hCostPhiMassCS_1pt2_2m_rec -> Write();
    hCostPhiMassCS_2pt3_2m_rec -> Write();
    hCostPhiMassCS_3pt4_2m_rec -> Write();
    hCostPhiMassCS_4pt5_2m_rec -> Write();
    hCostPhiMassCS_5pt6_2m_rec -> Write();
    hCostPhiMassCS_6pt7_2m_rec -> Write();
    hCostPhiMassCS_7pt8_2m_rec -> Write();
    hCostPhiMassCS_8pt9_2m_rec -> Write();
    hCostPhiMassCS_9pt10_2m_rec -> Write();
    hCostPhiMassCS_10pt11_2m_rec -> Write();
    hCostPhiMassCS_11pt12_2m_rec -> Write();
    hCostPhiMassCS_12ptinf_2m_rec -> Write();
}
