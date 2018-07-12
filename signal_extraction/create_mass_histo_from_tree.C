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

void create_mass_histo_from_tree(int pt_min, int pt_max){
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
  //string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";
  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + "_test.root";
  TFile *fileTree = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/TREE_FILTERED/TreeFiltered_full_statistics.root","READ");

  //============================================================================
  printf("---> Defining the output tree and histo ... \n");
  //============================================================================
  //string fileOutputName = dataset + ".root";
  string fileOutputName = dataset + "_test.root";
  TFile *fileOutput = new TFile(fileOutputName.c_str(),"RECREATE");
  Double_t DimuMass;

  if(pt_min == 0 && pt_max == 12){
    //============================================================================
    printf("---> Creating the integrated distributions ... \n");
    //============================================================================
    TTree *treeIntegrated = new TTree("treeIntegrated","treeIntegrated");
    treeIntegrated -> Branch("DimuMass",&DimuMass,"DimuMass/D");
    TH1D *histIntegrated = new TH1D("histIntegrated","histIntegrated",120,2,5);
    char treeName[100];

    double ptTmp = pt_min;
    while(ptTmp < pt_max){
      cout << ptTmp << " - " << ptTmp+1 << endl;

      sprintf(treeName,"treeHE_%ipt%i",(int) ptTmp,(int) ptTmp+1);
      Double_t DimuMass_unb;
      Double_t CostHE_unb, PhiHE_unb;

      TTree *tree = (TTree*) fileTree -> Get(treeName);
      tree -> SetBranchAddress("DimuMass_unb",&DimuMass_unb);
      tree -> SetBranchAddress("CostHE_unb",&CostHE_unb);
      tree -> SetBranchAddress("PhiHE_unb",&PhiHE_unb);

      for(int i = 0;i < tree -> GetEntries();i++){
        tree -> GetEntry(i);
        DimuMass = DimuMass_unb;
        treeIntegrated -> Fill();
        histIntegrated -> Fill(DimuMass);
      }
      ptTmp++;
    }

    fileOutput -> cd();
    treeIntegrated -> Write();
    histIntegrated -> Write();
    fileOutput -> Close();
    return;
  }

  //============================================================================
  printf("---> Reading the binning file ... \n");
  //============================================================================

  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  Binning *binning = (Binning*) fileBinning -> Get("Binning");

  vector <double> CostValues;
  vector <double> PhiValues;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  TTree *treeCost[NCostBins];
  TH1D *histCost[NCostBins];
  char treeCostName[100];
  for(int i = 0;i < NCostBins;i++){
    sprintf(treeCostName,"treeCost_%i",i);
    treeCost[i] = new TTree(treeCostName,treeCostName);
    treeCost[i] -> Branch("DimuMass",&DimuMass,"DimuMass/D");
    sprintf(treeCostName,"histCost_%i",i);
    histCost[i] = new TH1D(treeCostName,"",120,2,5);
  }

  TTree *treePhi[NPhiBins];
  TH1D *histPhi[NPhiBins];
  char treePhiName[100];
  for(int i = 0;i < NPhiBins;i++){
    sprintf(treePhiName,"treePhi_%i",i);
    treePhi[i] = new TTree(treePhiName,treePhiName);
    treePhi[i] -> Branch("DimuMass",&DimuMass,"DimuMass/D");
    sprintf(treePhiName,"histPhi_%i",i);
    histPhi[i] = new TH1D(treePhiName,"",120,2,5);
  }

  TTree *treeCostPhi[NCostBins][NPhiBins];
  TH1D *histCostPhi[NCostBins][NPhiBins];
  char treeCostPhiName[100];
  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
    sprintf(treeCostPhiName,"treeCost_%i_Phi_%i",i,j);
    treeCostPhi[i][j] = new TTree(treeCostPhiName,treeCostPhiName);
    treeCostPhi[i][j] -> Branch("DimuMass",&DimuMass,"DimuMass/D");
    sprintf(treeCostPhiName,"histCost_%i_Phi_%i",i,j);
    histCostPhi[i][j] = new TH1D(treeCostPhiName,"",120,2,5);
    }
  }

  //============================================================================
  printf("---> Reading the input tree and filling the output tree ... \n");
  //============================================================================
  double ptTmp = pt_min;
  while(ptTmp < pt_max){
    cout << ptTmp << " - " << ptTmp+1 << endl;

    char treeName[100];
    sprintf(treeName,"treeHE_%ipt%i",(int) ptTmp,(int) ptTmp+1);
    //TTree *tree = (TTree*) fileTree -> Get(Form("treeHE_%ipt%i",ptTmp,ptTmp+1));
    Double_t DimuMass_unb;
    Double_t CostHE_unb, PhiHE_unb;

    TTree *tree = (TTree*) fileTree -> Get(treeName);
    tree -> SetBranchAddress("DimuMass_unb",&DimuMass_unb);
    tree -> SetBranchAddress("CostHE_unb",&CostHE_unb);
    tree -> SetBranchAddress("PhiHE_unb",&PhiHE_unb);

    // Cost values
    int indexCost = 0;
    for(int i = 0;i < tree -> GetEntries();i++){
      tree -> GetEntry(i);
      while(CostHE_unb < CostValues[indexCost] || CostHE_unb >= CostValues[indexCost+1]) indexCost++;
      DimuMass = DimuMass_unb;
      if(PhiHE_unb > PhiValues[1] && PhiHE_unb < PhiValues[NPhiBins-1]){
        //treeCost[indexCost] -> Fill();
        histCost[indexCost] -> Fill(DimuMass);
      }
      indexCost = 0;
    }

    // Phi values
    int indexPhi = 0;
    for(int i = 0;i < tree -> GetEntries();i++){
      tree -> GetEntry(i);
      while(PhiHE_unb < PhiValues[indexPhi] || PhiHE_unb >= PhiValues[indexPhi+1]) indexPhi++;
      DimuMass = DimuMass_unb;
      if(CostHE_unb > CostValues[1] && CostHE_unb < CostValues[NCostBins-1]){
        //treePhi[indexPhi] -> Fill();
        histPhi[indexPhi] -> Fill(DimuMass);
      }
      indexPhi = 0;
    }

    // (Cost,Phi) values
    indexCost = 0;
    indexPhi = 0;
    for(int i = 0;i < tree -> GetEntries();i++){
      tree -> GetEntry(i);
      while(CostHE_unb < CostValues[indexCost] || CostHE_unb >= CostValues[indexCost+1]) indexCost++;
      while(PhiHE_unb < PhiValues[indexPhi] || PhiHE_unb >= PhiValues[indexPhi+1]) indexPhi++;
      DimuMass = DimuMass_unb;
      histCostPhi[indexCost][indexPhi] -> Fill(DimuMass);
      indexCost = 0;
      indexPhi = 0;
    }

    ptTmp++;
  }

  fileOutput -> cd();
  //for(int i = 0;i < NCostBins;i++){treeCost[i] -> Write(); histCost[i] -> Write();}
  //for(int i = 0;i < NPhiBins;i++){treePhi[i] -> Write(); histPhi[i] -> Write();}
  for(int i = 0;i < NCostBins;i++){histCost[i] -> Write();}
  for(int i = 0;i < NPhiBins;i++){histPhi[i] -> Write();}
  for(int i = 0;i < NCostBins;i++){
    for(int j = 0;j < NPhiBins;j++){
      histCostPhi[i][j] -> Write();
    }
  }
  fileOutput -> Close();

}
