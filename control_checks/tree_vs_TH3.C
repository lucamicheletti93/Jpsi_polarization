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
#endif

void tree_vs_TH3(){
  string ptRange = "3pt4";

  TFile *fileTH3 = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/Histos_full_statistics.root","READ");
  string histName = "hMassCostPhiHE_" + ptRange + "_2m";
  TH3D *hMassCostPhi = (TH3D*) fileTH3 -> Get(histName.c_str());
  TH1D *hTH3Mass = (TH1D*) hMassCostPhi -> ProjectionZ();

  TFile *fileTree = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/TREE_FILTERED/TreeFiltered_full_statistics.root","READ");
  string treeName = "treeHE_" + ptRange;
  TTree *tree = (TTree*) fileTree -> Get(treeName.c_str());

  Double_t DimuMass_unb;
  Double_t CostHE_unb, PhiHE_unb;

  tree -> SetBranchAddress("DimuMass_unb",&DimuMass_unb);
  tree -> SetBranchAddress("CostHE_unb",&CostHE_unb);
  tree -> SetBranchAddress("PhiHE_unb",&PhiHE_unb);

  TH1D *hTreeMass = new TH1D("hTreeMass","",120,2,5);
  hTreeMass -> SetLineColor(kRed);

  for(int i = 0;i < tree -> GetEntries();i++){
    tree -> GetEntry(i);
    hTreeMass -> Fill(DimuMass_unb);
  }

  hTH3Mass -> Draw();
  hTreeMass -> Draw("same");
}
