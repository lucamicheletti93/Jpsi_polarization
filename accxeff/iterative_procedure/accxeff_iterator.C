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
#include <TF2.h>
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

double Func_W(double *, double *);

void accxeff_iterator(int ptMin, int ptMax){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  double PI = TMath::Pi();

  ostringstream convertPtMin;
  convertPtMin << ptMin;
  string strPtMin =  convertPtMin.str();

  ostringstream convertPtMax;
  convertPtMax << ptMax;
  string strPtMax =  convertPtMax.str();

  string dataset = strPtMin + "pt" + strPtMax;
  string nameOption = "_test";
  string outputNameOption = "_extreme_value";

  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + nameOption + ".root";
  string fileTreeName = "~/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root";
  string fileAccxEffName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/NEW_GIT_OUTPUT/accxeff_" + dataset + nameOption + ".root";
  string fileOutputName = "~/GITHUB/Jpsi_polarization/2D_approach/data_analysis/accxeff/iterative_procedure/" + dataset + outputNameOption + ".root";

  TFile *fileBinning = new TFile(fileBinningName.c_str(),"READ");
  TFile *fileTree = new TFile(fileTreeName.c_str(),"READ");
  TFile *fileAccxEff = new TFile(fileAccxEffName.c_str(),"READ");

  //============================================================================
  printf("---> Setting histograms bins ... \n");
  //============================================================================
  vector <double> CostValues;
  vector <double> PhiValues;
  vector <int> CostBinsMin, CostBinsMax;
  vector <int> PhiBinsMin, PhiBinsMax;
  vector < vector <double> > CellAreaMatrix;

  Binning *binning = (Binning*) fileBinning -> Get("Binning");
  CostValues = binning -> GetCostValues();
  CostBinsMin = binning -> GetCostBinsMin();
  CostBinsMax = binning -> GetCostBinsMax();
  const int NCostBins = CostValues.size() - 1;
  PhiValues = binning -> GetPhiValues();
  PhiBinsMin = binning -> GetPhiBinsMin();
  PhiBinsMax = binning -> GetPhiBinsMax();
  const int NPhiBins = PhiValues.size() - 1;

  CellAreaMatrix = binning -> GetCellAreaMatrix();

  //============================================================================
  printf("---> Inizializing the tree ... \n");
  //============================================================================
  Int_t NDimu_gen;
  Double_t DimuPt_gen[3000], DimuY_gen[3000];
  Double_t CostHE_gen[3000], PhiHE_gen[3000], CostCS_gen[3000], PhiCS_gen[3000];

  Int_t NDimu_rec;
  Double_t DimuPt_rec[3000], DimuY_rec[3000];
  Double_t DimuMass_rec[3000];
  Int_t DimuMatch_rec[3000];
  Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];

  TTree *tree = (TTree*) fileTree -> Get("MCTree");
  tree -> SetBranchAddress("NDimu_gen",&NDimu_gen);
  tree -> SetBranchAddress("DimuPt_gen",DimuPt_gen);
  tree -> SetBranchAddress("DimuY_gen",DimuY_gen);
  tree -> SetBranchAddress("CostHE_gen",CostHE_gen);
  tree -> SetBranchAddress("PhiHE_gen",PhiHE_gen);
  tree -> SetBranchAddress("CostCS_gen",CostCS_gen);
  tree -> SetBranchAddress("PhiCS_gen",PhiCS_gen);

  tree -> SetBranchAddress("NDimu_rec",&NDimu_rec);
  tree -> SetBranchAddress("DimuPt_rec",DimuPt_rec);
  tree -> SetBranchAddress("DimuY_rec",DimuY_rec);
  tree -> SetBranchAddress("DimuMass_rec",DimuMass_rec);
  tree -> SetBranchAddress("DimuMatch_rec",DimuMatch_rec);
  tree -> SetBranchAddress("CostHE_rec",CostHE_rec);
  tree -> SetBranchAddress("PhiHE_rec",PhiHE_rec);
  tree -> SetBranchAddress("CostCS_rec",CostCS_rec);
  tree -> SetBranchAddress("PhiCS_rec",PhiCS_rec);

  //============================================================================
  printf("---> Weighting the input gen and rec distributions ... \n");
  //============================================================================
  int countDimuRec = 0;

  double weight2D = 0;
  double lambdaThInput;
  //if(ptMin == 6 && ptMax == 10){lambdaThInput = -0.201;} // stage 1
  //if(ptMin == 6 && ptMax == 10){lambdaThInput = -0.209;} // stage 2
  //if(ptMin == 6 && ptMax == 10){lambdaThInput = -0.379;} // stage 1 free sigma
  //if(ptMin == 6 && ptMax == 10){lambdaThInput = -0.398;} // stage 2 free sigma
  if(ptMin == 6 && ptMax == 10){lambdaThInput = -1.;} // stage 3 (extreme value)


  double nEntries = tree -> GetEntries();
  double percentage = 0;

  printf("Lambda Theta = %2.1f \n",lambdaThInput);

  TH2D *histCostPhiGen = new TH2D("histCostPhiGen","",100,-1,1,50,0,PI);
  TH2D *histCostPhiGenRebin = new TH2D("histCostPhiGenRebin","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  TH2D *histCostPhiRec = new TH2D("histCostPhiRec","",100,-1,1,50,0,PI);
  TH2D *histCostPhiRecRebin = new TH2D("histCostPhiRecRebin","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);

  TF2 *func2DVarPolWeight = new TF2("func2DVarPolWeight",Func_W,-1,1,0,PI,4);
  func2DVarPolWeight -> SetParameters(1000,lambdaThInput,0,0);

  for(int ii = 0;ii < nEntries;ii++){
    tree -> GetEntry(ii);
    percentage = ((double) ii)/((double) nEntries)*100;
    printf("work in progress %2.1f \r",percentage);

    for(int k = 0;k < NDimu_gen;k++){
      if(DimuPt_gen[k] > ptMin && DimuPt_gen[k] <= ptMax){
        weight2D = func2DVarPolWeight -> Eval(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]))/func2DVarPolWeight -> GetMaximum();
        histCostPhiGen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight2D);
        histCostPhiGenRebin -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight2D);
      }
    }

    for(int k = 0;k < NDimu_rec;k++){
      if(DimuPt_rec[k] > ptMin && DimuPt_rec[k] <= ptMax){
        if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
          if(DimuMatch_rec[k] == 2){
            if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
              histCostPhiRec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight2D);
              histCostPhiRecRebin -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight2D);
              //countDimuRec++;
            }
          }
        }
      }
    }
    //eventIndex++;
  //}
  //countDimuRec = 0;
  //eventIndex = 0;
  }
  printf("\n");

  //============================================================================
  printf("---> Computing AccxEff ... \n");
  //============================================================================
  histCostPhiRec -> Sumw2();
  histCostPhiGen -> Sumw2();
  histCostPhiRecRebin -> Sumw2();
  histCostPhiGenRebin -> Sumw2();

  TH2D *hist_accxeff_HE = new TH2D("hist_accxeff_HE","",100,-1,1,50,0,PI);
  hist_accxeff_HE -> Divide(histCostPhiRec,histCostPhiGen,1,1,"B");
  TH2D *hist_accxeff_HE_rebin = new TH2D("hist_accxeff_HE_rebin","",NCostBins,&CostValues[0],NPhiBins,&PhiValues[0]);
  hist_accxeff_HE_rebin -> Divide(histCostPhiRecRebin,histCostPhiGenRebin,1,1,"B");

  //============================================================================
  printf("---> Saving results ... \n");
  //============================================================================
  TFile *fileOutput = new TFile(fileOutputName.c_str(),"RECREATE");
  histCostPhiGen -> Write();
  histCostPhiGenRebin -> Write();
  histCostPhiRec -> Write();
  histCostPhiRecRebin -> Write();
  hist_accxeff_HE -> Write();
  hist_accxeff_HE_rebin -> Write();

  TCanvas *cRecPlots = new TCanvas("cRecPlots","cRecPlots",20,20,600,600);
  histCostPhiGen -> Draw("COLZ");

  fileOutput -> Close();

}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double Func_W(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];
  double L_thph = par[3];

  double costh = x[0];
  double phi = x[1];

  double cosph = TMath::Cos(phi);
  double cos2ph = TMath::Cos(2*phi);

  double W =  (N/(3 + L_th))*(1 + (L_th + L_ph*cos2ph)*costh*costh + 2*L_thph*costh*cosph*TMath::Sqrt(1 - costh*costh) + L_ph*cos2ph);
  return W;
}
