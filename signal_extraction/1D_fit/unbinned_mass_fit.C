//#if !defined(__CINT__) || defined(__MAKECINT__)
#ifndef __CINT__
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooCBShape.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooClassFactory.h"
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

#ifndef __CINT__
#include "CB2Pdf.h"
#include "/home/luca/GITHUB/Jpsi_polarization/2D_approach/data_analysis/Binning/Binning.h"
#endif

void fit_of_minv(int, string, TTree *);
using namespace RooFit;

vector <int> CostBinsMin, CostBinsMax;
vector <int> PhiBinsMin, PhiBinsMax;

//void fit_of_minv(TTree *, int, string, TH1D *, TH1D *);
//void fit_of_minv(TTree *, string, TH1D *, TH1D *);
//void fit_of_minv(TTree *);

////////////////////////////////////////////////////////////////////////////////
// To run in compiled mode :
//      - .L CB2Pdf.cxx+
//      - .x unbinned_mass_fit.C+
////////////////////////////////////////////////////////////////////////////////

void unbinned_mass_fit(int pt_min, int pt_max){

  //#ifdef __CINT__
    //gROOT -> ProcessLineSync(".x CB2Pdf.cxx+") ;
  //#endif

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
  string fileBinningName = "~/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/Binning/binning_" + dataset + ".root";

  string outputDir = "unbinned_1D_" + dataset;
  string plotsDirectory = "mkdir " + outputDir;
  gSystem -> Exec(plotsDirectory.c_str());

  const int NCostBins = 17;
  double CostValues[18] = {-1.00,-0.80,-0.60,-0.50,-0.40,-0.30,-0.20,-0.12,-0.04,0.04,0.12,0.20,0.30,0.40,0.50,0.60,0.80,1.00};
  const int NPhiBins = 10;
  double PhiValues[11] = {0.000000,0.502655,1.005310,1.256637,1.445133,1.570796,1.696460,1.884956,2.136283,2.638938,3.141593};

  //============================================================================
  printf("---> Fitting the tree dataset ... \n");
  //============================================================================
  TFile *fileInput = new TFile(Form("../%s.root",dataset.c_str()),"READ");
  TH1D *histNJpsiCost = new TH1D("histNJpsiCost","",NCostBins,CostValues);
  TH1D *histSigmaJpsiCost = new TH1D("histSigmaJpsiCost","",NCostBins,CostValues);

  /*
  for(int i = 1;i < NCostBins-1;i++){
      TTree* tree = (TTree*) fileInput -> Get(Form("treeCost_%i",i));
      fit_of_minv(tree,i,outputDir,histNJpsiCost,histSigmaJpsiCost);
  }
  //============================================================================
  printf("---> Saving results ... \n");
  //============================================================================
  string fileJpsiName = "unbinned_1D_" + dataset + "/" + dataset + ".root";
  TFile *fileJpsi = new TFile(fileJpsiName.c_str(),"RECREATE");
  histNJpsiCost -> Write();
  histSigmaJpsiCost -> Write();
  //histNJpsiPhi -> Write();
  //histSigmaJpsiPhi -> Write();
  fileJpsi -> Close();
  */
}
//==============================================================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//==============================================================================
//void fit_of_minv(TTree *tree, int indexCost, string outputDir, TH1D *histNJpsi, TH1D *histSigmaJpsi){
//void fit_of_minv(TTree *tree, string outputDir, TH1D *histNJpsi, TH1D *histSigmaJpsi){
void fit_of_minv(int index, string ciao, TTree *albero){

  /*double PI = TMath::Pi();
  RooRealVar DimuMass("DimuMass","DimuMass",2,5);
  //DimuMass.setRange("sig",2.9,3.3);
  //DimuMass.setRange("bck_left",2,2.9);
  //DimuMass.setRange("bck_right",3.3,5);
  RooRealVar CosthHE("CosthHE","CosthHE",-1,1);
  RooRealVar PhiHE("PhiHE","PhiHE",0,PI);
  RooDataSet datasetCost("datasetCost","datasetCost",RooArgSet(DimuMass),Import(*tree));
  RooRealVar mean_CB2("#it{m}_{J/#psi}","mean CB2",3.096,2.9,3.3);
  //RooRealVar width_CB2("#it{#sigma}_{J/#psi}","width CB2",0.07,0.05,0.09);
  RooRealVar width_CB2("#it{#sigma}_{J/#psi}","width CB2",0.07,0.01,0.2);
  RooRealVar alpha1_CB2("#alpha1_{CB2}","alpha1 CB2",0.98);
  RooRealVar n1_CB2("n1_{CB2}","n1 CB2",6.97);
  RooRealVar alpha2_CB2("#alpha2_{CB2}","alpha2 CB2",1.86);
  RooRealVar n2_CB2("n2_{CB2}","n2 CB2",14.99);
  CB2Pdf sig_CB2("sig_CB2","sig_CB2",DimuMass,mean_CB2,width_CB2,alpha1_CB2,n1_CB2,alpha2_CB2,n2_CB2);
  RooRealVar sig("signal","signal yield",1.70000e+03,500,3000);

  RooRealVar mean_VWG("mean_{VWG}","mean VWG",2,-10,10);
  RooRealVar alpha_VWG("#alpha_{VWG}","alpha VWG",0.5,-10,10);
  RooRealVar beta_VWG("#beta_{VWG}","beta VWG",0.2,-10,10);
  RooGenericPdf bck_VWG("bck_VWG","bck_VWG","exp(-(@0 - @1)*(@0 - @1)/(2*(@2 + @3*((@0 - @1)/@1))*(@2 + @3*((@0 - @1)/@1))))",RooArgList(DimuMass,mean_VWG,alpha_VWG,beta_VWG));
  RooRealVar bck("background","background yield",3.45165e+04,1000,500000);

  RooAddPdf model("model","CB2 + VWG",RooArgList(sig_CB2,bck_VWG),RooArgList(sig,bck));
  model.fitTo(datasetCost,Extended(),Save());

  histNJpsi -> SetBinContent(indexCost+1,sig.getVal());
  histSigmaJpsi -> SetBinContent(indexCost+1,width_CB2.getVal());

  TCanvas *cDimuMass = new TCanvas("cDimuMass","cDimuMass",20,20,600,600);
  RooPlot *DimuMass_frame = DimuMass.frame(Title("Dimuon mass distribution"));
  datasetCost.plotOn(DimuMass_frame);
  model.plotOn(DimuMass_frame,Components("bck_VWG"),LineStyle(kDashed),LineColor(kGray+1));
  model.plotOn(DimuMass_frame,Components("sig_CB2"),DrawOption("F"),FillColor(kRed+1),FillStyle(3344));
  model.plotOn(DimuMass_frame,Components("sig_CB2"),DrawOption("L"),LineColor(kRed+1));
  model.plotOn(DimuMass_frame);
  model.paramOn(DimuMass_frame);
  DimuMass_frame -> Draw();
  //cDimuMass -> SaveAs(Form("%s/treeCost_%i.png",plotsDirectory.c_str(),i));

  cDimuMass -> SaveAs(Form("%s/HE_%icost%i.png",outputDir.c_str(),CostBinsMin[indexCost],CostBinsMax[indexCost]));*/
}
