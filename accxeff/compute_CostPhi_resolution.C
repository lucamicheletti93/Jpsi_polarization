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

void compute_CostPhi_resolution(){
  //============================================================================
  printf("---> Setting main quantities ... \n");
  //============================================================================
  gStyle -> SetOptStat(0);
  string fileTreeName = "~/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root";
  TFile *fileTree = new TFile(fileTreeName.c_str(),"READ");

  TTree *tree = (TTree*) fileTree -> Get("MCTree");

  Int_t NDimu_gen;
  Double_t DimuPt_gen[3000], DimuY_gen[3000];
  Double_t CostHE_gen[3000], PhiHE_gen[3000], CostCS_gen[3000], PhiCS_gen[3000];

  tree -> SetBranchAddress("NDimu_gen",&NDimu_gen);
  tree -> SetBranchAddress("DimuPt_gen",DimuPt_gen);
  tree -> SetBranchAddress("DimuY_gen",DimuY_gen);
  tree -> SetBranchAddress("CostHE_gen",CostHE_gen);
  tree -> SetBranchAddress("PhiHE_gen",PhiHE_gen);
  tree -> SetBranchAddress("CostCS_gen",CostCS_gen);
  tree -> SetBranchAddress("PhiCS_gen",PhiCS_gen);

  Int_t NDimu_rec;
  Double_t DimuPt_rec[3000], DimuY_rec[3000];
  Double_t DimuMass_rec[3000];
  Int_t DimuMatch_rec[3000];
  Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];

  tree -> SetBranchAddress("NDimu_rec",&NDimu_rec);
  tree -> SetBranchAddress("DimuPt_rec",DimuPt_rec);
  tree -> SetBranchAddress("DimuY_rec",DimuY_rec);
  tree -> SetBranchAddress("DimuMass_rec",DimuMass_rec);
  tree -> SetBranchAddress("DimuMatch_rec",DimuMatch_rec);
  tree -> SetBranchAddress("CostHE_rec",CostHE_rec);
  tree -> SetBranchAddress("PhiHE_rec",PhiHE_rec);
  tree -> SetBranchAddress("CostCS_rec",CostCS_rec);
  tree -> SetBranchAddress("PhiCS_rec",PhiCS_rec);

  double PI = TMath::Pi();
  TH1D *histCostResolution = new TH1D("histCostResolution","cos#theta^{GEN} - cos#theta^{REC} (p_{T}>2 GeV/#it{c})",100,-0.5,0.5);
  TH1D *histPhiResolution = new TH1D("histPhiResolution","#varphi^{GEN} - #varphi^{REC} (p_{T}>2 GeV/#it{c})",100,-0.5,0.5);
  TH2D *histCostPhiResolution = new TH2D("histCostPhiResolution","p_{T}>2 GeV/#it{c}",100,-0.5,0.5,100,-0.5,0.5);

  TH1D *histCostResolution2pt4 = new TH1D("histCostResolution2pt4","cos#theta^{GEN} - cos#theta^{REC} (2<p_{T}<4 GeV/#it{c})",100,-0.5,0.5);
  TH1D *histPhiResolution2pt4 = new TH1D("histPhiResolution2pt4","#varphi^{GEN} - #varphi^{REC} (2<p_{T}<4 GeV/#it{c})",100,-0.5,0.5);
  TH2D *histCostPhiResolution2pt4 = new TH2D("histCostPhiResolution2pt4","2<p_{T}<4 GeV/#it{c}",100,-0.5,0.5,100,-0.5,0.5);

  TH1D *histCostResolution4pt6 = new TH1D("histCostResolution4pt6","cos#theta^{GEN} - cos#theta^{REC} (4<p_{T}<6 GeV/#it{c})",100,-0.5,0.5);
  TH1D *histPhiResolution4pt6 = new TH1D("histPhiResolution4pt6","#varphi^{GEN} - #varphi^{REC} (4<p_{T}<6 GeV/#it{c})",100,-0.5,0.5);
  TH2D *histCostPhiResolution4pt6 = new TH2D("histCostPhiResolution4pt6","4<p_{T}<6 GeV/#it{c}",100,-0.5,0.5,100,-0.5,0.5);

  TH1D *histCostResolution6pt10 = new TH1D("histCostResolution6pt10","cos#theta^{GEN} - cos#theta^{REC} (6<p_{T}<10 GeV/#it{c})",100,-0.5,0.5);
  TH1D *histPhiResolution6pt10 = new TH1D("histPhiResolution6pt10","#varphi^{GEN} - #varphi^{REC} (6<p_{T}<10 GeV/#it{c})",100,-0.5,0.5);
  TH2D *histCostPhiResolution6pt10 = new TH2D("histCostPhiResolution6pt10","6<p_{T}<10 GeV/#it{c}",100,-0.5,0.5,100,-0.5,0.5);

  double deltaCost = 0;
  double deltaPhi = 0;

  TH1D *histCostGen = new TH1D("histCostGen","",100,-1,1);
  TH1D *histCostRec = new TH1D("histCostRec","",100,-1,1);
  TH1D *histPhiGen = new TH1D("histPhiGen","",100,-PI,PI);
  TH1D *histPhiRec = new TH1D("histPhiRec","",100,-PI,PI);

  Int_t NEntries = tree -> GetEntries();
  for(int i = 0;i < NEntries;i++){
    printf("GENERATED %i -> %i : %3.2f%\r",NEntries,i,(double) i/NEntries*100);
    tree -> GetEntry(i);
      if(NDimu_gen == 1 && NDimu_rec == 1){
        if(DimuY_gen[0] > -4. && DimuY_gen[0] < -2.5 && DimuY_rec[0] > -4. && DimuY_rec[0] < -2.5){
          if(DimuMatch_rec[0] == 2){
            if(DimuMass_rec[0] > 2 && DimuMass_rec[0] < 5 && DimuPt_rec[0] > 2){
              deltaCost = CostHE_rec[0] - CostHE_gen[0];
              deltaPhi = PhiHE_rec[0] - PhiHE_gen[0];
              //printf("cost gen = %f , cost rec = %f -> delta = %f |phi gen = %f , phi rec = %f -> delta = %f \n",CostHE_gen[0],CostHE_rec[0],deltaCost,PhiHE_gen[0],PhiHE_rec[0],deltaPhi);
              //printf("phi gen = %f , phi rec = %f -> delta = %f \n",PhiHE_gen[0],PhiHE_rec[0],deltaPhi);

              if(DimuPt_rec[0] > 2 && DimuPt_rec[0] <= 4){
                histCostResolution2pt4 -> Fill(deltaCost);
                histPhiResolution2pt4 -> Fill(deltaPhi);
                histCostPhiResolution2pt4 -> Fill(deltaCost,deltaPhi);
              }
              if(DimuPt_rec[0] > 4 && DimuPt_rec[0] <= 6){
                histCostResolution4pt6 -> Fill(deltaCost);
                histPhiResolution4pt6 -> Fill(deltaPhi);
                histCostPhiResolution4pt6 -> Fill(deltaCost,deltaPhi);
              }
              if(DimuPt_rec[0] > 6 && DimuPt_rec[0] <= 10){
                histCostResolution6pt10 -> Fill(deltaCost);
                histPhiResolution6pt10 -> Fill(deltaPhi);
                histCostPhiResolution6pt10 -> Fill(deltaCost,deltaPhi);
              }

              histCostGen -> Fill(CostHE_gen[0]);
              histCostRec -> Fill(CostHE_rec[0]);
              histPhiGen -> Fill(PhiHE_gen[0]);
              histPhiRec -> Fill(PhiHE_rec[0]);

              histCostPhiResolution -> Fill(deltaCost,deltaPhi);
              histCostResolution -> Fill(deltaCost);
              histPhiResolution -> Fill(deltaPhi);
            }
          }
        }
      }
    }
    printf("\n");

    TCanvas *cHistCostPhi = new TCanvas("cHistCostPhi","cHistCostPhi",20,20,600,600);
    cHistCostPhi -> Divide(2,2);
    cHistCostPhi -> cd(1);
    histCostGen -> Draw();
    cHistCostPhi -> cd(2);
    histCostRec -> Draw();
    cHistCostPhi -> cd(3);
    histPhiGen -> Draw();
    cHistCostPhi -> cd(4);
    histPhiRec -> Draw();

    TCanvas *cHistCostResolution = new TCanvas("cHistCostResolution","cHistCostResolution",20,20,600,600);
    cHistCostResolution -> Divide(2,2);
    cHistCostResolution -> cd(1);
    histCostResolution2pt4 -> Draw();
    cHistCostResolution -> cd(2);
    histCostResolution4pt6 -> Draw();
    cHistCostResolution -> cd(3);
    histCostResolution6pt10 -> Draw();
    cHistCostResolution -> cd(4);
    histCostResolution -> Draw();
    //printf("Resolution Cost = %f \n",histCostResolution -> GetRMS());
    //histCostResolution -> Fit("gaus");

    printf("============================================ \n");
    printf("pt > 2  GeV/c RMS Cost = %f \n",histCostResolution -> GetRMS());
    printf("2 < pt < 4  GeV/c RMS Cost = %f \n",histCostResolution2pt4 -> GetRMS());
    printf("4 < pt < 6  GeV/c RMS Cost = %f \n",histCostResolution4pt6 -> GetRMS());
    printf("6 < pt < 10 GeV/c RMS Cost = %f \n",histCostResolution6pt10 -> GetRMS());

    TCanvas *cHistPhiResolution = new TCanvas("cHistPhiResolution","cHistPhiResolution",20,20,600,600);
    cHistPhiResolution -> Divide(2,2);
    cHistPhiResolution -> cd(1);
    histPhiResolution2pt4 -> Draw();
    cHistPhiResolution -> cd(2);
    histPhiResolution4pt6 -> Draw();
    cHistPhiResolution -> cd(3);
    histPhiResolution6pt10 -> Draw();
    cHistPhiResolution -> cd(4);
    histPhiResolution -> Draw();
    //printf("Resolution Phi = %f \n",histPhiResolution -> GetRMS());
    //histPhiResolution -> Fit("gaus");

    printf("============================================ \n");
    printf("pt > 2  GeV/c RMS Phi = %f \n",histPhiResolution -> GetRMS());
    printf("2 < pt < 4  GeV/c RMS Phi = %f \n",histPhiResolution2pt4 -> GetRMS());
    printf("4 < pt < 6  GeV/c RMS Phi = %f \n",histPhiResolution4pt6 -> GetRMS());
    printf("6 < pt < 10 GeV/c RMS Phi = %f \n",histPhiResolution6pt10 -> GetRMS());

    TCanvas *cHistCostPhiResolution = new TCanvas("cHistCostPhiResolution","cHistCostPhiResolution",20,20,600,600);
    cHistCostPhiResolution -> Divide(2,2);
    cHistCostPhiResolution -> cd(1);
    histCostPhiResolution2pt4 -> Draw("COLZ");
    cHistCostPhiResolution -> cd(2);
    histCostPhiResolution4pt6 -> Draw("COLZ");
    cHistCostPhiResolution -> cd(3);
    histCostPhiResolution6pt10 -> Draw("COLZ");
    cHistCostPhiResolution -> cd(4);
    histCostPhiResolution -> Draw("COLZ");

    cHistCostResolution -> SaveAs("/home/luca/Scrivania/CostResolution.png");
    cHistPhiResolution -> SaveAs("/home/luca/Scrivania/PhiResolution.png");
    cHistCostPhiResolution -> SaveAs("/home/luca/Scrivania/CostPhiResolution.png");


}
