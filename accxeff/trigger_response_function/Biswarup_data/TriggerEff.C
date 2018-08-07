#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TGraphErrors.h>
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <AliAODDimuon.h>
#endif   

void TriggerEff(char* period="LHC15o", char *type="JPsi"){

gStyle->SetPalette(1);
gStyle->SetCanvasColor(10);
gStyle->SetFrameFillColor(10);
gStyle->SetOptStat(111111);

//--------------------------------
// event summary
//--------------------------------

//TH1D *hNEv = new TH1D("hNEv","Muon events",3,0.,3.);
//TString namelabel[3]={"AllMuons","MuonsFromPsi","MuonsNotFromPsi"};
//for(int kk=0;kk<3;kk++) hNEv->GetXaxis()->SetBinLabel(kk+1,namelabel[kk]);
//======================================================================================================
//BEFORE
//======================================================================================================
TH1D *hMassBefore_25y4 = new TH1D("hMassBefore_25y4","hMassBefore_25y4",100,0.,5.);
TH1D *hPtBefore_25y4 = new TH1D("hPtBefore_25y4","hPtBefore_25y4",20,0.,20.);
TH1D *hYBefore_25y4 = new TH1D("hYBefore_25y4","hYBefore_25y4",6,-4,-2.5);

TH1D *hMassBefore_25y3 = new TH1D("hMassBefore_25y3","hMassBefore_25y3",100,0.,5.);
TH1D *hPtBefore_25y3 = new TH1D("hPtBefore_25y3","hPtBefore_25y3",20,0.,20.);
TH1D *hYBefore_25y3 = new TH1D("hYBefore_25y3","hYBefore_25y3",6,-4,-2.5);

TH1D *hMassBefore_3y35 = new TH1D("hMassBefore_3y35","hMassBefore_3y35",100,0.,5.);
TH1D *hPtBefore_3y35 = new TH1D("hPtBefore_3y35","hPtBefore_3y35",20,0.,20.);
TH1D *hYBefore_3y35 = new TH1D("hYBefore_3y35","hYBefore_3y35",6,-4,-2.5);

TH1D *hMassBefore_35y4 = new TH1D("hMassBefore_35y4","hMassBefore_35y4",100,0.,5.);
TH1D *hPtBefore_35y4 = new TH1D("hPtBefore_35y4","hPtBefore_35y4",20,0.,20.);
TH1D *hYBefore_35y4 = new TH1D("hYBefore_35y4","hYBefore_35y4",6,-4,-2.5);

TH1D *hMassBefore_c010 = new TH1D("hMassBefore_c010","hMassBefore_c010",100,0.,5.);
TH1D *hPtBefore_c010 = new TH1D("hPtBefore_c010","hPtBefore_c010",20,0.,20.);
TH1D *hYBefore_c010 = new TH1D("hYBefore_c010","hYBefore_c010",6,-4,-2.5);

TH1D *hMassBefore_c1020 = new TH1D("hMassBefore_c1020","hMassBefore_c1020",100,0.,5.);
TH1D *hPtBefore_c1020 = new TH1D("hPtBefore_c1020","hPtBefore_c1020",20,0.,20.);
TH1D *hYBefore_c1020 = new TH1D("hYBefore_c1020","hYBefore_c1020",6,-4,-2.5);

TH1D *hMassBefore_c2030 = new TH1D("hMassBefore_c2030","hMassBefore_c2030",100,0.,5.);
TH1D *hPtBefore_c2030 = new TH1D("hPtBefore_c2030","hPtBefore_c2030",20,0.,20.);
TH1D *hYBefore_c2030 = new TH1D("hYBefore_c2030","hYBefore_c2030",6,-4,-2.5);

TH1D *hMassBefore_c3040 = new TH1D("hMassBefore_c3040","hMassBefore_c3040",100,0.,5.);
TH1D *hPtBefore_c3040 = new TH1D("hPtBefore_c3040","hPtBefore_c3040",20,0.,20.);
TH1D *hYBefore_c3040 = new TH1D("hYBefore_c3040","hYBefore_c3040",6,-4,-2.5);

TH1D *hMassBefore_c4050 = new TH1D("hMassBefore_c4050","hMassBefore_c4050",100,0.,5.);
TH1D *hPtBefore_c4050 = new TH1D("hPtBefore_c4050","hPtBefore_c4050",20,0.,20.);
TH1D *hYBefore_c4050 = new TH1D("hYBefore_c4050","hYBefore_c4050",6,-4,-2.5);

TH1D *hMassBefore_c5060 = new TH1D("hMassBefore_c5060","hMassBefore_c5060",100,0.,5.);
TH1D *hPtBefore_c5060 = new TH1D("hPtBefore_c5060","hPtBefore_c5060",20,0.,20.);
TH1D *hYBefore_c5060 = new TH1D("hYBefore_c5060","hYBefore_c5060",6,-4,-2.5);

TH1D *hMassBefore_c6070 = new TH1D("hMassBefore_c6070","hMassBefore_c6070",100,0.,5.);
TH1D *hPtBefore_c6070 = new TH1D("hPtBefore_c6070","hPtBefore_c6070",20,0.,20.);
TH1D *hYBefore_c6070 = new TH1D("hYBefore_c6070","hYBefore_c6070",6,-4,-2.5);

TH1D *hMassBefore_c7080 = new TH1D("hMassBefore_c7080","hMassBefore_c7080",100,0.,5.);
TH1D *hPtBefore_c7080 = new TH1D("hPtBefore_c7080","hPtBefore_c7080",20,0.,20.);
TH1D *hYBefore_c7080 = new TH1D("hYBefore_c7080","hYBefore_c7080",6,-4,-2.5);

TH1D *hMassBefore_c8090 = new TH1D("hMassBefore_c8090","hMassBefore_c8090",100,0.,5.);
TH1D *hPtBefore_c8090 = new TH1D("hPtBefore_c8090","hPtBefore_c8090",20,0.,20.);
TH1D *hYBefore_c8090 = new TH1D("hYBefore_c8090","hYBefore_c8090",6,-4,-2.5);

TH1D *hMassBefore_c6090 = new TH1D("hMassBefore_c6090","hMassBefore_c6090",100,0.,5.);
TH1D *hPtBefore_c6090 = new TH1D("hPtBefore_c6090","hPtBefore_c6090",20,0.,20.);
TH1D *hYBefore_c6090 = new TH1D("hYBefore_c6090","hYBefore_c6090",6,-4,-2.5);
//======================================================================================================
//AFTER
//======================================================================================================
TH1D *hMassAfter_25y4 = new TH1D("hMassAfter_25y4","hMassAfter_25y4",100,0.,5.);
TH1D *hPtAfter_25y4 = new TH1D("hPtAfter_25y4","hPtAfter_25y4",20,0.,20.);
TH1D *hYAfter_25y4 = new TH1D("hYAfter_25y4","hYAfter_25y4",6,-4,-2.5);

TH1D *hMassAfter_25y3 = new TH1D("hMassAfter_25y3","hMassAfter_25y3",100,0.,5.);
TH1D *hPtAfter_25y3 = new TH1D("hPtAfter_25y3","hPtAfter_25y3",20,0.,20.);
TH1D *hYAfter_25y3 = new TH1D("hYAfter_25y3","hYAfter_25y3",6,-4,-2.5);

TH1D *hMassAfter_3y35 = new TH1D("hMassAfter_3y35","hMassAfter_3y35",100,0.,5.);
TH1D *hPtAfter_3y35 = new TH1D("hPtAfter_3y35","hPtAfter_3y35",20,0.,20.);
TH1D *hYAfter_3y35 = new TH1D("hYAfter_3y35","hYAfter_3y35",6,-4,-2.5);

TH1D *hMassAfter_35y4 = new TH1D("hMassAfter_35y4","hMassAfter_35y4",100,0.,5.);
TH1D *hPtAfter_35y4 = new TH1D("hPtAfter_35y4","hPtAfter_35y4",20,0.,20.);
TH1D *hYAfter_35y4 = new TH1D("hYAfter_35y4","hYAfter_35y4",6,-4,-2.5);

TH1D *hMassAfter_c010 = new TH1D("hMassAfter_c010","hMassAfter_c010",100,0.,5.);
TH1D *hPtAfter_c010 = new TH1D("hPtAfter_c010","hPtAfter_c010",20,0.,20.);
TH1D *hYAfter_c010 = new TH1D("hYAfter_c010","hYAfter_c010",6,-4,-2.5);

TH1D *hMassAfter_c1020 = new TH1D("hMassAfter_c1020","hMassAfter_c1020",100,0.,5.);
TH1D *hPtAfter_c1020 = new TH1D("hPtAfter_c1020","hPtAfter_c1020",20,0.,20.);
TH1D *hYAfter_c1020 = new TH1D("hYAfter_c1020","hYAfter_c1020",6,-4,-2.5);

TH1D *hMassAfter_c2030 = new TH1D("hMassAfter_c2030","hMassAfter_c2030",100,0.,5.);
TH1D *hPtAfter_c2030 = new TH1D("hPtAfter_c2030","hPtAfter_c2030",20,0.,20.);
TH1D *hYAfter_c2030 = new TH1D("hYAfter_c2030","hYAfter_c2030",6,-4,-2.5);

TH1D *hMassAfter_c3040 = new TH1D("hMassAfter_c3040","hMassAfter_c3040",100,0.,5.);
TH1D *hPtAfter_c3040 = new TH1D("hPtAfter_c3040","hPtAfter_c3040",20,0.,20.);
TH1D *hYAfter_c3040 = new TH1D("hYAfter_c3040","hYAfter_c3040",6,-4,-2.5);

TH1D *hMassAfter_c4050 = new TH1D("hMassAfter_c4050","hMassAfter_c4050",100,0.,5.);
TH1D *hPtAfter_c4050 = new TH1D("hPtAfter_c4050","hPtAfter_c4050",20,0.,20.);
TH1D *hYAfter_c4050 = new TH1D("hYAfter_c4050","hYAfter_c4050",6,-4,-2.5);

TH1D *hMassAfter_c5060 = new TH1D("hMassAfter_c5060","hMassAfter_c5060",100,0.,5.);
TH1D *hPtAfter_c5060 = new TH1D("hPtAfter_c5060","hPtAfter_c5060",20,0.,20.);
TH1D *hYAfter_c5060 = new TH1D("hYAfter_c5060","hYAfter_c5060",6,-4,-2.5);

TH1D *hMassAfter_c6070 = new TH1D("hMassAfter_c6070","hMassAfter_c6070",100,0.,5.);
TH1D *hPtAfter_c6070 = new TH1D("hPtAfter_c6070","hPtAfter_c6070",20,0.,20.);
TH1D *hYAfter_c6070 = new TH1D("hYAfter_c6070","hYAfter_c6070",6,-4,-2.5);

TH1D *hMassAfter_c7080 = new TH1D("hMassAfter_c7080","hMassAfter_c7080",100,0.,5.);
TH1D *hPtAfter_c7080 = new TH1D("hPtAfter_c7080","hPtAfter_c7080",20,0.,20.);
TH1D *hYAfter_c7080 = new TH1D("hYAfter_c7080","hYAfter_c7080",6,-4,-2.5);

TH1D *hMassAfter_c8090 = new TH1D("hMassAfter_c8090","hMassAfter_c8090",100,0.,5.);
TH1D *hPtAfter_c8090 = new TH1D("hPtAfter_c8090","hPtAfter_c8090",20,0.,20.);
TH1D *hYAfter_c8090 = new TH1D("hYAfter_c8090","hYAfter_c8090",6,-4,-2.5);

TH1D *hMassAfter_c6090 = new TH1D("hMassAfter_c6090","hMassAfter_c6090",100,0.,5.);
TH1D *hPtAfter_c6090 = new TH1D("hPtAfter_c6090","hPtAfter_c6090",20,0.,20.);
TH1D *hYAfter_c6090 = new TH1D("hYAfter_c6090","hYAfter_c6090",6,-4,-2.5);
//======================================================================================================
TFile *hData = new TFile("/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Data/Apt_Lpt/With_CINT7_MUFAST/Data_divide_Pt.root");
TH1D *hData25y4 =(TH1D*)hData->Get("hData25y4");
TH1D *hData25y3 =(TH1D*)hData->Get("hData25y3");
TH1D *hData3y35 =(TH1D*)hData->Get("hData3y35");
TH1D *hData35y4 =(TH1D*)hData->Get("hData35y4");

  TFile *fData_25y4 = new TFile("Data_divide_centrality.root");

  TH1D *hDatac010 = (TH1D*)fData_25y4->Get("hDatac010");
  TH1D *hDatac1020 = (TH1D*)fData_25y4->Get("hDatac1020");
  TH1D *hDatac2030 = (TH1D*)fData_25y4->Get("hDatac2030");
  TH1D *hDatac3040 = (TH1D*)fData_25y4->Get("hDatac3040");
  TH1D *hDatac4050 = (TH1D*)fData_25y4->Get("hDatac4050");
  TH1D *hDatac5060 = (TH1D*)fData_25y4->Get("hDatac5060");
  TH1D *hDatac6070 = (TH1D*)fData_25y4->Get("hDatac6070");
  TH1D *hDatac7080 = (TH1D*)fData_25y4->Get("hDatac7080");
  TH1D *hDatac8090 = (TH1D*)fData_25y4->Get("hDatac8090");
  TH1D *hDatac6090 = (TH1D*)fData_25y4->Get("hDatac6090");

TFile *hMC_25y4 = new TFile("/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Simulation/AllPt_Lpt/MC_divide_Pt.root");
TH1D *hMC25y4 = (TH1D*)hMC_25y4->Get("Divide_Pt");
TFile *fMC_25y3 = new TFile("/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Simulation/AllPt_Lpt/In_rapidity_bins/MC_25y3_divide_Pt_new.root");
TH1D *hMC25y3 = (TH1D*)fMC_25y3->Get("Divide_Pt");
TFile *fMC_3y35 = new TFile("/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Simulation/AllPt_Lpt/In_rapidity_bins/MC_3y35_divide_Pt_new.root");
TH1D *hMC3y35 = (TH1D*)fMC_3y35->Get("Divide_Pt");
TFile *fMC_35y4 = new TFile("/home/biswarup/VAF/PbPb5TeV/Trigger_efficiency/Simulation/AllPt_Lpt/In_rapidity_bins/MC_35y4_divide_Pt_new.root");
TH1D *hMC35y4 = (TH1D*)fMC_35y4->Get("Divide_Pt");

  TFile *fEmbed_25y4 = new TFile("/home/biswarup/VAF/PbPb5TeV/Reading_Embedding/weighting/Embedding_divide_centrality.root");
  TH1D *hMCc010 = (TH1D*)fEmbed_25y4->Get("hEmbedc010");
  TH1D *hMCc1020 = (TH1D*)fEmbed_25y4->Get("hEmbedc1020");
  TH1D *hMCc2030 = (TH1D*)fEmbed_25y4->Get("hEmbedc2030");
  TH1D *hMCc3040 = (TH1D*)fEmbed_25y4->Get("hEmbedc3040");
  TH1D *hMCc4050 = (TH1D*)fEmbed_25y4->Get("hEmbedc4050");
  TH1D *hMCc5060 = (TH1D*)fEmbed_25y4->Get("hEmbedc5060");
  TH1D *hMCc6070 = (TH1D*)fEmbed_25y4->Get("hEmbedc6070");
  TH1D *hMCc7080 = (TH1D*)fEmbed_25y4->Get("hEmbedc7080");
  TH1D *hMCc8090 = (TH1D*)fEmbed_25y4->Get("hEmbedc8090");
  TH1D *hMCc6090 = (TH1D*)fEmbed_25y4->Get("hEmbedc6090");

//
//funzione estratta dal fit e funzione shiftata. Prototipo della funzione:
//par[3]+par[0]*(1.+TMath::Erf((xx-par[1])/par[2]/TMath::Sqrt(2.)));
//nella funzione shiftata facciamo variare SOLO ptcut e sigma (par[1] e par [2])
//

TF1 *fbefore = new TF1("fbefore","0.05663+0.4656*(1.+TMath::Erf((x-1.135)/0.3357/TMath::Sqrt(2.)))",0,20); //MC
TF1 *fafter  = new TF1("fafter","0.1482+0.4151*(1.+TMath::Erf((x-1.049)/0.2867/TMath::Sqrt(2.)))",0,20); //Data

//TF1 *fafter  = new TF1("fafter","0.05663+0.4656*(1.+TMath::Erf((x-1.049)/0.2867/TMath::Sqrt(2.)))",0,20); //Data

//SumOfWeight
hMassBefore_25y4->Sumw2();
hPtBefore_25y4->Sumw2();
hYBefore_25y4->Sumw2();
//
hMassAfter_25y4->Sumw2();
hPtAfter_25y4->Sumw2(); 
hYAfter_25y4->Sumw2();

hMassBefore_25y3->Sumw2();
hPtBefore_25y3->Sumw2();
hYBefore_25y3->Sumw2();
//
hMassAfter_25y3->Sumw2();
hPtAfter_25y3->Sumw2(); 
hYAfter_25y3->Sumw2();

hMassBefore_3y35->Sumw2();
hPtBefore_3y35->Sumw2();
hYBefore_3y35->Sumw2();
//
hMassAfter_3y35->Sumw2();
hPtAfter_3y35->Sumw2(); 
hYAfter_3y35->Sumw2();

hMassBefore_35y4->Sumw2();
hPtBefore_35y4->Sumw2();
hYBefore_35y4->Sumw2();
//
hMassAfter_35y4->Sumw2();
hPtAfter_35y4->Sumw2(); 
hYAfter_35y4->Sumw2();

hMassBefore_c010->Sumw2();
hPtBefore_c010->Sumw2();
hYBefore_c010->Sumw2();
//
hMassAfter_c010->Sumw2();
hPtAfter_c010->Sumw2(); 
hYAfter_c010->Sumw2();

hMassBefore_c1020->Sumw2();
hPtBefore_c1020->Sumw2();
hYBefore_c1020->Sumw2();
//
hMassAfter_c1020->Sumw2();
hPtAfter_c1020->Sumw2(); 
hYAfter_c1020->Sumw2();

hMassBefore_c2030->Sumw2();
hPtBefore_c2030->Sumw2();
hYBefore_c2030->Sumw2();
//
hMassAfter_c2030->Sumw2();
hPtAfter_c2030->Sumw2(); 
hYAfter_c2030->Sumw2();

hMassBefore_c3040->Sumw2();
hPtBefore_c3040->Sumw2();
hYBefore_c3040->Sumw2();
//
hMassAfter_c3040->Sumw2();
hPtAfter_c3040->Sumw2(); 
hYAfter_c3040->Sumw2();
//
hMassBefore_c4050->Sumw2();
hPtBefore_c4050->Sumw2();
hYBefore_c4050->Sumw2();
//
hMassAfter_c4050->Sumw2();
hPtAfter_c4050->Sumw2(); 
hYAfter_c4050->Sumw2();

hMassBefore_c5060->Sumw2();
hPtBefore_c5060->Sumw2();
hYBefore_c5060->Sumw2();
//
hMassAfter_c5060->Sumw2();
hPtAfter_c5060->Sumw2(); 
hYAfter_c5060->Sumw2();

hMassBefore_c6070->Sumw2();
hPtBefore_c6070->Sumw2();
hYBefore_c6070->Sumw2();
//
hMassAfter_c6070->Sumw2();
hPtAfter_c6070->Sumw2(); 
hYAfter_c6070->Sumw2();

hMassBefore_c7080->Sumw2();
hPtBefore_c7080->Sumw2();
hYBefore_c7080->Sumw2();
//
hMassAfter_c7080->Sumw2();
hPtAfter_c7080->Sumw2(); 
hYAfter_c7080->Sumw2();

hMassBefore_c8090->Sumw2();
hPtBefore_c8090->Sumw2();
hYBefore_c8090->Sumw2();
//
hMassAfter_c8090->Sumw2();
hPtAfter_c8090->Sumw2(); 
hYAfter_c8090->Sumw2();

hMassBefore_c6090->Sumw2();
hPtBefore_c6090->Sumw2();
hYBefore_c6090->Sumw2();
//
hMassAfter_c6090->Sumw2();
hPtAfter_c6090->Sumw2(); 
hYAfter_c6090->Sumw2();
//--------------------------------
// Loop on runs
//--------------------------------
Int_t NTotRun=137;
int *Run = new int[NTotRun]; 
Double_t NAcc[300]; 
Double_t ErrNAcc[300]; 

Int_t firstchunk=1;
Int_t nchunk=200;
FILE *flist;
char runlist[300];
sprintf(runlist,"dummy.dat"); 
sprintf(runlist,"RunList_%s.dat",period); 
 
if((flist = fopen(runlist,"r")) == NULL){
   printf("Cannot open file %s \n", runlist);
   return;
} else printf("Reading List %s\n",runlist);

TChain *aodTree;
for(int irun=0; irun<NTotRun; irun++){
  fscanf(flist,"%d",&Run[irun]);
  if (feof(flist)!=0) {
    break;
  }
  printf("\nNumber = %d\n",irun);
  aodTree = new TChain("aodTree");
    
  char file[300];
  for(int ichunk=1;ichunk<nchunk;ichunk++){
    

    if(ichunk<10) sprintf(file,"/home/biswarup/VAF/PbPb5TeV/With_More_Histos/With_Many_More_Histos/28_12_2015/AXE/AXE/MC_Simulation/MC_files_v4/%d/00%d/AliAOD.Muons.root",Run[irun],ichunk);
     else if(ichunk>10 && ichunk<100)sprintf(file,"/home/biswarup/VAF/PbPb5TeV/With_More_Histos/With_Many_More_Histos/28_12_2015/AXE/AXE/MC_Simulation/MC_files_v4/%d/0%d/AliAOD.rootAliAOD.Muons.root",Run[irun],ichunk);	
      else if(ichunk>100)sprintf(file,"/home/biswarup/VAF/PbPb5TeV/With_More_Histos/With_Many_More_Histos/28_12_2015/AXE/AXE/MC_Simulation/MC_files_v4/%d/0%d/AliAOD.Muons.root",Run[irun],ichunk);   

    
        
      Long_t *dummy1 =0, *dummy2 =0, *dummy3 =0, *dummy4 =0;
      if (gSystem->GetPathInfo(file,dummy1,dummy2,dummy3, dummy4) == 0) {
        printf(" Opening %s\n",file);
        aodTree->Add(file);
      }
  } // end loop on chunks


  for(int nev=0;nev<aodTree->GetEntries();nev++){
   if((nev%10000)==0) printf("Reading ev %d\n",nev);
     
  aodTree->GetEvent(nev);
    
  //Read AOD event    
  AliAODEvent *event = new AliAODEvent();
  event->ReadFromTree(aodTree);
  
  //Read MC info   
  TClonesArray *mcarray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!mcarray) {
    printf("No MC info in AOD!\n");
    return;
  }  
    

   //**************************************************************************
   // Loop on reconstructed tracks    
   //**************************************************************************
	
	  if(event->GetNumberOfDimuons()==0 || event->GetNumberOfDimuons()>1) continue;
          for(int i=0;i<event->GetNumberOfDimuons();i++){
           AliAODDimuon *dimu = dynamic_cast<AliAODDimuon*>(event->GetDimuon(i));
           AliAODTrack *mu0 = dimu->GetMu(0); 
           AliAODTrack *mu1 = dimu->GetMu(1); 
                 
           Double_t DimuPt = dimu->Pt();
           Double_t DimuY = dimu->Y();
           Double_t DimuMass = dimu->Mass();
           Double_t RAbs_Mu0 = mu0->GetRAtAbsorberEnd();
           Double_t RAbs_Mu1 = mu1->GetRAtAbsorberEnd();
           Double_t Eta_Mu0 = mu0->Eta();
           Double_t Eta_Mu1 = mu1->Eta();
           Int_t Charge_Mu0 = mu0->Charge();
           Int_t Charge_Mu1 = mu1->Charge();
	   
	   if(DimuY<-4 || DimuY>-2.5) continue;
           if((RAbs_Mu0<17.6 || RAbs_Mu0>89.5) || (RAbs_Mu1<17.6 || RAbs_Mu1>89.5)) continue;
           if((Eta_Mu0<-4 || Eta_Mu0>-2.5) || (Eta_Mu1<-4 || Eta_Mu1>-2.5)) continue;
                 
           Int_t Match_Mu0=mu0->GetMatchTrigger();
           Int_t Match_Mu1=mu1->GetMatchTrigger();
           Double_t Pt_Mu0 = mu0->Pt(); 
           Double_t Pt_Mu1 = mu1->Pt();
	      	   
           if(Match_Mu0>=2 && Match_Mu1>=2){ 

	   //epsilon1 e epsilon2 valutati sulla funzione di fit iniziale
	   Double_t eps1_bef25y4 = hMC25y4->GetBinContent(hMC25y4->FindFixBin(Pt_Mu0));
           if(eps1_bef25y4 == 0)eps1_bef25y4 = 0.05;
           Double_t eps2_bef25y4 = hMC25y4->GetBinContent(hMC25y4->FindFixBin(Pt_Mu1));
           if(eps2_bef25y4 == 0)eps2_bef25y4 = 0.05;

	   Double_t weight_bef25y4 = eps1_bef25y4*eps2_bef25y4;
           //==========================================================================

           Double_t eps1_aft25y4 = hData25y4->GetBinContent(hData25y4->FindFixBin(Pt_Mu0));
           Double_t eps2_aft25y4 = hData25y4->GetBinContent(hData25y4->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aft25y4 = eps1_aft25y4*eps2_aft25y4;

           //==========================================================================

	   Double_t eps1_bef25y3 = hMC25y3->GetBinContent(hMC25y3->FindFixBin(Pt_Mu0));
           if(eps1_bef25y3 == 0)eps1_bef25y3 = 0.05;
           Double_t eps2_bef25y3 = hMC25y3->GetBinContent(hMC25y3->FindFixBin(Pt_Mu1));
           if(eps2_bef25y3 == 0)eps2_bef25y3 = 0.05;

	   Double_t weight_bef25y3 = eps1_bef25y3*eps2_bef25y3;
           //==========================================================================

           Double_t eps1_aft25y3 = hData25y3->GetBinContent(hData25y3->FindFixBin(Pt_Mu0));
           Double_t eps2_aft25y3 = hData25y3->GetBinContent(hData25y3->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aft25y3 = eps1_aft25y3*eps2_aft25y3;
           //==========================================================================

	   Double_t eps1_bef3y35 = hMC3y35->GetBinContent(hMC3y35->FindFixBin(Pt_Mu0));
           if(eps1_bef3y35 == 0)eps1_bef3y35 = 0.05;
           Double_t eps2_bef3y35 = hMC3y35->GetBinContent(hMC3y35->FindFixBin(Pt_Mu1));
           if(eps2_bef3y35 == 0)eps2_bef3y35 = 0.05;

	   Double_t weight_bef3y35 = eps1_bef3y35*eps2_bef3y35;
           //==========================================================================

           Double_t eps1_aft3y35 = hData3y35->GetBinContent(hData3y35->FindFixBin(Pt_Mu0));
           Double_t eps2_aft3y35 = hData3y35->GetBinContent(hData3y35->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aft3y35 = eps1_aft3y35*eps2_aft3y35;

           //==========================================================================

	   Double_t eps1_bef35y4 = hMC35y4->GetBinContent(hMC35y4->FindFixBin(Pt_Mu0));
           if(eps1_bef35y4 == 0)eps1_bef35y4 = 0.05;
           Double_t eps2_bef35y4 = hMC35y4->GetBinContent(hMC35y4->FindFixBin(Pt_Mu1));
           if(eps2_bef35y4 == 0)eps2_bef35y4 = 0.05;

	   Double_t weight_bef35y4 = eps1_bef35y4*eps2_bef35y4;

           //==========================================================================

           Double_t eps1_aft35y4 = hData35y4->GetBinContent(hData35y4->FindFixBin(Pt_Mu0));
           Double_t eps2_aft35y4 = hData35y4->GetBinContent(hData35y4->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aft35y4 = eps1_aft35y4*eps2_aft35y4;
           //==========================================================================
	   Double_t eps1_befc010 = hMCc010->GetBinContent(hMCc010->FindFixBin(Pt_Mu0));
           if(eps1_befc010 == 0)eps1_befc010 = 0.05;
           Double_t eps2_befc010 = hMCc010->GetBinContent(hMCc010->FindFixBin(Pt_Mu1));
           if(eps2_befc010 == 0)eps2_befc010 = 0.05;

	   Double_t weight_befc010 = eps1_befc010*eps2_befc010;
           //==========================================================================

           Double_t eps1_aftc010 = hDatac010->GetBinContent(hDatac010->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc010 = hDatac010->GetBinContent(hDatac010->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc010 = eps1_aftc010*eps2_aftc010;

           //==========================================================================

	   Double_t eps1_befc1020 = hMCc1020->GetBinContent(hMCc1020->FindFixBin(Pt_Mu0));
           if(eps1_befc1020 == 0)eps1_befc1020 = 0.05;
           Double_t eps2_befc1020 = hMCc1020->GetBinContent(hMCc1020->FindFixBin(Pt_Mu1));
           if(eps2_befc1020 == 0)eps2_befc1020 = 0.05;

	   Double_t weight_befc1020 = eps1_befc1020*eps2_befc1020;
           //==========================================================================

           Double_t eps1_aftc1020 = hDatac1020->GetBinContent(hDatac1020->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc1020 = hDatac1020->GetBinContent(hDatac1020->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc1020 = eps1_aftc1020*eps2_aftc1020;
           //==========================================================================

	   Double_t eps1_befc2030 = hMCc2030->GetBinContent(hMCc2030->FindFixBin(Pt_Mu0));
           if(eps1_befc2030 == 0)eps1_befc2030 = 0.05;
           Double_t eps2_befc2030 = hMCc2030->GetBinContent(hMCc2030->FindFixBin(Pt_Mu1));
           if(eps2_befc2030 == 0)eps2_befc2030 = 0.05;

	   Double_t weight_befc2030 = eps1_befc2030*eps2_befc2030;
           //==========================================================================

           Double_t eps1_aftc2030 = hDatac2030->GetBinContent(hDatac2030->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc2030 = hDatac2030->GetBinContent(hDatac2030->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc2030 = eps1_aftc2030*eps2_aftc2030;

           //==========================================================================

	   Double_t eps1_befc3040 = hMCc3040->GetBinContent(hMCc3040->FindFixBin(Pt_Mu0));
           if(eps1_befc3040 == 0)eps1_befc3040 = 0.05;
           Double_t eps2_befc3040 = hMCc3040->GetBinContent(hMCc3040->FindFixBin(Pt_Mu1));
           if(eps2_befc3040 == 0)eps2_befc3040 = 0.05;

	   Double_t weight_befc3040 = eps1_befc3040*eps2_befc3040;

           //==========================================================================

           Double_t eps1_aftc3040 = hDatac3040->GetBinContent(hDatac3040->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc3040 = hDatac3040->GetBinContent(hDatac3040->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc3040 = eps1_aftc3040*eps2_aftc3040;

           //==========================================================================

	   Double_t eps1_befc4050 = hMCc4050->GetBinContent(hMCc4050->FindFixBin(Pt_Mu0));
           if(eps1_befc4050 == 0)eps1_befc4050 = 0.05;
           Double_t eps2_befc4050 = hMCc4050->GetBinContent(hMCc4050->FindFixBin(Pt_Mu1));
           if(eps2_befc4050 == 0)eps2_befc4050 = 0.05;

	   Double_t weight_befc4050 = eps1_befc4050*eps2_befc4050;

           //==========================================================================

           Double_t eps1_aftc4050 = hDatac4050->GetBinContent(hDatac4050->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc4050 = hDatac4050->GetBinContent(hDatac4050->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc4050 = eps1_aftc4050*eps2_aftc4050;
           //==========================================================================
	   Double_t eps1_befc5060 = hMCc5060->GetBinContent(hMCc5060->FindFixBin(Pt_Mu0));
           if(eps1_befc5060 == 0)eps1_befc5060 = 0.05;
           Double_t eps2_befc5060 = hMCc5060->GetBinContent(hMCc5060->FindFixBin(Pt_Mu1));
           if(eps2_befc5060 == 0)eps2_befc5060 = 0.05;

	   Double_t weight_befc5060 = eps1_befc5060*eps2_befc5060;
           //==========================================================================

           Double_t eps1_aftc5060 = hDatac5060->GetBinContent(hDatac5060->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc5060 = hDatac5060->GetBinContent(hDatac5060->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc5060 = eps1_aftc5060*eps2_aftc5060;

           //==========================================================================

	   Double_t eps1_befc6070 = hMCc6070->GetBinContent(hMCc6070->FindFixBin(Pt_Mu0));
           if(eps1_befc6070 == 0)eps1_befc6070 = 0.05;
           Double_t eps2_befc6070 = hMCc6070->GetBinContent(hMCc6070->FindFixBin(Pt_Mu1));
           if(eps2_befc6070 == 0)eps2_befc6070 = 0.05;

	   Double_t weight_befc6070 = eps1_befc6070*eps2_befc6070;
           //==========================================================================

           Double_t eps1_aftc6070 = hDatac6070->GetBinContent(hDatac6070->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc6070 = hDatac6070->GetBinContent(hDatac6070->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc6070 = eps1_aftc6070*eps2_aftc6070;
           //==========================================================================

	   Double_t eps1_befc7080 = hMCc7080->GetBinContent(hMCc7080->FindFixBin(Pt_Mu0));
           if(eps1_befc7080 == 0)eps1_befc7080 = 0.05;
           Double_t eps2_befc7080 = hMCc7080->GetBinContent(hMCc7080->FindFixBin(Pt_Mu1));
           if(eps2_befc7080 == 0)eps2_befc7080 = 0.05;

	   Double_t weight_befc7080 = eps1_befc7080*eps2_befc7080;
           //==========================================================================

           Double_t eps1_aftc7080 = hDatac7080->GetBinContent(hDatac7080->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc7080 = hDatac7080->GetBinContent(hDatac7080->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc7080 = eps1_aftc7080*eps2_aftc7080;

           //==========================================================================

	   Double_t eps1_befc8090 = hMCc8090->GetBinContent(hMCc8090->FindFixBin(Pt_Mu0));
           if(eps1_befc8090 == 0)eps1_befc8090 = 0.05;
           Double_t eps2_befc8090 = hMCc8090->GetBinContent(hMCc8090->FindFixBin(Pt_Mu1));
           if(eps2_befc8090 == 0)eps2_befc8090 = 0.05;

	   Double_t weight_befc8090 = eps1_befc8090*eps2_befc8090;

           //==========================================================================

           Double_t eps1_aftc8090 = hDatac8090->GetBinContent(hDatac8090->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc8090 = hDatac8090->GetBinContent(hDatac8090->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc8090 = eps1_aftc8090*eps2_aftc8090;
           //==========================================================================

	   Double_t eps1_befc6090 = hMCc6090->GetBinContent(hMCc6090->FindFixBin(Pt_Mu0));
           if(eps1_befc6090 == 0)eps1_befc6090 = 0.05;
           Double_t eps2_befc6090 = hMCc6090->GetBinContent(hMCc6090->FindFixBin(Pt_Mu1));
           if(eps2_befc6090 == 0)eps2_befc6090 = 0.05;

	   Double_t weight_befc6090 = eps1_befc6090*eps2_befc6090;

           //==========================================================================

           Double_t eps1_aftc6090 = hDatac6090->GetBinContent(hDatac6090->FindFixBin(Pt_Mu0));
           Double_t eps2_aftc6090 = hDatac6090->GetBinContent(hDatac6090->FindFixBin(Pt_Mu1));
 
	   Double_t weight_aftc6090 = eps1_aftc6090*eps2_aftc6090;	   
	   //==========================================================================
	   //				   filling
	   //==========================================================================
   	   	   	   

	   if (DimuY>-4 && DimuY<-2.5){  

	   hMassBefore_25y4->Fill(DimuMass);
	   hPtBefore_25y4->Fill(DimuPt);
           hYBefore_25y4->Fill(DimuY);  

	   hMassAfter_25y4->Fill(DimuMass,weight_aft25y4/weight_bef25y4);
	   hPtAfter_25y4->Fill(DimuPt,weight_aft25y4/weight_bef25y4);
	   hYAfter_25y4->Fill(DimuY,weight_aft25y4/weight_bef25y4);

	   hMassBefore_25y3->Fill(DimuMass);
	   hPtBefore_25y3->Fill(DimuPt);
           hYBefore_25y3->Fill(DimuY);  

	   hMassAfter_25y3->Fill(DimuMass,weight_aft25y3/weight_bef25y3);
	   hPtAfter_25y3->Fill(DimuPt,weight_aft25y3/weight_bef25y3);
	   hYAfter_25y3->Fill(DimuY,weight_aft25y3/weight_bef25y3);

	   hMassBefore_3y35->Fill(DimuMass);
	   hPtBefore_3y35->Fill(DimuPt);
           hYBefore_3y35->Fill(DimuY);  

	   hMassAfter_3y35->Fill(DimuMass,weight_aft3y35/weight_bef3y35);
	   hPtAfter_3y35->Fill(DimuPt,weight_aft3y35/weight_bef3y35);
	   hYAfter_3y35->Fill(DimuY,weight_aft3y35/weight_bef3y35);

	   hMassBefore_35y4->Fill(DimuMass);
	   hPtBefore_35y4->Fill(DimuPt);
           hYBefore_35y4->Fill(DimuY);  

	   hMassAfter_35y4->Fill(DimuMass,weight_aft35y4/weight_bef35y4);
	   hPtAfter_35y4->Fill(DimuPt,weight_aft35y4/weight_bef35y4);
	   hYAfter_35y4->Fill(DimuY,weight_aft35y4/weight_bef35y4);

	   hMassBefore_c010->Fill(DimuMass);
	   hPtBefore_c010->Fill(DimuPt);
           hYBefore_c010->Fill(DimuY);  

	   hMassAfter_c010->Fill(DimuMass,weight_aftc010/weight_befc010);
	   hPtAfter_c010->Fill(DimuPt,weight_aftc010/weight_befc010);
	   hYAfter_c010->Fill(DimuY,weight_aftc010/weight_befc010);

	   hMassBefore_c1020->Fill(DimuMass);
	   hPtBefore_c1020->Fill(DimuPt);
           hYBefore_c1020->Fill(DimuY);  

	   hMassAfter_c1020->Fill(DimuMass,weight_aftc1020/weight_befc1020);
	   hPtAfter_c1020->Fill(DimuPt,weight_aftc1020/weight_befc1020);
	   hYAfter_c1020->Fill(DimuY,weight_aftc1020/weight_befc1020);

	   hMassBefore_c2030->Fill(DimuMass);
	   hPtBefore_c2030->Fill(DimuPt);
           hYBefore_c2030->Fill(DimuY);  

	   hMassAfter_c2030->Fill(DimuMass,weight_aftc2030/weight_befc2030);
	   hPtAfter_c2030->Fill(DimuPt,weight_aftc2030/weight_befc2030);
	   hYAfter_c2030->Fill(DimuY,weight_aftc2030/weight_befc2030);

	   hMassBefore_c3040->Fill(DimuMass);
	   hPtBefore_c3040->Fill(DimuPt);
           hYBefore_c3040->Fill(DimuY);  

	   hMassAfter_c3040->Fill(DimuMass,weight_aftc3040/weight_befc3040);
	   hPtAfter_c3040->Fill(DimuPt,weight_aftc3040/weight_befc3040);
	   hYAfter_c3040->Fill(DimuY,weight_aftc3040/weight_befc3040);

	   hMassBefore_c4050->Fill(DimuMass);
	   hPtBefore_c4050->Fill(DimuPt);
           hYBefore_c4050->Fill(DimuY);  

	   hMassAfter_c4050->Fill(DimuMass,weight_aftc4050/weight_befc4050);
	   hPtAfter_c4050->Fill(DimuPt,weight_aftc4050/weight_befc4050);
	   hYAfter_c4050->Fill(DimuY,weight_aftc4050/weight_befc4050);

	   hMassBefore_c5060->Fill(DimuMass);
	   hPtBefore_c5060->Fill(DimuPt);
           hYBefore_c5060->Fill(DimuY);  

	   hMassAfter_c5060->Fill(DimuMass,weight_aftc5060/weight_befc5060);
	   hPtAfter_c5060->Fill(DimuPt,weight_aftc5060/weight_befc5060);
	   hYAfter_c5060->Fill(DimuY,weight_aftc5060/weight_befc5060);

	   hMassBefore_c6070->Fill(DimuMass);
	   hPtBefore_c6070->Fill(DimuPt);
           hYBefore_c6070->Fill(DimuY);  

	   hMassAfter_c6070->Fill(DimuMass,weight_aftc6070/weight_befc6070);
	   hPtAfter_c6070->Fill(DimuPt,weight_aftc6070/weight_befc6070);
	   hYAfter_c6070->Fill(DimuY,weight_aftc6070/weight_befc6070);

	   hMassBefore_c7080->Fill(DimuMass);
	   hPtBefore_c7080->Fill(DimuPt);
           hYBefore_c7080->Fill(DimuY);  

	   hMassAfter_c7080->Fill(DimuMass,weight_aftc7080/weight_befc7080);
	   hPtAfter_c7080->Fill(DimuPt,weight_aftc7080/weight_befc7080);
	   hYAfter_c7080->Fill(DimuY,weight_aftc7080/weight_befc7080);

	   hMassBefore_c8090->Fill(DimuMass);
	   hPtBefore_c8090->Fill(DimuPt);
           hYBefore_c8090->Fill(DimuY);  

	   hMassAfter_c8090->Fill(DimuMass,weight_aftc8090/weight_befc8090);
	   hPtAfter_c8090->Fill(DimuPt,weight_aftc8090/weight_befc8090);
	   hYAfter_c8090->Fill(DimuY,weight_aftc8090/weight_befc8090);

	   hMassBefore_c6090->Fill(DimuMass);
	   hPtBefore_c6090->Fill(DimuPt);
           hYBefore_c6090->Fill(DimuY);  

	   hMassAfter_c6090->Fill(DimuMass,weight_aftc6090/weight_befc6090);
	   hPtAfter_c6090->Fill(DimuPt,weight_aftc6090/weight_befc6090);
	   hYAfter_c6090->Fill(DimuY,weight_aftc6090/weight_befc6090);
	   }//Rapidity loop
	   
           }//Match trig
          }//loop dimuone     
        }//event loop
      }//run loop   
  
      
   //output histogram
   char fnameout[100];
   sprintf(fnameout,"TriggerEff_%s.root",period);

   TFile *fout = new TFile(fnameout,"recreate");
   fout->cd();
   //
   hMassBefore_25y4->Write();
   hPtBefore_25y4->Write();
   hYBefore_25y4->Write();
   //
   hMassAfter_25y4->Write();
   hPtAfter_25y4->Write(); 
   hYAfter_25y4->Write();
   //
   hMassBefore_25y3->Write();
   hPtBefore_25y3->Write();
   hYBefore_25y3->Write();
   //
   hMassAfter_25y3->Write();
   hPtAfter_25y3->Write(); 
   hYAfter_25y3->Write();
   //
   hMassBefore_3y35->Write();
   hPtBefore_3y35->Write();
   hYBefore_3y35->Write();
   //
   hMassAfter_3y35->Write();
   hPtAfter_3y35->Write(); 
   hYAfter_3y35->Write();
   //
   hMassBefore_35y4->Write();
   hPtBefore_35y4->Write();
   hYBefore_35y4->Write();
   //
   hMassAfter_35y4->Write();
   hPtAfter_35y4->Write(); 
   hYAfter_35y4->Write();
   //
   hMassBefore_c010->Write();
   hPtBefore_c010->Write();
   hYBefore_c010->Write();
   //
   hMassAfter_c010->Write();
   hPtAfter_c010->Write(); 
   hYAfter_c010->Write();
   //
   hMassBefore_c1020->Write();
   hPtBefore_c1020->Write();
   hYBefore_c1020->Write();
   //
   hMassAfter_c1020->Write();
   hPtAfter_c1020->Write(); 
   hYAfter_c1020->Write();
   //
   hMassBefore_c2030->Write();
   hPtBefore_c2030->Write();
   hYBefore_c2030->Write();
   //
   hMassAfter_c2030->Write();
   hPtAfter_c2030->Write(); 
   hYAfter_c2030->Write();
   //
   hMassBefore_c3040->Write();
   hPtBefore_c3040->Write();
   hYBefore_c3040->Write();
   //
   hMassAfter_c3040->Write();
   hPtAfter_c3040->Write(); 
   hYAfter_c3040->Write();
   //
   hMassBefore_c4050->Write();
   hPtBefore_c4050->Write();
   hYBefore_c4050->Write();
   //
   hMassAfter_c4050->Write();
   hPtAfter_c4050->Write(); 
   hYAfter_c4050->Write();
   //
   hMassBefore_c5060->Write();
   hPtBefore_c5060->Write();
   hYBefore_c5060->Write();
   //
   hMassAfter_c5060->Write();
   hPtAfter_c5060->Write(); 
   hYAfter_c5060->Write();
   //
   hMassBefore_c6070->Write();
   hPtBefore_c6070->Write();
   hYBefore_c6070->Write();
   //
   hMassAfter_c6070->Write();
   hPtAfter_c6070->Write(); 
   hYAfter_c6070->Write();
   //
   hMassBefore_c7080->Write();
   hPtBefore_c7080->Write();
   hYBefore_c7080->Write();
   //
   hMassAfter_c7080->Write();
   hPtAfter_c7080->Write(); 
   hYAfter_c7080->Write();
   //
   hMassBefore_c8090->Write();
   hPtBefore_c8090->Write();
   hYBefore_c8090->Write();
   //
   hMassAfter_c8090->Write();
   hPtAfter_c8090->Write(); 
   hYAfter_c8090->Write();
   //
   hMassBefore_c6090->Write();
   hPtBefore_c6090->Write();
   hYBefore_c6090->Write();
   //
   hMassAfter_c6090->Write();
   hPtAfter_c6090->Write(); 
   hYAfter_c6090->Write();
   //
   fout->Close();

}
