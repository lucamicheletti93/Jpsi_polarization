/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskPbPbJPsi.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a list of histograms.
// The macro class can run on AOD or in the train with the ESD filter.
// R. Arnaldi
//
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskPbPbJPsi_CXX
//#define AliAnalysisTaskPbPbJPsi_CXX

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"  //for Centrality

#include "AliAnalysisTaskPbPbJPsi.h"

//ADD THE FOLLOWING 1 LINE!!!
#include "AliMuonTrackCuts.h" 

ClassImp(AliAnalysisTaskPbPbJPsi)
//__________________________________________________________________________
AliAnalysisTaskPbPbJPsi::AliAnalysisTaskPbPbJPsi() :
  AliAnalysisTaskSE(),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fAODEvent(0x0),
  fMuonTrackCuts(0x0),
  fOutput(0x0)  
{
  /// Default ctor.
}

//__________________________________________________________________________
AliAnalysisTaskPbPbJPsi::AliAnalysisTaskPbPbJPsi(const char *name) :
  AliAnalysisTaskSE(name),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fAODEvent(0x0),
  fMuonTrackCuts(0x0),
  fOutput(0x0)  
{
 //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskPbPbJPsi","Calling Constructor");
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "TestStandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);  //I select the Pdca cut in this way 
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  
  DefineOutput(1,TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskPbPbJPsi& AliAnalysisTaskPbPbJPsi::operator=(const AliAnalysisTaskPbPbJPsi& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskPbPbJPsi::AliAnalysisTaskPbPbJPsi(const AliAnalysisTaskPbPbJPsi& c) :
  AliAnalysisTaskSE(c),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fCountTotEv(c.fCountTotEv),
  fCountTrigger(c.fCountTrigger),
  fCountCINT7(c.fCountCINT7),
  fCountCMUL7(c.fCountCMUL7),
  fCountCMLL7(c.fCountCMLL7),
  fCountCMSL7(c.fCountCMSL7),
  fCountCMSH7(c.fCountCMSH7),
  fAODEvent(c.fAODEvent),
  //ADD THE FOLLOWING 1 LINE!!!
  fMuonTrackCuts(c.fMuonTrackCuts),
  fOutput(c.fOutput)
 {
  //
  // Copy Constructor									
  //
}

//___________________________________________________________________________
AliAnalysisTaskPbPbJPsi::~AliAnalysisTaskPbPbJPsi() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskPbPbJPsi","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutput;
}

//___________________________________________________________________________
//
////ADD THE FOLLOWING 6 LINES!!!
void AliAnalysisTaskPbPbJPsi::NotifyRun()
{
//fMuonTrackCuts->GetMuonTrackCuts().SetPassName("pass2_muon");
fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
fMuonTrackCuts->SetRun(fInputHandler);
}
//___________________________________________________________________________
void AliAnalysisTaskPbPbJPsi::UserCreateOutputObjects(){

 //
 // output objects creation
 //	 
 fOutput = new TList();
 fOutput->SetOwner(); 
 //
 // dimuon histos
 //
TH1D *hMassOS_notrig = new TH1D("hMassOS_notrig","hMassOS_notrig",600,0.,15.);
TH1D *hMassOS = new TH1D("hMassOS","hMassOS",600,0.,15.);
TH1D *hMassOS_CINT7 = new TH1D("hMassOS_CINT7","hMassOS_CINT7",600,0.,15.);
TH1D *hMassOS_CMUL7 = new TH1D("hMassOS_CMUL7","hMassOS_CMUL7",600,0.,15.);
TH1D *hMassOS_1m = new TH1D("hMassOS_1m","hMassOS_1m",600,0.,15.);
TH1D *hMassOS_CINT7_1m = new TH1D("hMassOS_CINT7_1m","hMassOS_CINT7_1m",600,0.,15.);
TH1D *hMassOS_CMUL7_1m = new TH1D("hMassOS_CMUL7_1m","hMassOS_CMUL7_1m",600,0.,15.);
TH1D *hMassOS_CMSL7_1m = new TH1D("hMassOS_CMSL7_1m","hMassOS_CMSL7_1m",600,0.,15.);
TH1D *hMassOS_2m = new TH1D("hMassOS_2m","hMassOS_2m",600,0.,15.);
TH1D *hMassOS_CINT7_2m = new TH1D("hMassOS_CINT7_2m","hMassOS_CINT7_2m",600,0.,15.);
TH1D *hMassOS_CMUL7_2m = new TH1D("hMassOS_CMUL7_2m","hMassOS_CMUL7_2m",600,0.,15.);
TH1D *hMassOS_CMSL7_2m = new TH1D("hMassOS_CMSL7_2m","hMassOS_CMSL7_2m",600,0.,15.);


TH1D *hMassPP_CMLL7_2m = new TH1D("hMassPP_CMLL7_2m","hMassPP_CMLL7_2m",600,0.,15.);
TH1D *hMassMM_CMLL7_2m = new TH1D("hMassMM_CMLL7_2m","hMassMM_CMLL7_2m",600,0.,15.);

TH1D *hPt_CMUL7_2m = new TH1D("hPt_CMUL7_2m","hPt_CMUL7_2m",150,0.,15.);
TH1D *hY_CMUL7_2m = new TH1D("hY_CMUL7_2m","hY_CMUL7_2m",100,-5,0);
 //
 // various histos
 //
TH1D *hNDimuons = new TH1D("hNDimuons","hNDimuons",50,0.,50.); 
TH1D *hNDimuons_2m = new TH1D("hNDimuons_2m","hNDimuons_2m",50,0.,50.); 
TH1D *hMu0MatchTrigger = new TH1D("hMu0MatchTrigger","hMu0MatchTrigger",4,0.,4.);
TH1D *hNContrib = new TH1D("hNContrib","hNContrib",2000,0.,2000.);
TH1D *hZVertex = new TH1D("hZVertex","hZVertex",2000,-100.,100.);
 //
 // single muon histos
 //
TH1D *hNMu = new TH1D("hNMu","hNMu",50,0.,50.);
TH1D *hNMuons_1m = new TH1D("hNMuons_1m","hNMuons_1m",50,0.,50.);
TH1D *hPtMu_Match = new TH1D("hPtMu_Match","hPtMu_Match",600,0.,15.);
 //
 // trigger summary
 //
TH1D *hNEv = new TH1D("hNEv","hNEv",6,0.,6.);
TString namelabel[6]={"TotEv","CINT7","CMUL7","CMLL7","CMSL7","CMSH7"};
for(int kk=0;kk<6;kk++) hNEv->GetXaxis()->SetBinLabel(kk+1,namelabel[kk]);
 // centrality
 TH1F *hPercentileV0M = new TH1F("hPercentileV0M","hPercentileV0M",201,-1,200);  //V0M
 TH1F *hPercentileV0M_CMUL7 = new TH1F("hPercentileV0M_CMUL7","hPercentileV0M_CMUL7",201,-1,200);
 TH1F *hPercentileV0M_CINT7 = new TH1F("hPercentileV0M_CINT7","hPercentileV0M_CINT7",201,-1,200);

TH1D *hAllpt_1M = new TH1D("hAllpt_1M","hAllpt_1M",500,0.,50.);      
TH1D *hLpt_1M = new TH1D("hLpt_1M","hLpt_1M",500,0.,50.);

TH1D *hAllpt_25y3 = new TH1D("hAllpt_25y3","hAllpt_25y3",500,0.,50.);   
TH1D *hLpt_25y3 = new TH1D("hLpt_25y3","hLpt_25y3",500,0.,50.);
TH1D *hAllpt_3y35 = new TH1D("hAllpt_3y35","hAllpt_3y35",500,0.,50.); 
TH1D *hLpt_3y35 = new TH1D("hLpt_3y35","hLpt_3y35",500,0.,50.);
TH1D *hAllpt_35y4 = new TH1D("hAllpt_35y4","hAllpt_35y4",500,0.,50.);
TH1D *hLpt_35y4 = new TH1D("hLpt_35y4","hLpt_35y4",500,0.,50.); 

TH1D *hAllpt_25y3_1M = new TH1D("hAllpt_25y3_1M","hAllpt_25y3_1M",500,0.,50.);   
TH1D *hLpt_25y3_1M = new TH1D("hLpt_25y3_1M","hLpt_25y3_1M",500,0.,50.);
TH1D *hAllpt_3y35_1M = new TH1D("hAllpt_3y35_1M","hAllpt_3y35_1M",500,0.,50.); 
TH1D *hLpt_3y35_1M = new TH1D("hLpt_3y35_1M","hLpt_3y35_1M",500,0.,50.);
TH1D *hAllpt_35y4_1M = new TH1D("hAllpt_35y4_1M","hAllpt_35y4_1M",500,0.,50.);
TH1D *hLpt_35y4_1M = new TH1D("hLpt_35y4_1M","hLpt_35y4_1M",500,0.,50.); 
  
TH1D *hAllpt_010_1M = new TH1D("hAllpt_010_1M","hAllpt_010_1M",500,0.,50.);
TH1D *hLpt_010_1M = new TH1D("hLpt_010_1M","hLpt_010_1M",500,0.,50.);
TH1D *hAllpt_1020_1M = new TH1D("hAllpt_1020_1M","hAllpt_1020_1M",500,0.,50.);
TH1D *hLpt_1020_1M = new TH1D("hLpt_1020_1M","hLpt_1020_1M",500,0.,50.);
TH1D *hAllpt_2030_1M = new TH1D("hAllpt_2030_1M","hAllpt_2030_1M",500,0.,50.);
TH1D *hLpt_2030_1M = new TH1D("hLpt_2030_1M","hLpt_2030_1M",500,0.,50.);
TH1D *hAllpt_3040_1M = new TH1D("hAllpt_3040_1M","hAllpt_3040_1M",500,0.,50.);
TH1D *hLpt_3040_1M = new TH1D("hLpt_3040_1M","hLpt_3040_1M",500,0.,50.);
TH1D *hAllpt_4050_1M = new TH1D("hAllpt_4050_1M","hAllpt_4050_1M",500,0.,50.);
TH1D *hLpt_4050_1M = new TH1D("hLpt_4050_1M","hLpt_4050_1M",500,0.,50.);
TH1D *hAllpt_5060_1M = new TH1D("hAllpt_5060_1M","hAllpt_5060_1M",500,0.,50.);
TH1D *hLpt_5060_1M = new TH1D("hLpt_5060_1M","hLpt_5060_1M",500,0.,50.);
TH1D *hAllpt_6070_1M = new TH1D("hAllpt_6070_1M","hAllpt_6070_1M",500,0.,50.);
TH1D *hLpt_6070_1M = new TH1D("hLpt_6070_1M","hLpt_6070_1M",500,0.,50.);
TH1D *hAllpt_7080_1M = new TH1D("hAllpt_7080_1M","hAllpt_7080_1M",500,0.,50.);
TH1D *hLpt_7080_1M = new TH1D("hLpt_7080_1M","hLpt_7080_1M",500,0.,50.);
TH1D *hAllpt_8090_1M = new TH1D("hAllpt_8090_1M","hAllpt_8090_1M",500,0.,50.);
TH1D *hLpt_8090_1M = new TH1D("hLpt_8090_1M","hLpt_8090_1M",500,0.,50.);
 

 fOutput->Add(hMassOS_notrig); 	
 fOutput->Add(hMassOS); 	
 fOutput->Add(hMassOS_CINT7); 	
 fOutput->Add(hMassOS_1m); 	
 fOutput->Add(hMassOS_CINT7_1m); 	
 fOutput->Add(hMassOS_2m); 	
 fOutput->Add(hMassOS_CINT7_2m); 	
 fOutput->Add(hMassOS_CMUL7); 	
 fOutput->Add(hMassOS_CMUL7_1m); 	
 fOutput->Add(hMassOS_CMUL7_2m); 	
 fOutput->Add(hMassOS_CMSL7_1m); 	
 fOutput->Add(hMassOS_CMSL7_2m); 	
 fOutput->Add(hMassPP_CMLL7_2m); 	
 fOutput->Add(hMassMM_CMLL7_2m); 	
 fOutput->Add(hPt_CMUL7_2m);
 fOutput->Add(hY_CMUL7_2m);
 fOutput->Add(hMu0MatchTrigger);
 fOutput->Add(hNContrib);
 fOutput->Add(hZVertex);
 fOutput->Add(hPercentileV0M);
 fOutput->Add(hPercentileV0M_CMUL7);
 fOutput->Add(hPercentileV0M_CINT7);

 fOutput->Add(hNMu); 	
 fOutput->Add(hNMuons_1m); 	
 fOutput->Add(hPtMu_Match); 	
  
 fOutput->Add(hNEv);


 fOutput->Add(hNDimuons);	
 fOutput->Add(hNDimuons_2m);

 fOutput->Add(hAllpt_1M);	 	
 fOutput->Add(hLpt_1M);

 fOutput->Add(hAllpt_25y3);
 fOutput->Add(hLpt_25y3); 
 fOutput->Add(hAllpt_3y35); 
 fOutput->Add(hLpt_3y35);
 fOutput->Add(hAllpt_35y4); 
 fOutput->Add(hLpt_35y4);

 	 	
 fOutput->Add(hAllpt_25y3_1M);
 fOutput->Add(hLpt_25y3_1M);
 fOutput->Add(hAllpt_3y35_1M);
 fOutput->Add(hLpt_3y35_1M); 	 	
 fOutput->Add(hAllpt_35y4_1M);	 	
 fOutput->Add(hLpt_35y4_1M); 
 	 	
 fOutput->Add(hAllpt_010_1M);
 fOutput->Add(hLpt_010_1M);
 fOutput->Add(hAllpt_1020_1M);
 fOutput->Add(hLpt_1020_1M);
 fOutput->Add(hAllpt_2030_1M);
 fOutput->Add(hLpt_2030_1M);
 fOutput->Add(hAllpt_3040_1M);
 fOutput->Add(hLpt_3040_1M);
 fOutput->Add(hAllpt_4050_1M);
 fOutput->Add(hLpt_4050_1M);
 fOutput->Add(hAllpt_5060_1M);
 fOutput->Add(hLpt_5060_1M);
 fOutput->Add(hAllpt_6070_1M);
 fOutput->Add(hLpt_6070_1M);
 fOutput->Add(hAllpt_7080_1M);
 fOutput->Add(hLpt_7080_1M);
 fOutput->Add(hAllpt_8090_1M);
 fOutput->Add(hLpt_8090_1M);
	
 fOutput->ls(); 
 
 PostData(1,fOutput); 
 
} 

//_________________________________________________
void AliAnalysisTaskPbPbJPsi::UserExec(Option_t *)
{

//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }
  
  char hname[200];
  AliAODHeader *aodheader=dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
  TString firedtrigger = aodheader->GetFiredTriggerClasses();
  
  Bool_t TriggerSelected=kFALSE;
  Bool_t TriggerSelected_CINT7=kFALSE;
  Bool_t TriggerSelected_CMUL7=kFALSE;
  Bool_t TriggerSelected_CMLL7=kFALSE;
  Bool_t TriggerSelected_CMSL7=kFALSE;
  Bool_t TriggerSelected_CMSH7=kFALSE;

   if(firedtrigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE; 
    else TriggerSelected = kFALSE; 
//    if(firedtrigger.Contains("CINT7-B-NOPF-CENT") || firedtrigger.Contains("CV0L7-B-NOPF-CENT")) TriggerSelected_CINT7 = kTRUE; 
    if(firedtrigger.Contains("CINT7-B-NOPF-MUFAST")) TriggerSelected_CINT7 = kTRUE;
//     if(firedtrigger.Contains("CINT7-B-NOPF-MUFAST) || firedtrigger.Contains("CINT7-B-NOPF-CENT") || firedtrigger.Contains("CINT5-B-NOPF-CENT")
//       || firedtrigger.Contains("CINT10-B-NOPF-CENTNOTRD")|| firedtrigger.Contains("CTRUE-B-NOPF-CENT")|| firedtrigger.Contains("C0TVX-B-NOPF-CENT")
//       || firedtrigger.Contains("C0TVX-B-NOPF-CENTNOTRD")|| firedtrigger.Contains("CINT7EJ1-B-NOPF-CENTNOPMD")) TriggerSelected_CINT7 = kTRUE; 
    else TriggerSelected_CINT7 = kFALSE; 
    if(firedtrigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected_CMUL7 = kTRUE; 
    else TriggerSelected_CMUL7 = kFALSE; 
    if(firedtrigger.Contains("CMLL7-B-NOPF-MUFAST")) TriggerSelected_CMLL7 = kTRUE; 
    else TriggerSelected_CMLL7 = kFALSE; 
    if(firedtrigger.Contains("CMSL7-B-NOPF-MUFAST")) TriggerSelected_CMSL7 = kTRUE; 
    else TriggerSelected_CMSL7 = kFALSE; 
    if(firedtrigger.Contains("CMSH7-B-NOPF-MUFAST")) TriggerSelected_CMSH7 = kTRUE; 
    else TriggerSelected_CMSH7 = kFALSE; 
  
  fCountTotEv++;
  if(TriggerSelected) fCountTrigger++;
  if(TriggerSelected_CINT7) fCountCINT7++;
  if(TriggerSelected_CMUL7) fCountCMUL7++;
  if(TriggerSelected_CMLL7) fCountCMLL7++;
  if(TriggerSelected_CMSL7) fCountCMSL7++;
  if(TriggerSelected_CMSH7) fCountCMSH7++;
  
  // centrality 
  // https://twiki.cern.ch/twiki/bin/view/ALICE/CentralityCodeSnippets#Tables_with_centrality_bins_for
  AliMultSelection *MultSelection = (AliMultSelection * ) fAODEvent->FindListObject("MultSelection");
  Float_t PercV0M = MultSelection->GetMultiplicityPercentile("V0M");  //VOM
  ((TH1F*)(fOutput->FindObject("hPercentileV0M")))->Fill(PercV0M);
  if(TriggerSelected_CINT7) ((TH1F*)(fOutput->FindObject("hPercentileV0M_CINT7")))->Fill(PercV0M);
  if(TriggerSelected_CMUL7) ((TH1F*)(fOutput->FindObject("hPercentileV0M_CMUL7")))->Fill(PercV0M);
  
  // general infos
  Double_t ZVertex;
  if(TriggerSelected){
    AliAODVertex *PrimVertex =  fAODEvent->GetPrimaryVertex();
    ZVertex = PrimVertex->GetZ();
    ((TH1D*)(fOutput->FindObject("hZVertex")))->Fill(PrimVertex->GetZ());
    ((TH1D*)(fOutput->FindObject("hNContrib")))->Fill(PrimVertex->GetNContributors());
    ((TH1D*)(fOutput->FindObject("hNDimuons")))->Fill(fAODEvent->GetNumberOfDimuons());
    ((TH1D*)(fOutput->FindObject("hNMu")))->Fill(fAODEvent->GetNumberOfMuonTracks());
  }  
  
//   // loop on single muons
   Int_t mumatch=0;
   for(int i=0;i<fAODEvent->GetNumberOfTracks();i++){
     AliAODTrack *track=dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
     if( ! fMuonTrackCuts->IsSelected(track) ) continue;              
     if(track->IsMuonTrack()){
      if(TriggerSelected_CINT7){
     if(track->Eta()>-4 && track->Eta()<-2.5){
       if(track->GetRAtAbsorberEnd()>17.6 && track->GetRAtAbsorberEnd()){

     if(track->Eta()>-3 && track->Eta()<-2.5){
		       if(track->GetMatchTrigger()>=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_25y3")))->Fill(track->Pt());                                               
                                  }
		       if(track->GetMatchTrigger()>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_25y3")))->Fill(track->Pt());                                               
                                  }
                                 }//2.5<eta<3

     if(track->Eta()>-3.5 && track->Eta()<-3){
		       if(track->GetMatchTrigger()>=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_3y35")))->Fill(track->Pt());                                               
                                  }
		       if(track->GetMatchTrigger()>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_3y35")))->Fill(track->Pt());                                               
                                  }
                                 }//3<eta<3.5

     if(track->Eta()>-4 && track->Eta()<-3.5){
		       if(track->GetMatchTrigger()>=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_35y4")))->Fill(track->Pt());                                               
                                  }
		       if(track->GetMatchTrigger()>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_35y4")))->Fill(track->Pt());                                               
                                  }
                                 }//3.5<eta<4

                    }//Rabs
                  }//eta cut
         if(track->GetMatchTrigger()>1) { 
 	   mumatch++;
	   printf("mumatch in loop =%d\n",mumatch);
 	  ((TH1D*)(fOutput->FindObject("hPtMu_Match")))->Fill(track->Pt());
        }
      }//loop on CINT7
     }//muon track   
   }//loop on tracks   
   if(mumatch!=0)printf("mumatch = %d\n",mumatch);
    ((TH1D*)(fOutput->FindObject("hNMuons_1m")))->Fill((Double_t)mumatch);       
 
  Int_t ndimu=0;  
  // loop on dimuons

  for(int nd=0;nd<fAODEvent->GetNumberOfDimuons();nd++){
    AliAODDimuon *dimu = dynamic_cast<AliAODDimuon*>(fAODEvent->GetDimuon(nd));
    AliAODTrack *mu0 = dimu->GetMu(0); 
    AliAODTrack *mu1 = dimu->GetMu(1); 

  //ADD THE FOLLOWING 2 LINES!!!
    if (!fMuonTrackCuts->IsSelected(mu0)) continue;
    if (!fMuonTrackCuts->IsSelected(mu1)) continue;  
    
    Double_t DimuMass=999;
    Double_t DimuPt=-999;
    Double_t DimuY=-999;
    Double_t Match_Mu0=-999;
    Double_t Match_Mu1=-999;
    Double_t Pt_Mu0=-999;
    Double_t Pt_Mu1=-999;
    Double_t Eta_Mu0=-999;
    Double_t Eta_Mu1=-999;
    Double_t RAbs_Mu0=-999;
    Double_t RAbs_Mu1=-999;	  
    
    Match_Mu0=mu0->GetMatchTrigger();
    Match_Mu1=mu1->GetMatchTrigger();
    Eta_Mu0=mu0->Eta();
    Eta_Mu1=mu1->Eta();
    RAbs_Mu0=mu0->GetRAtAbsorberEnd();
    RAbs_Mu1=mu1->GetRAtAbsorberEnd();       
    Pt_Mu0=mu0->Pt();
    Pt_Mu1=mu1->Pt();
    
    ((TH1D*)(fOutput->FindObject("hMu0MatchTrigger")))->Fill(Match_Mu0);
    
     if((Eta_Mu0>-4 && Eta_Mu0<-2.5) && (Eta_Mu1>-4 && Eta_Mu1<-2.5)){
       if((RAbs_Mu0>17.6 && RAbs_Mu0<89.5) && (RAbs_Mu1>17.6 && RAbs_Mu1<89.5)){
         DimuY=dimu->Y();
         if(DimuY>-4 && DimuY<-2.5){
            DimuMass=dimu->Mass();
            if(Match_Mu0>1 && Match_Mu1>1){  
              if(TriggerSelected_CMLL7){
	       if(dimu->Charge()>0)((TH1D*)(fOutput->FindObject("hMassPP_CMLL7_2m")))->Fill(DimuMass);    
 	       if(dimu->Charge()<0)((TH1D*)(fOutput->FindObject("hMassMM_CMLL7_2m")))->Fill(DimuMass);  
// 	       ndimu++;
	      }
               ndimu++;
 	   }  
 	  
            if(dimu->Charge()==0) {
              DimuPt=dimu->Pt();
 	       //if(DimuPt>15) continue;
		
    		((TH1D*)(fOutput->FindObject("hMassOS_notrig")))->Fill(DimuMass);
    		if(TriggerSelected)((TH1D*)(fOutput->FindObject("hMassOS")))->Fill(DimuMass);
    		if(TriggerSelected_CINT7) ((TH1D*)(fOutput->FindObject("hMassOS_CINT7")))->Fill(DimuMass);
    		if(TriggerSelected_CMUL7) ((TH1D*)(fOutput->FindObject("hMassOS_CMUL7")))->Fill(DimuMass);
        	
    		if(Match_Mu0>1 && Match_Mu1>1){  
    		  if(TriggerSelected){
		    ((TH1D*)(fOutput->FindObject("hMassOS_2m")))->Fill(DimuMass);
		  }  
    		  if(TriggerSelected_CINT7) {
		    ((TH1D*)(fOutput->FindObject("hMassOS_CINT7_2m")))->Fill(DimuMass);
		  }  
		  if(TriggerSelected_CMSL7) ((TH1D*)(fOutput->FindObject("hMassOS_CMSL7_2m")))->Fill(DimuMass);
     		  if(TriggerSelected_CMUL7) {
		  

		  ((TH1D*)(fOutput->FindObject("hMassOS_CMUL7_2m")))->Fill(DimuMass);
		  ((TH1D*)(fOutput->FindObject("hPt_CMUL7_2m")))->Fill(DimuPt);
		  ((TH1D*)(fOutput->FindObject("hY_CMUL7_2m")))->Fill(DimuY);
 		  } //CMUL7 loop	   
 		} //2 mu matching loop
  
   	       if (Match_Mu0>1 || Match_Mu1>1){  
    		  if(TriggerSelected)((TH1D*)(fOutput->FindObject("hMassOS_1m")))->Fill(DimuMass);
    		  if(TriggerSelected_CINT7) ((TH1D*)(fOutput->FindObject("hMassOS_CINT7_1m")))->Fill(DimuMass);
    		  if(TriggerSelected_CMUL7) ((TH1D*)(fOutput->FindObject("hMassOS_CMUL7_1m")))->Fill(DimuMass);
		  if(TriggerSelected_CMSL7) ((TH1D*)(fOutput->FindObject("hMassOS_CMSL7_1m")))->Fill(DimuMass);
               }

    		  if(TriggerSelected_CINT7) {

		    if(DimuPt>0 && DimuPt<50){ 

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_1M")))->Fill(Pt_Mu0);
                                        }
		       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_1M")))->Fill(Pt_Mu0);
                                     }
		       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_1M")))->Fill(Pt_Mu1);                                               
                                  }

         if(DimuY>-3 && DimuY<-2.5){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_25y3_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_25y3_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_25y3_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_25y3_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//2.5<y<3

         if(DimuY>-3.5 && DimuY<-3){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_3y35_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_3y35_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_3y35_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_3y35_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//3<y<3.5

         if(DimuY>-4 && DimuY<-3.5){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_35y4_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_35y4_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_35y4_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_35y4_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//3.5<y<4

         if(PercV0M<90){

         if(PercV0M>=0 && PercV0M<10){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_010_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_010_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_010_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_010_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//010

         if(PercV0M>=10 && PercV0M<20){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_1020_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_1020_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_1020_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_1020_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//1020

         if(PercV0M>=20 && PercV0M<30){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_2030_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_2030_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_2030_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_2030_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//2030

         if(PercV0M>=30 && PercV0M<40){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_3040_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_3040_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_3040_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_3040_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//3040

         if(PercV0M>=40 && PercV0M<50){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_4050_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_4050_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_4050_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_4050_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//4050

         if(PercV0M>=50 && PercV0M<60){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_5060_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_5060_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_5060_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_5060_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//5060

         if(PercV0M>=60 && PercV0M<70){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_6070_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_6070_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_6070_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_6070_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//6070

         if(PercV0M>=70 && PercV0M<80){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_7080_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_7080_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_7080_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_7080_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//7080

         if(PercV0M>=80 && PercV0M<90){

		       if(Match_Mu0 >=1){   
                            ((TH1D*)(fOutput->FindObject("hAllpt_8090_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=1){
                            ((TH1D*)(fOutput->FindObject("hAllpt_8090_1M")))->Fill(Pt_Mu1);                                               
                                  }
		       if(Match_Mu0>=2){   
                            ((TH1D*)(fOutput->FindObject("hLpt_8090_1M")))->Fill(Pt_Mu0);
                                    }
                       if(Match_Mu1>=2){
                            ((TH1D*)(fOutput->FindObject("hLpt_8090_1M")))->Fill(Pt_Mu1);                                               
                                  }

                                 }//8090
                               }//<90
                              }//DimuPt cut
                          }//MB trigger
              } //charge
            } //Y
          } //Rabs
        } //eta
      } //dimuon loop

   if(ndimu!=0) printf("ndimu=%d\n",ndimu);
   ((TH1D*)(fOutput->FindObject("hNDimuons_2m")))->Fill((Double_t)ndimu);

((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(1,(Double_t)fCountTotEv);
((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(2,(Double_t)fCountCINT7);
((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(3,(Double_t)fCountCMUL7);
((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(4,(Double_t)fCountCMLL7);
((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(5,(Double_t)fCountCMSL7);
((TH1D*)(fOutput->FindObject("hNEv")))->SetBinContent(6,(Double_t)fCountCMSH7);

PostData(1,fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskPbPbJPsi::Terminate(Option_t *) 
{
//
// Draw histos
//
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  Int_t xmin=20; 
  Int_t ymin=20;
  
  printf("Using beam Energy=%f \n",fBeamEnergy);

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
//  fOutput = static_cast<TList*> (GetOutputData(1));
  TH1D *hMassOS_1m = static_cast<TH1D*> (fOutput->FindObject("hMassOS_1m"));    
  TH1D *hMassOS_2m = static_cast<TH1D*> (fOutput->FindObject("hMassOS_2m"));    
  TH1D *hNEv = static_cast<TH1D*> (fOutput->FindObject("hNEv"));  
  
  xmin+=20; ymin+=20;
  TCanvas *c2_MuonDistributions = new TCanvas("c2_MuonDistributions","OS Mass",xmin,ymin,600,600);
  c2_MuonDistributions->Divide(2,2);
  c2_MuonDistributions->cd(1);
  gPad->SetLogy(1);
  hMassOS_1m->Draw("e");
  c2_MuonDistributions->cd(2);
  gPad->SetLogy(1);
  hMassOS_2m->Draw("e");
  c2_MuonDistributions->cd(3);
  hNEv->Draw();
  hNEv->SetFillColor(kOrange);
 }


