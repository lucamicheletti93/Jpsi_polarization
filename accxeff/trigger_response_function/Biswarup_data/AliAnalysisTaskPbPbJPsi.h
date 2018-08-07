#ifndef AliAnalysisTaskPbPbJPsi_H
#define AliAnalysisTaskPbPbJPsi_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TObjArray;
class AliVParticle;
class AliAODEvent;
class TLorentzVector;
//ADD THE FOLLOWING LINE!!!
class AliMuonTrackCuts;
//class TList;

class AliAnalysisTaskPbPbJPsi : public AliAnalysisTaskSE {
  public:

  virtual void NotifyRun();          // Implement the Notify run to search for the new parameters at each new runs

  AliAnalysisTaskPbPbJPsi();
  AliAnalysisTaskPbPbJPsi(const char *name);
  virtual ~AliAnalysisTaskPbPbJPsi();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void Terminate(Option_t *);
  //virtual void NotifyRun(); 
  
  void SetBeamEnergy(Double_t en) {fBeamEnergy=en;}
  void SetAnalysisType(const char* type) {fkAnalysisType=type;}
  void SetPeriod(TString period) {fPeriod=period;}
    
 private:
  AliAnalysisTaskPbPbJPsi(const AliAnalysisTaskPbPbJPsi&);
  AliAnalysisTaskPbPbJPsi& operator=(const AliAnalysisTaskPbPbJPsi&);
     
 //protected:
     
  Double_t fBeamEnergy;   // Energy of the beam (required for the CS angle)    
  const char* fkAnalysisType; //ESD or AOD based analysis
  TString fPeriod; //period
  Int_t fCountTotEv;   // counter
  Int_t fCountTrigger;   // counter
  Int_t fCountCINT7;   // counter
  Int_t fCountCMUL7; //counter
  Int_t fCountCMLL7; //counter
  Int_t fCountCMSL7; //counter
  Int_t fCountCMSH7; //counter
  AliAODEvent* fAODEvent;      //! AOD event  //tolgo !
  AliMuonTrackCuts* fMuonTrackCuts;    // Use the class as a data member
  TList *fOutput;  //!< List of histograms for data
  
 ClassDef(AliAnalysisTaskPbPbJPsi,1);
};

#endif
