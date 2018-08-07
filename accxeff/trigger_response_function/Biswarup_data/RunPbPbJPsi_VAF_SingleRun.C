void RunPbPbJPsi_VAF_SingleRun(
             Int_t runNumber = 195949, 
	     TString period="LHC15o",
             TString datasetpath1="BasePath=/alice/data/2015/",
             TString datasetpath2="/%09d/muon_calo_pass1/AOD/%/;",
             TString filename="aod_archive.zip;",
             TString anchor="AliAOD.Muons.root;",
             TString tree="aodTree;",
             TString mode="remote",
             Bool_t usePhysicsSelection = kTRUE,
             Int_t numEvents = 999999999,
             Int_t firstEvent = 0)
{
  TString datasetWithRun;
  datasetWithRun.Form(datasetpath2.Data(), runNumber);

  TString dataset = "Find;"+datasetpath1+period+datasetWithRun+"FileName="+filename+"Anchor="+anchor
  +"Tree=/"+tree+"Mode="+mode;
  printf("Dataset= %s\n",dataset.Data());
     
  TString extraLibs = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice:ANALYSISaliceBase:CORRFW:OADB"; // extraLibs = "ANALYSIS:OADB:ANALYSISalice:CORRFW:OADB:PWGmuon";

  TList *list = new TList(); 
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS", extraLibs.Data()));
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));  // important: creates token on every PROOF worker

  TProof::Open("pod://");

  // Check the dataset before running the analysis!
  gProof->ShowDataSet( dataset.Data() );
  //return;  // <-- uncomment this to test search before running the analysis!

  // A single AliRoot package for *all* AliRoot versions: new on VAF
  TFile::Cp("http://alibrary.web.cern.ch/alibrary/vaf/AliceVaf.par", "AliceVaf.par");
  gProof->UploadPackage("AliceVaf.par");
  gProof->EnablePackage("AliceVaf.par", list);  // this "list" is the same as always

  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train");

  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

//-----------------------------------------------------------------------------------------------------------------------
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"); //added for centrality
  AliMultSelectionTask *task = AddTaskMultSelection(kFALSE);
  task->SetUseDefaultCalib(kTRUE);
//-----------------------------------------------------------------------------------------------------------------------

  gROOT->LoadMacro("AddTaskOnTheFlyAliAODDimuon.C");   //to create Dimuon branch    

  gProof->Load("AliAnalysisTaskPbPbJPsi.cxx+");  // DON'T use double '+' when running multiple times: it uselessly recompiles everything!
  gROOT->LoadMacro("AddTaskPbPbJPsi.C");

  AliAnalysisTaskPbPbJPsi *PbPbJPsiTask = AddTaskPbPbJPsi(runNumber,period,usePhysicsSelection);
  if (usePhysicsSelection) {
//    PbPbJPsiTask->SelectCollisionCandidates(AliVEvent::kINT7);
//PbPbJPsiTask->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kINT7 | AliVEvent::kMuonLikePB | AliVEvent::kMuonSingleLowPt7); 
//PbPbJPsiTask->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kINT7 | AliVEvent::kMuonLikePB); 
//PbPbJPsiTask->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kINT7 | AliVEvent::kMuonLikePB | AliVEvent::kMuonSingleLowPt7 | AliVEvent::kMUS7); 
//PbPbJPsiTask->SelectCollisionCandidates(AliVEvent::kMuonUnlikePB | AliVEvent::kMuonLikePB | AliVEvent::kMuonSingleLowPt7 | AliVEvent::kMUS7); 
}
 
  if (!mgr->InitAnalysis()) return;

  mgr->StartAnalysis("proof", dataset, numEvents, firstEvent);

}
