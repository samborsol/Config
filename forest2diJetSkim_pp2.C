#include <ctime>
#include <TMath.h>
#include "cutsAndBin_bgk.h"
#include <TLorentzVector.h>
#include "commonUtility.h"

#include <iostream>
#include <algorithm>
#include <vector>

using std::stringstream;
using std::vector;
using std::abs;
     
static const long MAXTREESIZE = 10000000000;

bool arrange(float i, float j) {return (i>j);}
TString getDayAndTime();

void forest2diJetSkim_pp2(
                           TString fname = "/u/user/bekim/758p3_UPC/src/Muon_test/4_UPCTriggers_pp_reco_170427.root",
			   TString outputFname = "upcDiJetSkim170919", 
			   TString trig = "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1",  // "HLT_HIUPCL1SingleEG5NotHF2_v1",    // "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1" //  "HLT_HIUPCL1SingleEG5NotHF2_v1"
			   TString jetCollection = "ak4PFJetAnalyzer", // "akPu5PFJetAnalyzer",
			   float minjPt = 0,
			   int nevt=-1
			   ) {  
  
  using namespace std;
  TFile *f1 = new TFile(fname.Data());

  TFile *f2 = new TFile("h_eff.root");
  TH1F *heff = (TH1F*)f2->Get("h4");
  
  TTree *HltTree = (TTree*)f1->Get("hltanalysis/HltTree");
  TTree *hiTree  = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");
//  TTree *trackTree = (TTree*)f1->Get("anaTrack/trackTree");
  TTree *t = (TTree*)f1->Get(Form("%s/t",jetCollection.Data()));
  TTree *trackTree  = (TTree*)f1->Get("ppTrack/trackTree");
  TTree *hbhe = (TTree*)f1->Get("rechitanalyzer/hbhe");
  TTree *hf = (TTree*)f1->Get("rechitanalyzer/hf");
  TTree *ee = (TTree*)f1->Get("rechitanalyzer/ee");
  TTree *eb = (TTree*)f1->Get("rechitanalyzer/eb");
  //  TTree *akpu5pf = (TTree*)f1->Get("t");
  //  TTree *pfTree  = (TTree*)f1->Get("pfcandAnalyzer/pfTree");
  TString dayTime = getDayAndTime();
  TFile* newfile = new TFile(Form("skimmedFiles/TrkCutsppreco_norechitcut_%s_trig%s_jetCollection%s_minJetPt%d_%s.root",outputFname.Data(), trig.Data(), jetCollection.Data(), (int)minjPt, dayTime.Data()  ),"recreate");
//  TFile* newfile = new TFile("skimmedFiles/test.root", "recreate");  

  t->AddFriend(HltTree);
  t->AddFriend(hiTree);
  t->AddFriend(trackTree);
  t->AddFriend(hbhe);
  t->AddFriend(hf);
  t->AddFriend(ee);
  t->AddFriend(eb);

//  hf->AddFriend(t);
//  hf->AddFriend(hiTree);
   // import the tree to the RooDataSet
   Int_t           trigBit;
   TBranch        *b_trigBit;
   if (trig == "" ) {  cout << " No Trigger selection! " << endl ;}     
   else { 
     HltTree->SetBranchAddress(trig.Data(), &trigBit, &b_trigBit);
   }

  Int_t           Run;
  TBranch        *b_Run;   //!
  HltTree->SetBranchAddress("Run", &Run, &b_Run);
  Int_t           LumiBlock;
  TBranch        *b_LumiBlock;   //!
  HltTree->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
  ULong64_t       Event;
  TBranch        *b_Event;   //!
  HltTree->SetBranchAddress("Event", &Event, &b_Event);

  UInt_t          run;
  ULong64_t       evt;
  UInt_t          lumi;
  TBranch        *b_run;
  TBranch        *b_evt;
  TBranch        *b_lumi;
  hiTree->SetBranchAddress("run", &run, &b_run);
  hiTree->SetBranchAddress("evt", &evt, &b_evt);
  hiTree->SetBranchAddress("lumi", &lumi, &b_lumi);

  /// Photon inputs : 
/*  Int_t           nPFpart;
  vector<int>     *pfId;
  vector<float>   *pfPt;
  vector<float>   *pfEnergy;
  vector<float>   *pfVsPtInitial;
  vector<float>   *pfEta;
  vector<float>   *pfPhi;
  nPFpart = 0;
  pfId = 0;
  pfPt = 0;
  pfEnergy = 0;
  pfVsPtInitial = 0;
  pfEta = 0;
  pfPhi = 0;

  TBranch        *b_nPFpart;   //!
  TBranch        *b_pfId;   //!
  TBranch        *b_pfPt;   //!
  TBranch        *b_pfEnergy;   //!
  TBranch        *b_pfVsPtInitial;   //!
  TBranch        *b_pfEta;   //!
  TBranch        *b_pfPhi;   //!
  
   pfTree->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
   pfTree->SetBranchAddress("pfId", &pfId, &b_pfId);
   pfTree->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
   pfTree->SetBranchAddress("pfEnergy", &pfEnergy, &b_pfEnergy);
   pfTree->SetBranchAddress("pfVsPtInitial", &pfVsPtInitial, &b_pfVsPtInitial);
   pfTree->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
   pfTree->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);
*/
  // Jet inputs :
  Int_t nref;
  Float_t jtpt[200] = {0};
  Float_t jteta[200] = {0};
  Float_t jty[200] = {0};
  Float_t jtphi[200] = {0};
  Float_t jtm[200] = {0};
  Float_t jtmass[200] = {0};
  Float_t pt[200] = {0};
  Float_t eta[200] = {0};
  Float_t y[200] = {0};
  Float_t phi[200] = {0};
  Float_t m[200] = {0};
//  Float_t e[200] = {0};
  Float_t mass[200] = {0}; 
  Float_t jtpt2[200] = {0};
  Float_t jteta2[200] = {0};
  Float_t jty2[200] = {0};
  Float_t jtphi2[200] = {0};
  Float_t jtm2[200] = {0};
  Float_t jtmass2[200] = {0};
  Float_t pt2[200] = {0};
  Float_t eta2[200] = {0};
  Float_t y2[200] = {0};
  Float_t phi2[200] = {0};
  Float_t m2[200] = {0};
  Float_t e2[200] = {0};
  Float_t mass2[200] = {0};

  TBranch *b_nref;
  TBranch *b_jtpt;
  TBranch *b_jteta;
  TBranch *b_jty;
  TBranch *b_jtphi;
  TBranch *b_jtm;

  t->SetBranchAddress("nref", &nref, &b_nref);
  t->SetBranchAddress("jtpt", jtpt, &b_jtpt);
  t->SetBranchAddress("jteta", jteta, &b_jteta);
  t->SetBranchAddress("jty", jty, &b_jty);
  t->SetBranchAddress("jtphi", jtphi, &b_jtphi);
  t->SetBranchAddress("jtm", jtm, &b_jtm);

   // HiTree inputs : 
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiHFplusEta4;   //!
   TBranch        *b_hiHFminusEta4;   //!
   hiTree->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   hiTree->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   hiTree->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   hiTree->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   hiTree->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   hiTree->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);

   Int_t           nTrk;
   Float_t         trkEta[33885];
   Float_t         trkPhi[33885];
   Float_t         trkPt[33885];
   Float_t         trkPtError[33885];
   Bool_t          highPurity[33885];
   Float_t         trkDxy1[33885];   //[nTrk]
   Float_t         trkDxyError1[33885];   //[nTrk]
   Float_t         trkDz1[33885];   //[nTrk]
   Float_t         trkDzError1[33885];   //[nTrk]
   UChar_t         trkNHit[33885];
   Float_t         trkChi2[33885];
   UChar_t         trkNdof[33885];
   UChar_t         trkNlayer[33885];
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkEta;
   TBranch        *b_trkPhi;
   TBranch        *b_trkPt;
   TBranch        *b_trkPtError;
   TBranch        *b_highPurity;
   TBranch        *b_trkDxy1;
   TBranch        *b_trkDxyError1;
   TBranch        *b_trkDz1;
   TBranch        *b_trkDzError1;
   TBranch        *b_trkNHit;
   TBranch        *b_trkChi2;
   TBranch        *b_trkNdof;
   TBranch        *b_trkNlayer;
   trackTree->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   trackTree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   trackTree->SetBranchAddress("trkPhi", trkPhi, &b_trkPhi);
   trackTree->SetBranchAddress("trkPt", trkPt, &b_trkPt);
   trackTree->SetBranchAddress("trkPtError", trkPtError, &b_trkPtError);
   trackTree->SetBranchAddress("highPurity", highPurity, &b_highPurity);
   trackTree->SetBranchAddress("trkDxy1", trkDxy1, &b_trkDxy1);
   trackTree->SetBranchAddress("trkDxyError1", trkDxyError1, &b_trkDxyError1);
   trackTree->SetBranchAddress("trkDz1", trkDz1, &b_trkDz1);
   trackTree->SetBranchAddress("trkDzError1", trkDzError1, &b_trkDzError1);
   trackTree->SetBranchAddress("trkNHit", trkNHit, &b_trkNHit);
   trackTree->SetBranchAddress("trkChi2", trkChi2, &b_trkChi2);
   trackTree->SetBranchAddress("trkNdof", trkNdof, &b_trkNdof);
   trackTree->SetBranchAddress("trkNlayer", trkNlayer, &b_trkNlayer);

   Int_t           hben;
   Float_t         hbee[1000];   
   Float_t         hbeet[1000];   
   Float_t         hbeeta[1000];   
   Float_t         hbephi[1000];   
   Float_t         hbeperp[1000];   
   Bool_t          hbeisjet[1000];   
   Int_t           hbedepth[1000];   
   TBranch        *b_hben;   
   TBranch        *b_hbee;   
   TBranch        *b_hbeet;   
   TBranch        *b_hbeeta;   
   TBranch        *b_hbephi;   
   TBranch        *b_hbeperp;   
   TBranch        *b_hbeisjet;   
   TBranch        *b_hbedepth;   
   hbhe->SetBranchAddress("n", &hben, &b_hben);
   hbhe->SetBranchAddress("e", hbee, &b_hbee);
   hbhe->SetBranchAddress("et", hbeet, &b_hbeet);
   hbhe->SetBranchAddress("eta", hbeeta, &b_hbeeta);
   hbhe->SetBranchAddress("phi", hbephi, &b_hbephi);
   hbhe->SetBranchAddress("perp", hbeperp, &b_hbeperp);
   hbhe->SetBranchAddress("isjet", hbeisjet, &b_hbeisjet);
   hbhe->SetBranchAddress("depth", hbedepth, &b_hbedepth); 

   Int_t           hfn;
   Float_t         hfe[1000];
   Float_t         hfet[1000];
   Float_t         hfeta[1000];
   Float_t         hfphi[1000];
   Float_t         hfperp[1000];
   Bool_t          hfisjet[1000];
   Int_t           hfdepth[1000];
   Float_t         totalEtHFplus;
   Float_t         totalEtHFminus;
   Float_t         totalEtHF;
   TBranch        *b_hfn;
   TBranch        *b_hfe;
   TBranch        *b_hfet;
   TBranch        *b_hfeta;
   TBranch        *b_hfphi;
   TBranch        *b_hfperp;
   TBranch        *b_hfisjet;
   TBranch        *b_hfdepth;
   TBranch        *b_totalEtHFplus;
   TBranch        *b_totalEtHFminus;
   TBranch        *b_totalEtHF;
   hf->SetBranchAddress("n", &hfn, &b_hfn);
   hf->SetBranchAddress("e", hfe, &b_hfe);
   hf->SetBranchAddress("et", hfet, &b_hfet);
   hf->SetBranchAddress("eta", hfeta, &b_hfeta);
   hf->SetBranchAddress("phi", hfphi, &b_hfphi);
   hf->SetBranchAddress("perp", hfperp, &b_hfperp);
   hf->SetBranchAddress("isjet", hfisjet, &b_hfisjet);
   hf->SetBranchAddress("depth", hfdepth, &b_hfdepth);
   hf->SetBranchAddress("totalEtHFplus",&totalEtHFplus,&b_totalEtHFplus);
   hf->SetBranchAddress("totalEtHFminus",&totalEtHFminus,&b_totalEtHFminus);
   hf->SetBranchAddress("totalEtHF",&totalEtHF,&b_totalEtHF);

   Int_t           een;
   Float_t         eee[1000];
   Float_t         eeet[1000];
   Float_t         eeeta[1000];
   Float_t         eephi[1000];
   Float_t         eeperp[1000];
   Bool_t          eeisjet[1000];
   TBranch        *b_een;
   TBranch        *b_eee;
   TBranch        *b_eeet;
   TBranch        *b_eeeta;
   TBranch        *b_eephi;
   TBranch        *b_eeperp;
   TBranch        *b_eeisjet;
   ee->SetBranchAddress("n", &een, &b_een);
   ee->SetBranchAddress("e", eee, &b_eee);
   ee->SetBranchAddress("et", eeet, &b_eeet);
   ee->SetBranchAddress("eta", eeeta, &b_eeeta);
   ee->SetBranchAddress("phi", eephi, &b_eephi);
   ee->SetBranchAddress("perp", eeperp, &b_eeperp);
   ee->SetBranchAddress("isjet", eeisjet, &b_eeisjet);

   Int_t           ebn;
   Float_t         ebe[1000];
   Float_t         ebet[1000];
   Float_t         ebeta[1000];
   Float_t         ebphi[1000];
   Float_t         ebperp[1000];
   Bool_t          ebisjet[1000];
   TBranch        *b_ebn;
   TBranch        *b_ebe;
   TBranch        *b_ebet;
   TBranch        *b_ebeta;
   TBranch        *b_ebphi;
   TBranch        *b_ebperp;
   TBranch        *b_ebisjet;
   eb->SetBranchAddress("n", &ebn, &b_ebn);
   eb->SetBranchAddress("e", ebe, &b_ebe);
   eb->SetBranchAddress("et", ebet, &b_ebet);
   eb->SetBranchAddress("eta", ebeta, &b_ebeta);
   eb->SetBranchAddress("phi", ebphi, &b_ebphi);
   eb->SetBranchAddress("perp", ebperp, &b_ebperp);
   eb->SetBranchAddress("isjet", ebisjet, &b_ebisjet);

  ////////////////////////////////////////////////////////////////////////
  //////////////////  dijet tree 
  ////////////////////////////////////////////////////////////////////////
   UPCEvent event;
   TTree *eventTree = new TTree("evt","Event Tree");
   eventTree->SetMaxTreeSize(MAXTREESIZE);
   eventTree->Branch("event",&event.run,eventBranchString.Data());
//   eventTree->Branch("

/*   DiPhoPF dp;
   TTree *dpTree = new TTree("diphoton","Photon Tree");
   dpTree->SetMaxTreeSize(MAXTREESIZE);
   dpTree->Branch("dp",&dp.mass,dpBranchStringPF);
*/
   UPCdiJet dj;
   TTree *djTree = new TTree("dijet","Jet Tree");
   djTree->SetMaxTreeSize(MAXTREESIZE);
   djTree->Branch("dj",&dj.nJet,djBranchString);

   UPCnTrk Track;
   TTree *trkTree = new TTree("Track","Track Tree");
   trkTree->SetMaxTreeSize(MAXTREESIZE);
   trkTree->Branch("trkVar",&Track.nTrack,nTrkString);

   int ntrk;
   float floatntrk;
   const int MAXtrk = 50000; // to accomodate 100 smeared trks, need to be careful with ram
   float trkpt[MAXtrk];
   float trketa[MAXtrk];
   float trkphi[MAXtrk];
   float trkHFplus;
   float trkHFminus;

   TTree *newtrkTree = new TTree("fullTrkTree","Track Tree 2");
   newtrkTree->SetMaxTreeSize(MAXTREESIZE);
   newtrkTree->Branch("ntrk",&ntrk,"ntrk/I");
   newtrkTree->Branch("floatntrk",&floatntrk,"floatntrk/F");
   newtrkTree->Branch("pT",trkpt,"pT[ntrk]/F");
   newtrkTree->Branch("Eta",trketa,"Eta[ntrk]/F");
   newtrkTree->Branch("Phi",trkphi,"Phi[ntrk]/F");
   newtrkTree->Branch("trkHFplus",&trkHFplus,"trkHFplus/F");
   newtrkTree->Branch("trkHFminus",&trkHFminus,"trkHFminus/F");

   float floatnTrkabsEtaover1p5;
   float hemax;

   TTree *trkheTree = new TTree("trkandHE","TrackVsHE"); 
   trkheTree->SetMaxTreeSize(MAXTREESIZE);
   trkheTree->Branch("floatnTrkabsEtaover1p5",&floatnTrkabsEtaover1p5,"floatnTrkabsEtaover1p5/F");
   trkheTree->Branch("hemax",&hemax,"hemax/F");

   Int_t           HBHEn;
   Float_t         HBHEe[1000];   
   Float_t         HBHEet[1000];   
   Float_t         HBHEeta[1000];   
   Float_t         HBHEphi[1000];   
   Float_t         HBHEperp[1000];   
   Bool_t          HBHEisjet[1000];   
   Int_t           HBHEdepth[1000];   
   Float_t         HBmax;
   Float_t         HEmax;
   Int_t           HFn;
   Float_t         HFe[1000];
   Float_t         HFet[1000];
   Float_t         HFeta[1000];
   Float_t         HFphi[1000];
   Float_t         HFperp[1000];
   Bool_t          HFisjet[1000];
   Int_t           HFdepth[1000];
   Float_t         HFtotal;
   Float_t         HFplus;
   Float_t         HFminus;
   Float_t         HFmax;
   Int_t           EEn;
   Float_t         EEe[1000];
   Float_t         EEet[1000];
   Float_t         EEeta[1000];
   Float_t         EEphi[1000];
   Float_t         EEperp[1000];
   Bool_t          EEisjet[1000];
   Float_t         EEmax;
   Int_t           EBn;
   Float_t         EBe[1000];
   Float_t         EBet[1000];
   Float_t         EBeta[1000];
   Float_t         EBphi[1000];
   Float_t         EBperp[1000];
   Bool_t          EBisjet[1000];
   Float_t         EBmax;

   TTree *calTree = new TTree("Cal","Cal Tree");
   calTree->SetMaxTreeSize(MAXTREESIZE);
   calTree->Branch("HBHEn",&HBHEn,"HBHEn/I");
   calTree->Branch("HBHEe",HBHEe,"HBHEe/F");
   calTree->Branch("HBHEeta",HBHEeta,"HBHEeta/F");
   calTree->Branch("HBHEphi",HBHEphi,"HBHEphi/F");
   calTree->Branch("HBHEperp",HBHEperp,"HBHEperp/F");
   calTree->Branch("HBHEisjet",HBHEisjet,"HBHEisjet/F");
   calTree->Branch("HBHEdepth",HBHEdepth,"HBHEdepth/F"); 
   calTree->Branch("HBmax",&HBmax,"HBmax/F");
   calTree->Branch("HEmax",&HEmax,"HEmax/F");
   calTree->Branch("HFn",&HFn,"HFn/I");
   calTree->Branch("HFe",HFe,"HFe/F");
   calTree->Branch("HFeta",HFeta,"HFeta/F");
   calTree->Branch("HFphi",HFphi,"HFphi/F");
   calTree->Branch("HFperp",HFperp,"HFperp/F");
   calTree->Branch("HFisjet",HFisjet,"HFisjet/F");
   calTree->Branch("HFdepth",HFdepth,"HFdepth/F");
   calTree->Branch("HFtotal",&HFtotal,"HFtotal/F");
   calTree->Branch("HFplus",&HFplus,"HFplus/F");
   calTree->Branch("HFminus",&HFminus,"HFminus/F");
   calTree->Branch("HFmax",&HFmax,"HFmax/F");
   calTree->Branch("EEn",&EEn,"EEn/I");
   calTree->Branch("EEe",EEe,"EEe/F");
   calTree->Branch("EEeta",EEeta,"EEeta/F");
   calTree->Branch("EEphi",EEphi,"EEphi/F");
   calTree->Branch("EEperp",EEperp,"EEperp/F");
   calTree->Branch("EEisjet",EEisjet,"EEisjet/F");
   calTree->Branch("EEmax",&EEmax,"EEmax/F");
   calTree->Branch("EBn",&EBn,"EBn/I");
   calTree->Branch("EBe",EBe,"EBe/F");
   calTree->Branch("EBeta",EBeta,"EBeta/F");
   calTree->Branch("EBphi",EBphi,"EBphi/F");
   calTree->Branch("EBperp",EBperp,"EBperp/F");
   calTree->Branch("EBisjet",EBisjet,"EBisjet/F");
   calTree->Branch("EBmax",&EBmax,"EBmax/F");

   Int_t bin1;
   Int_t bin2;
   Float_t pT1;
   Float_t pT2;
   Float_t eff1;
   Float_t eff2;
   Float_t w1;
   Float_t w2;
   Float_t DjpT;

   TTree *CSTree = new TTree("CS","CS Tree");
   CSTree->SetMaxTreeSize(MAXTREESIZE);
   CSTree->Branch("bin1",&bin1,"bin1/I");
   CSTree->Branch("bin2",&bin2,"bin2/I");
   CSTree->Branch("pT1",&pT1,"pT1/F");
   CSTree->Branch("pT2",&pT2,"pT2/F");
   CSTree->Branch("eff1",&eff1,"eff1/F");
   CSTree->Branch("eff2",&eff2,"eff2/F");
   CSTree->Branch("w1",&w1,"w1/F");
   CSTree->Branch("w2",&w2,"w2/F");
   CSTree->Branch("DjpT",&DjpT,"DjpT/F");

   Int_t d_phi;

   vector<float> jetpt;
   vector<float> jeteta;
   vector<float> jetphi;
   vector<float> jetmass;
   vector<float>::iterator iter_jetpt;

   vector<float> HBe_cand;
   vector<float> HEe_cand;
   vector<float> HFe_cand;
   vector<float> EBe_cand;
   vector<float> EEe_cand;
   vector<float>::iterator iter_HBe_cand;
   vector<float>::iterator iter_HEe_cand;
   vector<float>::iterator iter_HFe_cand;
   vector<float>::iterator iter_EBe_cand;
   vector<float>::iterator iter_EEe_cand;

   Float_t HB_max;
   Float_t HE_max;
   Float_t HF_max;
   Float_t EB_max;
   Float_t EE_max;

   Int_t numjt;
   // event loop start
   //   if(nevt == -1) nevt = pfTree->GetEntries();   

   Int_t numevt = 0;
   if(nevt == -1) nevt = hf->GetEntries(); 
   for (int iev=0; iev<nevt; iev++) {
     if(iev%100000==0)
     { 
       cout << ">>>>> EVENT " << iev << " / " << hf->GetEntries() <<  " ("<<(int)(100.*iev/hf->GetEntries()) << "%)" << endl; 
     }
     HltTree->GetEntry(iev);

   
     if ( (trig != "" ) && (!trigBit) ) {
       continue;
     }
     
     t->GetEntry(iev);
     //hiTree->GetEntry(iev);
     // trigger selection 

     ///// Call the values /////
     pT1 = 0;
     pT2 = 0;
     w1 = 0;
     w2 = 0;
     eff1 = 0;
     eff2 = 0;
     bin1 = 0;
     bin2 = 0;
     DjpT = 0;   
     
     HFtotal = 0;
     HFplus = 0;
     HFminus = 0;    
 
     HB_max = 0;
     HE_max = 0;
     HF_max = 0;
     EB_max = 0;
     EE_max = 0;

     HBe_cand.clear();
     HEe_cand.clear();
     HFe_cand.clear();
     EBe_cand.clear();
     EEe_cand.clear();

     HBHEn = hben;
     for(Int_t d = 0; d != hben; d++)
     {
       HBHEe[d] = hbee[d];
       HBHEet[d] = hbeet[d];
       HBHEeta[d] = hbeeta[d];
       HBHEphi[d] = hbephi[d];
       HBHEperp[d] = hbeperp[d];
       HBHEisjet[d] = hbeisjet[d];
       HBHEdepth[d] = hbedepth[d];
       if(hbeeta[d] < 1.5)
       {
         HBe_cand.push_back(hbee[d]);
       }
       if(hbeeta[d] > 1.5 && hbeeta[d] < 3.0)
       {
         HEe_cand.push_back(hbee[d]);
       }
     }
     HFn = hfn;
     for(Int_t e = 0; e != hfn; e++)
     {
       HFe[e] = hfe[e];
       HFet[e] = hfet[e];
       HFeta[e] = hfeta[e];
       HFphi[e] = hfphi[e];
       HFperp[e] = hfperp[e];
       HFisjet[e] = hfisjet[e];
       HFdepth[e] = hfdepth[e];
       HFe_cand.push_back(hfe[e]);
     }
     HFtotal = totalEtHF;
     HFminus = totalEtHFminus;
     HFplus = totalEtHFplus;
     //cout << "evt : " << iev << " HF : " << HFtotal << endl;
     EEn = een;
     for(Int_t f = 0; f != een; f++)
     {
       EEe[f] = eee[f];
       EEet[f] = eeet[f];
       EEeta[f] = eeeta[f];
       EEphi[f] = eephi[f];
       EEperp[f] = eeperp[f];
       EEisjet[f] = eeisjet[f];
       EEe_cand.push_back(eee[f]);
     }
     EBn = ebn;
     for(Int_t g = 0; g != ebn; g++)
     {
       EBe[g] = ebe[g];
       EBet[g] = ebet[g];
       EBeta[g] = ebeta[g];
       EBphi[g] = ebphi[g];
       EBperp[g] = ebperp[g];
       EBisjet[g] = ebisjet[g];
       EBe_cand.push_back(ebe[g]);
     }
     
     iter_HBe_cand = max_element(HBe_cand.begin(), HBe_cand.end());
     iter_HEe_cand = max_element(HEe_cand.begin(), HEe_cand.end());
     iter_HFe_cand = max_element(HFe_cand.begin(), HFe_cand.end());
     iter_EEe_cand = max_element(EEe_cand.begin(), EEe_cand.end());
     iter_EBe_cand = max_element(EBe_cand.begin(), EBe_cand.end()); 
     if(HBe_cand.size() > 0)
     {
       HB_max = *iter_HBe_cand;
     }
     else
     {
       HB_max = 0;
     }
     if(HEe_cand.size() > 0)
     {
       HE_max = *iter_HEe_cand;
     }
     else
     {
       HE_max = 0;
     }
     if(HFe_cand.size() > 0)
     {
       HF_max = *iter_HFe_cand;
     }
     else
     {
       HF_max = 0;
     }
     if(EEe_cand.size() > 0)
     {
       EE_max = *iter_EEe_cand;
     }
     else
     {
       EE_max = 0;
     }
     if(EBe_cand.size() > 0)
     {
       EB_max = *iter_EBe_cand;
     }
     else
     {
       EB_max = 0;
     }

     HBmax = HB_max;
     HEmax = HE_max;
     hemax = HE_max;
     HFmax = HF_max;
     EEmax = EE_max;
     EBmax = EB_max;
     //cout << "HF_max : " << HF_max << " HE_max : " << HE_max << " HB_max : " << HB_max << endl;

     /*if(HF_max >= 3.0 || HE_max >= 1.95 || HB_max >= 1.18)
     {
       continue;
     }*/

     event.clear();
     event.run = run;
     event.lumi = lumi;
     event.event = evt; 
     event.vz = -99;
     event.hfsum = HFtotal; 
     event.hfplus = HFplus; 
     event.hfminus = HFminus;    

     if(Run==262988 && LumiBlock==663 && Event==89487086)
     {
       cout << "262988 : 663 : 89487086  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89469683)
     {
       cout << "262988 : 663 : 89469683  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89482692)
     {
       cout << "262988 : 663 : 89482692  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89472230)
     {
       cout << "262988 : 663 : 89472230  hfsum: " << HFtotal << endl;
     }
     if(Run==262988 && LumiBlock==663 && Event==89553231)
     {
       cout << "262988 : 663 : 89553231  hfsum: " << HFtotal << endl;
     }
     
     dj.clear();

     Track.clear();
     Track.nTrack = 0;
     ntrk = 0;
     floatntrk = 0;
     trkHFplus = HFplus;
     trkHFminus = HFminus;
     floatnTrkabsEtaover1p5 = 0;
     for(Int_t c = 0; c != nTrk; c++)
     {
       // Yongsun:  Here is where the quality cuts of trak will enter.
       //  pt cutk, pTerr/pt cut, etc etc in near fugure.  ask Daniel and Sasha which track quality cut you have you have to use.  
       if ( fabs(trkEta[c])> 2.4 ) 
	     continue;   // eta cut 
       //
       //
       //
       // From now you should not make any cut below :

       if(highPurity[c] == 1 && trkPtError[c]/trkPt[c] < 0.1 && fabs(trkDz1[c]/trkDzError1[c]) < 3 && fabs(trkDxy1[c]/trkDxyError1[c]) < 3 && trkNHit[c] >= 11 && trkChi2[c]/(float)trkNdof[c]/(float)trkNlayer[c] < 0.15 && trkPt[c] > 0.2)
       {
       trketa[ntrk] = trkEta[c];
       trkphi[ntrk] = trkPhi[c];
       trkpt[ntrk] = trkPt[c];
       ntrk = ntrk + 1;
       floatntrk = (float)ntrk;

       Track.nTrack = Track.nTrack + 1;

       if(fabs(trkEta[c]) > 1.5)
       {
         Track.nTrkabsEtaover1p5 = Track.nTrkabsEtaover1p5 + 1;
         floatnTrkabsEtaover1p5 = floatnTrkabsEtaover1p5 + 1;
       }
       if(fabs(trkEta[c]) < 1.5)
       {
         Track.nTrkabsEtaunder1p5 = Track.nTrkabsEtaunder1p5 + 1;
       }
       if(trkEta[c] >= -2.5 && trkEta[c] < -2.0)
       {
         Track.nTrketam2to2p5 = Track.nTrketam2to2p5 + 1;
       }
       else if(trkEta[c] >= -2 && trkEta[c] < -1.5)
       {
         Track.nTrketam1p5to2 = Track.nTrketam1p5to2 + 1; 
       }
       else if(trkEta[c] >= -1.5 && trkEta[c] < -1.0)
       {
         Track.nTrketam1to1p5 = Track.nTrketam1to1p5 + 1 ;
       }
       else if(trkEta[c] >= -1.0 && trkEta[c] < -0.5)   
       {
         Track.nTrketam0p5to1 = Track.nTrketam0p5to1 + 1;    
       }
       else if(trkEta[c] >= -0.5 && trkEta[c] < 0.0)
       {
         Track.nTrketam0to0p5 = Track.nTrketam0to0p5 + 1;
       }
       else if(trkEta[c] >= 0.0 && trkEta[c] < 0.5)
       {
         Track.nTrketa0to0p5 = Track.nTrketa0to0p5 + 1;
       }
       else if(trkEta[c] >= 0.5 && trkEta[c] < 1.0)
       {
         Track.nTrketa0p5to1 = Track.nTrketa0p5to1 + 1 ;
       }
       else if(trkEta[c] >= 1.0 && trkEta[c] < 1.5)
       {
         Track.nTrketa1to1p5 = Track.nTrketa1to1p5 +1 ;
       }
       else if(trkEta[c] >= 1.5 && trkEta[c] < 2.0)
       {
         Track.nTrketa1p5to2 = Track.nTrketa1p5to2 + 1 ;
       }
       else if(trkEta[c] >= 2.0 && trkEta[c] < 2.5)
       {
         Track.nTrketa2to2p5 = Track.nTrketa2to2p5 + 1 ;
       }
       }
       else
       {
         continue;
       }
     }

     numjt = 0;

     jetpt.clear();
     jeteta.clear();
     jetphi.clear();
     jetmass.clear();

     for(Int_t a = 0; a != nref; a++)      {
       TLorentzVector jt;
       jt.SetPtEtaPhiM( jtpt[a], jteta[a], jtphi[a], jtm[a] );
       // some jet cuts 
       //if(TMath::Abs(jteta[a]) < 2.4 && jtpt[a] > minjPt)
       //
       if(TMath::Abs(jteta[a]) > 2.4) // Should still have eta cut 
    	 continue ;
       
       // From here now, do not introduce any other cuts until the loop is over. 
       jetpt.push_back(jtpt[a]);
       jeteta.push_back(jteta[a]);
       jetphi.push_back(jtphi[a]);
       jetmass.push_back(jtmass[a]);
       numjt += 1;
     }
     
     if ( numjt == 0 )  { 
       jtpt2[0] = -1 ;
       jteta2[0] = -1 ;
       jtphi2[0] = -1 ;
       jtm2[0] = -1 ;
       jtpt2[1] = -1 ;
       jteta2[1] = -1 ;
       jtphi2[1] = -1 ;
       jtm2[1] = -1 ;

     }
     else if ( numjt == 1 ) { 
       jtpt2[0] = jetpt[0];
       jteta2[0] = jeteta[0];
       jtphi2[0] = jetphi[0];
       jtm2[0]   = jetmass[0];

       jtpt2[1] = -1 ;
       jteta2[1] = -1 ;
       jtphi2[1] = -1 ;
       jtm2[1] = -1 ;

     }
     else {     // if numjt > 1 
       sort(jetpt.begin(), jetpt.end(), arrange);
       //cout << numjt << " " << jetpt.size() << " " << jetpt[0] << " " << jetpt[1] << " " << jetpt[2] << " " << jetpt[3] << endl;
       //if(numjt > 0 && jetpt.size() > 0 && numjt == jetpt.size())
       //{
       //  cout << "DAMN!" << endl;
       //}
       for(Int_t b = 0; b != nref; b++) {
    	 if(jtpt[b] == jetpt[0])
  	   {
	       jtpt2[0] = jtpt[b];
	       jteta2[0] = jteta[b];
	       jtphi2[0] = jtphi[b];
  	     jtm2[0] = jtm[b];
	     }
    	 else if(jtpt[b] == jetpt[1])	
       {
  	     jtpt2[1] = jtpt[b];
	       jteta2[1] = jteta[b];
	       jtphi2[1] = jtphi[b];
         jtm2[1] = jtm[b];
	     }
     }
   }
//     jetpt.clear();
//     jeteta.clear();
//     jetphi.clear();
//     jetmass.clear();
      
     d_phi = 0;
     d_phi = TMath::Abs(getDPHI(jtphi2[0], jtphi2[1]));
      if(d_phi > TMath::Pi())
	{
	  d_phi = 2 * 3.141592653589 - d_phi;
	}
      //cout << "d_phi : " << d_phi << endl;
      if(d_phi > 3.141592653589)
	{
	  cout << "Not good" << endl;
	}
      //cout << Run << " " << LumiBlock << " " << Event << endl;
      pt[0] = jtpt2[0];
      eta[0] = jteta2[0];
      phi[0] = jtphi2[0];
      mass[0] = jtm2[0];
      pt[1] = jtpt2[1];
      eta[1] = jteta2[1];
      phi[1] = jtphi2[1];
      mass[1] = jtm2[1];
 
      //cout << pt[0] << "  " << eta[0] << "  " << pt[1] << "  " << eta[1] << endl;     
      TLorentzVector jt1vec, jt2vec, djvec;
      jt1vec.SetPtEtaPhiM( pt[0], eta[0], phi[0], mass[0] );
      jt2vec.SetPtEtaPhiM( pt[1], eta[1], phi[1], mass[1] );
      djvec = jt1vec + jt2vec ;
      
      dj.nJet = nref;
      dj.mass = djvec.M();
      dj.pt   = djvec.Pt();
      dj.y   = djvec.Rapidity();
      dj.phi   = djvec.Phi();
      dj.eta  = djvec.Eta();
      dj.dphi = TMath::Abs(getDPHI( phi[0], phi[1] ));
      dj.dpt = TMath::Abs(pt[0] - pt[1]);
      dj.deta = TMath::Abs(eta[0] - eta[1]);
      dj.aj = TMath::Abs(pt[0] - pt[1])/(pt[0] + pt[1]);
      dj.pt1 = pt[0];
      dj.eta1 = eta[0];
      dj.phi1 = phi[0];
      dj.e1 = jt1vec.E();
      dj.pt2 = pt[1];
      dj.eta2 = eta[1];
      dj.phi2 = phi[1];
      dj.e2 = jt2vec.E();
    
      pT1 = jtpt2[0];
      pT2 = jtpt2[1]; 
      bin1 = heff->FindBin(pt[0]);
      eff1 = heff->GetBinContent(bin1);
      w1 = 1/eff1;
      bin2 = heff->FindBin(pt[1]);
      eff2 = heff->GetBinContent(bin2);
      w2 = 1/eff2;
      DjpT = djvec.Pt();
   
      eventTree->Fill();
      djTree->Fill();
      trkTree->Fill();   
      newtrkTree->Fill();
      trkheTree->Fill();
      calTree->Fill();
      CSTree->Fill(); 
      jetpt.clear();
      jeteta.clear();
      jetphi.clear();
      jetmass.clear(); 
   } //end of event loop
   // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
   //   cout << numevt << endl;
   CSTree->Write();
   calTree->Write();
   trkheTree->Write();
   newtrkTree->Write();
   trkTree->Write();
   djTree->Write();
   eventTree->Write();
   
   newfile->Close();
} 



TString getDayAndTime() { 
  time_t currentTime;
  struct tm *localTime;
  
  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time
  
  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  //  int Sec    = localTime->tm_sec;
  return Form("%dy_%dm_%dd_%dh_%dm",Year,Month,Day,Hour,Min);
}
