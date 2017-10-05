#include <TSystem.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>
#include <TMath.h>
#include <TF1.h>
#include <TArrow.h>
#include <TColor.h>
#include "cutsAndBin_bgk.h"
#include "commonUtility.h"
#include "plotStyle/SONGKYO.h"
#include "plotStyle/tdrstyle.C"
#include "plotStyle/CMS_lumi_raaCent.C"

using namespace std;
void UPC_trigger_plots_dijet_HF()
{ 
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 3; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  TH1::SetDefaultSumw2();
  
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);

  TTree *fChain;
  TTree *fChain2;
  Int_t fCurrent;
  TFile *data = TFile::Open("/u/user/bekim/758p3_UPC/src/Muon_test/skimmedFiles/TrkCutsppreco_norechitcut_upcDiJetSkim170919_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak4PFJetAnalyzer_minJetPt0_2017y_9m_19d_13h_28m.root");

//  TFile *f2 = TFile::Open("h_eff_v2.root");
//  TH1F *heff = (TH1F*)f2->Get("h4");

  TTree *trkTree = (TTree*) data -> Get("Track");
  TTree *ftrkTree = (TTree*) data -> Get("fullTrkTree");
  TTree *jetTree = (TTree*) data -> Get("dijet");
  TTree *evtTree = (TTree*) data -> Get("evt");
  TTree *csTree = (TTree*) data -> Get("CS");
  TTree *calTree = (TTree*) data -> Get("Cal");

/*
  TTree *trkTree2 = (TTree*) data2 -> Get("Track");
  TTree *ftrkTree2 = (TTree*) data2 -> Get("fullTrkTree");
  TTree *jetTree2 = (TTree*) data2 -> Get("dijet");
  TTree *evtTree2 = (TTree*) data2 -> Get("evt");
  TTree *csTree2 = (TTree*) data2 -> Get("CS");
  TTree *calTree2 = (TTree*) data2 -> Get("Cal");
  TTree *trkheTree2 = (TTree*) data2 -> Get("trkandHE");
*/

  TCanvas* c1 = new TCanvas("c1","",600,600);
  TCanvas* c2 = new TCanvas("c2","",600,600);
  TCanvas* c3 = new TCanvas("c3","",600,600);
  TCanvas* c4 = new TCanvas("c4","",600,600);
  TCanvas* c5 = new TCanvas("c5","",600,600);
  TCanvas* c6 = new TCanvas("c6","",600,600);
  TCanvas* c7 = new TCanvas("c7","",600,600);
  TCanvas* c8 = new TCanvas("c8","",600,600);
  TCanvas* c9 = new TCanvas("c9","",600,600);
  TCanvas* c10 = new TCanvas("c10","",600,600);
  TCanvas* c11 = new TCanvas("c11","",600,600);
  TCanvas* c12 = new TCanvas("c12","",600,600);
  TCanvas* c13 = new TCanvas("c13","",600,600);
  TCanvas* c14 = new TCanvas("c14","",600,600);
  TCanvas* c15 = new TCanvas("c15","",600,600);
  TCanvas* c16 = new TCanvas("c16","",600,600);
  TCanvas* c17 = new TCanvas("c17","",600,600);
  TCanvas* c18 = new TCanvas("c18","",600,600);
  TCanvas* cv1 = new TCanvas("cv1","",600,600);
  TCanvas* cv2 = new TCanvas("cv2","",600,600);
  TCanvas* cv3 = new TCanvas("cv3","",600,600);
  TCanvas* cv4 = new TCanvas("cv4","",600,600);
  TCanvas* cv5 = new TCanvas("cv5","",600,600);
  TCanvas* cv6 = new TCanvas("cv6","",600,600);
  TCanvas* cv7 = new TCanvas("cv7","",600,600);
  TCanvas* cv8 = new TCanvas("cv8","",600,600);
  TCanvas* cv9 = new TCanvas("cv9","",600,600);
  TCanvas* cv10 = new TCanvas("cv10","",600,600);
  TCanvas* cv11 = new TCanvas("cv11","",600,600);
  TCanvas* cv12 = new TCanvas("cv12","",600,600);
/*  TCanvas* cv13 = new TCanvas("cv13","",600,600);
  TCanvas* cv14 = new TCanvas("cv14","",600,600);
  TCanvas* cv15 = new TCanvas("cv15","",600,600);
  TCanvas* cv16 = new TCanvas("cv16","",600,600);
  TCanvas* cv17 = new TCanvas("cv17","",600,600);
  TCanvas* cv18 = new TCanvas("cv18","",600,600);*/

  TH1F *h1 = new TH1F("h1", " ", 20, 0, 1);
  TH1F *h2 = new TH1F("h2", " ", 20, 0, 1);
  TH1F *h3 = new TH1F("h3", " ", 20, 0, 1);
  TH1F *h4 = new TH1F("h4", " ", 20, 0, 1);
  TH1F *h5 = new TH1F("h5", " ", 20, 0, 0.07);
  TH1F *h6 = new TH1F("h6", " ", 20, 0, 0.07);
  TH1F *h7 = new TH1F("h7", " ", 20, 0, 0.07);
  TH1F *h8 = new TH1F("h8", " ", 20, 0, 0.07);
  TH1F *h9 = new TH1F("h9", " ", 10, 0, 1);
  TH1F *h10 = new TH1F("h10", " ", 10, 0, 1);
  TH1F *h11 = new TH1F("h11", " ", 10, 0, 1);
  TH1F *h12 = new TH1F("h12", " ", 10, 0, 1);
  TH1F *h13 = new TH1F("h13", " ", 20, 0, 3.15);
  TH1F *h14 = new TH1F("h14", " ", 20, 0, 3.15);
  TH1F *h15 = new TH1F("h15", " ", 20, 0, 3.15);
  TH1F *h16 = new TH1F("h16", " ", 20, 0, 3.15);
  TH1F *h17 = new TH1F("h17", " ", 20, 0, 50);
  TH1F *h18 = new TH1F("h18", " ", 20, 0, 50);
  TH1F *h19 = new TH1F("h19", " ", 20, 0, 50);
  TH1F *h20 = new TH1F("h20", " ", 20, 0, 50);
  TH1F *h21 = new TH1F("h21", " ", 20, 0, 20);
  TH1F *h22 = new TH1F("h22", " ", 20, 0, 20);
  TH1F *h23 = new TH1F("h23", " ", 20, 0, 20);
  TH1F *h24 = new TH1F("h24", " ", 20, 0, 20);
  TH1F *h25 = new TH1F("h25", " ", 20, 0, 20);
  TH1F *h26 = new TH1F("h26", " ", 20, 0, 20);
  TH1F *h27 = new TH1F("h27", " ", 20, 0, 20);
  TH1F *h28 = new TH1F("h28", " ", 20, 0, 20);
  TH1F *h29 = new TH1F("h29", " ", 20, 0, 20);
  TH1F *h30 = new TH1F("h30", " ", 20, 0, 20);
  TH1F *h31 = new TH1F("h31", " ", 20, 0, 20);
  TH1F *h32 = new TH1F("h32", " ", 20, 0, 20);
  TH1F *h33 = new TH1F("h33", " ", 20, 0, 20);
  TH1F *h34 = new TH1F("h34", " ", 20, 0, 20);
  TH1F *h35 = new TH1F("h35", " ", 20, 0, 20);
  TH1F *h36 = new TH1F("h36", " ", 20, 0, 100);
  TH1F *h37 = new TH1F("h37", " ", 20, 0, 10);

  TH2F *th1 = new TH2F("th1", " ", 25, 0, 25, 25, 0, 50);
  TH2F *th2 = new TH2F("th2", " ", 25, 0, 3.15, 25, 0, 5);
  TH2F *th3 = new TH2F("th3", " ", 25, 0, 0.07, 25, 0, 50);
  TH2F *th4 = new TH2F("th4", " ", 25, 0, 0.07, 25, 0, 50);
  TH2F *th5 = new TH2F("th5", " ", 25, 0, 100, 25, 0, 7);

struct TRACK
{
  Int_t           trkVar_nTrk;
  Int_t           trkVar_nTrkabsEtaover1p5;
  Int_t           trkVar_nTrkabsEtaunder1p5;
  Int_t           trkVar_nTrketa0to0p5;
  Int_t           trkVar_nTrketa0p5to1;
  Int_t           trkVar_nTrketa1to1p5;
  Int_t           trkVar_nTrketa1p5to2;
  Int_t           trkVar_nTrketa2to2p5;
  Int_t           trkVar_nTrketam0to0p5;
  Int_t           trkVar_nTrketam0p5to1;
  Int_t           trkVar_nTrketam1to1p5;
  Int_t           trkVar_nTrketam1p5to2;
  Int_t           trkVar_nTrketam2to2p5;
};

  TRACK           Trk;

  trkTree->SetBranchAddress("trkVar", &Trk);

  Int_t           HBHEn;
  Float_t         HBHEe;
  Float_t         HBHEeta;
  Float_t         HBHEphi;
  Float_t         HBHEperp;
  Float_t         HBHEisjet;
  Float_t         HBHEdepth;
  Float_t         HBmax;
  Float_t         HEmax;
  Int_t           HFn;
  Float_t         HFe;
  Float_t         HFeta;
  Float_t         HFphi;
  Float_t         HFperp;
  Float_t         HFisjet;
  Float_t         HFdepth;
  Float_t         HFtotal;
  Float_t         HFplus;
  Float_t         HFminus;
  Float_t         HFmax;
  Int_t           EEn;
  Float_t         EEe;
  Float_t         EEeta;
  Float_t         EEphi;
  Float_t         EEperp;
  Float_t         EEisjet;
  Float_t         EEmax;
  Int_t           EBn;
  Float_t         EBe;
  Float_t         EBeta;
  Float_t         EBphi;
  Float_t         EBperp;
  Float_t         EBisjet;
  Float_t         EBmax;

  TBranch        *b_HBHEn;   //!
  TBranch        *b_HBHEe;   //!
  TBranch        *b_HBHEeta;   //!
  TBranch        *b_HBHEphi;   //!
  TBranch        *b_HBHEperp;   //!
  TBranch        *b_HBHEisjet;   //!
  TBranch        *b_HBHEdepth;   //!
  TBranch        *b_HBmax;   //!
  TBranch        *b_HEmax;   //!
  TBranch        *b_HFn;   //!
  TBranch        *b_HFe;   //!
  TBranch        *b_HFeta;   //!
  TBranch        *b_HFphi;   //!
  TBranch        *b_HFperp;   //!
  TBranch        *b_HFisjet;   //!
  TBranch        *b_HFdepth;   //!
  TBranch        *b_HFtotal;   //!
  TBranch        *b_HFplus;   //!
  TBranch        *b_HFminus;   //!
  TBranch        *b_HFmax;   //!
  TBranch        *b_EEn;   //!
  TBranch        *b_EEe;   //!
  TBranch        *b_EEeta;   //!
  TBranch        *b_EEphi;   //!
  TBranch        *b_EEperp;   //!
  TBranch        *b_EEisjet;   //!
  TBranch        *b_EEmax;   //!
  TBranch        *b_EBn;   //!
  TBranch        *b_EBe;   //!
  TBranch        *b_EBeta;   //!
  TBranch        *b_EBphi;   //!
  TBranch        *b_EBperp;   //!
  TBranch        *b_EBisjet;   //!
  TBranch        *b_EBmax;   //!

  calTree->SetBranchAddress("HBHEn", &HBHEn, &b_HBHEn);
  calTree->SetBranchAddress("HBHEe", &HBHEe, &b_HBHEe);
  calTree->SetBranchAddress("HBHEeta", &HBHEeta, &b_HBHEeta);
  calTree->SetBranchAddress("HBHEphi", &HBHEphi, &b_HBHEphi);
  calTree->SetBranchAddress("HBHEperp", &HBHEperp, &b_HBHEperp);
  calTree->SetBranchAddress("HBHEisjet", &HBHEisjet, &b_HBHEisjet);
  calTree->SetBranchAddress("HBHEdepth", &HBHEdepth, &b_HBHEdepth);
  calTree->SetBranchAddress("HBmax", &HBmax, &b_HBmax);
  calTree->SetBranchAddress("HEmax", &HEmax, &b_HEmax);
  calTree->SetBranchAddress("HFn", &HFn, &b_HFn);
  calTree->SetBranchAddress("HFe", &HFe, &b_HFe);
  calTree->SetBranchAddress("HFeta", &HFeta, &b_HFeta);
  calTree->SetBranchAddress("HFphi", &HFphi, &b_HFphi);
  calTree->SetBranchAddress("HFperp", &HFperp, &b_HFperp);
  calTree->SetBranchAddress("HFisjet", &HFisjet, &b_HFisjet);
  calTree->SetBranchAddress("HFdepth", &HFdepth, &b_HFdepth);
  calTree->SetBranchAddress("HFtotal", &HFtotal, &b_HFtotal);
  calTree->SetBranchAddress("HFplus", &HFplus, &b_HFplus);
  calTree->SetBranchAddress("HFminus", &HFminus, &b_HFminus);
  calTree->SetBranchAddress("HFmax", &HFmax, &b_HFmax);
  calTree->SetBranchAddress("EEn", &EEn, &b_EEn);
  calTree->SetBranchAddress("EEe", &EEe, &b_EEe);
  calTree->SetBranchAddress("EEeta", &EEeta, &b_EEeta);
  calTree->SetBranchAddress("EEphi", &EEphi, &b_EEphi);
  calTree->SetBranchAddress("EEperp", &EEperp, &b_EEperp);
  calTree->SetBranchAddress("EEisjet", &EEisjet, &b_EEisjet);
  calTree->SetBranchAddress("EEmax", &EEmax, &b_EEmax);
  calTree->SetBranchAddress("EBn", &EBn, &b_EBn);
  calTree->SetBranchAddress("EBe", &EBe, &b_EBe);
  calTree->SetBranchAddress("EBeta", &EBeta, &b_EBeta);
  calTree->SetBranchAddress("EBphi", &EBphi, &b_EBphi);
  calTree->SetBranchAddress("EBperp", &EBperp, &b_EBperp);
  calTree->SetBranchAddress("EBisjet", &EBisjet, &b_EBisjet);
  calTree->SetBranchAddress("EBmax", &EBmax, &b_EBmax);

struct DIJET
{
  Int_t           dj_nJet;
  Float_t         dj_mass;
  Float_t         dj_pt;
  Float_t         dj_y;
  Float_t         dj_phi;
  Float_t         dj_eta;
  Float_t         dj_dphi;
  Float_t         dj_dpt;
  Float_t         dj_deta;
  Float_t         dj_aj;
  Float_t         dj_pt1;
  Float_t         dj_eta1;
  Float_t         dj_phi1;
  Float_t         dj_e1;
  Float_t         dj_pt2;
  Float_t         dj_eta2;
  Float_t         dj_phi2;
  Float_t         dj_e2;
};

  DIJET           djt;

  jetTree->SetBranchAddress("dj", &djt);

  Int_t           bin1;
  Int_t           bin2;
  Float_t         pT1;
  Float_t         pT2;
  Float_t         eff1;
  Float_t         eff2;
  Float_t         w1;
  Float_t         w2;
  Float_t         DjpT;
 
  TBranch        *b_bin1;   //!
  TBranch        *b_bin2;   //!
  TBranch        *b_pT1;   //!
  TBranch        *b_pT2;   //!
  TBranch        *b_eff1;   //!
  TBranch        *b_eff2;   //!
  TBranch        *b_w1;   //!
  TBranch        *b_w2;   //!
  TBranch        *b_DjpT;   //!

  csTree->SetBranchAddress("bin1", &bin1, &b_bin1);
  csTree->SetBranchAddress("bin2", &bin2, &b_bin2);
  csTree->SetBranchAddress("pT1", &pT1, &b_pT1);
  csTree->SetBranchAddress("pT2", &pT2, &b_pT2);
  csTree->SetBranchAddress("eff1", &eff1, &b_eff1);
  csTree->SetBranchAddress("eff2", &eff2, &b_eff2);
  csTree->SetBranchAddress("w1", &w1, &b_w1);
  csTree->SetBranchAddress("w2", &w2, &b_w2);
  csTree->SetBranchAddress("DjpT", &DjpT, &b_DjpT);

struct EVENT
{
  Int_t           event_run;
  Int_t           event_lumi;
  Int_t           event_event;
  Int_t           event_nPho;
  Int_t           event_nTrk;
  Float_t         event_hfsum;
  Float_t         event_hfplus;
  Float_t         event_hfminus;
  Float_t         event_vz;
};

  EVENT           event;

  evtTree->SetBranchAddress("event", &event);

  Int_t           ntrk;
  Float_t         floatntrk;
  Float_t         pT[52];   //[ntrk]
  Float_t         Eta[52];   //[ntrk]
  Float_t         Phi[52];   //[ntrk]
  Float_t         trkHFplus;
  Float_t         trkHFminus;

  TBranch        *b_ntrk;   //!
  TBranch        *b_floatntrk;   //!
  TBranch        *b_pT;   //!
  TBranch        *b_Eta;   //!
  TBranch        *b_Phi;   //!
  TBranch        *b_trkHFplus;   //!
  TBranch        *b_trkHFminus;

  ftrkTree->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
  ftrkTree->SetBranchAddress("floatntrk", &floatntrk, &b_floatntrk);
  ftrkTree->SetBranchAddress("pT", pT, &b_pT);
  ftrkTree->SetBranchAddress("Eta", Eta, &b_Eta);
  ftrkTree->SetBranchAddress("Phi", Phi, &b_Phi);
  ftrkTree->SetBranchAddress("trkHFplus", &trkHFplus, &b_trkHFplus);
  ftrkTree->SetBranchAddress("trkHFminus", &trkHFminus, &b_trkHFminus);

  jetTree->AddFriend(evtTree);
  jetTree->AddFriend(trkTree);
  trkTree->AddFriend(evtTree);
  trkTree->AddFriend(jetTree);
  evtTree->AddFriend(jetTree);
  evtTree->AddFriend(trkTree);
  csTree->AddFriend(jetTree);
  csTree->AddFriend(trkTree);
  csTree->AddFriend(evtTree);
  calTree->AddFriend(jetTree);
  calTree->AddFriend(trkTree);
  calTree->AddFriend(evtTree);

  jetTree->AddFriend(ftrkTree);
  evtTree->AddFriend(ftrkTree);
  trkTree->AddFriend(ftrkTree);

  Int_t eta1pluseta2plus = 0;
  Int_t eta1minuseta2minus = 0;
  Int_t eta1pluseta2minus = 0;
  Int_t eta1minuseta2plus = 0;

  TH2F *pe[500];
  TH2F *pe2[500];
  for(Int_t a = 0; a != 500; a++)
  {
    pe[a] = new TH2F(Form("phiVSeta_jets_%d",a), "pe", 50, -2.5, 2.5, 50, -3.15, 3.15);
    pe2[a] = new TH2F(Form("phiVSeta_Tracks_%d",a), "pe2", 50, -2.5, 2.5, 50, -3.15, 3.15);
  }

  Int_t tot = 0;
  tot = jetTree->GetEntries();
  //tot = 10000;
  cout << "total events : " << tot << endl;
  Int_t te = 0;
  for(Int_t i = 0; i != tot; i++)
  {
    if(i%10000==0)
    {
      cout << ">>>>> EVENT " << i << " / " << jetTree->GetEntries() <<  " ("<<(int)(100.*i/jetTree->GetEntries()) << "%)" << endl;
    }
    /*dj_dphi = 0;
    dj_eta = 0;
    dj_eta1 = 0;
    dj_eta2 = 0;
    dj_pt1 = 0;
    dj_pt2 = 0;
    trkVar_nTrkabsEtaover1p5 = 0;
    event_hfplus = 0;
    event_hfminus = 0;*/
  
    jetTree->GetEntry(i);

    //cout << "jet phi1 : " << djt.dj_phi1 << "  hfplus : " << event.event_hfplus << "  hfminus : " << event.event_hfminus << "  pt1 : " << djt.dj_pt1 << "  pt2 : " << djt.dj_pt2 << "  eta1 : " << djt.dj_eta1 << "  eta2 : " << djt.dj_eta2 << "  nTrk over 1.5 :" << Trk.trkVar_nTrkabsEtaover1p5 << endl;
//    if(djt.dj_dphi > 2 && event.event_hfplus < 5 && event.event_hfminus < 5 && djt.dj_pt1 > 20 && djt.dj_pt2 > 20 && fabs(djt.dj_eta1) < 1.5 && fabs(djt.dj_eta2) < 1.5 && Trk.trkVar_nTrkabsEtaover1p5 == 0)
    if(djt.dj_pt1 > 20 && djt.dj_pt2 > 20 && fabs(djt.dj_eta1) < 1.5 && fabs(djt.dj_eta2) < 1.5)
    {     
      /*pe[te]->Fill(djt.dj_eta1, djt.dj_phi1,4);
      pe[te]->Fill(djt.dj_eta2, djt.dj_phi2,4);
      //cout << "phi1 : " << djt.dj_phi1 << " eta1 : " << djt.dj_eta1 << " phi2 : " << djt.dj_phi2 << " eta2 : " << djt.dj_eta2 << endl;
      for(Int_t b = 0; b != ntrk; b++)
      {
        pe2[te]->Fill(Eta[b], Phi[b]);
        //cout << "trkPhi[" << b << "] : " << Phi[b] << " trkEta[" << b << "] : " << Eta[b] << endl;
      }

      c1->Clear();
      c1->cd();
      gPad->SetBottomMargin(0.14);
      gPad->SetTopMargin(0.067);
      //Double_t sc1 = 1/h1->Integral();
      //h1->Scale(sc1);
      pe[te]->SetMarkerSize(1.2);
      pe[te]->GetXaxis()->SetTitle("#eta");      
      //pe[te]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      //pe[te]->GetXaxis()->SetTitle("X_{#gamma}");
      //pe[te]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
      //pe[te]->GetXaxis()->SetTitle("A_{j}");
      //pe[te]->GetXaxis()->SetTitle("#Delta#phi");
      //pe[te]->GetXaxis()->SetTitle("p_{T,2}/p_{T,1}");
      //pe[te]->GetXaxis()->SetTitle("nTracks");
      //pe[te]->GetXaxis()->SetTitleOffset(1.3);
      pe[te]->GetXaxis()->CenterTitle();
      pe[te]->GetYaxis()->SetTitle("#phi");
      //pe[te]->GetYaxis()->SetTitle("Cross section (nb/GeV)");
      //pe[te]->GetYaxis()->SetTitle("Entries/N");
      //pe[te]->GetYaxis()->SetTitle("Entries");
      //pe[te]->GetYaxis()->SetTitleOffset(1.3);
      pe[te]->GetYaxis()->CenterTitle();
      pe[te]->SetFillColor(2);
      pe[te]->Draw("BOX");
      pe2[te]->SetFillColor(4);
      pe2[te]->Draw("BOX same");

      TLegend *leg = new TLegend(0.8,0.70,0.9,0.85);
      //  leg->SetFillStyle(0);
      leg->SetTextSize(0.020);
      SetLegendStyle(leg);
      leg->AddEntry(pe[0],"Jets","F");
      leg->AddEntry(pe2[0],"Tracks","F");
      TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.020);
      leg->Draw("same");*/

      //c1->SaveAs(Form("/u/user/bekim/758p3_UPC/src/Muon_test/phietaplots/phiVSeta_%d.png",te));
 
      if(djt.dj_eta1 > 0 && djt.dj_eta2 > 0)
      {
        h1->Fill(djt.dj_aj);
        h5->Fill((djt.dj_mass/5020)*TMath::Exp(-djt.dj_y));
        h9->Fill(djt.dj_pt2/djt.dj_pt1);
        h13->Fill(djt.dj_dphi);
        h17->Fill(djt.dj_pt);
        th3->Fill((djt.dj_mass/5020)*TMath::Exp(-djt.dj_y), djt.dj_pt);
        h24->Fill(event.event_hfplus);
        h28->Fill(event.event_hfminus);
        h32->Fill(event.event_hfsum);
        eta1pluseta2plus += 1;
      }
      if(djt.dj_eta1 < 0 && djt.dj_eta2 < 0)
      {
        h2->Fill(djt.dj_aj);
        h6->Fill((djt.dj_mass/5020)*TMath::Exp(-djt.dj_y));
        h10->Fill(djt.dj_pt2/djt.dj_pt1);
        h14->Fill(djt.dj_dphi);
        h18->Fill(djt.dj_pt);
        th4->Fill((djt.dj_mass/5020)*TMath::Exp(-djt.dj_y), djt.dj_pt);
        h25->Fill(event.event_hfplus);
        h29->Fill(event.event_hfminus);
        h33->Fill(event.event_hfsum);
        eta1minuseta2minus += 1;
      }
      if(djt.dj_eta1 > 0 && djt.dj_eta2 < 0)
      {
        h3->Fill(djt.dj_aj);
        h7->Fill((djt.dj_mass/5020)*TMath::Exp(-djt.dj_y));
        h11->Fill(djt.dj_pt2/djt.dj_pt1);
        h15->Fill(djt.dj_dphi);
        h19->Fill(djt.dj_pt);
        h26->Fill(event.event_hfplus);
        h30->Fill(event.event_hfminus);
        h34->Fill(event.event_hfsum);
        eta1pluseta2minus += 1;
      }
      if(djt.dj_eta1 < 0 && djt.dj_eta2 > 0)
      {
        h4->Fill(djt.dj_aj);
        h8->Fill((djt.dj_mass/5020)*TMath::Exp(-djt.dj_y));
        h12->Fill(djt.dj_pt2/djt.dj_pt1);
        h16->Fill(djt.dj_dphi);
        h20->Fill(djt.dj_pt);
        h27->Fill(event.event_hfplus);
        h31->Fill(event.event_hfminus);
        h35->Fill(event.event_hfsum);
        eta1minuseta2plus += 1;
      }
      h21->Fill(event.event_hfplus);
      h22->Fill(event.event_hfminus);
      h23->Fill(event.event_hfsum);
      h36->Fill((djt.dj_pt)*(djt.dj_pt));
      h37->Fill(0.5*djt.dj_dpt);
      th5->Fill((djt.dj_pt)*(djt.dj_pt), 0.5*djt.dj_dpt);
      te += 1;
    }
    else
    {
      continue;
    }
  }
  cout << "total event : " << te << endl;
  cout << "eta1 > 0 & eta2 > 0 : " << eta1pluseta2plus << endl;
  cout << "eta1 < 0 & eta2 < 0 : " << eta1minuseta2minus << endl;
  cout << "eta1 > 0 & eta2 < 0 : " << eta1pluseta2minus << endl;
  cout << "eta1 < 0 & eta2 > 0 : " << eta1minuseta2plus << endl;

  c2->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc1 = 1/h1->Integral();
  h1->Scale(sc1);
  h1->SetMarkerSize(1.2);
  //h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //  h1->GetXaxis()->SetTitle("X_{#gamma}");
  //  h1->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  h1->GetXaxis()->SetTitle("A_{j}");
  //  h1->GetXaxis()->SetTitle("#Delta#phi");
  //  h1->GetXaxis()->SetTitle("p_{T,2}/p_{T,1}");
  //  h1->GetXaxis()->SetTitle("nTracks");
  //  h1->GetXaxis()->SetTitleOffset(1.3);
  h1->GetXaxis()->CenterTitle();
  //h1->GetYaxis()->SetTitle("Cross section (nb/GeV)");
  h1->GetYaxis()->SetTitle("Entries/N");
  //h1->GetYaxis()->SetTitle("Entries");
  //  h3->GetYaxis()->SetTitleOffset(1.3);
  h1->GetYaxis()->CenterTitle();
  h1->SetLineColor(1);
  h1->SetMarkerColor(1);
  h1->SetMarkerStyle(kOpenCircle);
  //TH1ScaleByWidth(h1);
  h1->Sumw2();
  //h1->Scale(1/401.); //Only for Mass, X_gamma, pT
  h1->SetMaximum(0.68);
  h1->SetMinimum(0);
  h1->Draw("PE");
  Double_t sc2 = 1/h2->Integral();
  h2->Scale(sc2);
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h2->SetMarkerStyle(kOpenDiamond);
  h2->Draw("PE same");
  Double_t sc3 = 1/h3->Integral();
  h3->Scale(sc3);
  h3->SetLineColor(4);
  h3->SetMarkerColor(4);
  h3->SetMarkerStyle(kOpenTriangleUp);
  h3->Draw("PE same");
  Double_t sc4 = 1/h4->Integral();
  h4->Scale(sc4);
  h4->SetLineColor(6);
  h4->SetMarkerColor(6);
  h4->SetMarkerStyle(kOpenTriangleDown);
  h4->Draw("PE same"); 
  
  TLegend *leg = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg);
  leg->AddEntry(h1,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg->AddEntry(h2,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg->AddEntry(h3,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg->AddEntry(h4,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  //TLatex *lat = new TLatex(); lat->SetNDC(); lat->SetTextSize(0.025);
  leg->Draw("same");

  c3->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc5 = 1/h5->Integral();
  h5->Scale(sc5);
  h5->SetMarkerSize(1.2);
  h5->GetXaxis()->SetTitle("X_{#gamma}");
  h5->GetXaxis()->CenterTitle();
  h5->GetYaxis()->SetTitle("Entries/N");
  h5->GetYaxis()->CenterTitle();
  h5->SetLineColor(1);
  h5->SetMarkerColor(1);
  h5->SetMarkerStyle(kOpenCircle);
  h5->Sumw2();
  h5->SetMaximum(1);
  h5->SetMinimum(0);
  h5->Draw("PE");
  Double_t sc6 = 1/h6->Integral();
  h6->Scale(sc6);
  h6->SetLineColor(2);
  h6->SetMarkerColor(2);
  h6->SetMarkerStyle(kOpenDiamond);
  h6->Draw("PE same");
  Double_t sc7 = 1/h7->Integral();
  h7->Scale(sc7);
  h7->SetLineColor(4);
  h7->SetMarkerColor(4);
  h7->SetMarkerStyle(kOpenTriangleUp);
  h7->Draw("PE same");
  Double_t sc8 = 1/h8->Integral();
  h8->Scale(sc8);
  h8->SetLineColor(6);
  h8->SetMarkerColor(6);
  h8->SetMarkerStyle(kOpenTriangleDown);
  h8->Draw("PE same");

  TLegend *leg2 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg2);
  leg2->AddEntry(h5,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg2->AddEntry(h6,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg2->AddEntry(h7,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg2->AddEntry(h8,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg2->Draw("same");

  c4->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc9 = 1/h9->Integral();
  h9->Scale(sc9);
  h9->SetMarkerSize(1.2);
  h9->GetXaxis()->SetTitle("p_{T,2}/p_{T,1}");
  h9->GetXaxis()->CenterTitle();
  h9->GetYaxis()->SetTitle("Entries/N");
  h9->GetYaxis()->CenterTitle();
  h9->SetLineColor(1);
  h9->SetMarkerColor(1);
  h9->SetMarkerStyle(kOpenCircle);
  h9->Sumw2();
  h9->SetMaximum(0.7);
  h9->SetMinimum(0);
  h9->Draw("PE");
  Double_t sc10 = 1/h10->Integral();
  h10->Scale(sc10);
  h10->SetLineColor(2);
  h10->SetMarkerColor(2);
  h10->SetMarkerStyle(kOpenDiamond);
  h10->Draw("PE same");
  Double_t sc11 = 1/h11->Integral();
  h11->Scale(sc11);
  h11->SetLineColor(4);
  h11->SetMarkerColor(4);
  h11->SetMarkerStyle(kOpenTriangleUp);
  h11->Draw("PE same");
  Double_t sc12 = 1/h12->Integral();
  h12->Scale(sc12);
  h12->SetLineColor(6);
  h12->SetMarkerColor(6);
  h12->SetMarkerStyle(kOpenTriangleDown);
  h12->Draw("PE same");

  TLegend *leg3 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg3);
  leg3->AddEntry(h9,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg3->AddEntry(h10,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg3->AddEntry(h11,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg3->AddEntry(h12,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg3->Draw("same");

  c5->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc13 = 1/h13->Integral();
  h13->Scale(sc13);
  h13->SetMarkerSize(1.2);
  h13->GetXaxis()->SetTitle("#Delta#phi");
  h13->GetXaxis()->CenterTitle();
  h13->GetYaxis()->SetTitle("Entries/N");
  h13->GetYaxis()->CenterTitle();
  h13->SetLineColor(1);
  h13->SetMarkerColor(1);
  h13->SetMarkerStyle(kOpenCircle);
  h13->Sumw2();
  h13->SetMaximum(0.92);
  h13->SetMinimum(0);
  h13->Draw("PE");
  Double_t sc14 = 1/h14->Integral();
  h14->Scale(sc14);
  h14->SetLineColor(2);
  h14->SetMarkerColor(2);
  h14->SetMarkerStyle(kOpenDiamond);
  h14->Draw("PE same");
  Double_t sc15 = 1/h15->Integral();
  h15->Scale(sc15);
  h15->SetLineColor(4);
  h15->SetMarkerColor(4);
  h15->SetMarkerStyle(kOpenTriangleUp);
  h15->Draw("PE same");
  Double_t sc16 = 1/h16->Integral();
  h16->Scale(sc16);
  h16->SetLineColor(6);
  h16->SetMarkerColor(6);
  h16->SetMarkerStyle(kOpenTriangleDown);
  h16->Draw("PE same");

  TLegend *leg4 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg4);
  leg4->AddEntry(h13,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg4->AddEntry(h14,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg4->AddEntry(h15,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg4->AddEntry(h16,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg4->Draw("same");

  c6->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc17 = 1/h17->Integral();
  h17->Scale(sc17);
  h17->SetMarkerSize(1.2);
  h17->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h17->GetXaxis()->CenterTitle();
  h17->GetYaxis()->SetTitle("Entries/N");
  h17->GetYaxis()->CenterTitle();
  h17->SetLineColor(1);
  h17->SetMarkerColor(1);
  h17->SetMarkerStyle(kOpenCircle);
  h17->Sumw2();
  h17->SetMaximum(1);
  //h17->SetMinimum(0);
  h17->Draw("PE");
  Double_t sc18 = 1/h18->Integral();
  h18->Scale(sc18);
  h18->SetLineColor(2);
  h18->SetMarkerColor(2);
  h18->SetMarkerStyle(kOpenDiamond);
  h18->Draw("PE same");
  Double_t sc19 = 1/h19->Integral();
  h19->Scale(sc19);
  h19->SetLineColor(4);
  h19->SetMarkerColor(4);
  h19->SetMarkerStyle(kOpenTriangleUp);
  h19->Draw("PE same");
  Double_t sc20 = 1/h20->Integral();
  h20->Scale(sc20);
  h20->SetLineColor(6);
  h20->SetMarkerColor(6);
  h20->SetMarkerStyle(kOpenTriangleDown);
  h20->Draw("PE same");
  c6->SetLogy();

  TLegend *leg5 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg5);
  leg5->AddEntry(h17,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg5->AddEntry(h18,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg5->AddEntry(h19,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg5->AddEntry(h20,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg5->Draw("same");

  c7->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  th3->GetXaxis()->SetTitle("X_{#gamma}");
  th3->GetXaxis()->CenterTitle();
  th3->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  th3->GetYaxis()->CenterTitle();
  th3->Draw("COLZ");

  c8->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  th4->GetXaxis()->SetTitle("X_{#gamma}");
  th4->GetXaxis()->CenterTitle();
  th4->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  th4->GetYaxis()->CenterTitle();
  th4->Draw("COLZ");

  c9->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc21 = 1/h21->Integral();
  h21->Scale(sc21);
  h21->SetMarkerSize(1.2);
  h21->GetXaxis()->SetTitle("HF+ (GeV)");
  h21->GetXaxis()->CenterTitle();
  h21->GetYaxis()->SetTitle("Entries/N");
  h21->GetYaxis()->CenterTitle();
  h21->SetLineColor(1);
  h21->SetMarkerColor(1);
  h21->SetMarkerStyle(kOpenCircle);
  h21->Sumw2();
  h21->SetMaximum(0.4);
  //h21->SetMinimum(0);
  h21->Draw("PE");

  c10->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc22 = 1/h22->Integral();
  h22->Scale(sc22);
  h22->SetMarkerSize(1.2);
  h22->GetXaxis()->SetTitle("HF- (GeV)");
  h22->GetXaxis()->CenterTitle();
  h22->GetYaxis()->SetTitle("Entries/N");
  h22->GetYaxis()->CenterTitle();
  h22->SetLineColor(1);
  h22->SetMarkerColor(1);
  h22->SetMarkerStyle(kOpenCircle);
  h22->Sumw2();
  h22->SetMaximum(0.4);
  //h22->SetMinimum(0);
  h22->Draw("PE"); 

  c11->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc23 = 1/h23->Integral();
  h23->Scale(sc23);
  h23->SetMarkerSize(1.2);
  h23->GetXaxis()->SetTitle("HFsum (GeV)");
  h23->GetXaxis()->CenterTitle();
  h23->GetYaxis()->SetTitle("Entries/N");
  h23->GetYaxis()->CenterTitle();
  h23->SetLineColor(1);
  h23->SetMarkerColor(1);
  h23->SetMarkerStyle(kOpenCircle);
  h23->Sumw2();
  h23->SetMaximum(0.4);
  //h23->SetMinimum(0);
  h23->Draw("PE");

  c12->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc24 = 1/h24->Integral();
  h24->Scale(sc24);
  h24->SetMarkerSize(1.2);
  h24->GetXaxis()->SetTitle("HF+ (GeV)");
  h24->GetXaxis()->CenterTitle();
  h24->GetYaxis()->SetTitle("Entries/N");
  h24->GetYaxis()->CenterTitle();
  h24->SetLineColor(1);
  h24->SetMarkerColor(1);
  h24->SetMarkerStyle(kOpenCircle);
  h24->Sumw2();
  h24->SetMaximum(0.4);
  h24->Draw("PE");
  Double_t sc25 = 1/h25->Integral();
  h25->Scale(sc25);
  h25->SetLineColor(2);
  h25->SetMarkerColor(2);
  h25->SetMarkerStyle(kOpenDiamond);
  h25->Draw("PE same");
  Double_t sc26 = 1/h26->Integral();
  h26->Scale(sc26);
  h26->SetLineColor(4);
  h26->SetMarkerColor(4);
  h26->SetMarkerStyle(kOpenTriangleUp);
  h26->Draw("PE same");
  Double_t sc27 = 1/h27->Integral();
  h27->Scale(sc27);
  h27->SetLineColor(6);
  h27->SetMarkerColor(6);
  h27->SetMarkerStyle(kOpenTriangleDown);
  h27->Draw("PE same");  

  TLegend *leg6 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg6);
  leg6->AddEntry(h24,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg6->AddEntry(h25,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg6->AddEntry(h26,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg6->AddEntry(h27,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg6->Draw("same");

  c13->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc28 = 1/h28->Integral();
  h28->Scale(sc28);
  h28->SetMarkerSize(1.2);
  h28->GetXaxis()->SetTitle("HF- (GeV)");
  h28->GetXaxis()->CenterTitle();
  h28->GetYaxis()->SetTitle("Entries/N");
  h28->GetYaxis()->CenterTitle();
  h28->SetLineColor(1);
  h28->SetMarkerColor(1);
  h28->SetMarkerStyle(kOpenCircle);
  h28->Sumw2();
  h28->SetMaximum(0.4);
  h28->Draw("PE");
  Double_t sc29 = 1/h29->Integral();
  h29->Scale(sc29);
  h29->SetLineColor(2);
  h29->SetMarkerColor(2);
  h29->SetMarkerStyle(kOpenDiamond);
  h29->Draw("PE same");
  Double_t sc30 = 1/h30->Integral();
  h30->Scale(sc30);
  h30->SetLineColor(4);
  h30->SetMarkerColor(4);
  h30->SetMarkerStyle(kOpenTriangleUp);
  h30->Draw("PE same");
  Double_t sc31 = 1/h31->Integral();
  h31->Scale(sc31);
  h31->SetLineColor(6);
  h31->SetMarkerColor(6);
  h31->SetMarkerStyle(kOpenTriangleDown);
  h31->Draw("PE same");

  TLegend *leg7 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg7);
  leg7->AddEntry(h28,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg7->AddEntry(h29,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg7->AddEntry(h30,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg7->AddEntry(h31,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg7->Draw("same");

  c14->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc32 = 1/h32->Integral();
  h32->Scale(sc32);
  h32->SetMarkerSize(1.2);
  h32->GetXaxis()->SetTitle("HFsum (GeV)");
  h32->GetXaxis()->CenterTitle();
  h32->GetYaxis()->SetTitle("Entries/N");
  h32->GetYaxis()->CenterTitle();
  h32->SetLineColor(1);
  h32->SetMarkerColor(1);
  h32->SetMarkerStyle(kOpenCircle);
  h32->Sumw2();
  h32->SetMaximum(0.4);
  h32->Draw("PE");
  Double_t sc33 = 1/h33->Integral();
  h33->Scale(sc33);
  h33->SetLineColor(2);
  h33->SetMarkerColor(2);
  h33->SetMarkerStyle(kOpenDiamond);
  h33->Draw("PE same");
  Double_t sc34 = 1/h34->Integral();
  h34->Scale(sc34);
  h34->SetLineColor(4);
  h34->SetMarkerColor(4);
  h34->SetMarkerStyle(kOpenTriangleUp);
  h34->Draw("PE same");
  Double_t sc35 = 1/h35->Integral();
  h35->Scale(sc35);
  h35->SetLineColor(6);
  h35->SetMarkerColor(6);
  h35->SetMarkerStyle(kOpenTriangleDown);
  h35->Draw("PE same");

  TLegend *leg8 = new TLegend(0.35,0.70,0.85,0.85);
  SetLegendStyle(leg8);
  leg8->AddEntry(h32,"#eta_{jet,1} > 0 & #eta_{jet,2} > 0","lp");
  leg8->AddEntry(h33,"#eta_{jet,1} < 0 & #eta_{jet,2} < 0","lp");
  leg8->AddEntry(h34,"#eta_{jet,1} > 0 & #eta_{jet,2} < 0","lp");
  leg8->AddEntry(h35,"#eta_{jet,1} < 0 & #eta_{jet,2} > 0","lp");
  leg8->Draw("same");

  c15->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc36 = 1/h36->Integral();
  h36->Scale(sc36);
  h36->SetMarkerSize(1.2);
  h36->GetXaxis()->SetTitle("#Delta^{2} (GeV^{2})");
  h36->GetXaxis()->CenterTitle();
  h36->GetYaxis()->SetTitle("Entries/N");
  h36->GetYaxis()->CenterTitle();
  h36->SetLineColor(1);
  h36->SetMarkerColor(1);
  h36->SetMarkerStyle(kOpenCircle);
  h36->Sumw2();
  h36->SetMaximum(0.2);
  h36->Draw("PE");

  c16->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  Double_t sc37 = 1/h37->Integral();
  h37->Scale(sc37);
  h37->SetMarkerSize(1.2);
  h37->GetXaxis()->SetTitle("P ((GeV/c)^{2})");
  h37->GetXaxis()->CenterTitle();
  h37->GetYaxis()->SetTitle("Entries/N");
  h37->GetYaxis()->CenterTitle();
  h37->SetLineColor(1);
  h37->SetMarkerColor(1);
  h37->SetMarkerStyle(kOpenCircle);
  h37->Sumw2();
  h37->SetMaximum(0.2);
  h37->Draw("PE");

  c17->cd();
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  th5->GetXaxis()->SetTitle("#Delta (GeV^{2})");
  th5->GetXaxis()->CenterTitle();
  th5->GetYaxis()->SetTitle("P ((GeV/c)^{2})");
  th5->GetYaxis()->CenterTitle();
  th5->Draw("COLZ");

  c2->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/aj_norechitcut_forHF.png");
  c3->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/Xgamma_norechitcut_forHF.png");
  c4->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/ratept_norechitcut_forHF.png");
  c5->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/deltaphi_norechitcut_forHF.png");
  c6->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/djpt_norechitcut_forHF.png");
  c7->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/djptVSxgammabothplus_norechitcut_forHF.png");
  c8->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/djptVSxgammabothminus_norechitcut_forHF.png"); 
  c9->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/HFplus_forHF.png");
  c10->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/HFminus_forHF.png");
  c11->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/HFsum_forHF.png");
  c12->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/HFplus_forHF_4jtetaselection.png");
  c13->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/HFminus_forHF_4jtetaselection.png");
  c14->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/HFsum_forHF_4jtetaselection.png"); 
  c15->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/proton_recoil_momentumsquare.png");
  c16->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/dijet_relative_momentum.png");  
  c17->SaveAs("/u/user/bekim/758p3_UPC/src/Muon_test/plotsforetajet/prm_Vs_drm.png");

  return;
}
