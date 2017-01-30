#include <ctime>
#include <TMath.h>
#include <TLorentzVector.h>
#include "commonUtility.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <TCut.h>

using std::stringstream;
using std::vector;
using std::abs;
  
int kAll = 0;
int kSDfwd = 1;  // HFp > 5,   HFm < 5
int kSDbwd = 2;  // HFp < 5,   HFm > 5


void analysis_midRapidity() {
  using namespace std;
  TH1::SetDefaultSumw2();

  TFile *f1 = new TFile("upcDiJetSkim170124_trigHLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1_jetCollectionak5PFJetAnalyzer_minJetPt0_2017y_1m_24d_12h_30m.root");
  TTree* dijet = (TTree*)f1->Get("dijet");
  TTree* Track = (TTree*)f1->Get("Track");
  TTree* fullTrkTree = (TTree*)f1->Get("fullTrkTree");
  TTree* evt = (TTree*)f1->Get("evt");
  Track->AddFriend("evt");
  Track->AddFriend("fullTrkTree");
  Track->AddFriend("dijet");
  
  TCut fwdHFcut = "(evt.hfplus > 5) && (evt.hfminus < 5) ";
  TCut bwdHFcut = "(evt.hfplus < 5) && (evt.hfminus > 5) ";
  TCut UPCcut = "(evt.hfplus < 5) && (evt.hfminus < 5) ";

  Track->SetAlias("nTrkAllEta", "nTrketam2to2p5 + nTrketam1p5to2 + nTrketam1to1p5 + nTrketam0p5to1 + nTrketam0to0p5 + nTrketa0to0p5 + nTrketa0p5to1 + nTrketa1to1p5 + nTrketa1p5to2 + nTrketa2to2p5");
  Track->SetAlias("nTrkForward", "nTrketa1p5to2 + nTrketa2to2p5 + nTrketam2to2p5 + nTrketam1p5to2");

  TCanvas* c3 = new TCanvas("c3","",400,400);
  TH1D* hdphi = new TH1D("hdphi",";#Delta #phi",20,0,3.14);

  TCut jetPtCut = "dijet.pt1>30 && dijet.pt2>30 && abs(dijet.eta1)<1.5 && abs(dijet.eta2)<1.5";
  Track->Draw(Form("dijet.dphi>>%s",hdphi->GetName()), UPCcut && jetPtCut && "nTrkForward==0");
  handsomeTH1(hdphi,1);
  hdphi->Draw();

  TCanvas* c5 = new TCanvas("c5","",400,400);
  TH1D* jeteta= new TH1D("jeteta",";#eta;Entries",40,-2.4,2.4);
  Track->Draw(Form("dijet.y>>%s",jeteta->GetName()), UPCcut && jetPtCut && "nTrkForward==0");
  handsomeTH1(jeteta,1);
  jeteta->Draw();
  
  TCanvas* c4 = new TCanvas("c4","",400,400);
  TH1D* aj15= new TH1D("ajNoTrk15",";A_{J};Entries",20,0,0.5);
  Track->Draw(Form("dijet.aj>>%s",aj15->GetName()), UPCcut && jetPtCut && "nTrkForward == 0" && "dijet.dphi>2*3.141592/3");
  handsomeTH1(aj15,1);
  cout << aj15->GetEntries() << endl;
  aj15->Draw();

  
}
