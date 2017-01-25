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


void analysis_rapidityGap( ) {   
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

  Track->SetAlias("nTrkBelow15", "nTrketam2to2p5 + nTrketam1p5to2 + nTrketam1to1p5 + nTrketam0p5to1 + nTrketam0to0p5 + nTrketa0to0p5 + nTrketa0p5to1 + nTrketa1to1p5");
  Track->SetAlias("nTrkBelow10", "nTrketam2to2p5 + nTrketam1p5to2 + nTrketam1to1p5 + nTrketam0p5to1 + nTrketam0to0p5 + nTrketa0to0p5 + nTrketa0p5to1");
  Track->SetAlias("nTrkBelow05", "nTrketam2to2p5 + nTrketam1p5to2 + nTrketam1to1p5 + nTrketam0p5to1 + nTrketam0to0p5 + nTrketa0to0p5");
  Track->SetAlias("nTrkBelow00", "nTrketam2to2p5 + nTrketam1p5to2 + nTrketam1to1p5 + nTrketam0p5to1 + nTrketam0to0p5");
  TH1D* nTrkBkwMid[4] ; 
  

  nTrkBkwMid[kAll] = new TH1D("nTrkBkwMid_all",";N_{trk}^{-2.5 < eta< 1.5};Entries",100,-.5,99.5);
  nTrkBkwMid[kSDfwd] = (TH1D*)nTrkBkwMid[kAll]->Clone("nTrkBkwMid_sdfwd");
  nTrkBkwMid[kSDbwd] = (TH1D*)nTrkBkwMid[kAll]->Clone("nTrkBkwMid_sdbwd");
  Track->Draw(Form("nTrkBelow15>>%s",nTrkBkwMid[kAll]->GetName() ) );
  Track->Draw(Form("nTrkBelow15>>%s",nTrkBkwMid[kSDfwd]->GetName() ), fwdHFcut );
  Track->Draw(Form("nTrkBelow15>>%s",nTrkBkwMid[kSDbwd]->GetName() ), bwdHFcut );
  TCanvas* c1 = new TCanvas("c1","",400,400);
  handsomeTH1(nTrkBkwMid[kAll],1);
  handsomeTH1(nTrkBkwMid[kSDfwd],2);
  handsomeTH1(nTrkBkwMid[kSDbwd],4);
  scaleInt(nTrkBkwMid[kAll]);
  scaleInt(nTrkBkwMid[kSDfwd]);
  scaleInt(nTrkBkwMid[kSDbwd]);
  nTrkBkwMid[kAll]->Draw("hist");
  nTrkBkwMid[kSDfwd]->Draw("same");
  nTrkBkwMid[kSDbwd]->Draw("same");
  gPad->SetLogy();
  TLegend* leg1 = new TLegend(0.5046176,0.65,0.9,0.9,NULL,"brNDC");
  easyLeg(leg1,"Trgi: EG5, NotHF2, SingleTrack");
  leg1->AddEntry(nTrkBkwMid[kAll],"All");
  leg1->AddEntry(nTrkBkwMid[kSDfwd],"HF+ > 5GeV,  HF- < 5GeV","pl");
  leg1->AddEntry(nTrkBkwMid[kSDbwd],"HF+ < 5GeV,  HF- > 5GeV","pl");
  leg1->Draw();

  TCanvas* c2 = new TCanvas("c2","",400,400);
  TH1D* nTrkUpTo15Dijet = (TH1D*)nTrkBkwMid[kAll]->Clone("nTrkUpTo15Dijet");
  TH1D* nTrkUpTo10Dijet = (TH1D*)nTrkBkwMid[kAll]->Clone("nTrkUpTo10Dijet");
  TH1D* nTrkUpTo05Dijet = (TH1D*)nTrkBkwMid[kAll]->Clone("nTrkUpTo05Dijet");
  TH1D* nTrkUpTo00Dijet = (TH1D*)nTrkBkwMid[kAll]->Clone("nTrkUpTo00Dijet");
  TCut fwdDijetCut = "dijet.pt1>20 && dijet.pt2>20 && dijet.y > 1.5";

  nTrkUpTo15Dijet->Reset();
  nTrkUpTo10Dijet->Reset();
  nTrkUpTo05Dijet->Reset();
  nTrkUpTo00Dijet->Reset();
  Track->Draw(Form("nTrkBelow15>>%s",nTrkUpTo15Dijet->GetName()), fwdHFcut && fwdDijetCut );
  Track->Draw(Form("nTrkBelow10>>%s",nTrkUpTo10Dijet->GetName()), fwdHFcut && fwdDijetCut );
  Track->Draw(Form("nTrkBelow05>>%s",nTrkUpTo05Dijet->GetName()), fwdHFcut && fwdDijetCut );
  Track->Draw(Form("nTrkBelow00>>%s",nTrkUpTo00Dijet->GetName()), fwdHFcut && fwdDijetCut );
  handsomeTH1(nTrkUpTo15Dijet,1);
  handsomeTH1(nTrkUpTo10Dijet,2);
  handsomeTH1(nTrkUpTo05Dijet,4);
  handsomeTH1(nTrkUpTo00Dijet,7);
  nTrkUpTo00Dijet->SetXTitle("Number of tracks");
  nTrkUpTo00Dijet->Draw("");
  nTrkUpTo15Dijet->Draw("same");
  nTrkUpTo10Dijet->Draw("same");
  nTrkUpTo05Dijet->Draw("same");
  gPad->SetLogy();
  TLegend* leg2 = new TLegend(0.5046176,0.65,0.9,0.9,NULL,"brNDC");
  easyLeg(leg2,"N_{trk} in ...");
  leg2->AddEntry(nTrkUpTo15Dijet,"#eta < 1.5","pl");
  leg2->AddEntry(nTrkUpTo10Dijet,"#eta < 1.0","pl");
  leg2->AddEntry(nTrkUpTo05Dijet,"#eta < 0.5","pl");
  leg2->AddEntry(nTrkUpTo00Dijet,"#eta < 0","pl");
  leg2->Draw();


  TCanvas* c3 = new TCanvas("c3","",400,400);
  TH1D* aj15= new TH1D("ajNoTrk15",";A_{J};Entries",20,0,0.5);
  TH1D* aj10= (TH1D*)aj15->Clone("ajNoTrk10");
  TH1D* aj05= (TH1D*)aj15->Clone("ajNoTrk05");
  TH1D* aj00= (TH1D*)aj15->Clone("ajNoTrk00");

  Track->Draw(Form("dijet.aj>>%s",aj15->GetName()), fwdHFcut && fwdDijetCut && "nTrkBelow15==0");
  Track->Draw(Form("dijet.aj>>%s",aj10->GetName()), fwdHFcut && fwdDijetCut && "nTrkBelow10==0");
  Track->Draw(Form("dijet.aj>>%s",aj05->GetName()), fwdHFcut && fwdDijetCut && "nTrkBelow05==0");
  Track->Draw(Form("dijet.aj>>%s",aj00->GetName()), fwdHFcut && fwdDijetCut && "nTrkBelow00==0");
  scaleInt(aj15);   scaleInt(aj10);   scaleInt(aj05);   scaleInt(aj00); 
  handsomeTH1(aj15,1);
  handsomeTH1(aj10,2);
  handsomeTH1(aj05,4);
  handsomeTH1(aj00,7);
  aj05->Draw();
  aj00->Draw("same");
  aj15->Draw("same");
  aj10->Draw("same");
  cout << "aj15 ave = " << aj15->GetMean() << " +/- " << aj15->GetMeanError() << endl;
  cout << "aj10 ave = " << aj10->GetMean() << " +/- " << aj10->GetMeanError() << endl;
  cout << "aj05 ave = " << aj05->GetMean() << " +/- " << aj05->GetMeanError() << endl;
  cout << "aj00 ave = " << aj00->GetMean() << " +/- " << aj00->GetMeanError() << endl;
}
