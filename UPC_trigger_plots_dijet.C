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
#include "JpsiFunc.h"

using namespace std;
void UPC_trigger_plots_dijet()
{
  TTree *fChain;
  TTree *fChain2;
  Int_t fCurrent;
  TFile *data = TFile::Open("/home/kbg777/CMSwork/skimmedFiles/upcDiJetSkim2_HLT_HIUPCL1SingleEG5NotHF2_v1.root","read");
 
  TTree *jetTree = (TTree*) data -> Get("dijet");
  TTree *evtTree = (TTree*) data -> Get("evt");

  TCanvas* c1 = new TCanvas("c1","",700,700);
  TCanvas* c2 = new TCanvas("c2","",700,700);
  TCanvas* c3 = new TCanvas("c3","",700,700);
  TCanvas* c4 = new TCanvas("c4","",700,700);
  TCanvas* c5 = new TCanvas("c5","",700,700);
  TCanvas* c6 = new TCanvas("c6","",700,700);
  TCanvas* c7 = new TCanvas("c7","",700,700);
  TCanvas* c8 = new TCanvas("c8","",700,700);
  TCanvas* c9 = new TCanvas("c9","",700,700);
  TCanvas* c10 = new TCanvas("c10","",700,700);
  TCanvas* c11 = new TCanvas("c11","",700,700);
  TCanvas* c12 = new TCanvas("c12","",700,700);
  TCanvas* c13 = new TCanvas("c13","",700,700);
  TCanvas* c14 = new TCanvas("c14","",700,700);
  TCanvas* c15 = new TCanvas("c15","",700,700);
  TCanvas* c16 = new TCanvas("c16","",700,700);
  TCanvas* c17 = new TCanvas("c17","",700,700);
  TCanvas* c18 = new TCanvas("c18","",700,700);
  TCanvas* cv1 = new TCanvas("cv1","",700,700);
  TCanvas* cv2 = new TCanvas("cv2","",700,700);
  TCanvas* cv3 = new TCanvas("cv3","",700,700);
  TCanvas* cv4 = new TCanvas("cv4","",700,700);
  TCanvas* cv5 = new TCanvas("cv5","",700,700);
  TCanvas* cv6 = new TCanvas("cv6","",700,700);
  TCanvas* cv7 = new TCanvas("cv7","",700,700);
  TCanvas* cv8 = new TCanvas("cv8","",700,700);
  TCanvas* cv9 = new TCanvas("cv9","",700,700);
  TCanvas* cv10 = new TCanvas("cv10","",700,700);
  TCanvas* cv11 = new TCanvas("cv11","",700,700);
  TCanvas* cv12 = new TCanvas("cv12","",700,700);
  TCanvas* cv13 = new TCanvas("cv13","",700,700);
  TCanvas* cv14 = new TCanvas("cv14","",700,700);
  TCanvas* cv15 = new TCanvas("cv15","",700,700);
  TCanvas* cv16 = new TCanvas("cv16","",700,700);
  TCanvas* cv17 = new TCanvas("cv17","",700,700);
  TCanvas* cv18 = new TCanvas("cv18","",700,700);

  TH2F *h1 = new TH2F("h1", " ", 50, 0, 100, 50, 0, 100);
  h1->GetXaxis()->SetTitle("p_{T,1} (GeV/c)");
  h1->GetXaxis()->SetTitleOffset(1.3);
  h1->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->Sumw2();
  Double_t sc1 = 1/h1->Integral();
  h1->Scale(sc1);
  jetTree->AddFriend(evtTree);
  cv1->cd();
  jetTree->Draw("pt2:pt1 >> h1", "hfplus < 5 && hfminus < 5");
  c1->cd();
  h1->Draw("colz");
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h1 entries : " <<  h1->GetEntries() << endl;
  c1->SaveAs("pt1_vs_pt2_bothHF_under_5GeV.C");

  TH2F *h2 = new TH2F("h2", " ", 50, 0, 100, 50, 0, 100);
  h2->GetXaxis()->SetTitle("p_{T,1} (GeV/c)");
  h2->GetXaxis()->SetTitleOffset(1.3);
  h2->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  h2->GetYaxis()->SetTitleOffset(1.3);
  h2->Sumw2();
  Double_t sc2 = 1/h2->Integral();
  h2->Scale(sc2);
  cv2->cd();
  jetTree->Draw("pt2:pt1 >> h2", "hfplus < 5 && hfminus > 5");
  c2->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h2 entries : " <<  h2->GetEntries() << endl;
  h2->Draw("colz");
  c2->SaveAs("pt1_vs_pt2_HFminus_over_5GeV.C");
 
  TH2F *h3 = new TH2F("h3", " ", 50, 0, 100, 50, 0, 100);
  h3->GetXaxis()->SetTitle("p_{T,1} (GeV/c)");
  h3->GetXaxis()->SetTitleOffset(1.3);
  h3->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  h3->GetYaxis()->SetTitleOffset(1.3);
  h3->Sumw2();
  Double_t sc3 = 1/h3->Integral();
  h3->Scale(sc3);
  cv3->cd();
  jetTree->Draw("pt2:pt1 >> h3", "hfplus > 5 && hfminus < 5");
  c3->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h3 entries : " <<  h3->GetEntries() << endl;
  h3->Draw("colz");
  c3->SaveAs("pt1_vs_pt2_HFplus_over_5GeV.C");
 
  TH2F *h4 = new TH2F("h4", " ", 50, -3.5, 3.5, 50, -3.5, 3.5);
  h4->GetXaxis()->SetTitle("#eta_{1}");
  h4->GetXaxis()->SetTitleOffset(1.3);
  h4->GetYaxis()->SetTitle("#eta_{2}");
  h4->GetYaxis()->SetTitleOffset(1.3);
  h4->Sumw2();
  Double_t sc4 = 1/h4->Integral();
  h4->Scale(sc4);
  cv4->cd();
  jetTree->Draw("eta2:eta1 >> h4", "hfplus < 5 && hfminus < 5");
  c4->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h4 entries : " <<  h4->GetEntries() << endl;
  h4->Draw("colz");
  c4->SaveAs("eta1_vs_eta2_both_under_5GeV.C");
  
  TH2F *h5 = new TH2F("h5", " ", 50, -3.5, 3.5, 50, -3.5, 3.5);
  h5->GetXaxis()->SetTitle("#eta_{1}");
  h5->GetXaxis()->SetTitleOffset(1.3);
  h5->GetYaxis()->SetTitle("#eta_{2}");
  h5->GetYaxis()->SetTitleOffset(1.3);
  h5->Sumw2();
  Double_t sc5 = 1/h5->Integral();
  h5->Scale(sc5);
  cv5->cd();
  jetTree->Draw("eta2:eta1 >> h5", "hfplus < 5 && hfminus > 5");
  c5->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h5 entries : " <<  h5->GetEntries() << endl;
  h5->Draw("colz");
  c5->SaveAs("eta1_vs_eta2_HFminus_over_5GeV.C");

  TH2F *h6 = new TH2F("h6", " ", 50, -3.5, 3.5, 50, -3.5, 3.5);
  h6->GetXaxis()->SetTitle("#eta_{1}");
  h6->GetXaxis()->SetTitleOffset(1.3);
  h6->GetYaxis()->SetTitle("#eta_{2}");
  h6->GetYaxis()->SetTitleOffset(1.3);
  h6->Sumw2();
  Double_t sc6 = 1/h6->Integral();
  h6->Scale(sc6);
  cv6->cd();
  jetTree->Draw("eta2:eta1 >> h6", "hfplus > 5 && hfminus < 5");
  c6->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h6 entries : " <<  h6->GetEntries() << endl;
  h6->Draw("colz");
  c6->SaveAs("eta1_vs_eta2_HFplus_over_5GeV.C");

  TH1F *h7 = new TH1F("h7", " ", 50, 0, 100);
  h7->SetLineColor(2);
  h7->SetMarkerColor(2);
  h7->SetMarkerStyle(kOpenCircle);
  h7->GetXaxis()->SetTitle("#Delta p_{T}");
  h7->GetXaxis()->SetTitleOffset(1.3);
  h7->Sumw2();
  Double_t sc7 = 1/h7->Integral();
  h7->Scale(sc7);
  cv7->cd();
  jetTree->Draw("dpt >> h7", "hfplus < 5 && hfminus < 5");
  c7->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h7 entries : " <<  h7->GetEntries() << endl;
  h7->Draw("PE");
  c7->SetLogy();
  c7->SaveAs("delta_pt_bothHF_under_5GeV.C");

  TH1F *h8 = new TH1F("h8", " ", 50, 0, 100);
  h8->SetLineColor(2);
  h8->SetMarkerColor(2);
  h8->SetMarkerStyle(kOpenCircle);
  h8->GetXaxis()->SetTitle("#Delta p_{T}");
  h8->GetXaxis()->SetTitleOffset(1.3);
  h8->Sumw2();
  Double_t sc8 = 1/h8->Integral();
  h8->Scale(sc8);
  cv8->cd();
  jetTree->Draw("dpt >> h8", "hfplus < 5 && hfminus > 5");
  c8->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h8 entries : " <<  h8->GetEntries() << endl;
  h8->Draw("PE");
  c8->SetLogy();
  c8->SaveAs("delta_pt_HFminus_over_5GeV.C");

  TH1F *h9 = new TH1F("h9", " ", 50, 0, 100);
  h9->SetLineColor(2);
  h9->SetMarkerColor(2);
  h9->SetMarkerStyle(kOpenCircle);
  h9->GetXaxis()->SetTitle("#Delta p_{T}");
  h9->GetXaxis()->SetTitleOffset(1.3);
  h9->Sumw2();
  Double_t sc9 = 1/h9->Integral();
  h9->Scale(sc9);
  cv9->cd();
  jetTree->Draw("dpt >> h9", "hfplus > 5 && hfminus < 5");
  c9->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h9 entries : " <<  h9->GetEntries() << endl;
  h9->Draw("PE");
  c9->SetLogy();
  c9->SaveAs("delta_pt_HFplus_over_5GeV.C");

  TH1F *h10 = new TH1F("h10", " ", 50, -3.5, 3.5);
  h10->SetLineColor(2);
  h10->SetMarkerColor(2);
  h10->SetMarkerStyle(kOpenCircle);
  h10->GetXaxis()->SetTitle("#eta");
  h10->GetXaxis()->SetTitleOffset(1.3);
  h10->Sumw2();
  Double_t sc10 = 1/h10->Integral();
  h10->Scale(sc10);
  cv10->cd();
  jetTree->Draw("eta >> h10", "hfplus < 5 && hfminus < 5");
  c10->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h10 entries : " <<  h10->GetEntries() << endl;
  h10->Draw("PE");
  c10->SetLogy();
  c10->SaveAs("dijet_eta_bothHF_under_5GeV.C");

  TH1F *h11 = new TH1F("h11", " ", 50, -3.5, 3.5);
  h11->SetLineColor(2);
  h11->SetMarkerColor(2);
  h11->SetMarkerStyle(kOpenCircle);
  h11->GetXaxis()->SetTitle("#eta");
  h11->GetXaxis()->SetTitleOffset(1.3);
  h11->Sumw2();
  Double_t sc11 = 1/h11->Integral();
  h11->Scale(sc11);
  cv11->cd();
  jetTree->Draw("eta >> h11", "hfplus < 5 && hfminus > 5");
  c11->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h11 entries : " <<  h11->GetEntries() << endl;
  h11->Draw("PE");
  c11->SetLogy();
  c11->SaveAs("dijet_eta_HFminus_over_5GeV.C");
  
  TH1F *h12 = new TH1F("h12", " ", 50, -3.5, 3.5);
  h12->SetLineColor(2);
  h12->SetMarkerColor(2);
  h12->SetMarkerStyle(kOpenCircle);
  h12->GetXaxis()->SetTitle("#eta");
  h12->GetXaxis()->SetTitleOffset(1.3);
  h12->Sumw2();
  Double_t sc12 = 1/h12->Integral();
  h12->Scale(sc12);
  cv12->cd();
  jetTree->Draw("eta >> h12", "hfplus > 5 && hfminus < 5");
  c12->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("PE");
  cout <<"h12 entries : " <<  h12->GetEntries() << endl;
  h12->Draw();
  c12->SetLogy();
  c12->SaveAs("dijet_eta_HFplus_over_5GeV.C");

  TH1F *h13 = new TH1F("h13", " ", 50, 0, 1);
  h13->SetLineColor(2);
  h13->SetMarkerColor(2);
  h13->SetMarkerStyle(kOpenCircle);
  h13->GetXaxis()->SetTitle("A_{j}");
  h13->GetXaxis()->SetTitleOffset(1.3);
  h13->Sumw2();
  Double_t sc13 = 1/h13->Integral();
  h13->Scale(sc13);
  cv13->cd();
  jetTree->Draw("aj >> h13", "hfplus < 5 && hfminus < 5");
  c13->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("PE");
  cout <<"h13 entries : " <<  h13->GetEntries() << endl;
  h13->Draw("PE");
  c13->SaveAs("dijet_Aj_bothHF_under_5GeV.C");

  TH1F *h14 = new TH1F("h14", " ", 50, 0, 1);
  h14->SetLineColor(2);
  h14->SetMarkerColor(2);
  h14->SetMarkerStyle(kOpenCircle);
  h14->GetXaxis()->SetTitle("A_{j}");
  h14->GetXaxis()->SetTitleOffset(1.3);
  h14->Sumw2();
  Double_t sc14 = 1/h14->Integral();
  h14->Scale(sc14);
  cv14->cd();
  jetTree->Draw("aj >> h14", "hfplus < 5 && hfminus > 5");
  c14->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("PE");
  cout <<"h14 entries : " <<  h14->GetEntries() << endl;
  h14->Draw("PE");
  c14->SaveAs("dijet_Aj_HFminus_over_5GeV.C");

  TH1F *h15 = new TH1F("h15", " ", 50, 0, 1);
  h15->SetLineColor(2);
  h15->SetMarkerColor(2);
  h15->SetMarkerStyle(kOpenCircle);
  h15->GetXaxis()->SetTitle("A_{j}");
  h15->GetXaxis()->SetTitleOffset(1.3);
  h15->Sumw2();
  Double_t sc15 = 1/h15->Integral();
  h15->Scale(sc15);
  cv15->cd();
  jetTree->Draw("aj >> h15", "hfplus > 5 && hfminus < 5");
  c15->cd();
  gStyle->SetOptStat(0);
  gStyle->SetTitle("PE");
  cout <<"h15 entries : " <<  h15->GetEntries() << endl;
  h15->Draw("PE");
  c15->SaveAs("dijet_Aj_HFplus_over_5GeV.C");

  TH2F *h16 = new TH2F("h16", " ", 50, 0, 3.5, 50, 0, 5);
  h16->GetXaxis()->SetTitle("#Delta#phi");
  h16->GetXaxis()->SetTitleOffset(1.3);
  h16->GetYaxis()->SetTitle("#Delta#eta");
  h16->GetYaxis()->SetTitleOffset(1.3);
  h16->Sumw2();
  Double_t sc16 = 1/h16->Integral();
  h16->Scale(sc16);
  cv16->cd();
  jetTree->Draw("deta:dphi >> h16", "hfplus < 5 && hfminus < 5");
  c16->cd();
  h16->Draw("colz");
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h16 entries : " <<  h16->GetEntries() << endl;
  c16->SaveAs("delta_phi_vs_delta_eta_bothHF_under_5GeV.C");

  TH2F *h17 = new TH2F("h17", " ", 50, 0, 3.5, 50, 0, 5);
  h17->GetXaxis()->SetTitle("#Delta#phi");
  h17->GetXaxis()->SetTitleOffset(1.3);
  h17->GetYaxis()->SetTitle("#Delta#eta");
  h17->GetYaxis()->SetTitleOffset(1.3);
  h17->Sumw2();
  Double_t sc17 = 1/h17->Integral();
  h17->Scale(sc17);
  cv17->cd();
  jetTree->Draw("deta:dphi >> h17", "hfplus < 5 && hfminus > 5");
  c17->cd();
  h17->Draw("colz");
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h17 entries : " <<  h17->GetEntries() << endl;
  c17->SaveAs("delta_phi_vs_delta_eta_HFminus_over_5GeV.C");

  TH2F *h18 = new TH2F("h18", " ", 50, 0, 3.5, 50, 0, 5);
  h18->GetXaxis()->SetTitle("#Delta#phi");
  h18->GetXaxis()->SetTitleOffset(1.3);
  h18->GetYaxis()->SetTitle("#Delta#eta");
  h18->GetYaxis()->SetTitleOffset(1.3);
  h18->Sumw2();
  Double_t sc18 = 1/h18->Integral();
  h18->Scale(sc18);
  cv18->cd();
  jetTree->Draw("deta:dphi >> h18", "hfplus > 5 && hfminus < 5");
  c18->cd();
  h18->Draw("colz");
  gStyle->SetOptStat(0);
  gStyle->SetTitle("");
  cout <<"h18 entries : " <<  h18->GetEntries() << endl;
  c18->SaveAs("delta_phi_vs_delta_eta_HFplus_over_5GeV.C");


/*
  c1->SaveAs("pt1_vs_pt2_bothHF_under_5GeV.C");
  c2->SaveAs("pt1_vs_pt2_HFminus_over_5GeV.C");
  c3->SaveAs("pt1_vs_pt2_HFplus_over_5GeV.C");
  c4->SaveAs("eta1_vs_eta2_both_under_5GeV.C");
  c5->SaveAs("eta1_vs_eta2_HFminus_over_5GeV.C");
  c6->SaveAs("eta1_vs_eta2_HFplus_over_5GeV.C");
  c7->SaveAs("delta_pt_bothHF_under_5GeV.C");
  c8->SaveAs("delta_pt_HFminus_over_5GeV.C");
  c9->SaveAs("delta_pt_HFplus_over_5GeV.C");
*/
  /*
  TH1F *h1 = new TH1F("h1"," ", 25, 0, 100);
  TH1F *h2 = new TH1F("h2"," ", 25, 0, 100);
  TH1F *h3 = new TH1F("h3"," ", 25, 0, 1.5);
  TH1F *h4 = new TH1F("h4"," ", 25, 0, 1.5);
  TH1F *h5 = new TH1F("h5"," ", 25, 0, 100);
  TH1F *h6 = new TH1F("h6"," ", 25, 0, 100);
  TH1F *h7 = new TH1F("h7"," ", 25, 0, 100);
  TH1F *h8 = new TH1F("h8"," ", 25, 0, 100);
  TH1F *h9 = new TH1F("h9"," ", 10, 0, 1);
  h9->SetMaximum(0.5);
  TH1F *h10 = new TH1F("h10"," ", 10, 0, 1);
  //h10->SetMaximum(0.5);
  TH1F *h11 = new TH1F("h11"," ", 25, 0, 3.2);
  TH1F *h12 = new TH1F("h12"," ", 25, 0, 3.2);
  TH1F *h13 = new TH1F("h13"," ", 25, 0, 100);
  TH1F *h14 = new TH1F("h14"," ", 25, 0, 100);
  TH1F *h15 = new TH1F("h15"," ", 25, 0, 5);
  TH1F *h16 = new TH1F("h16"," ", 25, 0, 5);
  TH2F *h17 = new TH2F("h17"," ", 25, 0, 3.5, 25, 0, 5);
  TH1F *h18 = new TH1F("h18"," ", 50, -3.5, 3.5);
  TH1F *h19 = new TH1F("h19"," ", 50, -3.5, 3.5);
  TH1F *h20 = new TH1F("h20"," ", 25, 0, 1);
  TH1F *h21 = new TH1F("h21"," ", 25, 0, 1);
  TH1F *h22 = new TH1F("h22"," ", 10, 0, 1);
  TH1F *h23 = new TH1F("h23"," ", 10, 0, 1);
*/
  return;
}
