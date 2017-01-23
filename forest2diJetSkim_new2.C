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

void forest2diJetSkim_new2(
		      TString fname = "/u/user/bekim/758p3_UPC/src/Muon_test/4_UPCTriggers_170119.root",
		      TString outputFname = "upcDiJetSkim170123", 
		      TString trig = "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1",  // "HLT_HIUPCL1SingleEG5NotHF2_v1",    // "HLT_HIUPCSingleEG5NotHF2Pixel_SingleTrack_v1" //  "HLT_HIUPCL1SingleEG5NotHF2_v1"
		      TString jetCollection = "ak5PFJetAnalyzer", // "akPu5PFJetAnalyzer",
		      float minjPt = 0,
		      int nevt=-1
		      ) {  
  
  using namespace std;
  TFile *f1 = new TFile(fname.Data());
  
  TTree *HltTree = (TTree*)f1->Get("hltanalysis/HltTree");
  TTree *hiTree  = (TTree*)f1->Get("hiEvtAnalyzer/HiTree");
  TTree *trackTree = (TTree*)f1->Get("anaTrack/trackTree");
  TTree *t = (TTree*)f1->Get(Form("%s/t",jetCollection.Data()));
  //  TTree *trackTree  = (TTree*)f1->Get("ppTrack/trackTree");
  //  TTree *akpu5pf = (TTree*)f1->Get("t");
  //  TTree *pfTree  = (TTree*)f1->Get("pfcandAnalyzer/pfTree");
  TString dayTime = getDayAndTime();
  TFile* newfile = new TFile(Form("skimmedFiles/%s_trig%s_jetCollection%s_minJetPt%d_%s.root",outputFname.Data(), trig.Data(), jetCollection.Data(), (int)minjPt, dayTime.Data()  ),"recreate");
   
  t->AddFriend(HltTree);
  t->AddFriend(hiTree);
  t->AddFriend(trackTree);
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
  Float_t jtenergy[200] = {0};
  Float_t jtmass[200] = {0};
  Float_t pt[200] = {0};
  Float_t eta[200] = {0};
  Float_t y[200] = {0};
  Float_t phi[200] = {0};
  Float_t m[200] = {0};
  Float_t e[200] = {0};
  Float_t mass[200] = {0}; 
  Float_t jtpt2[200] = {0};
  Float_t jteta2[200] = {0};
  Float_t jty2[200] = {0};
  Float_t jtphi2[200] = {0};
  Float_t jtm2[200] = {0};
  Float_t jtenergy2[200] = {0};
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
  TBranch *b_jtenergy;

  t->SetBranchAddress("nref", &nref, &b_nref);
  t->SetBranchAddress("jtpt", jtpt, &b_jtpt);
  t->SetBranchAddress("jteta", jteta, &b_jteta);
  t->SetBranchAddress("jty", jty, &b_jty);
  t->SetBranchAddress("jtphi", jtphi, &b_jtphi);
  t->SetBranchAddress("jtm", jtm, &b_jtm);
  t->SetBranchAddress("jtenergy", jtenergy, &b_jtenergy);

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
   Float_t         trkPt[33885];
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkEta;
   TBranch        *b_trkPt;
   trackTree->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   trackTree->SetBranchAddress("trkEta", trkEta, &b_trkEta);
   trackTree->SetBranchAddress("trkPt", trkPt, &b_trkPt);

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
   trkTree->Branch("track",&Track.nTrack,nTrkString);

   int ntrk;
   const int MAXtrk = 50000; // to accomodate 100 smeared trks, need to be careful with ram
   float trkpt[MAXtrk];
   float trketa[MAXtrk];

   TTree *newtrkTree = new TTree("Track2","Track Tree 2");
   newtrkTree->SetMaxTreeSize(MAXTREESIZE);
   newtrkTree->Branch("ntrk",&ntrk,"ntrk/I");
   newtrkTree->Branch("pt",trkpt,"pt[ntrk]/F");
   newtrkTree->Branch("eta",trketa,"eta[ntrk]/F");

   Int_t d_phi;

   vector<float> jetpt;
   vector<float>::iterator iter_jetpt;

   Int_t numjt;
   // event loop start
   //   if(nevt == -1) nevt = pfTree->GetEntries();   

   Int_t numevt = 0;
   if(nevt == -1) nevt = t->GetEntries(); 
   for (int iev=0; iev<500000; iev++) {
     if(iev%10000==0)
     { 
       cout << ">>>>> EVENT " << iev << " / " << t->GetEntries() <<  " ("<<(int)(100.*iev/t->GetEntries()) << "%)" << endl; 
     }
     HltTree->GetEntry(iev);
     jetpt.clear();

     if ( (trig != "" ) && (!trigBit) ) {
       continue;
     }
     
     t->GetEntry(iev);
     
     // trigger selection 

     ///// Call the values /////
     
     event.clear();
     event.run = Run;
     event.lumi = LumiBlock;
     event.event = Event; 
     event.vz = -99;
     event.hfsum = hiHF; 
     event.hfplus = hiHFplus; 
     event.hfminus = hiHFminus; 
     
     dj.clear();

     Track.clear();
     Track.nTrack = nTrk;

     //Track2->clear();
     ntrk = nTrk;

     for(Int_t c = 0; c != nTrk; c++)
     {
       trkpt[c] = trkPt[c];
       trketa[c] = trkEta[c];
       if(trkEta[c] >= -2.5 && trkEta[c] < -2.0)
       {
         Track.nTrketam2to2p5 = nTrk;
       }
       else if(trkEta[c] >= -2 && trkEta[c] < -1.5)
       {
         Track.nTrketam1p5to2 = nTrk;
       }
       else if(trkEta[c] >= -1.5 && trkEta[c] < -1.0)
       {
         Track.nTrketam1to1p5 = nTrk;
       }
       else if(trkEta[c] >= -1.0 && trkEta[c] < -0.5)   
       {
         Track.nTrketam0p5to1 = nTrk;     
       }
       else if(trkEta[c] >= -0.5 && trkEta[c] < 0.0)
       {
         Track.nTrketam0to0p5 = nTrk;
       }
       else if(trkEta[c] >= 0.0 && trkEta[c] < 0.5)
       {
         Track.nTrketa0to0p5 = nTrk;
       }
       else if(trkEta[c] >= 0.5 && trkEta[c] < 1.0)
       {
         Track.nTrketa0p5to1 = nTrk;
       }
       else if(trkEta[c] >= 1.0 && trkEta[c] < 1.5)
       {
         Track.nTrketa1to1p5 = nTrk;
       }
       else if(trkEta[c] >= 1.5 && trkEta[c] < 2.0)
       {
         Track.nTrketa1p5to2 = nTrk;
       }
       else if(trkEta[c] >= 2.0 && trkEta[c] < 2.5)
       {
         Track.nTrketa2to2p5 = nTrk;
       }
       else
       {
         continue;
       }
     } 

    Int_t num[200] = {0};
    numjt = 0;
    int ljInd = -1;
    //if(nref >= 2)
    {
      for(Int_t a = 0; a != nref; a++)
      {
        TLorentzVector jt;
        jt.SetPtEtaPhiM( jtpt[a], jteta[a], jtphi[a], jtm[a] );
        //if(TMath::Abs(jteta[a]) < 2.4 && jtpt[a] > minjPt)
        {
          num[a] = 1;
          numjt += 1;
          jetpt.push_back(jtpt[a]);
        }
        //else
        //{
        //  num[a] = 0;
        //}
        //jetpt.push_back(jtpt[a]);
      }
      //if(numjt == 2)
      {
      sort(jetpt.begin(), jetpt.end(), arrange);
      //cout << nref << " " << jetpt.size() << " " << jetpt[0] << " " << jetpt[1] << " " << jetpt[2] << " " << jetpt[3] << endl;
      for(Int_t b = 0; b != nref; b++)
      {
        if(jtpt[b] == jetpt[0])
        {
          jtpt2[0] = jtpt[b];
          jteta2[0] = jteta[b];
          jtphi2[0] = jtphi[b];
          //jtenergy2[0] = jtenergy[b];
          jtm2[0] = jtm[b];
        }
        else if(jtpt[b] == jetpt[1])
        {
          jtpt2[1] = jtpt[b];
          jteta2[1] = jteta[b];
          jtphi2[1] = jtphi[b];
          //jtenergy2[1] = jtenergy[b];
          jtm2[1] = jtm[b];
        }
        else
        {
          continue;
        }
      }
      //cout << jtpt2[0] << " " << jtpt2[1] << endl;
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
      //e[0] = jtenergy2[0];
      mass[0] = jtm2[0];
      pt[1] = jtpt2[1];
      eta[1] = jteta2[1];
      phi[1] = jtphi2[1];
      //e[1] = jtenergy2[1];
      mass[1] = jtm2[1];
         
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
      //dj.e1 = e[0];
      dj.pt2 = pt[1];
      dj.eta2 = eta[1];
      dj.phi2 = phi[1];
      //dj.e2 = e[1];
    }
    //else
    //{
    //  continue;
    //}
    }
    //else
    //{
    // continue;
    //}
    eventTree->Fill();
    djTree->Fill();
    trkTree->Fill();   
    newtrkTree->Fill(); 
    jetpt.clear(); 
   } //end of event loop
   // *==*==*==*==*==*==* Output file  *==*==*==*==*==*==* //
//   cout << numevt << endl;
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
