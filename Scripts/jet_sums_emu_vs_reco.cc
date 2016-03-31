#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMultiGraph.h"

#include <iostream>
#include <sstream>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetFilterDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

void print_canvas(TCanvas * c1, string str1, string str2){
  printName1 = str1 + ".svg";
  printName2 = "svg " + str2;
  c1->Print(printName1.c_str(), printName2.c_str());

  printName1 = str1 + ".png";
  printName2 = "png " + str2;
  c1->Print(printName1.c_str(), printName2.c_str());

  printName1 = str1 + ".pdf";
  printName2 = "pdf " + str2;
  c1->Print(printName1.c_str(), printName2.c_str());
}

void format_dist(TH2F * hist, const char name [], const char title []){
  hist->SetName(name);
  hist->SetTitle(title);
  hist->SetBins(100,0.,500.,100,0.,500.);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetTitleOffset(1.55);
}

void main_func(Int_t setPU, Int_t passMuonFil){
  //string sample = "DM_ZtoLL";
  string sample = "WJetsToLNu";
  Int_t nevents=3.e5;

  Bool_t setLogZAxis = false;
  //Int_t setPU = 3; // 0 = ALL, 1 = Low, 2 = Med, 3 = High
  //Int_t passMuonFil = 2; // 0 = no Muon filter, 1 = Pass, 2 = Fail

  string setPUTitle [] = {"", ", Low PU", ", Med PU", ", High PU"};
  string setPUName [] = {"", "_LowPU", "_MedPU", "_HighPU"};
  string passMuonFilTitle [] = {"", " and >=2 Mu's (pT>=20 GeV, isLooseMuon==1.)", " and <2 Mu's (pT>=20 GeV, isLooseMuon==1.)"};
  string passMuonFilName [] = {"", "_PassMuFilter", "_FailMuFilter"};
  string addTitle = ", EtCorr, MHT recalc (muMult=0.0)"+passMuonFilTitle[passMuonFil] + setPUTitle[setPU];
  //" (muMult=0.0)"
  string addName  = "_EtCorr_mhtRecalc_muMult0"+passMuonFilName[passMuonFil] + setPUName[setPU];
  //"_muMult0"

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.13);

  // make trees
  // /L1Ntuple.root
  // const char * filename = "/afs/cern.ch/work/s/sbreeze/public/L1Ntuple.root";
  // const char * filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/L1Ntuple.root";
  //const char * filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/DM_MonoZToLL/L1Ntuple_DM_ZtoLL.root";
  const char * filename = "/afs/cern.ch/work/s/sbreeze/public/jets_and_sums/WJetsToLNu/L1NTuple_WJetsToLNu.root";
  // const char * filename = "srm://gfe02.grid.hep.ph.ic.ac.uk:8443/srm/managerv2?SFN=/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sbreeze/ttHTobb_M125_13TeV_powheg_pythia8/crab_l1Ntuples-mc_ttHtobb_RunIIFall15DR76-EMU_RAW2DIGI/160311_131852/0000/L1Ntuple_5.root";
  TFile * file = new TFile(filename);
  //TFile * file = new TFile("/afs/cern.ch/user/a/ashtipli/private/work/run-data.root");
  
  file->ls();
  TTree * treeL1Up  = (TTree*) file->Get("l1JetRecoTree/JetRecoTree");
  if (! treeL1Up){
    cout << "ERROR: could not open jet reco tree\n";
    return;
  }
  TTree * treeL1UpEmu  = (TTree*) file->Get("l1UpgradeEmuTree/L1UpgradeTree");
  if (! treeL1UpEmu){
    cout << "ERROR: could not open tree\n";
    return;
  }
  TTree * treeL1Reco = (TTree*) file->Get("l1RecoTree/RecoTree");
  if (! treeL1Reco){
    cout << "ERROR: could not open reco tree\n";
    return;
  }
  TTree * treeL1RecoMetFilters = (TTree*) file->Get("l1MetFilterRecoTree/MetFilterRecoTree");
  if (! treeL1RecoMetFilters){
    cout << "ERROR: could not open reco met filter tree\n";
    return;
  }
  TTree * treeL1RecoMuons = (TTree*) file->Get("l1MuonRecoTree/Muon2RecoTree");
  if (! treeL1RecoMuons){
    cout << "ERROR: could not open reco muon tree\n";
    return;
  }
  TTree * treeL1RecoJets = (TTree*) file->Get("l1JetRecoTree/JetRecoTree");

  // set branch addresses
  // L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisL1UpgradeDataFormat     *upgradeEmu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMetDataFormat       *reco_       = new L1Analysis::L1AnalysisRecoMetDataFormat();
  L1Analysis::L1AnalysisRecoVertexDataFormat    *vertex_     = new L1Analysis::L1AnalysisRecoVertexDataFormat();
  L1Analysis::L1AnalysisRecoMetFilterDataFormat *metFilters_ = new L1Analysis::L1AnalysisRecoMetFilterDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat     *muons_      = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  L1Analysis::L1AnalysisRecoJetDataFormat       *jets_       = new L1Analysis::L1AnalysisRecoJetDataFormat();
  
  //treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);
  treeL1Up->SetBranchAddress("Sums", &reco_);
  treeL1UpEmu->SetBranchAddress("L1Upgrade", &upgradeEmu_);
  treeL1Reco->SetBranchAddress("Vertex", &vertex_);
  treeL1RecoMetFilters->SetBranchAddress("MetFilters", &metFilters_);
  treeL1RecoMuons->SetBranchAddress("Muon", &muons_);
  treeL1RecoJets->SetBranchAddress("Jet", &jets_);

  TH2F * mhtHists = new TH2F();
  stringstream ss1, ss2;
  ss1 << fixed << setprecision(0) << "mhtEmuRecoCorr" << addName;
  ss2 << fixed << setprecision(0) << "DM_Ztoll, #sqrt{s} = 13 TeV" << addTitle << ";MHT emu (GeV);MHT reco(GeV)";
  format_dist(mhtHists, ss1.str().c_str(), ss2.str().c_str());
  ss1.str("");
  ss2.str("");
  
  TH2F * metHists = new TH2F();
  ss1 << "metHist" << addName;
  ss2 << "DM_Ztoll, #sqrt{s} = 13 TeV" << addTitle << ";MET emu (GeV);MET reco (GeV)";
  format_dist(metHists, ss1.str().c_str(), ss2.str().c_str());
  ss1.str("");
  ss2.str("");
  

  // get entries
  Long64_t nentries = treeL1UpEmu->GetEntriesFast();
  if (nevents>nentries) nevents=nentries;

  std::cout << "Running over " << nevents << ", nentries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nevents;jentry++){

    if(((jentry+1)%1000)==0){
      std::cout << "Done " << jentry+1  << " events...\r" << std::flush;
    }

    treeL1Up->GetEntry(jentry);
    treeL1UpEmu->GetEntry(jentry);
    treeL1Reco->GetEntry(jentry);
    treeL1RecoMetFilters->GetEntry(jentry);
    treeL1RecoMuons->GetEntry(jentry);
    treeL1RecoJets->GetEntry(jentry);

    // Muon variables
    Int_t isPassMuonCount = 0;
    Bool_t isPassMuonFilter = false;

    // Muon Loop
    for(unsigned int iter=0; iter<muons_->nMuons; ++iter){
      double pt = muons_->pt[iter];
      double isMedMu = muons_->isMediumMuon[iter];
      if( pt >= 20. && isMedMu == 1. ) ++isPassMuonCount;
    }
    if( isPassMuonCount >= 2 ) isPassMuonFilter = true;

    // Jet Loop
    double mht_ex = 0.0;
    double mht_ey = 0.0;
    for(unsigned int iter=0; iter<jets_->nJets; ++iter){
      double et = jets_->etCorr[iter];
      double phi = jets_->phi[iter];
      double muMult = jets_->muMult[iter];
      if(muMult == 0.0){
        mht_ex += et*cos(phi);
        mht_ey += et*sin(phi);
      }
    }
    double mht_recalc = sqrt(mht_ex*mht_ex + mht_ey*mht_ey);
	
    // Sums variables 
    double etSumEmu  = -1.0;
    double htSumEmu  = -1.0;
    double metSumEmu = -1.0;
    double mhtSumEmu = -1.0;
	
    // Sums Loop
    for(unsigned int iter=0; iter<upgradeEmu_->nSums; ++iter){
      double et = upgradeEmu_->sumEt[iter];
      if (upgradeEmu_->sumType[iter] == L1Analysis::kTotalEt)   etSumEmu  = et;
      if (upgradeEmu_->sumType[iter] == L1Analysis::kTotalHt)   htSumEmu  = et;
      if (upgradeEmu_->sumType[iter] == L1Analysis::kMissingEt) metSumEmu = et;
      if (upgradeEmu_->sumType[iter] == L1Analysis::kMissingHt) mhtSumEmu = et;
    }

    Int_t nVert = vertex_->nVtx;
    Bool_t PUsetting [] = {true, nVert<=15, nVert>15 && nVert<=30, nVert>30};
    Bool_t failPassMuFil[] = {true, isPassMuonFilter, !isPassMuonFilter};

    if( PUsetting[setPU] ) mhtHists->Fill(mhtSumEmu, mht_recalc);//reco_->mHt);
    if( metFilters_->hbheNoiseFilter == 1.0 && failPassMuFil[passMuonFil] && PUsetting[setPU] ) metHists->Fill(metSumEmu, reco_->met);
  }

  // Create canvases for MHT distributions
  TCanvas * t1 = new TCanvas("comparison_canvas1", "Comparison Canvas", 800, 750);
  mhtHists->DrawClone("colz");
  string printName1 = "mht_emu_reco_corr" + addName;
  string printName2 = "Title:MHT Emu vs Reco Comparison" + addTitle;
  if( setLogZAxis ) t1->SetLogz();
  print_canvas(t1, printName1, printName2);

  // Create canvases for MET distributions
  TCanvas * t2 = new TCanvas("Comparison_canvas2", "Comparison Canvas", 800, 750);
  metHists->DrawClone("colz");
  printName1 = "met_emu_cut_reco_corr" + addName;
  printName2 = "Title:MET Emu vs Reco Comparison" + addTitle;
  if( setLogZAxis ) t2->SetLogz();
  print_canvas(t2, printName1, printName2);
}

void jet_sums_emu_vs_reco_2(){
  //for(Int_t passMuonFil=0; passMuonFil<3; ++passMuonFil){
  Int_t passMuonFil=0;
    for(Int_t setPU=0; setPU<4; ++setPU){
      main_func(setPU, passMuonFil);
    }
  //}
}
