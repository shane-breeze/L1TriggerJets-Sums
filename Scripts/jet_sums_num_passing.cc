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

void format_dist(TH1F & hist, const char name [], const char title [], Color_t col){
  hist.SetName(name);
  hist.SetTitle(title);
  hist.SetBins(200,0.,2000.);
  hist.GetXaxis()->SetTitleOffset(1.3);
  hist.GetYaxis()->SetTitleOffset(1.55);
  hist.SetLineColor(col);
}

void format_efficiency(TGraphAsymmErrors * graph, const char name [], const char title [], Color_t col){
  graph->SetName(name);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleOffset(1.3);
  graph->GetYaxis()->SetTitleOffset(1.55);
  graph->SetMaximum(1.1);
  graph->SetMinimum(0.0);
  graph->SetLineColor(col);
}

void jet_sums_turnon(){
  Int_t nevents=5.e4;

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
  const char * filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/DM_MonoZToLL/L1Ntuple_DM_ZtoLL.root";
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

  // set branch addresses
  // L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisL1UpgradeDataFormat     *upgradeEmu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMetDataFormat       *reco_       = new L1Analysis::L1AnalysisRecoMetDataFormat();
  L1Analysis::L1AnalysisRecoVertexDataFormat    *vertex_     = new L1Analysis::L1AnalysisRecoVertexDataFormat();
  L1Analysis::L1AnalysisRecoMetFilterDataFormat *metFilters_ = new L1Analysis::L1AnalysisRecoMetFilterDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat     *muons_      = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  
  //treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);
  treeL1Up->SetBranchAddress("Sums", &reco_);
  treeL1UpEmu->SetBranchAddress("L1Upgrade", &upgradeEmu_);
  treeL1Reco->SetBranchAddress("Vertex", &vertex_);
  treeL1RecoMetFilters->SetBranchAddress("MetFilters", &metFilters_);
  treeL1RecoMuons->SetBranchAddress("Muon", &muons_);

  Int_t NmhtCuts = 20;
  Double_t mhtCuts [] = {0.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.};
  TH1F * mhtHists = new TH1F[NmhtCuts];
  stringstream ss1, ss2;
  for (Int_t it=0; it<NmhtCuts; ++it){
    ss1 << fixed << setprecision(0) << "mhtHist" << mhtCuts[it] << "_MuFilter";
    ss2 << fixed << setprecision(0) << "DM_Ztoll, emu_mht cuts, #sqrt{s} = 13 TeV, Mu Filter;MHT reco (GeV); Number of events";
    format_dist(mhtHists[it], ss1.str().c_str(), ss2.str().c_str(), kBlue);
    ss1.str("");
    ss2.str("");
  }

  Int_t NmetCuts = 20;
  Double_t metCuts [] = {0.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.};
  TH1F * metHists = new TH1F[NmetCuts];
  for (Int_t it=0; it<NmetCuts; ++it){
    ss1 << "metHist" << metCuts[it] << "_MuFilter";
    ss2 << "DM_Ztoll, emu_met cuts, #sqrt{s} = 13 TeV, Mu Filter;MET reco (GeV); Number of events";
    format_dist(metHists[it], ss1.str().c_str(), ss2.str().c_str(), kBlue);
    ss1.str("");
    ss2.str("");
  }

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
    for(Int_t it=0; it<NmhtCuts; ++it){
      if( mhtSumEmu >= mhtCuts[it] ) mhtHists[it].Fill(reco_->mHt);
    }
    for(Int_t it=0; it<NmetCuts; ++it){
      if( metSumEmu >= metCuts[it] && metFilters_->hbheNoiseFilter == 1.0 && isPassMuonFilter ) metHists[it].Fill(reco_->met);
    }
  }

  // Create canvases for MHT distributions
  for(Int_t it=0; it<NmhtCuts; ++it){
    ss1 << "emu_mht >= " << mhtCuts[it] << " GeV";
    if( it == 0 ) ss1.str("Total");
    ss1.str("");
  }
  string printName1 = "mht_emu_cut_reco_distribution_muonFilter";
  string printName2 = "Title:MHT Reco Distribution with Emulator Cuts and Muon Filter";

  // Create canvases for MET distributions
  for(Int_t it=0; it<NmetCuts; ++it){
    ss1 << "emu_met >= " << metCuts[it] << " GeV";
    if( it == 0 ) ss1.str("Total");
    ss1.str("");
  }
  printName1 = "met_emu_cut_reco_distribution_muonFilter";
  printName2 = "Title:MET Reco Distribution with Emulator Cuts and Muon Filter";
}
