// ROOT classes
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

// C++ classes
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

// L1Trigger classes (data structures to be exact)
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoVertexDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetFilterDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"


// Print the canvas to svg, png and pdf files
void print_canvas(TCanvas * c1, string str1, string str2, string saveDir){
  printName1 = saveDir + "/" + str1 + ".svg";
  printName2 = "svg " + str2;
  c1->Print(printName1.c_str(), printName2.c_str());

  printName1 = saveDir + "/" + str1 + ".png";
  printName2 = "png " + str2;
  c1->Print(printName1.c_str(), printName2.c_str());

  printName1 = saveDir + "/" + str1 + ".pdf";
  printName2 = "pdf " + str2;
  c1->Print(printName1.c_str(), printName2.c_str());
}

// Format histograms
void format_dist(TH2F * hist, const char name [], const char title []){
  hist->SetName(name);
  hist->SetTitle(title);
  hist->SetBins(50,0.,500.,50,0.,500.);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetTitleOffset(1.55);
}

// Format the TGraphAsymmErrors (the turn-on curve graph)
void format_efficiency(TGraphAsymmErrors * graph, const char name [], const char title [], Color_t col){
  graph->SetName(name);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleOffset(1.3);
  graph->GetYaxis()->SetTitleOffset(1.55);
  graph->SetMaximum(1.1);
  graph->SetMinimum(0.0);
  graph->SetLineColor(col);
}

// Muon Filter
Bool_t muonFilter(L1Analysis::L1AnalysisRecoMuon2DataFormat * muons_, Double_t ptCut, Int_t nMuonsToPassFilter){
  Int_t isPassMuonCount = 0;
  Bool_t isPassMuonFilter = false;
  for(unsigned int iter=0; iter<muons_->nMuons; ++iter){
    double pt = muons_->pt[iter];
    double isLooseMu = muons_->isLooseMuon[iter];
    if( pt >= ptCut && isLooseMu == 1. ) ++isPassMuonCount;
  }
  if( isPassMuonCount >= nMuonsToPassFilter ) isPassMuonFilter = true;
  return isPassMuonFilter;
}

// Recalculate MHT using jets
Double_t recalculateMHT(L1Analysis::L1AnalysisRecoJetDataFormat * jets_){
  Double_t mht_ex = 0.0;
  Double_t mht_ey = 0.0;
  for(unsigned int iter=0; iter<jets_->nJets; ++iter){
    Double_t et = jets_->etCorr[iter];
    Double_t phi = jets_->phi[iter];
    Double_t muMult = jets_->muMult[iter];
    if(muMult == 0.0){
      mht_ex += et*cos(phi);
      mht_ey += et*sin(phi);
    }
  }
  return sqrt(mht_ex*mht_ex + mht_ey*mht_ey);
}

// Recalculate MET using Calo Towers
Double_t recalculateMET(L1Analysis::L1AnalysisL1CaloTowerDataFormat * caloTowers_, Int_t iEtaMax){
  Double_t metX = 0.0;
  Double_t metY = 0.0;
  for(unsigned int jTower=0; jTower<caloTowers_->nTower; ++jTower){
    Int_t ieta = caloTowers_->ieta[jTower];
    Int_t iphi = caloTowers_->iphi[jTower];
    Int_t iet = caloTowers_->iet[jTower];
    if( abs(ieta) < iEtaMax){
      Double_t phi = (Double_t)iphi * TMath::Pi()/36.;
      Double_t et = 0.5 * (Double_t)iet;
      metX += et * TMath::Cos(phi);
      metY += et * TMath::Sin(phi);
    }
  }
  return TMath::Sqrt(metX*metX + metY*metY);
}

// Where most things happen
void main_func(string puCat, string muonFilterCat, Bool_t doMHT){
  // Set some variables. Change as required
  // string sample = "DM_ZToLL";
  string sample = "WprimeToENu";
  // string sample = "ttHTobb";
  //string sample = "WprimeToMuNu";
  string location = "/afs/cern.ch/work/s/sbreeze/CMSSW_8_0_2/src/"+sample+"/Plots";
  string filename = "/afs/cern.ch/work/s/sbreeze/CMSSW_8_0_2/src/"+sample+"/Ntuples/L1NTuple_"+sample+".root";
  Int_t nevents=1.e6;
  Int_t nMuonsToPassFilter = 1;
  Int_t iEtaMax = 29;
  Bool_t setLogXAxis = false;
  Double_t ptCut = 10.0;

  // Create maps for the title and names of histograms and such
  // For PU:
  std::map<string, string> puTitleMap;

  puTitleMap["All"]  = "";
  puTitleMap["Low"]  = ", Low PU";
  puTitleMap["Med"]  = ", Med PU";
  puTitleMap["High"] = ", High PU";

  std::map<string, string> puNameMap;
  puNameMap["All"]  = "";
  puNameMap["Low"]  = "_LowPU";
  puNameMap["Med"]  = "_MedPU";
  puNameMap["High"] = "_HighPU";

  // For Muon Filter:
  stringstream sst1, sst2, ssn1, ssn2;
  sst1 << fixed << setprecision(0) << " and >=" << nMuonsToPassFilter << " Mu (pT>=" << ptCut << " GeV, isLooseMuon==1.)";
  sst2 << fixed << setprecision(0) << " and <" << nMuonsToPassFilter << " Mu (pT>="  << ptCut << " GeV, isLooseMuon==1.)";
  ssn1 << fixed << setprecision(0) << "_PassMuFilterLoose" << ptCut;
  ssn2 << fixed << setprecision(0) << "_FailMuFilterLoose" << ptCut;
  std::map<string, string> muonFilterTitleMap;
  muonFilterTitleMap["None"] = "";
  muonFilterTitleMap["Pass"] = sst1.str(); 
  muonFilterTitleMap["Fail"] = sst2.str();
  sst1.str("");
  sst2.str("");

  std::map<string, string> muonFilterNameMap;
  muonFilterNameMap["None"] = "";
  muonFilterNameMap["Pass"] = ssn1.str();
  muonFilterNameMap["Fail"] = ssn2.str();
  ssn1.str("");
  ssn2.str("");

  // Strings to be added on the end of histogram names and titles:
  string addTitle = muonFilterTitleMap[muonFilterCat] + puTitleMap[puCat];
  //string addTitle = ", EtCorr, MHT recalc (muMult=0.0)"+...
  string addName  = muonFilterNameMap[muonFilterCat] + puNameMap[puCat];
  //string addName  = "_etCorr_mhtRecalc_muMult0"+...


  // make trees
  // const char * filename = "/afs/cern.ch/work/s/sbreeze/public/L1Ntuple.root";
  // const char * filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/L1Ntuple.root";
  // const string filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/DM_MonoZToLL/L1Ntuple_DM_ZtoLL.root";
  // const char * filename = "/afs/cern.ch/work/s/sbreeze/public/jets_and_sums/WJetsToLNu/L1NTuple_WJetsToLNu.root";
  // const char * filename = "srm://gfe02.grid.hep.ph.ic.ac.uk:8443/srm/managerv2?SFN=/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sbreeze/ttHTobb_M125_13TeV_powheg_pythia8/crab_l1Ntuples-mc_ttHtobb_RunIIFall15DR76-EMU_RAW2DIGI/160311_131852/0000/L1Ntuple_5.root";
  TFile * file = new TFile(filename.c_str());
  //file->ls();
  
  // Map for Tree Category -> Tree path
  map<string, string> treeCat_treePath_map;
  treeCat_treePath_map["Reco Jet"]        = "l1JetRecoTree/JetRecoTree";
  treeCat_treePath_map["Emu Upgrade"]     = "l1UpgradeEmuTree/L1UpgradeTree";
  treeCat_treePath_map["Reco"]            = "l1RecoTree/RecoTree";
  treeCat_treePath_map["Reco MET Filter"] = "l1MetFilterRecoTree/MetFilterRecoTree";
  treeCat_treePath_map["Reco Muon"]       = "l1MuonRecoTree/Muon2RecoTree";
  treeCat_treePath_map["L1 Calo Tower"]   = "l1CaloTowerEmuTree/L1CaloTowerTree";

  // Get the Trees
  map<string, TTree*> treeCat_tree_map;
  for(auto it = treeCat_treePath_map.begin(); it != treeCat_treePath_map.end(); ++it){
    treeCat_tree_map[it->first] = (TTree*) file->Get(it->second.c_str());
    if( !treeCat_tree_map[it->first] ){
      cout << "ERROR: could not open path to " << it->first << ": " << it->second << endl;
      return;
    }
  }

  // set branch addresses
  // L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisL1UpgradeDataFormat      *upgradeEmu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMetDataFormat        *reco_       = new L1Analysis::L1AnalysisRecoMetDataFormat();
  L1Analysis::L1AnalysisRecoVertexDataFormat     *vertex_     = new L1Analysis::L1AnalysisRecoVertexDataFormat();
  L1Analysis::L1AnalysisRecoMetFilterDataFormat  *metFilters_ = new L1Analysis::L1AnalysisRecoMetFilterDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat      *muons_      = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  L1Analysis::L1AnalysisRecoJetDataFormat        *jets_       = new L1Analysis::L1AnalysisRecoJetDataFormat();
  L1Analysis::L1AnalysisL1CaloTowerDataFormat    *caloTowers_ = new L1Analysis::L1AnalysisL1CaloTowerDataFormat();
  
  //treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);
  treeCat_tree_map["Reco Jet"       ]->SetBranchAddress("Sums"       , &reco_      );
  treeCat_tree_map["Emu Upgrade"    ]->SetBranchAddress("L1Upgrade"  , &upgradeEmu_);
  treeCat_tree_map["Reco"           ]->SetBranchAddress("Vertex"     , &vertex_    );
  treeCat_tree_map["Reco MET Filter"]->SetBranchAddress("MetFilters" , &metFilters_);
  treeCat_tree_map["Reco Muon"      ]->SetBranchAddress("Muon"       , &muons_     );
  treeCat_tree_map["Reco Jet"       ]->SetBranchAddress("Jet"        , &jets_      );
  treeCat_tree_map["L1 Calo Tower"  ]->SetBranchAddress("L1CaloTower", &caloTowers_);

  // Create the histograms for mht
  TH2F * metHists = new TH2F();
  stringstream ss1, ss2;
  ss1 << fixed << setprecision(0) << "metHist" << addName;
  ss2 << fixed << setprecision(0) << sample << ", #sqrt{s} = 13 TeV" << addTitle << ";L1 MET (GeV); Recalculated L1 MET (GeV)";
  format_dist(metHists, ss1.str().c_str(), ss2.str().c_str());
  ss1.str("");
  ss2.str("");

  // get entries
  Long64_t nentries = treeCat_tree_map["Reco Jet"]->GetEntriesFast();
  if( nevents > nentries) nevents=nentries;

  // Which data sample and settings are used (printed out for you)
  std::cout << "Running over the dataset " << sample << addTitle << std::endl;
  std::cout << "Looping over " << nevents << ", out of a total of " << nentries << std::endl;
  for (int jentry=0; jentry<nevents;jentry++){
    // Print out every 1000th event
    if(((jentry+1)%1000)==0){
      std::cout << "Done " << jentry+1  << " events...\r" << std::flush;
    }

    // Fill the trees with the jentry-th event 
    for(auto it = treeCat_tree_map.begin(); it != treeCat_tree_map.end(); ++it){
      it->second->GetEntry(jentry);
    }

    Bool_t isPassMuonFilter = muonFilter(muons_, ptCut, nMuonsToPassFilter);
    Double_t mht_recalc = recalculateMHT(jets_);
    Double_t metRecalc = recalculateMET(caloTowers_, iEtaMax);
	
    // Sums Loop
    double etSumEmu  = -1.0;
    double htSumEmu  = -1.0;
    double metSumEmu = -1.0;
    double mhtSumEmu = -1.0;
    for(unsigned int iter=0; iter<upgradeEmu_->nSums; ++iter){
      double et = upgradeEmu_->sumEt[iter];
      if (upgradeEmu_->sumType[iter] == L1Analysis::kTotalEt)   etSumEmu  = et;
      if (upgradeEmu_->sumType[iter] == L1Analysis::kTotalHt)   htSumEmu  = et;
      if (upgradeEmu_->sumType[iter] == L1Analysis::kMissingEt) metSumEmu = et;
      if (upgradeEmu_->sumType[iter] == L1Analysis::kMissingHt) mhtSumEmu = et;
    }

    // PU splitting
    Int_t nVert = vertex_->nVtx;
    map<string, Bool_t> puCat_puBool_map;
    puCat_puBool_map["All" ] = true;
    puCat_puBool_map["Low" ] = nVert <= 15;
    puCat_puBool_map["Med" ] = nVert > 15 && nVert <= 30;
    puCat_puBool_map["High"] = nVert > 30;

    // muon Filter splitting
    map<string, Bool_t> muFilter_muFilterBool_map;
    muFilter_muFilterBool_map["None"] = true;
    muFilter_muFilterBool_map["Pass"] = isPassMuonFilter;
    muFilter_muFilterBool_map["Fail"] = !isPassMuonFilter;

    // Fill the histograms
    metHists->Fill(metSumEmu, metRecalc);//mht_recalc
  }

  // Create canvases for MHT distributions
  ss1 << "Dist_canvas_mht" << addName;
  TCanvas * t1 = new TCanvas(ss1.str().c_str(), "Distributions Canvas", 800, 750);
  ss1.str("");
  metHists->DrawClone("colz");
  string printName1 = "L1MetvsRecalcTowerMetDistribution" + addName;
  string printName2 = "Title:L1 MET vs Recalculated Tower MET" + addTitle;
  if( setLogXAxis ) t1->SetLogx();
  print_canvas(t1, printName1, printName2, location);
}

// Function called at the start (so I conveniently placed it at the end)
void offlineMetVsTriggerTowerMet(){
  // Global style setting for ROOT
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPadRightMargin(0.13);

  // Set which combinations of muon filter and PU category to loop through
  // - Muon filter can be off for MHT (MHT doesn't include the muon in it's calculation
  std::map<string, vector<string> > muFilter_puCat_selection;
  vector<string> noneMuFilterVector = {"All"};//, "Low", "Med", "High"};
  //vector<string> passMuFilterVector = {"All", "Low", "Med", "High"};
  //vector<string> failMuFilterVector = {"All", "Low", "Med", "High"};
  muFilter_puCat_selection["None"] = noneMuFilterVector;
  //muFilter_puCat_selection["Pass"] = passMuFilterVector;
  //muFilter_puCat_selection["Fail"] = failMuFilterVector;

  for(auto itMu=muFilter_puCat_selection.begin(); itMu != muFilter_puCat_selection.end(); ++itMu){
    Bool_t doMHT = true; // Muon filter is needed for MET only.
    if( itMu->first != "None" ) doMHT = false;
    for(auto itPu=itMu->second.begin(); itPu != itMu->second.end(); ++itPu){
      main_func(*itPu, itMu->first, doMHT);
    }
  }
}
