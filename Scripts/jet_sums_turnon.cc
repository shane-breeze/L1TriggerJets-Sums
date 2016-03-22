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

// Format histograms
void format_dist(TH1F & hist, const char name [], const char title [], Color_t col){
  hist.SetName(name);
  hist.SetTitle(title);
  hist.SetBins(50,0.,500.);
  hist.GetXaxis()->SetTitleOffset(1.3);
  hist.GetYaxis()->SetTitleOffset(1.55);
  hist.SetLineColor(col);
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

// Where most things happen
void main_func(string puCat, string muonFilterCat, Bool_t doMHT){
  // Set some variables. Change as required
  string sample = "DM_ZtoLL";
  //string sample = "WJetsToLNu";
  Int_t nevents=3.e5;
  Bool_t setLogXAxis = false;

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
  std::map<string, string> muonFilterTitleMap;
  muonFilterTitleMap["None"] = "";
  muonFilterTitleMap["Pass"] = " and >=2 Mu's (pT>=20 GeV, isLooseMuon==1.)";
  muonFilterTitleMap["Fail"] = " and <2 Mu's (pT>=20 GeV, isLooseMuon==1.)";

  std::map<string, string> muonFilterNameMap;
  muonFilterNameMap["None"] = "";
  muonFilterNameMap["Pass"] = "_PassMuFilter";
  muonFilterNameMap["Fail"] = "_FailMuFilter";

  // Strings to be added on the end of histogram names and titles:
  string addTitle = "" + muonFilterTitleMap[muonFilterCat] + puTitleMap[puCat];
  //string addTitle = ", EtCorr, MHT recalc (muMult=0.0)"+...
  string addName  = "" + muonFilterNameMap[muonFilterCat] + puNameMap[puCat];
  //string addName  = "_etCorr_mhtRecalc_muMult0"+...


  // make trees
  // const char * filename = "/afs/cern.ch/work/s/sbreeze/public/L1Ntuple.root";
  // const char * filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/L1Ntuple.root";
  const string filename = "/afs/cern.ch/user/s/sbreeze/public/jets_and_sums/DM_MonoZToLL/L1Ntuple_DM_ZtoLL.root";
  // const char * filename = "/afs/cern.ch/work/s/sbreeze/public/jets_and_sums/WJetsToLNu/L1NTuple_WJetsToLNu.root";
  // const char * filename = "srm://gfe02.grid.hep.ph.ic.ac.uk:8443/srm/managerv2?SFN=/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/sbreeze/ttHTobb_M125_13TeV_powheg_pythia8/crab_l1Ntuples-mc_ttHtobb_RunIIFall15DR76-EMU_RAW2DIGI/160311_131852/0000/L1Ntuple_5.root";
  TFile * file = new TFile(filename.c_str());
  file->ls();
  //TFile * file = new TFile("/afs/cern.ch/user/a/ashtipli/private/work/run-data.root");
  
  // Map for Tree Category -> Tree path
  map<string, string> treeCat_treePath_map;
  treeCat_treePath_map["Reco Jet"]        = "l1JetRecoTree/JetRecoTree";
  treeCat_treePath_map["Emu Upgrade"]     = "l1UpgradeEmuTree/L1UpgradeTree";
  treeCat_treePath_map["Reco"]            = "l1RecoTree/RecoTree";
  treeCat_treePath_map["Reco MET Filter"] = "l1MetFilterRecoTree/MetFilterRecoTree";
  treeCat_treePath_map["Reco Muon"]       = "l1MuonRecoTree/Muon2RecoTree";

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
  L1Analysis::L1AnalysisL1UpgradeDataFormat     *upgradeEmu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMetDataFormat       *reco_       = new L1Analysis::L1AnalysisRecoMetDataFormat();
  L1Analysis::L1AnalysisRecoVertexDataFormat    *vertex_     = new L1Analysis::L1AnalysisRecoVertexDataFormat();
  L1Analysis::L1AnalysisRecoMetFilterDataFormat *metFilters_ = new L1Analysis::L1AnalysisRecoMetFilterDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat     *muons_      = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  L1Analysis::L1AnalysisRecoJetDataFormat       *jets_       = new L1Analysis::L1AnalysisRecoJetDataFormat();
  
  //treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);
  treeCat_tree_map["Reco Jet"       ]->SetBranchAddress("Sums"      , &reco_      );
  treeCat_tree_map["Emu Upgrade"    ]->SetBranchAddress("L1Upgrade" , &upgradeEmu_);
  treeCat_tree_map["Reco"           ]->SetBranchAddress("Vertex"    , &vertex_    );
  treeCat_tree_map["Reco MET Filter"]->SetBranchAddress("MetFilters", &metFilters_);
  treeCat_tree_map["Reco Muon"      ]->SetBranchAddress("Muon"      , &muons_     );
  treeCat_tree_map["Reco Jet"       ]->SetBranchAddress("Jet"       , &jets_      );


  // The mhtCuts used for the turn-ons (and the colours used)
  Int_t NmhtCuts = 20;
  Double_t mhtCuts [] = {0.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.};
  Color_t myCol [] = {kRed,kRed-9,kOrange+7,kOrange,
      kYellow,kYellow-9,kSpring+10,kSpring-9,
      kGreen,kGreen-9,kTeal+7,kTeal-4,
      kCyan,kAzure+5,kBlue,kBlue-9,
      kViolet+5,kViolet,kMagenta,kPink};

  // Create the histograms for mht
  TH1F * mhtHists = new TH1F[NmhtCuts];
  stringstream ss1, ss2;
  for (Int_t it=0; it<NmhtCuts; ++it){
    ss1 << fixed << setprecision(0) << "mhtHist" << mhtCuts[it] << addName;
    ss2 << fixed << setprecision(0) << sample << ", emu_mht cuts, #sqrt{s} = 13 TeV" << addTitle << ";MHT reco (GeV); Number of events";
    format_dist(mhtHists[it], ss1.str().c_str(), ss2.str().c_str(), myCol[it]);
    ss1.str("");
    ss2.str("");
  }

  // The metCuts used for the turn-ons
  Int_t NmetCuts = 20;
  Double_t metCuts [] = {0.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.};

  // Create the histograms for met
  TH1F * metHists = new TH1F[NmetCuts];
  for (Int_t it=0; it<NmetCuts; ++it){
    ss1 << "metHist" << metCuts[it] << addName;
    ss2 << sample << ", emu_met cuts, #sqrt{s} = 13 TeV" << addTitle << ";MET reco (GeV); Number of events";
    format_dist(metHists[it], ss1.str().c_str(), ss2.str().c_str(), myCol[it]);
    ss1.str("");
    ss2.str("");
  }

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

    // Muon Loop
    Int_t isPassMuonCount = 0;
    Bool_t isPassMuonFilter = false;
    for(unsigned int iter=0; iter<muons_->nMuons; ++iter){
      double pt = muons_->pt[iter];
      double isLooseMu = muons_->isLooseMuon[iter];
      if( pt >= 20. && isLooseMu == 1. ) ++isPassMuonCount;
    }
    if( isPassMuonCount >= 2 ) isPassMuonFilter = true;

    // Jet Loop - i.e. recalculate MHT
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
    for(Int_t it=0; it<NmhtCuts; ++it){
      if( mhtSumEmu >= mhtCuts[it] && puCat_puBool_map[puCat] ) mhtHists[it].Fill(mht_recalc);//reco_->mHt);
    }
    for(Int_t it=0; it<NmetCuts; ++it){
      if( metSumEmu >= metCuts[it] && metFilters_->hbheNoiseFilter == 1.0 && puCat_puBool_map[puCat] && muFilter_muFilterBool_map[muonFilterCat] ) metHists[it].Fill(reco_->met);
    }
  }

  // Create canvases for MHT distributions
  ss1 << "Dist_canvas_mht" << addName;
  TCanvas * t1 = new TCanvas(ss1.str().c_str(), "Distributions Canvas", 800, 750);
  ss1.str("");
  TLegend * tleg1 = new TLegend(0.6,0.4,0.8,0.9);
  for(Int_t it=0; it<NmhtCuts; ++it){
    if( it==0 ) mhtHists[it].DrawClone("same");
    else mhtHists[it].DrawClone("same");
    ss1 << "emu_mht >= " << mhtCuts[it] << " GeV";
    if( it == 0 ) ss1.str("Total");
    tleg1->AddEntry(&mhtHists[it],ss1.str().c_str());
    ss1.str("");
  }
  tleg1->Draw();
  string printName1 = "mht_emu_cut_reco_distribution" + addName;
  string printName2 = "Title:MHT Reco Distribution with Emulator Cuts" + addTitle;
  if( setLogXAxis ) t1->SetLogx();
  t1->SetLogy();
  if( muonFilterCat == "None" && doMHT ) print_canvas(t1, printName1, printName2);

  // Create canvases for MET distributions
  ss1 << "Dist_canvas_met" << addName;
  TCanvas * t2 = new TCanvas(ss1.str().c_str(), "Distributions Canvas", 800, 750);
  ss1.str("");
  TLegend * tleg2 = new TLegend(0.6,0.4,0.8,0.9);
  for(Int_t it=0; it<NmetCuts; ++it){
    if( it==0 ) metHists[it].DrawClone();
    else metHists[it].DrawClone("same");
    ss1 << "emu_met >= " << metCuts[it] << " GeV";
    if( it == 0 ) ss1.str("Total");
    tleg2->AddEntry(&mhtHists[it],ss1.str().c_str());
    ss1.str("");
  }
  tleg2->Draw();
  printName1 = "met_emu_cut_reco_distribution" + addName;
  printName2 = "Title:MET Reco Distribution with Emulator Cuts" + addTitle;
  if( setLogXAxis ) t2->SetLogx();
  t2->SetLogy();
  print_canvas(t2, printName1, printName2);

  // Create the TGraphAsymmetricErrors for the efficiencies
  // MHT
  ss1 << "turn_on_mht" << addName;
  TCanvas * t3 = new TCanvas(ss1.str().c_str(),"Turn-on Curves MHT", 800, 750);
  ss1.str("");
  TLegend * tleg3 = new TLegend(0.78,0.4,0.98,0.85);
  for(Int_t it=1; it<NmhtCuts; ++it){
    TGraphAsymmErrors * tgae1 = new TGraphAsymmErrors();
    ss1 << "mhtEfficiency" << mhtCuts[it] << addName;
    ss2 << sample << ", emu_mht cuts, #sqrt{s} = 13 TeV" << addTitle << ";MHT reco (GeV);Efficiency";
    format_efficiency(tgae1, ss1.str().c_str(), ss2.str().c_str(), myCol[it]);

    const TH1F * hPass = (const TH1F*)mhtHists[it].Clone();
    const TH1F * hTotal = (const TH1F*)mhtHists[0].Clone();
    tgae1->Divide(hPass,hTotal);

    if( it==1 ) tgae1->DrawClone();
    else tgae1->DrawClone("same");
    ss1.str("");
    ss1 << "emu_mht >= " << mhtCuts[it] << " GeV";
    tleg3->AddEntry(tgae1,ss1.str().c_str());
    ss1.str("");
    ss2.str("");
  }
  tleg3->Draw();
  printName1 = sample + "_mht_efficiency" + addName;
  printName2 = "Title:" + sample + " MHT Efficiency" + addTitle;
  if( setLogXAxis ) t3->SetLogx();
  if( muonFilterCat == "None" && doMHT ) print_canvas(t3, printName1, printName2);

  // MET
  ss1 << "turn_on_met" << addName;
  TCanvas * t4 = new TCanvas(ss1.str().c_str(),"MET Turn-on Curves", 800, 750);
  ss1.str("");
  TLegend * tleg4 = new TLegend(0.78,0.4,0.98,0.85);
  for(Int_t it=1; it<NmetCuts; ++it){
    TGraphAsymmErrors * tgae1 = new TGraphAsymmErrors();
    ss1 << "metEfficiency" << mhtCuts[it] << addName;
    ss2 << sample << ", emu_met cuts, #sqrt{s} = 13 TeV" << addTitle << ";MET reco (GeV);Efficiency";
    format_efficiency(tgae1, ss1.str().c_str(), ss2.str().c_str(), myCol[it]);

    const TH1F * hPass = (const TH1F*)metHists[it].Clone();
    const TH1F * hTotal = (const TH1F*)metHists[0].Clone();
    tgae1->Divide(hPass,hTotal);

    if( it==1 ) tgae1->DrawClone();
    else tgae1->DrawClone("same");
    ss1.str("");
    ss1 << "emu_met >= " << metCuts[it] << " GeV";
    tleg4->AddEntry(tgae1,ss1.str().c_str());
    ss1.str("");
    ss2.str("");
  }
  tleg4->Draw();
  printName1 = sample + "_met_efficiency" + addName;
  printName2 = "Title:" + sample + " MET Efficiency" + addTitle;
  if( setLogXAxis ) t4->SetLogx();
  print_canvas(t4, printName1, printName2);
}

// Function called at the start (so I conveniently placed it at the end)
void jet_sums_turnon(){
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
  vector<string> noneMuFilterVector = {"All", "Low", "Med", "High"};
  vector<string> passMuFilterVector = {"All", "Low", "Med", "High"};
  vector<string> failMuFilterVector = {"All", "Low", "Med", "High"};
  muFilter_puCat_selection["None"] = noneMuFilterVector;
  muFilter_puCat_selection["Pass"] = passMuFilterVector;
  muFilter_puCat_selection["Fail"] = failMuFilterVector;

  for(auto itMu=muFilter_puCat_selection.begin(); itMu != muFilter_puCat_selection.end(); ++itMu){
    Bool_t doMHT = true; // Muon filter is needed for MET only.
    if( itMu->first != "None" ) doMHT = false;
    for(auto itPu=itMu->second.begin(); itPu != itMu->second.end(); ++itPu){
      main_func(*itPu, itMu->first, doMHT);
    }
  }
}
