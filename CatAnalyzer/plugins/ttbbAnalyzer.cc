// -*- C++ -*-
//
// Package:    ttbbAnalyzer
// Class:      ttbbAnalyzer
//
/**\class ttbbAnalyzer ttbbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Javier Brochero Cifuentes
//
//

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenWeights.h"

#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Kinematic Reconstruction
#include "CATTools/CatAnalyzer/interface/LepJetsFitter.h"

// JES Components
#include "CondFormats/JetMETObjects/interface/Utilities.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"


#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//
using namespace cat;

class ttbbAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ttbbAnalyzer(const edm::ParameterSet&);
  ~ttbbAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  //----------------------------------------------------------------
  bool IsSelectMuon    (const cat::Muon     & i_muon_candidate);
  bool IsVetoMuon      (const cat::Muon     & i_muon_candidate);
  bool IsSelectElectron(const cat::Electron & i_electron_candidate);
  bool IsVetoElectron  (const cat::Electron & i_electron_candidate);

  bool isMC_ ;

  int TTbarMC_; // 0->No ttbar, 1->ttbar Signal, 2->ttbar Background
  int TTbarCatMC_;
  unsigned int SkimNJets_;
  bool KFUsebtag_;
  bool CSVPosConKF_;
  // Trigger Names
  std::vector<string> triggerNameDataEl_;
  std::vector<string> triggerNameDataMu_;
  std::vector<string> triggerNameMCEl_;
  std::vector<string> triggerNameMCMu_;
  // JES File
  string JESUncFile_;
  
  // Event Weights
  edm::EDGetTokenT<float>                        genWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           pdfWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           scaleUpWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           scaleDownWeightToken_;
  edm::EDGetTokenT<cat::GenWeights>              hdampWeightToken_;
  edm::EDGetTokenT<float>                        puWeightToken_;
  edm::EDGetTokenT<float>                        puUpWeightToken_;
  edm::EDGetTokenT<float>                        puDownWeightToken_;
  edm::EDGetTokenT<float>                        CSVWeightToken_;
  edm::EDGetTokenT<std::vector<float>>           CSVSysWeightToken_;
  // Object Collections
  edm::EDGetTokenT<reco::GenParticleCollection>  genToken_;
  edm::EDGetTokenT<reco::GenJetCollection>       genJetToken_;
  edm::EDGetTokenT<int>                          genttbarHiggsCatToken_;
  edm::EDGetTokenT<cat::GenTopCollection>        genttbarCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
  edm::EDGetTokenT<cat::JetCollection>           jetToken_;
  edm::EDGetTokenT<cat::METCollection>           metToken_;
  edm::EDGetTokenT<int>                          pvToken_;
  edm::EDGetTokenT<std::vector<vector<int>>>     JetMotherToken_;
  // Trigger
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits2_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

// ----------member data ---------------------------

  TTree *tree;
  TTree *gentree;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Tree Branches
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Event info
  int b_Event, b_Run, b_Lumi_Number;
  float b_GenWeight;
  std::vector<float> *b_ScaleWeight;
  std::vector<float> *b_PDFWeight;
  std::vector<float> *b_hdampWeight;
  // PU/Vertices
  std::vector<float> *b_PUWeight;
  int b_nGoodPV;
  // Channel
  int b_GenChannel;
  int b_Channel;
  //  Event categorization
  std::vector<int>   *b_GenConeCatID;
  int b_GenHiggsCatID;
  // Jet info from event categorization
  std::vector<float> *b_GenCone_gJet_pT;
  std::vector<float> *b_GenCone_gJet_eta;
  std::vector<float> *b_GenCone_gJet_phi;
  std::vector<float> *b_GenCone_gJet_E;
  std::vector<int>   *b_GenCone_gJetFlavW;
  std::vector<int>   *b_GenCone_gJetIndex;
  int b_GenCone_NgJetsW;
  // MET
  float b_MET, b_MET_phi;
  // GEN Leptons
  float b_GenLepton_pT;
  float b_GenLepton_eta;
  float b_GenLepton_phi;
  float b_GenLepton_E;
  // GEN Jets
  std::vector<float> *b_GenJet_pT;
  std::vector<float> *b_GenJet_eta;
  std::vector<float> *b_GenJet_phi;
  std::vector<float> *b_GenJet_E;
  // GEN Neutrino
  float b_GenNu_pT;
  float b_GenNu_eta;
  float b_GenNu_phi;
  float b_GenNu_E;
  // Leptons
  float b_Lepton_pT;
  float b_Lepton_eta;
  float b_Lepton_phi;
  float b_Lepton_E;
  float b_Lepton_LES;
  // First additional b-jets 
  float b_addbjet1_pt; 
  float b_addbjet1_eta; 
  float b_addbjet1_phi; 
  float b_addbjet1_e; 
  // Second additional b-jets 
  float b_addbjet2_pt;
  float b_addbjet2_eta;
  float b_addbjet2_phi;
  float b_addbjet2_e; 
  // Jet Mother (MC Studies)
  std::vector<int> *b_GenJet_mom, *b_GenJet_GenConeMom;
  // Jets
  std::vector<float> *b_Jet_pT;
  std::vector<float> *b_Jet_eta;
  std::vector<float> *b_Jet_phi;
  std::vector<float> *b_Jet_E;
  std::vector<int>   *b_Jet_Index, *b_Jet_GenConeMom;
  std::vector<int>   *b_Jet_MatchedGenJetIndex;
  float               b_DRAddJets;
  // Jet Flavour
  std::vector<int> *b_Jet_partonFlavour;
  std::vector<int> *b_Jet_hadronFlavour;
  // JES and JER
  std::vector< std::vector<float> > *b_Jet_JESCom_Up, *b_Jet_JESCom_Down;
  std::vector<float> *b_Jet_JES_Up, *b_Jet_JES_Down;
  std::vector<float> *b_Jet_JER_Up, *b_Jet_JER_Nom, *b_Jet_JER_Down;
  // b-Jet discriminant
  std::vector<float> *b_Jet_CSV;
  std::vector<float> *b_Jet_SF_CSV;
  std::vector<float> *b_Jet_SF_CSVcsv;
  // c-Jet discriminant
  std::vector<float> *b_Jet_CvsL, *b_Jet_CvsB;
  // Kinematic Reconstruction
  float b_Kin_Chi2;
  // KR Leptons
  float b_KinNu_pT;
  float b_KinNu_eta;
  float b_KinNu_phi;
  float b_KinNu_E;
  // KR Jets
  std::vector<float> *b_KinJet_pT;
  std::vector<float> *b_KinJet_eta;
  std::vector<float> *b_KinJet_phi;
  std::vector<float> *b_KinJet_E;
  std::vector<int>   *b_KinJet_Index;
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Histograms
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  
  // Histograms: Number of Events and Weights
  TH1D *EventInfo, *ScaleWeights;
  // Scale factor evaluators
  BTagWeightEvaluator SF_CSV_;
};

//
// constants, enums and typedefs
//
// JES components 
const int nsrc = 25;
const char* srcnames[nsrc] =
  {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
   "Fragmentation", 
   "SinglePionECAL", "SinglePionHCAL", 
   "FlavorQCD",
   "TimePtEta", 
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
   "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
   "RelativeBal",
   "RelativeFSR",
   "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF",
   "PileUpDataMC",
   "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtHF"};


// static data member definitions
//

// constructors and destructor
//
ttbbAnalyzer::ttbbAnalyzer(const edm::ParameterSet& iConfig):
  TTbarMC_    (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0)),
  TTbarCatMC_ (iConfig.getUntrackedParameter<int>("TTbarCatLabel", 0)),
  SkimNJets_  (iConfig.getUntrackedParameter<unsigned int>("Skim_N_Jets", 0)),
  KFUsebtag_  (iConfig.getUntrackedParameter<bool>("KFUsebtagLabel", true)),
  CSVPosConKF_(iConfig.getUntrackedParameter<bool>("CSVPosConKFLabel", true)),
  triggerNameDataEl_(iConfig.getUntrackedParameter<std::vector<string>>("triggerNameDataEl")),
  triggerNameDataMu_(iConfig.getUntrackedParameter<std::vector<string>>("triggerNameDataMu")),
  triggerNameMCEl_  (iConfig.getUntrackedParameter<std::vector<string>>("triggerNameMCEl")),
  triggerNameMCMu_  (iConfig.getUntrackedParameter<std::vector<string>>("triggerNameMCMu")),
  JESUncFile_       (iConfig.getUntrackedParameter<string>("JESUncFile"))
{

  // CSV Scale Factors   
  SF_CSV_.initCSVWeight(false, "csvv2");
  
  // Weights
  auto genWeightLabel = iConfig.getParameter<edm::InputTag>("genWeightLabel");
  // aMC@NLO
  genWeightToken_        = consumes<float>             (edm::InputTag(genWeightLabel.label()));
  // PDF
  pdfWeightToken_       = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "pdf"));
  // Scale
  scaleUpWeightToken_   = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaleup"));
  scaleDownWeightToken_ = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaledown"));
  // hdamp
  hdampWeightToken_       = consumes<cat::GenWeights>(edm::InputTag("genWeight"));
  // PileUp
  auto puWeightLabel = iConfig.getParameter<edm::InputTag>("puWeightLabel");
  puWeightToken_         = consumes<float>             (edm::InputTag(puWeightLabel.label()));
  puUpWeightToken_       = consumes<float>             (edm::InputTag(puWeightLabel.label(),"up"));
  puDownWeightToken_     = consumes<float>             (edm::InputTag(puWeightLabel.label(),"dn"));
  // CSV Weights
  auto csvWeightLabel = iConfig.getParameter<edm::InputTag>("csvWeightLabel");
  CSVWeightToken_        = consumes<float>             (edm::InputTag(csvWeightLabel.label()));
  CSVSysWeightToken_     = consumes<std::vector<float>>(edm::InputTag(csvWeightLabel.label(), "syst"));
  // GEN level
  genToken_              = consumes<reco::GenParticleCollection>  (iConfig.getParameter<edm::InputTag>("genLabel"));
  genJetToken_           = consumes<reco::GenJetCollection>       (iConfig.getParameter<edm::InputTag>("genJetLabel"));  
  //ttbar Categorization
  genttbarCatToken_      = consumes<cat::GenTopCollection>        (iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));
  genttbarHiggsCatToken_ = consumes<int>                          (iConfig.getParameter<edm::InputTag>("genHiggsCatLabel"));
  // Object Collections
  muonToken_         = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_     = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_          = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_          = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_           = consumes<int>                          (iConfig.getParameter<edm::InputTag>("pvLabel"));
  JetMotherToken_    = consumes<vector<vector<int>>>          (iConfig.getParameter<edm::InputTag>("JetMother"));
  // Trigger  
  auto triggerLabel  = iConfig.getParameter<edm::InputTag>("triggerBits");
  triggerBits_       = consumes<edm::TriggerResults>                    (edm::InputTag(triggerLabel.label(),"","HLT"));
  triggerBits2_      = consumes<edm::TriggerResults>                    (edm::InputTag(triggerLabel.label(),"","HLT2"));
  triggerObjects_    = consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerObjects"));
  
  // Vector initialization
  // Weights
  b_PUWeight     = new std::vector<float>;
  b_PDFWeight    = new std::vector<float>;  
  b_ScaleWeight  = new std::vector<float>;  
  b_hdampWeight  = new std::vector<float>;  
  // ttbar Categorization
  b_GenConeCatID      = new std::vector<int>;

  b_GenCone_gJet_pT   = new std::vector<float>;
  b_GenCone_gJet_eta  = new std::vector<float>;
  b_GenCone_gJet_phi  = new std::vector<float>;
  b_GenCone_gJet_E    = new std::vector<float>;
  b_GenCone_gJetIndex = new std::vector<int>;
  b_GenCone_gJetFlavW = new std::vector<int>;
  b_GenJet_pT = new std::vector<float>;
  b_GenJet_eta = new std::vector<float>;
  b_GenJet_phi = new std::vector<float>;
  b_GenJet_E  = new std::vector<float>;
  b_GenJet_mom = new std::vector<int>;
  b_GenJet_GenConeMom = new std::vector<int>;
  // Jets
  b_Jet_pT   = new std::vector<float>;
  b_Jet_eta   = new std::vector<float>;
  b_Jet_phi   = new std::vector<float>;
  b_Jet_E    = new std::vector<float>;
  b_Jet_Index= new std::vector<int>;
  b_Jet_MatchedGenJetIndex = new std::vector<int>;
  b_Jet_GenConeMom=new std::vector<int>;
  b_Jet_CSV  = new std::vector<float>;
  b_Jet_SF_CSV     = new std::vector<float>;
  b_Jet_SF_CSVcsv  = new std::vector<float>;
  b_Jet_CvsL    = new std::vector<float>;  
  b_Jet_CvsB    = new std::vector<float>;  
  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;
  b_Jet_JESCom_Up   = new std::vector< std::vector<float> >;
  b_Jet_JESCom_Down = new std::vector< std::vector<float> >;
  b_Jet_JES_Up   = new std::vector<float>;
  b_Jet_JES_Down = new std::vector<float>;
  b_Jet_JER_Up   = new std::vector<float>;
  b_Jet_JER_Nom  = new std::vector<float>;
  b_Jet_JER_Down = new std::vector<float>;
  // Kinematic Reconstruction
  b_KinJet_pT  = new std::vector<float>;
  b_KinJet_eta = new std::vector<float>;
  b_KinJet_phi = new std::vector<float>;
  b_KinJet_E   = new std::vector<float>;

  b_KinJet_Index = new std::vector<int>;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "event/I");
  tree->Branch("run",        &b_Run,         "run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "luminumber/I");

  tree->Branch("genweight",  &b_GenWeight,   "genweight/F");
  tree->Branch("PUWeight",   "std::vector<float>", &b_PUWeight);

  tree->Branch("GoodPV",     &b_nGoodPV,     "GoodPV/I");

  tree->Branch("channel",    &b_Channel,     "channel/I");

  tree->Branch("MET",     &b_MET,     "MET/F");
  tree->Branch("MET_phi", &b_MET_phi, "MET_phi/F");

  tree->Branch("lepton_pT",  &b_Lepton_pT,  "lepton_pT/F");
  tree->Branch("lepton_eta", &b_Lepton_eta, "lepton_eta/F");
  tree->Branch("lepton_phi", &b_Lepton_phi, "lepton_phi/F");
  tree->Branch("lepton_E" ,  &b_Lepton_E,   "lepton_E/F" );


  tree->Branch("jet_pT",           "std::vector<float>", &b_Jet_pT);
  tree->Branch("jet_eta",          "std::vector<float>", &b_Jet_eta);
  tree->Branch("jet_phi",          "std::vector<float>", &b_Jet_phi);
  tree->Branch("jet_E" ,           "std::vector<float>", &b_Jet_E );
  tree->Branch("jet_index" ,       "std::vector<int>",   &b_Jet_Index );
  tree->Branch("jet_gencone_mom" , "std::vector<int>",   &b_Jet_GenConeMom );
  tree->Branch("jet_CSV" ,         "std::vector<float>", &b_Jet_CSV );
  tree->Branch("jet_SF_CSV",       "std::vector<float>", &b_Jet_SF_CSV );
  tree->Branch("jet_SF_CSVcsv",    "std::vector<float>", &b_Jet_SF_CSVcsv );
  tree->Branch("jet_CvsL",         "std::vector<float>", &b_Jet_CvsL );
  tree->Branch("jet_CvsB",         "std::vector<float>", &b_Jet_CvsB );

  tree->Branch("jet_partonFlavour", "std::vector<int>",   &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>",   &b_Jet_hadronFlavour);

  // Kinematic Reconstruction
  tree->Branch("kin_chi2",   &b_Kin_Chi2,  "kin_chi2/F");
  tree->Branch("kinnu_pT",   &b_KinNu_pT,  "kinnu_pT/F");
  tree->Branch("kinnu_eta",  &b_KinNu_eta, "kinnu_eta/F");
  tree->Branch("kinnu_phi",  &b_KinNu_phi, "kinnu_phi/F");
  tree->Branch("kinnu_E",    &b_KinNu_E,   "kinnu_E/F");
  
  tree->Branch("kinjet_pT",    "std::vector<float>", &b_KinJet_pT);
  tree->Branch("kinjet_eta",   "std::vector<float>", &b_KinJet_eta);
  tree->Branch("kinjet_phi",   "std::vector<float>", &b_KinJet_phi);
  tree->Branch("kinjet_E",     "std::vector<float>", &b_KinJet_E);
  tree->Branch("kinjet_index", "std::vector<int>",   &b_KinJet_Index);

  // Systematic variations
  tree->Branch("lepton_LES", &b_Lepton_LES, "lepton_LES/F" );

  tree->Branch("jet_JESCom_Up",     "std::vector< std::vector<float> >", &b_Jet_JESCom_Up );
  tree->Branch("jet_JESCom_Down",   "std::vector< std::vector<float> >", &b_Jet_JESCom_Down );
  tree->Branch("jet_JES_Up",        "std::vector<float>", &b_Jet_JES_Up );
  tree->Branch("jet_JES_Down",      "std::vector<float>", &b_Jet_JES_Down );

  tree->Branch("jet_JER_Up",        "std::vector<float>", &b_Jet_JER_Up );
  tree->Branch("jet_JER_Nom",       "std::vector<float>", &b_Jet_JER_Nom );
  tree->Branch("jet_JER_Down",      "std::vector<float>", &b_Jet_JER_Down );

 
  // GEN Variables (only ttbarSignal)
  if(TTbarMC_ == 1){
   
    tree->Branch("jet_MatchedGenJetIndex", "std::vector<int>",  &b_Jet_MatchedGenJetIndex);

    tree->Branch("genconecatid" , "std::vector<int>",   &b_GenConeCatID);
    tree->Branch("genhiggscatid",  &b_GenHiggsCatID,   "genhiggscatid/I");

    tree->Branch("gencone_gjet_pT" ,   "std::vector<float>", &b_GenCone_gJet_pT);
    tree->Branch("gencone_gjet_eta" ,  "std::vector<float>", &b_GenCone_gJet_eta);
    tree->Branch("gencone_gjet_phi" ,  "std::vector<float>", &b_GenCone_gJet_phi);
    tree->Branch("gencone_gjet_E" ,    "std::vector<float>", &b_GenCone_gJet_E);
    tree->Branch("gencone_gjetIndex" , "std::vector<int>",   &b_GenCone_gJetIndex);
    tree->Branch("gencone_gJetFlavW" , "std::vector<int>",   &b_GenCone_gJetFlavW);
    tree->Branch("gencone_NgjetsW",    &b_GenCone_NgJetsW,   "gencone_NgjetsW/I");
    tree->Branch("draddjets",          &b_DRAddJets,         "draddjets/F");

    tree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
    tree->Branch("genlepton_pT",  &b_GenLepton_pT,  "genlepton_pT/F");
    tree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");
    tree->Branch("genlepton_phi", &b_GenLepton_phi, "genlepton_phi/F");
    tree->Branch("genlepton_E",   &b_GenLepton_E,   "genlepton_E/F");
    tree->Branch("gennu_pT",      &b_GenNu_pT,      "gennu_pT/F");
    tree->Branch("gennu_eta",     &b_GenNu_eta,     "gennu_eta/F");
    tree->Branch("gennu_phi",     &b_GenNu_phi,     "gennu_phi/F");
    tree->Branch("gennu_E",       &b_GenNu_E,       "gennu_E/F");

    tree->Branch("genjet_pT",          "std::vector<float>", &b_GenJet_pT);
    tree->Branch("genjet_eta",         "std::vector<float>", &b_GenJet_eta);
    tree->Branch("genjet_phi",         "std::vector<float>", &b_GenJet_phi);
    tree->Branch("genjet_E",           "std::vector<float>", &b_GenJet_E);
    tree->Branch("genjet_mom",         "std::vector<int>",   &b_GenJet_mom);
    tree->Branch("genjet_gencone_mom", "std::vector<int>",   &b_GenJet_GenConeMom);

    tree->Branch("addbjet1_pt",  &b_addbjet1_pt,  "addbjet1_pt/F"); 
    tree->Branch("addbjet1_eta", &b_addbjet1_eta, "addbjet1_eta/F"); 
    tree->Branch("addbjet1_phi", &b_addbjet1_phi, "addbjet1_phi/F"); 
    tree->Branch("addbjet1_e",   &b_addbjet1_e,   "addbjet1_e/F"); 

    tree->Branch("addbjet2_pt",  &b_addbjet2_pt,  "addbjet2_pt/F");
    tree->Branch("addbjet2_eta", &b_addbjet2_eta, "addbjet2_eta/F");
    tree->Branch("addbjet2_phi", &b_addbjet2_phi, "addbjet2_phi/F");
    tree->Branch("addbjet2_e",   &b_addbjet2_e,   "addbjet2_e/F");  

    // Systematic variations for ttbar    
    tree->Branch("pdfweight",   "std::vector<float>", &b_PDFWeight );
    tree->Branch("scaleweight", "std::vector<float>", &b_ScaleWeight );
    tree->Branch("hdampweight", "std::vector<float>", &b_hdampWeight );


    //GEN TREE
    gentree = fs->make<TTree>("gentree", "TopGENTree");
    gentree->Branch("genweight",     &b_GenWeight,     "genweight/F");
    gentree->Branch("genchannel",    &b_GenChannel,    "genchannel/I");
 
    gentree->Branch("genconecatid" ,  "std::vector<int>",   &b_GenConeCatID);
    gentree->Branch("genhiggscatid",  &b_GenHiggsCatID,     "genhiggscatid/I");
    
    gentree->Branch("draddjets",     &b_DRAddJets,     "draddjets/F");
    gentree->Branch("genlepton_pT",  &b_GenLepton_pT,  "genlepton_pT/F");
    gentree->Branch("genlepton_eta", &b_GenLepton_eta, "genlepton_eta/F");

    gentree->Branch("scaleweight",        "std::vector<float>", &b_ScaleWeight );

    gentree->Branch("genjet_pT",          "std::vector<float>", &b_GenJet_pT);
    gentree->Branch("genjet_eta",         "std::vector<float>", &b_GenJet_eta);
    gentree->Branch("genjet_phi",         "std::vector<float>", &b_GenJet_phi);
    gentree->Branch("genjet_E",           "std::vector<float>", &b_GenJet_E);
    gentree->Branch("genjet_mom",         "std::vector<int>",   &b_GenJet_mom);
    gentree->Branch("genjet_gencone_mom", "std::vector<int>",   &b_GenJet_GenConeMom);
  }

  EventInfo = fs->make<TH1D>("EventInfo","Event Information",9,0,9);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
  EventInfo->GetXaxis()->SetBinLabel(3,"ttbb Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(4,"ttbb Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(5,"ttbb Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(6,"ttjj Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(7,"ttjj Lep Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(8,"ttjj Lep (tau lep decay) Events pT>20");
  EventInfo->GetXaxis()->SetBinLabel(9,"Sum of PU Weights");

  ScaleWeights = fs->make<TH1D>("ScaleWeights","Event Weights:",6,0,6);
  ScaleWeights->GetXaxis()->SetBinLabel(1,"muR=Nom  muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(2,"muR=Nom  muF=Down");
  ScaleWeights->GetXaxis()->SetBinLabel(3,"muR=Up   muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(4,"muR=Up   muF=Up");
  ScaleWeights->GetXaxis()->SetBinLabel(5,"muR=Down muF=Nom");
  ScaleWeights->GetXaxis()->SetBinLabel(6,"muR=Down muF=Down");


}


ttbbAnalyzer::~ttbbAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete b_PUWeight;
  delete b_PDFWeight;
  delete b_hdampWeight;
  delete b_ScaleWeight;

  delete b_GenConeCatID;
  delete b_GenCone_gJet_pT;
  delete b_GenCone_gJet_eta;
  delete b_GenCone_gJet_phi;
  delete b_GenCone_gJet_E;
  delete b_GenCone_gJetIndex;
  delete b_GenCone_gJetFlavW;

  delete b_GenJet_pT;
  delete b_GenJet_eta;
  delete b_GenJet_phi;
  delete b_GenJet_E;
  delete b_GenJet_mom;
  delete b_GenJet_GenConeMom;

  delete b_Jet_pT;
  delete b_Jet_eta;
  delete b_Jet_phi;
  delete b_Jet_E;
  delete b_Jet_Index;
  delete b_Jet_GenConeMom;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_MatchedGenJetIndex;

  delete b_Jet_JESCom_Up;
  delete b_Jet_JESCom_Down;
  delete b_Jet_JES_Up;
  delete b_Jet_JES_Down;
  delete b_Jet_JER_Up;
  delete b_Jet_JER_Nom;
  delete b_Jet_JER_Down;

  delete b_Jet_SF_CSV;
  delete b_Jet_SF_CSVcsv;

  delete b_Jet_CSV;
  delete b_Jet_CvsL;
  delete b_Jet_CvsB;

  delete b_KinJet_pT;
  delete b_KinJet_eta;
  delete b_KinJet_phi;
  delete b_KinJet_E;
  
  delete b_KinJet_Index;

}

//
// member functions
//

// ------------ method called for each event  ------------
void ttbbAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_PUWeight   ->clear();
  b_ScaleWeight->clear();
  b_PDFWeight  ->clear();
  b_hdampWeight->clear();

  b_GenConeCatID->clear();
  b_GenCone_gJet_pT->clear();
  b_GenCone_gJet_eta->clear();
  b_GenCone_gJet_phi->clear();
  b_GenCone_gJet_E->clear();
  b_GenCone_gJetIndex->clear();
  b_GenCone_gJetFlavW->clear();


  b_GenJet_pT  ->clear();
  b_GenJet_eta ->clear();
  b_GenJet_phi ->clear();
  b_GenJet_E   ->clear();
  b_GenJet_mom->clear();
  b_GenJet_GenConeMom->clear();

  b_Jet_pT    ->clear();
  b_Jet_eta   ->clear();
  b_Jet_phi   ->clear();
  b_Jet_E     ->clear();
  b_Jet_Index->clear();
  b_Jet_GenConeMom->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_MatchedGenJetIndex->clear();

  b_Jet_JESCom_Up  ->clear();
  b_Jet_JESCom_Down->clear();
  b_Jet_JES_Up  ->clear();
  b_Jet_JES_Down->clear();
  b_Jet_JER_Up  ->clear();
  b_Jet_JER_Nom ->clear();
  b_Jet_JER_Down->clear();

  b_Jet_CSV      ->clear();
  b_Jet_SF_CSV   ->clear();
  b_Jet_SF_CSVcsv->clear();
  b_Jet_CvsL     ->clear();
  b_Jet_CvsB     ->clear();

  b_KinJet_pT ->clear();
  b_KinJet_eta->clear();
  b_KinJet_phi->clear();
  b_KinJet_E  ->clear();
  
  b_KinJet_Index->clear();

  b_addbjet1_pt  = -1.0;
  b_addbjet1_eta = -1.0;
  b_addbjet1_phi = -1.0;
  b_addbjet1_e   = -1.0;

  b_addbjet2_pt  = -1.0;
  b_addbjet2_eta = -1.0;
  b_addbjet2_phi = -1.0;
  b_addbjet2_e   = -1.0;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  isMC_ = !iEvent.isRealData();

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();

  EventInfo->Fill(0.5, 1.0);         // Number of Events

  if(isMC_) {

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // PU Info
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    edm::Handle<float> PUWeight;
    iEvent.getByToken(puWeightToken_, PUWeight);
    b_PUWeight->push_back(*PUWeight); // Central

    EventInfo->Fill(8.5, *PUWeight); // Sum of PUWeights

    edm::Handle<float> PUWeight_Up;
    iEvent.getByToken(puUpWeightToken_, PUWeight_Up);
    b_PUWeight->push_back(*PUWeight_Up); //Syst. Up

    edm::Handle<float> PUWeight_Down;
    iEvent.getByToken(puDownWeightToken_, PUWeight_Down);
    b_PUWeight->push_back(*PUWeight_Down); //Syst. Down

    //---------------------------------------------------------------------------
    // Weights at Generation Level: aMC@NLO
    //---------------------------------------------------------------------------

    edm::Handle<float> genWeight;

    iEvent.getByToken(genWeightToken_, genWeight);
    b_GenWeight = *genWeight;

    EventInfo->Fill(1.5, b_GenWeight); // Sum of aMC@NLO Weights

  }

  else{
    b_PUWeight->push_back(1.0);
    b_GenWeight = 1.0;
  }

  //---------------------------------------------------------------------------
  // Weights for Syst. Scale and PDF: ttbar
  //---------------------------------------------------------------------------
  if(TTbarMC_ == 1 ) {
    edm::Handle<std::vector<float>> scaleUpWeightsHandle, scaleDownWeightsHandle;
    iEvent.getByToken(scaleUpWeightToken_,   scaleUpWeightsHandle);
    iEvent.getByToken(scaleDownWeightToken_, scaleDownWeightsHandle);

    // muR/muF Scale Weights
    if(scaleUpWeightsHandle->size()!=0 && scaleDownWeightsHandle->size()!=0){
      b_ScaleWeight->push_back(scaleUpWeightsHandle  ->at(0)); // muR=Nom  muF=Up
      b_ScaleWeight->push_back(scaleDownWeightsHandle->at(0)); // muR=Nom  muF=Down
      b_ScaleWeight->push_back(scaleUpWeightsHandle  ->at(1)); // muR=Up   muF=Nom
      b_ScaleWeight->push_back(scaleUpWeightsHandle  ->at(2)); // muR=Up   muF=Up
      b_ScaleWeight->push_back(scaleDownWeightsHandle->at(1)); // muR=Down muF=Nom
      b_ScaleWeight->push_back(scaleDownWeightsHandle->at(2)); // muR=Down muF=Down
    }
    
    // Sum of muR/muF Scale Weights
    for(unsigned int iscale = 0; iscale< b_ScaleWeight->size(); iscale++)
      ScaleWeights->Fill(iscale, b_ScaleWeight->at(iscale)); 
    
    // PDF Weights
    edm::Handle<std::vector<float>> PDFWeightsHandle;
    iEvent.getByToken(pdfWeightToken_,   PDFWeightsHandle);

    for ( auto& w : *PDFWeightsHandle ) b_PDFWeight->push_back(w);

    // hdamp: Only for ttbar. hdamp weights start at 227
    edm::Handle<cat::GenWeights> hdampWeightsHandle;
    iEvent.getByToken(hdampWeightToken_,   hdampWeightsHandle);
    if (hdampWeightsHandle->weights().size()>=254){
      for(int ihdamp=227; ihdamp<254; ihdamp++) b_hdampWeight->push_back(hdampWeightsHandle->weights().at(ihdamp)); 
    }
    
  }
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Generated Particles (For Pythia8)
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  int nGenLep  = 999;

  bool IsCat = false;

  if(TTbarMC_ > 0) {

    //---------------------------------------------------------------------------
    // Event Categorization Using Higgs Code
    // // Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GenHFHadronMatcher
    //---------------------------------------------------------------------------
    edm::Handle<int> genttbarHiggsCatHandle;
    iEvent.getByToken( genttbarHiggsCatToken_, genttbarHiggsCatHandle );
    b_GenHiggsCatID = *genttbarHiggsCatHandle;

    //---------------------------------------------------------------------------
    // Event Categorization Using Cone
    //---------------------------------------------------------------------------
    edm::Handle<cat::GenTopCollection> genttbarConeCat;
    iEvent.getByToken( genttbarCatToken_, genttbarConeCat );

    // [0]: Decay mode. Just as a Cross Check!
    b_GenConeCatID->push_back(genttbarConeCat->begin()->semiLeptonic(0));
    // [1]: Number of Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NJets20());
    // [2]: Number of b-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NbJets20());
    // [3]: Number of c-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NcJets20());
    // [4]: Number of b-Jets Not comming from the top
    b_GenConeCatID->push_back(genttbarConeCat->begin()->NbJets20NoTop());
    // [5]: Number of add Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()-> NaddJets20());
    // [6]: Number of add b-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()-> NaddbJets20());
    // [7]: Number of add c-Jets
    b_GenConeCatID->push_back(genttbarConeCat->begin()-> NaddcJets20());

    // additional Jets Information
    // p4
    math::XYZTLorentzVector gJetGenCone[6];
    gJetGenCone[0] = genttbarConeCat->begin()->bJetsFromTop1();
    gJetGenCone[1] = genttbarConeCat->begin()->bJetsFromTop2();
    gJetGenCone[2] = genttbarConeCat->begin()->JetsFromW1();
    gJetGenCone[3] = genttbarConeCat->begin()->JetsFromW2();
    gJetGenCone[4] = genttbarConeCat->begin()->addJets1();
    gJetGenCone[5] = genttbarConeCat->begin()->addJets2();

    b_GenCone_NgJetsW = genttbarConeCat->begin()-> NWJets();
    
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()-> Wquarkflav1());
    b_GenCone_gJetFlavW->push_back(genttbarConeCat->begin()-> Wquarkflav2());

    // All Jets info from GenCone
    for (int ijGT=0; ijGT<6; ijGT++){
      b_GenCone_gJet_pT  ->push_back(gJetGenCone[ijGT].Pt());
      b_GenCone_gJet_eta ->push_back(gJetGenCone[ijGT].Eta());
      b_GenCone_gJet_phi ->push_back(gJetGenCone[ijGT].Phi());
      b_GenCone_gJet_E   ->push_back(gJetGenCone[ijGT].E());
    }

    // adding additional b jet four-momentum 
    b_addbjet1_pt  = genttbarConeCat->begin()->addbJets1().Pt();
    b_addbjet1_eta = genttbarConeCat->begin()->addbJets1().Eta();
    b_addbjet1_phi = genttbarConeCat->begin()->addbJets1().Phi();
    b_addbjet1_e   = genttbarConeCat->begin()->addbJets1().E();

    b_addbjet2_pt  = genttbarConeCat->begin()->addbJets2().Pt();
    b_addbjet2_eta = genttbarConeCat->begin()->addbJets2().Eta();
    b_addbjet2_phi = genttbarConeCat->begin()->addbJets2().Phi();
    b_addbjet2_e   = genttbarConeCat->begin()->addbJets2().E();

    // DR 
    b_DRAddJets = genttbarConeCat->begin()->dRaddJets();

    if(genttbarConeCat->begin()-> NaddbJets20() > 1){
      EventInfo->Fill(2.5, 1.0); // Number of ttbb Events	
      if(genttbarConeCat->begin()->semiLeptonic(-1)) EventInfo->Fill(3.5, 1.0); // Number of ttbb Events (Includes tau) 
      if(genttbarConeCat->begin()->semiLeptonic(0))  EventInfo->Fill(4.5, 1.0); // Number of ttbb Events (Includes tau leptonic decay) 
    }
    if(genttbarConeCat->begin()-> NaddJets20() > 1){
      EventInfo->Fill(5.5, 1.0); // Number of ttjj Events	
      if(genttbarConeCat->begin()->semiLeptonic(-1)) EventInfo->Fill(6.5, 1.0); // Number of ttjj Events (Includes tau) 
      if(genttbarConeCat->begin()->semiLeptonic(0))  EventInfo->Fill(7.5, 1.0); // Number of ttjj Events (Includes tau leptonic decay) 
    }
    
    
    //---------------------------------------------------------------------------
    // Using the GenChannel from GenTop categorization
    //---------------------------------------------------------------------------
    nGenLep = genttbarConeCat->begin()->semiLeptonic(0); // semiLeptonic(0) includes tau leptonic decay
    //---------------------------------------------------------------------------

    // Category
    bool Isttjj = false;
    bool Isttbb = false;
    bool Isttb  = false;
    bool Isttcc = false;
    bool IsttLF = false;
    bool Istt   = false;

    // Categorization based in the Full Ph-Sp
    // Requires ttjj events to be categorized
    if(genttbarConeCat->begin()->NaddJets20() > 1) Isttjj = true;

    if      (Isttjj && genttbarConeCat->begin()->NaddbJets20() > 1) Isttbb = true;
    else if (Isttjj && genttbarConeCat->begin()->NaddbJets20() > 0) Isttb  = true;
    else if (Isttjj && genttbarConeCat->begin()->NaddcJets20() > 1) Isttcc = true;
    else if (Isttjj && genttbarConeCat->begin()->NaddJets20()  > 1) IsttLF = true;
    else Istt = true;

    // Categorization based in the Visible Ph-Sp
    // if(genttbarConeCat->begin()->NbJets20() > 1 && 
    //    genttbarConeCat->begin()->NJets20()  > 5) Isttjj = true;

    // if      (genttbarConeCat->begin()->NbJets20() > 3  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5) Isttbb = true;
    // else if (genttbarConeCat->begin()->NbJets20() > 2  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5) Isttb  = true;
    // else if (genttbarConeCat->begin()->NbJets20() > 1  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5  && 
    // 	     genttbarConeCat->begin()->NcJets20() > 1) Isttcc = true;
    // else if (genttbarConeCat->begin()->NbJets20() > 1  && 
    // 	     genttbarConeCat->begin()->NJets20()  > 5) IsttLF = true;
    // else Istt = true;


    if (TTbarCatMC_ == 0)           IsCat = true;
    if (Isttbb && TTbarCatMC_ == 1) IsCat = true;
    if (Isttb  && TTbarCatMC_ == 2) IsCat = true;
    if (Isttcc && TTbarCatMC_ == 3) IsCat = true;
    if (IsttLF && TTbarCatMC_ == 4) IsCat = true;
    if (Istt   && TTbarCatMC_ == 5) IsCat = true;
    if (Isttjj && TTbarCatMC_ == 6) IsCat = true;


    if(isMC_ && TTbarMC_== 1 && 
       IsCat && 
       nGenLep == 1){

      if (genttbarConeCat->begin()->semiLeptonicMuo())      b_GenChannel = 0;
      else if (genttbarConeCat->begin()->semiLeptonicEle()) b_GenChannel = 1;
	
      if(genttbarConeCat->begin()->lepton1().pt() != 0.){
	b_GenLepton_pT  = genttbarConeCat->begin()->lepton1().pt();
	b_GenLepton_eta = genttbarConeCat->begin()->lepton1().eta();
      }
      else{
	b_GenLepton_pT  = genttbarConeCat->begin()->lepton2().pt();
	b_GenLepton_eta = genttbarConeCat->begin()->lepton2().eta();
      }

      //---------------------------------------------------------------------------
      // Generated particles: Lepton, neutrino and jets
      //---------------------------------------------------------------------------
	
      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByToken(genToken_, genParticles);
	
      TLorentzVector Genlepton(0,0,0,0); 
      TLorentzVector Gennu(0,0,0,0);

      // Gen Status: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
      // abs(id) == 6  // t-Quark
      // abs(id) == 5  // b-Quark
      // abs(id) == 24 // W-Boson
      // abs(id) == 15 // Tau
      // abs(id) == 11 // Electron
      // abs(id) == 13 // Muon

      for (unsigned int i = 0; i < genParticles->size(); i++){

	const reco::GenParticle & gp = (*genParticles)[i];
	const int id = gp.pdgId();
	// GenFlags
	reco::GenStatusFlags statusflags = gp.statusFlags();
	
	bool CorrGenMuEl = (abs(id) == 11 || abs(id) == 13) && gp.status() == 1 && statusflags.isPrompt();
	bool CorrGenNu   = (abs(id) == 12 || abs(id) == 14 || abs(id) == 16) && gp.status() == 1 && statusflags.isPrompt();
	bool CorrGenTau  = abs(id) == 15 && gp.status() == 2 && statusflags.isPrompt();

	// Only leptons
	if(CorrGenMuEl || CorrGenTau || CorrGenNu){

	  const reco::Candidate *it     = 0;
	  const reco::Candidate *mom    = 0;
	  const reco::Candidate *itmom  = 0;
	  const reco::Candidate *mommom = 0;

	  int momid = id;
	  it = (&gp);
	  // This loop searches the particle's mother
	  while(momid == id){
	    if(it != 0){
	      mom = it->mother();
	      if(mom != 0) momid = mom->pdgId();
	      else momid = 0;
	      if(momid == id) it = mom;
	    } // if(it != 0)
	    else momid = 0;
	  } // while(momid == id)

	  int mommomid = momid;

	  if(mom != 0){
	    itmom = mom;
	    // This loop searches the mother's mom of the particle
	    while (mommomid == momid){
	      if(itmom !=0){
		mommom = itmom->mother();
		if(mommom != 0) mommomid = mommom->pdgId();
		else mommomid = 0;
		if(mommomid == momid) itmom = mommom->mother();
	      }
	      else mommomid = 0;
	    } // if(mom != 0)
	  } // while(mommomid == momid)
	    
	  if (abs(momid) == 24 && abs(mommomid) == 6){

	    if (abs(id) == 13 || abs(id) == 11) Genlepton.SetPtEtaPhiE(gp.pt(), gp.eta(), gp.phi(), gp.energy());  

	    if (abs(id) == 12 || abs(id) == 14 || abs(id) == 16) Gennu.SetPtEtaPhiE(gp.pt(), gp.eta(), gp.phi(), gp.energy());  
	      
	    if (abs(id) == 15){ // Taus 
	      for(unsigned int h = 0; h <  gp.numberOfDaughters(); h++) {
		const reco::Candidate *gd = gp.daughter(h);
		const int taudauid = gd->pdgId();
		if (abs(taudauid) == 13 || abs(taudauid) == 11) Genlepton.SetPtEtaPhiE(gd->pt(), gd->eta(), gd->phi(), gd->energy());  
	      } // for(taus' daughters)
	    } // if(taus)
	        
	  } // if(t->W)
	    
	}// if (mu || e || tau)
      } // for(genParticles)

	// Full lepton and neutrino information
      b_GenNu_pT  = Gennu.Pt();
      b_GenNu_eta = Gennu.Eta();
      b_GenNu_phi = Gennu.Phi();
      b_GenNu_E   = Gennu.E();
	
      b_GenLepton_pT  = Genlepton.Pt();
      b_GenLepton_eta = Genlepton.Eta();
      b_GenLepton_phi = Genlepton.Phi();
      b_GenLepton_E   = Genlepton.E();

      edm::Handle<reco::GenJetCollection> genJets;
      iEvent.getByToken(genJetToken_, genJets);
	
      edm::Handle<vector<vector<int>>> JetMom;
      iEvent.getByToken(JetMotherToken_, JetMom);	
	
      for (unsigned int j = 0; j < genJets->size(); j++){
	const cat::GenJet & gjet = (*genJets)[j];
	if(std::abs(gjet.eta())< 2.5 &&
	   gjet.pt()> 20){
	    
	  b_GenJet_pT ->push_back(gjet.pt());
	  b_GenJet_eta->push_back(gjet.eta());
	  b_GenJet_phi->push_back(gjet.phi());
	  b_GenJet_E  ->push_back(gjet.energy());

	  b_GenJet_pT->push_back(gjet.pt());

	  std::vector<int> moms = (*JetMom)[j];	    
	  int WIDJet = 0, TopIDJet = 0; 
	  for(unsigned int nmom = 0; nmom < moms.size() ; nmom++){
	    int momID = moms.at(nmom);
	    if(std::abs(momID) == 24) WIDJet   = momID;
	    if(std::abs(momID) == 6 ) TopIDJet = momID;
	  }

	  if (WIDJet != 0) b_GenJet_mom->push_back(WIDJet);
	  else if (TopIDJet != 0)  b_GenJet_mom->push_back(TopIDJet);
	  else b_GenJet_mom->push_back(0);
 
	  b_GenJet_GenConeMom->push_back(0);
	}// if(Good genJet)
      }// for(genJet)	

      // GenTop Top, W and Add GEN-Jets 
      int GenConeMomID[6] = {6, 6, 24, 24, 0, 0};
      for (int ijGT=0; ijGT<6; ijGT++){
	int GenConeGenJetIndex = -999;
	float mGenConegJet_E = (*b_GenCone_gJet_E)[ijGT];
	auto itGenConeGenJet = std::find((*b_GenJet_E).begin(), (*b_GenJet_E).end(), mGenConegJet_E);
	if(itGenConeGenJet != (*b_GenJet_E).end()) GenConeGenJetIndex = std::distance((*b_GenJet_E).begin(), itGenConeGenJet);
	b_GenCone_gJetIndex->push_back(GenConeGenJetIndex);

	if(GenConeGenJetIndex != -999) (*b_GenJet_GenConeMom)[GenConeGenJetIndex] = GenConeMomID[ijGT];

      }// for(ijGT)

      if(b_GenJet_pT->size() >= SkimNJets_) gentree->Fill();

    }// if (MC && ttbarSignal && IsCat && nGenLep == 1)

  } // if(TTbarMC>0)


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Primary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  edm::Handle<int> pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );

  b_nGoodPV = *pvHandle;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Missing E_T
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  Handle<cat::METCollection> MET;
  iEvent.getByToken(metToken_, MET);

  // MET-PF
  b_MET     = MET->begin()->pt();
  b_MET_phi = MET->begin()->phi();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  cat::ElectronCollection selectedElectrons;
  cat::ElectronCollection vetoElectrons;

  Handle<cat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  for (unsigned int i = 0; i < electrons->size() ; i++) {
    const cat::Electron & electron = electrons->at(i);
    
    if( IsSelectElectron( electron ) ) selectedElectrons.push_back( electron );
    else if( IsVetoElectron( electron ) ) vetoElectrons.push_back( electron ); // does not Include selected electrons
  }
  
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  cat::MuonCollection selectedMuons;
  cat::MuonCollection vetoMuons;

  Handle<cat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (unsigned int i = 0; i < muons->size() ; i++) {
    const cat::Muon & muon = muons->at(i);
    
    if( IsSelectMuon( muon) ) selectedMuons.push_back( muon);
      else if( IsVetoMuon( muon) ) vetoMuons.push_back( muon); // does not Include selected muons
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  TLorentzVector lepton;

  float lepton_LES = 0.0;
  int ch_tag  = 999;

  if(selectedMuons.size()     == 1 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 0 &&
     vetoElectrons.size()     == 0){
    lepton.SetPtEtaPhiE(selectedMuons[0].pt(), selectedMuons[0].eta(), selectedMuons[0].phi(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets
    //LES
    if(isMC_) lepton_LES = selectedMuons[0].shiftedEn();
  }

  if(selectedMuons.size()     == 0 &&
     vetoMuons.size()         == 0 &&
     selectedElectrons.size() == 1 &&
     vetoElectrons.size()     == 0){
    ch_tag = 1; //electron + jets
    lepton.SetPtEtaPhiE(selectedElectrons[0].pt(), selectedElectrons[0].eta(), selectedElectrons[0].phi(), selectedElectrons[0].energy());
    //LES 
    if(isMC_) lepton_LES = selectedElectrons[0].shiftedEn();
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // HLTrigger
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  std::vector<string> triggerNameEl_, triggerNameMu_;
  if(isMC_){
    triggerNameEl_ = triggerNameMCEl_;
    triggerNameMu_ = triggerNameMCMu_;
  } 
  else{
    triggerNameEl_ = triggerNameDataEl_;
    triggerNameMu_ = triggerNameDataMu_;  
  }  
  
  bool EvTrigger = false; 
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerBits_, triggerBits);
  if( !iEvent.getByToken(triggerBits_, triggerBits) ){
    iEvent.getByToken(triggerBits2_, triggerBits);
  } 
  iEvent.getByToken(triggerObjects_, triggerObjects);
  const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
  AnalysisHelper trigHelper = AnalysisHelper(triggerNames, triggerBits, triggerObjects);

  bool IsTriggerMu = false;
  bool IsTriggerEl = false;
  
  // DEBUG: Print all triggers
  // trigHelper.listFiredTriggers();

  for(std::vector<string>::iterator TrMu_it = triggerNameMu_.begin(); TrMu_it != triggerNameMu_.end(); TrMu_it++){
    IsTriggerMu = trigHelper.triggerFired(*TrMu_it);
    // No trigger
    if(*TrMu_it == "notrigger") IsTriggerMu = true; 
    if (IsTriggerMu) break;
  }
  
  for(std::vector<string>::iterator TrEl_it = triggerNameEl_.begin(); TrEl_it != triggerNameEl_.end(); TrEl_it++){
    IsTriggerEl = trigHelper.triggerFired(*TrEl_it);
    // No trigger
    if(*TrEl_it == "notrigger") IsTriggerEl = true;
    if (IsTriggerEl) break;
  }

  if ( (ch_tag == 0 && IsTriggerMu) ||
       (ch_tag == 1 && IsTriggerEl) 
       ) {
    EvTrigger = true;
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Fill Tree with events that have ONLY one lepton
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Check Gen Level for ttbar sample
  if (TTbarMC_ > 0){
    // Signal ttbar event
    if(TTbarMC_ == 1 && nGenLep != 1) ch_tag = 999;
    // Background ttbar event
    if(TTbarMC_ == 2 && nGenLep == 1) ch_tag = 999;
    // ttbar Category
    if(!IsCat) ch_tag = 999;
  } // if(TTbarMC_ >0)

  if (ch_tag<2 && EvTrigger){ // Single lepton event

    b_Channel  = ch_tag;

    b_Lepton_pT  = lepton.Pt();
    b_Lepton_eta = lepton.Eta();
    b_Lepton_phi = lepton.Phi();
    b_Lepton_E   = lepton.E();

    b_Lepton_LES = lepton_LES;

    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Jets
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    
    Handle<cat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    // Arrange Jets by CSV discriminator
    std::vector<int>  JetIndex;
    for (unsigned int i = 0; i < jets->size() ; i++) JetIndex.push_back(i);
    for(unsigned int ijet=0; ijet < JetIndex.size(); ijet++){
      for(unsigned int jjet=ijet+1; jjet < JetIndex.size(); jjet++){
        const cat::Jet & jet_i = jets->at(JetIndex[ijet]);
        const cat::Jet & jet_j = jets->at(JetIndex[jjet]);

        float iJetCSV = jet_i.bDiscriminator(BTAG_CSVv2);
        float jJetCSV = jet_j.bDiscriminator(BTAG_CSVv2);

        if(jJetCSV > iJetCSV){
          float tempIndex = JetIndex[ijet];
          JetIndex[ijet] = JetIndex[jjet];
          JetIndex[jjet] = tempIndex;
        }//if(jJetCSV > iJetCSV)
      }
    }

    // Initialize SF_btag
    // Jet_SF_CSV[SystVariations];
    float Jet_SF_CSV[19], Jet_SF_CSVcsv[19];
    for (unsigned int iu=0; iu<19; iu++){ Jet_SF_CSV[iu] = 1.0; Jet_SF_CSVcsv[iu] = 1.0;}

    // Central Value
    Handle<float> rSF_CSV;
    iEvent.getByToken(CSVWeightToken_, rSF_CSV);
    Jet_SF_CSV[0] = *rSF_CSV;

    // Full set of systematics variations estimated for jets with pT>30
    Handle<std::vector<float>> rsysSF_CSV;
    iEvent.getByToken(CSVSysWeightToken_, rsysSF_CSV);
    for (unsigned int icsv = 0; icsv < rsysSF_CSV->size() ; icsv++)  Jet_SF_CSV[icsv+1] = rsysSF_CSV->at(icsv); 

    std::vector<cat::Jet> selectedJets;
    // Run again over all Jets (CSV order)
    for (unsigned int i = 0; i < JetIndex.size() ; i++) {

      const cat::Jet & jet = jets->at(JetIndex[i]);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection (pT>20GeV to take into account SYST Variations)
      if(std::abs(jet.eta()) < 2.4 && jet.pt() > 20. && jet.LooseId()) goodJet = true;
      // Jet Cleaning
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if(dr_LepJet > 0.4) cleanJet = true;

      if(goodJet && cleanJet){

        selectedJets.push_back( jet );
        // Basic variables
        b_Jet_pT ->push_back(jet.pt());
        b_Jet_eta->push_back(jet.eta());
        b_Jet_phi->push_back(jet.phi());
        b_Jet_E  ->push_back(jet.energy());
        b_Jet_Index ->push_back(JetIndex[i]);	

        // Parton Flavour
        b_Jet_partonFlavour->push_back(jet.partonFlavour());
        b_Jet_hadronFlavour->push_back(jet.hadronFlavour());

        // b-tag discriminant
        float jet_btagDis_CSV = jet.bDiscriminator(BTAG_CSVv2);
        b_Jet_CSV ->push_back(jet_btagDis_CSV);
        // c-tag discriminant
        float jet_btagDis_CvsL = jet.bDiscriminator(CTAG_CvsL);
        b_Jet_CvsL ->push_back(jet_btagDis_CvsL);
        float jet_btagDis_CvsB = jet.bDiscriminator(CTAG_CvsB);
	b_Jet_CvsB ->push_back(jet_btagDis_CvsB);

        if(isMC_) {
          // JES
          b_Jet_JES_Up  ->push_back(jet.shiftedEnUp());
          b_Jet_JES_Down->push_back(jet.shiftedEnDown());

	  // Full set of components
	  std::vector<float> jesCom_Up, jesCom_Down;
	  for (int isrc = 0; isrc < nsrc; isrc++) {
	    const char *name = srcnames[isrc];
	    JetCorrectorParameters *p = new JetCorrectorParameters(edm::FileInPath("CATTools/CatProducer/data/JEC/" + JESUncFile_).fullPath(), name);
	    JetCorrectionUncertainty *jetunc = new JetCorrectionUncertainty(*p);
	    jetunc->setJetPt(jet.pt());
	    jetunc->setJetEta(jet.eta());
	    float jesUp = jetunc->getUncertainty(true);
	    jesCom_Up.push_back(jesUp);
	    jetunc->setJetPt(jet.pt());
	    jetunc->setJetEta(jet.eta());
	    float jesDown = jetunc->getUncertainty(false);
	    jesCom_Down.push_back(jesDown);

	    delete p;
	    delete jetunc;
	  } // for isrc

	  b_Jet_JESCom_Up->push_back(jesCom_Up);
	  b_Jet_JESCom_Down->push_back(jesCom_Down);
	  
          // JER
          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
          b_Jet_JER_Up   ->push_back(jet.smearedResUp());
          b_Jet_JER_Nom  ->push_back(jet.smearedRes());
          b_Jet_JER_Down ->push_back(jet.smearedResDown());

          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
          // Saving the central SF and the 18 syst. unc. for:
	  // pT_Jets > 30 GeV
          if(jet.pt() > 30.) for (unsigned int iu=0; iu<19; iu++) Jet_SF_CSVcsv[iu] *= SF_CSV_.getSF(jet, iu);

	  // GEN-Jets matched with RECO-Jets
	  // Reference to gen object
	  if(TTbarMC_== 1){
	    int MatchedGenJetIndex = -999;
	    int rJetGenConeMom = -999;
	    if(jet.genJet()){ 
	      float mGenJet_E = jet.genJet()->energy();
	      auto itGenJet = std::find((*b_GenJet_E).begin(), (*b_GenJet_E).end(), mGenJet_E);
	      if(itGenJet != (*b_GenJet_E).end()) MatchedGenJetIndex = std::distance((*b_GenJet_E).begin(), itGenJet);

	      if(MatchedGenJetIndex != -999){
		if      (MatchedGenJetIndex == (*b_GenCone_gJetIndex)[0] ||
		         MatchedGenJetIndex == (*b_GenCone_gJetIndex)[1]) rJetGenConeMom = 6;
		else if (MatchedGenJetIndex == (*b_GenCone_gJetIndex)[2] ||
			 MatchedGenJetIndex == (*b_GenCone_gJetIndex)[3]) rJetGenConeMom = 24;
		else if (MatchedGenJetIndex == (*b_GenCone_gJetIndex)[4] ||
			 MatchedGenJetIndex == (*b_GenCone_gJetIndex)[5]) rJetGenConeMom = 0;
	      } // if(MatchedGenJetIndex != -999) 
	    } // if(jet.genJet())

	    b_Jet_MatchedGenJetIndex->push_back(MatchedGenJetIndex);

	    b_Jet_GenConeMom->push_back(rJetGenConeMom);
	    
	  } // if(TTbarMC_== 1)
        } // if(isMC_)	
      }// if(GoodJets)
    }// for(AllJets)

    // CSV SF and sys unc
    for (unsigned int iu=0; iu<19; iu++) {b_Jet_SF_CSV->push_back(1.0); b_Jet_SF_CSVcsv->push_back(1.0);}
    (*b_Jet_SF_CSV)[0] = Jet_SF_CSV[0]; //Central
    (*b_Jet_SF_CSVcsv)[0] = Jet_SF_CSVcsv[0]; //Central

    // To save only the error
    for (unsigned int iu=1; iu<19; iu++){
      (*b_Jet_SF_CSV)[iu] = std::abs(Jet_SF_CSV[iu] - Jet_SF_CSV[0]) ; // Syst. Unc.
      (*b_Jet_SF_CSVcsv)[iu] = std::abs(Jet_SF_CSVcsv[iu] - Jet_SF_CSVcsv[0]) ; // Syst. Unc.
    }
    //---------------------------------------------------------------------------
    // Kinematic Reconstruction
    //---------------------------------------------------------------------------
    TLorentzVector Kinnu, Kinblrefit, Kinbjrefit, Kinj1refit, Kinj2refit;
    Kinnu.SetPtEtaPhiE(0,0,0,0);
    Kinblrefit.SetPtEtaPhiE(0,0,0,0);
    Kinbjrefit.SetPtEtaPhiE(0,0,0,0);
    Kinj1refit.SetPtEtaPhiE(0,0,0,0);
    Kinj2refit.SetPtEtaPhiE(0,0,0,0);
    
    std::vector<int> KinBestIndices;
    KinBestIndices.push_back(-999);
    KinBestIndices.push_back(-999);
    KinBestIndices.push_back(-999);
    KinBestIndices.push_back(-999);
    float bestchi2 = 0;
    
    if(b_Jet_pT->size() > 3){
      
      TLorentzVector KinLep;
      KinLep = lepton;
      
      TLorentzVector KinMET;
      KinMET.SetPtEtaPhiE(b_MET, 0.0, b_MET_phi, b_MET);
      
      std::vector<ComJet> KinJets;
      for (unsigned int kj=0; kj<b_Jet_pT->size(); kj++){
	ComJet kjet;
	kjet.SetPtEtaPhiE((*b_Jet_pT)[kj],(*b_Jet_eta)[kj],(*b_Jet_phi)[kj],(*b_Jet_E)[kj]);
	kjet.CSV = (*b_Jet_CSV)[kj];
	
	KinJets.push_back(kjet);
      }
      
      ttbb::FindHadronicTop(KinLep, KinJets, KinMET, KFUsebtag_, CSVPosConKF_, KinBestIndices, bestchi2, Kinnu, Kinblrefit, Kinbjrefit, Kinj1refit, Kinj2refit);
      
    } //if(NJets > 3)   
    
    // LorentzVector for Jets. Same order as KinBestIndices
    b_KinJet_pT->push_back(Kinbjrefit.Pt());
    b_KinJet_pT->push_back(Kinj1refit.Pt());
    b_KinJet_pT->push_back(Kinj2refit.Pt()); 
    b_KinJet_pT->push_back(Kinblrefit.Pt());
    
    b_KinJet_eta->push_back(Kinbjrefit.Eta());
    b_KinJet_eta->push_back(Kinj1refit.Eta());
    b_KinJet_eta->push_back(Kinj2refit.Eta()); 
    b_KinJet_eta->push_back(Kinblrefit.Eta());
    
    b_KinJet_phi->push_back(Kinbjrefit.Phi());
    b_KinJet_phi->push_back(Kinj1refit.Phi());
    b_KinJet_phi->push_back(Kinj2refit.Phi()); 
    b_KinJet_phi->push_back(Kinblrefit.Phi());
    
    b_KinJet_E->push_back(Kinbjrefit.E());
    b_KinJet_E->push_back(Kinj1refit.E());
    b_KinJet_E->push_back(Kinj2refit.E()); 
    b_KinJet_E->push_back(Kinblrefit.E());
    
    b_KinNu_pT  = Kinnu.Pt();
    b_KinNu_eta = Kinnu.Eta();
    b_KinNu_phi = Kinnu.Phi();
    b_KinNu_E   = Kinnu.E();
    
    for(unsigned int iki=0; iki<KinBestIndices.size(); iki++) b_KinJet_Index->push_back(KinBestIndices.at(iki));
    
    b_Kin_Chi2 = bestchi2;
        
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    // Fill Tree with event at 1 lepton cut level
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------
    
    if(b_Jet_pT->size() >= SkimNJets_) tree->Fill();
        
  }// if(ch_tag)
  
}

//------------- Good Muon Selection -----------------------
bool ttbbAnalyzer::IsSelectMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Tight selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.isTightMuon());

  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 20);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.1); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.15 );

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodMuon;
}
//------------- Loose Muon Selection -----------------------
bool ttbbAnalyzer::IsVetoMuon(const cat::Muon & i_muon_candidate)
{
  bool GoodMuon=true;

  // Loose selection already defined into CAT::Muon
  GoodMuon &= (i_muon_candidate.isLooseMuon());

  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 15);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes PU subtraction
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - 0.5*PUIso) )/ ecalpt;

  GoodMuon &=( i_muon_candidate.relIso( 0.4 ) < 0.25 );

  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool ttbbAnalyzer::IsSelectElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.1);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // Electron cut based selection
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Summer16-80X-V1-medium") > 0.0;
  
  // Electron MVA selection (Tight: WP80)
  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_for_7_4_12_Spring15_MVA
  // if (std::abs(i_electron_candidate.eta()) < 0.8)        GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.988153; // Barrel
  // else if (std::abs(i_electron_candidate.eta()) > 0.8 &&
  // 	   std::abs(i_electron_candidate.eta()) < 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.967910; // Barrel
  // else if (std::abs(i_electron_candidate.eta()) > 1.479) GoodElectron &= i_electron_candidate.electronID("mvaEleID-Spring15-25ns-Trig-V1-wp80") > 0.841729; // EndCaps

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;

  // Isolation is already included in the cut-based cuts, it is not needed 
  // if ( !doLooseLepton_ ) {
  //   if ( std::abs(i_electron_candidate.scEta()) <= 1.479)   GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.0695 );
  //   else GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.0821 );
  // }

  // Effective Area Parametrization can be found in:
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}
//------------- Loose Electron Selection -----------------------
bool ttbbAnalyzer::IsVetoElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Summer16-80X-V1-veto") > 0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // if ( std::abs(i_electron_candidate.scEta()) <= 1.479) GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.175 );
  // else GoodElectron &= ( i_electron_candidate.relIso( 0.3 ) < 0.159 );
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbbAnalyzer);

