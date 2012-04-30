#ifndef NtupleWriter_h
#define NtupleWriter_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "UHHAnalysis/NtupleWriter/interface/Objects.h"

//
// class declaration
//



//allow Nmax different electron, muon, jet collections
const int Nmax=12;

class NtupleWriter : public edm::EDAnalyzer {
   public:
      explicit NtupleWriter(const edm::ParameterSet&);
      ~NtupleWriter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      TFile *outfile;
      TTree *tr;
      TString fileName;

      bool doElectrons;
      bool doMuons;
      bool doTaus;
      bool doJets;
      bool doTopJets;
      bool doGenTopJets;
      bool doMET;
      bool doPhotons;
      bool doGenInfo;
      bool doAllGenParticles;
      bool doLumiInfo;
      bool doPV;
      bool doTrigger;

      unsigned int run;
      unsigned int luminosityBlock;
      int event;
      bool isRealData;
      bool HBHENoiseFilterResult;

      float intgDelLumi;
      float intgRecLumi;
      float totalDelLumi;
      float totalRecLumi;

      std::vector<std::string> electron_sources;
      std::vector<Electron> eles[Nmax];

      std::vector<std::string> muon_sources;
      std::vector<Muon> mus[Nmax];

      std::vector<std::string> tau_sources;
      std::vector<Tau> taus[Nmax];
      double tau_ptmin;
      double tau_etamax;

      std::vector<std::string> jet_sources;
      std::vector<Jet> jets[Nmax];
      double jet_ptmin;
      double jet_etamax;

      std::vector<std::string> topjet_sources;
      std::vector<TopJet> topjets[Nmax];
      double topjet_ptmin;
      double topjet_etamax;

      std::vector<std::string> gentopjet_sources;
      std::vector<TopJet> gentopjets[Nmax];
      double gentopjet_ptmin;
      double gentopjet_etamax;

      std::vector<std::string> photon_sources;
      std::vector<Photon> phs[Nmax];

      std::vector<std::string> met_sources;
      MET met[Nmax];

      std::vector<std::string> pv_sources;
      std::vector<PrimaryVertex> pvs[Nmax];
      
      float beamspot_x0;
      float beamspot_y0;
      float beamspot_z0;

      GenInfo genInfo;
      std::vector<GenParticle> genps;

      std::vector<std::string> trigger_prefixes;
      //std::map<std::string, bool> triggerResults;
      std::vector<std::string> triggerNames;
      std::vector<bool> triggerResults;
/*       std::vector<int> L1_prescale; */
/*       std::vector<int> HLT_prescale; */
      
      HLTConfigProvider hlt_cfg;
      bool newrun;
      bool previouslumiblockwasfilled;

      JetCorrectionUncertainty *jecUnc;
};


#endif
