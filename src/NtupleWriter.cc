// -*- C++ -*-
//
// Package:    NtupleWriter
// Class:      NtupleWriter
// 
/**\class NtupleWriter NtupleWriter.cc NtupleWriter/src/NtupleWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Peiffer,,,Uni Hamburg
//         Created:  Tue Mar 13 08:43:34 CET 2012
// $Id: NtupleWriter.cc,v 1.6 2012/04/11 15:15:44 peiffer Exp $
//
//

#include "UHHAnalysis//NtupleWriter/interface/NtupleWriter.h"



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtupleWriter::NtupleWriter(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  //edm::Service<TFileService> fs;
  //tr = fs->make<TTree>("AnalysisTree","AnalysisTree");

  fileName =  iConfig.getParameter<std::string>("fileName");
  outfile = new TFile(fileName,"RECREATE");
  outfile->cd();
  tr = new TTree("AnalysisTree","AnalysisTree");

  std::string name;

  doElectrons = iConfig.getParameter<bool>("doElectrons");
  doMuons = iConfig.getParameter<bool>("doMuons");
  doTaus = iConfig.getParameter<bool>("doTaus"); 
  doJets = iConfig.getParameter<bool>("doJets");
  doPhotons = iConfig.getParameter<bool>("doPhotons");
  doMET = iConfig.getParameter<bool>("doMET");
  doGenInfo = iConfig.getParameter<bool>("doGenInfo");
  doPV = iConfig.getParameter<bool>("doPV");
  doTopJets = iConfig.getParameter<bool>("doTopJets");
  doTrigger = iConfig.getParameter<bool>("doTrigger");

  // initialization of tree variables

  tr->Branch("run",&run);
  tr->Branch("event",&event);
  tr->Branch("luminosityBlock",&luminosityBlock);
  tr->Branch("isRealData",&isRealData);
  tr->Branch("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  tr->Branch("beamspot_x0",&beamspot_x0);
  tr->Branch("beamspot_y0",&beamspot_y0);
  tr->Branch("beamspot_z0",&beamspot_z0);

  if(doElectrons){
    electron_sources = iConfig.getParameter<std::vector<std::string> >("electron_sources");
    for(size_t j=0; j< electron_sources.size(); ++j){  
      tr->Branch( electron_sources[j].c_str(), "std::vector<Electron>", &eles[j]); 
    }
  }
  if(doMuons){
    muon_sources = iConfig.getParameter<std::vector<std::string> >("muon_sources");
    for(size_t j=0; j< muon_sources.size(); ++j){  
      tr->Branch( muon_sources[j].c_str(), "std::vector<Muon>", &mus[j]);
    }
  }
  if(doTaus){
    tau_sources = iConfig.getParameter<std::vector<std::string> >("tau_sources");
    tau_ptmin = iConfig.getParameter<double> ("tau_ptmin");
    tau_etamax = iConfig.getParameter<double> ("tau_etamax");
    for(size_t j=0; j< tau_sources.size(); ++j){  
      tr->Branch( tau_sources[j].c_str(), "std::vector<Tau>", &taus[j]);
    }
  }
  if(doJets){
    jet_sources = iConfig.getParameter<std::vector<std::string> >("jet_sources");
    jet_ptmin = iConfig.getParameter<double> ("jet_ptmin");
    jet_etamax = iConfig.getParameter<double> ("jet_etamax");
    for(size_t j=0; j< jet_sources.size(); ++j){  
      tr->Branch( jet_sources[j].c_str(), "std::vector<Jet>", &jets[j]);
    }
  }
  if(doTopJets){
    topjet_sources = iConfig.getParameter<std::vector<std::string> >("topjet_sources");
    topjet_ptmin = iConfig.getParameter<double> ("topjet_ptmin");
    topjet_etamax = iConfig.getParameter<double> ("topjet_etamax");
    for(size_t j=0; j< topjet_sources.size(); ++j){  
      tr->Branch( topjet_sources[j].c_str(), "std::vector<TopJet>", &topjets[j]);
    }
  }
  if(doPhotons){
    photon_sources = iConfig.getParameter<std::vector<std::string> >("photon_sources");
    for(size_t j=0; j< photon_sources.size(); ++j){  
      tr->Branch( photon_sources[j].c_str(), "std::vector<Photon>", &phs[j]);
    }
  }
  if(doMET){
    met_sources = iConfig.getParameter<std::vector<std::string> >("met_sources");
    for(size_t j=0; j< met_sources.size(); ++j){  
      tr->Branch(met_sources[j].c_str(), "MET", &met[j]);
    }
  }
  if(doPV){
    pv_sources = iConfig.getParameter<std::vector<std::string> >("pv_sources");
    for(size_t j=0; j< pv_sources.size(); ++j){  
      tr->Branch( pv_sources[j].c_str(), "std::vector<PrimaryVertex>", &pvs[j]);
    }
  }
  if(doGenInfo){
    tr->Branch("genInfo","GenInfo",&genInfo);
    tr->Branch("GenParticles","std::vector<GenParticle>", &genps);
  }
  if(doTrigger){
    trigger_prefixes = iConfig.getParameter<std::vector<std::string> >("trigger_prefixes");
    //tr->Branch("triggerResults","std::map<std::string, bool>",&triggerResults);
    tr->Branch("triggerNames", "std::vector<std::string>", &triggerNames);  
    tr->Branch("triggerResults", "std::vector<bool>", &triggerResults);
    tr->Branch("L1_prescale", "std::vector<int>", &L1_prescale);
    tr->Branch("HLT_prescale", "std::vector<int>", &HLT_prescale);
  }
  newrun = true;
}


NtupleWriter::~NtupleWriter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
NtupleWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //   using namespace edm;


   // ------------- common variables ------------- 
   
   run = iEvent.id().run();
   event = iEvent.id().event();
   luminosityBlock = iEvent.luminosityBlock();
   isRealData      = iEvent.isRealData();

   if(isRealData){
     edm::Handle<bool> bool_handle;
     iEvent.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),bool_handle);
     HBHENoiseFilterResult = *(bool_handle.product());
   }
   else HBHENoiseFilterResult = false;

   // ------------- primary vertices and beamspot  -------------

   if(doPV){
     for(size_t j=0; j< pv_sources.size(); ++j){
       pvs[j].clear();
       
       edm::Handle< std::vector<reco::Vertex> > pv_handle;
       iEvent.getByLabel(pv_sources[j], pv_handle);
       const std::vector<reco::Vertex>& reco_pvs = *(pv_handle.product());
       
       for (unsigned int i = 0; i <  reco_pvs.size(); ++i) {
	 reco::Vertex reco_pv = reco_pvs[i];

	 PrimaryVertex pv;
	 pv.x =  reco_pv.x();
	 pv.y =  reco_pv.y();
	 pv.z =  reco_pv.z();
	 pv.nTracks =  reco_pv.nTracks();
	 //pv.isValid =  reco_pv.isValid();
	 pv.chi2 =  reco_pv.chi2();
	 pv.ndof =  reco_pv.ndof();	 

	 pvs[j].push_back(pv);
       }
     }
   }
   
   edm::Handle<reco::BeamSpot> beamSpot;
   iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpot);
   const reco::BeamSpot & bsp = *beamSpot;
   
   beamspot_x0 = bsp.x0();
   beamspot_y0 = bsp.y0();
   beamspot_z0 = bsp.z0();

   // ------------- electrons -------------   
   if(doElectrons){
     for(size_t j=0; j< electron_sources.size(); ++j){
       eles[j].clear();
       edm::Handle< std::vector<pat::Electron> > ele_handle;
       iEvent.getByLabel(electron_sources[j], ele_handle);
       const std::vector<pat::Electron>& pat_electrons = *(ele_handle.product());
       
       for (unsigned int i = 0; i < pat_electrons.size(); ++i) {
	 pat::Electron pat_ele = pat_electrons[i];
	 Electron ele;
	 
	 ele.charge =  pat_ele.charge();
	 ele.pt =  pat_ele.pt();
	 ele.eta =  pat_ele.eta();
	 ele.phi =  pat_ele.phi();
	 ele.energy =  pat_ele.energy();
	 ele.vertex_x = pat_ele.vertex().x();
	 ele.vertex_y = pat_ele.vertex().y();
	 ele.vertex_z = pat_ele.vertex().z();
	 ele.supercluster_eta = pat_ele.superCluster()->eta();
	 ele.supercluster_phi = pat_ele.superCluster()->phi();
	 ele.dB = pat_ele.dB();
	 //ele.particleIso = pat_ele.particleIso();
	 ele.neutralHadronIso = pat_ele.neutralHadronIso();
	 ele.chargedHadronIso = pat_ele.chargedHadronIso();
	 ele.trackIso = pat_ele.trackIso();
	 ele.puChargedHadronIso = pat_ele.puChargedHadronIso();
	 ele.gsfTrack_trackerExpectedHitsInner_numberOfLostHits = pat_ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
	 ele.gsfTrack_px= pat_ele.gsfTrack()->px();
	 ele.gsfTrack_py= pat_ele.gsfTrack()->py();
	 ele.gsfTrack_pz= pat_ele.gsfTrack()->pz();
	 ele.gsfTrack_vx= pat_ele.gsfTrack()->vx();
	 ele.gsfTrack_vy= pat_ele.gsfTrack()->vy();
	 ele.gsfTrack_vz= pat_ele.gsfTrack()->vz();
	 eles[j].push_back(ele);
       }
     }
   }

   // ------------- muons ------------- 
   if(doMuons){
     for(size_t j=0; j< muon_sources.size(); ++j){
       mus[j].clear();
       
       edm::Handle< std::vector<pat::Muon> > mu_handle;
       iEvent.getByLabel(muon_sources[j], mu_handle);
       const std::vector<pat::Muon>& pat_muons = *(mu_handle.product());
       
       for (unsigned int i = 0; i <pat_muons.size() ; ++i) {
	 pat::Muon pat_mu = pat_muons[i];

	 Muon mu;
	 mu.charge =  pat_mu.charge();
	 mu.pt =  pat_mu.pt();
	 mu.eta =  pat_mu.eta();
	 mu.phi =  pat_mu.phi();
	 mu.energy =  pat_mu.energy();
	 mu.vertex_x = pat_mu.vertex().x();
	 mu.vertex_y = pat_mu.vertex().y();
	 mu.vertex_z = pat_mu.vertex().z();
	 mu.dB = pat_mu.dB();
	 //mu.particleIso = pat_mu.particleIso();
	 mu.neutralHadronIso = pat_mu.neutralHadronIso();
	 mu.chargedHadronIso = pat_mu.chargedHadronIso();
	 mu.trackIso = pat_mu.trackIso();
	 mu.puChargedHadronIso = pat_mu.puChargedHadronIso();
	 mu.isGlobalMuon = pat_mu.isGlobalMuon();
	 mu.isStandAloneMuon = pat_mu.isStandAloneMuon();
	 mu.isTrackerMuon = pat_mu.isTrackerMuon();
	 mu.numberOfMatchedStations = pat_mu.numberOfMatchedStations();
	 reco::TrackRef globalTrack = pat_mu.globalTrack();
	 if(!globalTrack.isNull()){
	   mu.globalTrack_chi2 = globalTrack->chi2();
	   mu.globalTrack_ndof = globalTrack->ndof();
	   mu.globalTrack_d0 = globalTrack->d0();	 
	   mu.globalTrack_d0Error = globalTrack->d0Error();
	   mu.globalTrack_numberOfValidHits = globalTrack->numberOfValidHits();
	   mu.globalTrack_numberOfLostHits = globalTrack->numberOfLostHits();
	 }
	 else{
	   mu.globalTrack_chi2 = 0;
	   mu.globalTrack_ndof = 0;
	   mu.globalTrack_d0 = 0;
	   mu.globalTrack_d0Error = 0;
	   mu.globalTrack_numberOfValidHits = 0;
	   mu.globalTrack_numberOfLostHits = 0;
	 }
	 reco::TrackRef innerTrack = pat_mu.innerTrack();
	 if(!innerTrack.isNull()){
	   mu.innerTrack_chi2 = innerTrack->chi2();
	   mu.innerTrack_ndof = innerTrack->ndof();
	   mu.innerTrack_d0 = innerTrack->d0();	 
	   mu.innerTrack_d0Error = innerTrack->d0Error();
	   mu.innerTrack_numberOfValidHits = innerTrack->numberOfValidHits();
	   mu.innerTrack_numberOfLostHits = innerTrack->numberOfLostHits();
	 }
	 else{
	   mu.innerTrack_chi2 = 0;
	   mu.innerTrack_ndof = 0;
	   mu.innerTrack_d0 = 0;
	   mu.innerTrack_d0Error = 0;
	   mu.innerTrack_numberOfValidHits = 0;
	   mu.innerTrack_numberOfLostHits = 0;
	 }
	 reco::TrackRef outerTrack = pat_mu.outerTrack();
	 if(!outerTrack.isNull()){
	   mu.outerTrack_chi2 = outerTrack->chi2();
	   mu.outerTrack_ndof = outerTrack->ndof();
	   mu.outerTrack_d0 = outerTrack->d0();	 
	   mu.outerTrack_d0Error = outerTrack->d0Error();
	   mu.outerTrack_numberOfValidHits = outerTrack->numberOfValidHits();
	   mu.outerTrack_numberOfLostHits = outerTrack->numberOfLostHits();
	 } 
	 else{
	   mu.outerTrack_chi2 = 0;
	   mu.outerTrack_ndof = 0;
	   mu.outerTrack_d0 = 0;
	   mu.outerTrack_d0Error = 0;
	   mu.outerTrack_numberOfValidHits = 0;
	   mu.outerTrack_numberOfLostHits = 0;
	 }

	 mus[j].push_back(mu);
       }
     }
   }
   // ------------- taus ------------- 

   if(doTaus){
     for(size_t j=0; j< tau_sources.size(); ++j){
       taus[j].clear();

       
       edm::Handle< std::vector<pat::Tau> > tau_handle;
       iEvent.getByLabel(tau_sources[j], tau_handle);
       const std::vector<pat::Tau>& pat_taus = *(tau_handle.product());
       
       for (unsigned int i = 0; i < pat_taus.size(); ++i) {
	 pat::Tau pat_tau = pat_taus[i];
	 if(pat_tau.pt() < tau_ptmin) continue;
	 if(fabs(pat_tau.eta()) > tau_etamax) continue;

	 Tau tau;
	 tau.charge =  pat_tau.charge();
	 tau.pt =  pat_tau.pt();
	 tau.eta =  pat_tau.eta();
	 tau.phi =  pat_tau.phi();
	 tau.energy =  pat_tau.energy();

	 taus[j].push_back(tau);
       }
     }
   }

   // ------------- jets -------------
   if(doJets){
     for(size_t j=0; j< jet_sources.size(); ++j){
             
       jets[j].clear();

       edm::Handle< std::vector<pat::Jet> > jet_handle;
       iEvent.getByLabel(jet_sources[j], jet_handle);
       const std::vector<pat::Jet>& pat_jets = *(jet_handle.product());
  
       for (unsigned int i = 0; i < pat_jets.size(); ++i) {
	 pat::Jet pat_jet = pat_jets[i];
	 if(pat_jet.pt() < jet_ptmin) continue;
	 if(fabs(pat_jet.eta()) > jet_etamax) continue;
// 	 std::cout << "available btag discriminators: " << std::endl;
// 	 const std::vector<std::pair<std::string, float> > & dis = pat_jets[i].getPairDiscri();
// 	 for(size_t k=0; k<dis.size(); ++k){
// 	   std::cout << dis[k].first << std::endl;
// 	 }

	 Jet jet;
	 jet.charge = pat_jet.charge();
	 jet.pt = pat_jet.pt();
	 jet.eta = pat_jet.eta();
	 jet.phi = pat_jet.phi();
	 jet.energy = pat_jet.energy();
	 jet.numberOfDaughters =pat_jet.numberOfDaughters();
	 const reco::TrackRefVector&  jettracks = pat_jet.associatedTracks();
	 jet.nTracks = jettracks.size();
	 jet.jetArea = pat_jet.jetArea();
	 jet.pileup = pat_jet.pileup();
	 if(pat_jet.isPFJet()){
	   jet.neutralEmEnergyFraction =pat_jet.neutralEmEnergyFraction();
	   jet.neutralHadronEnergyFraction =pat_jet.neutralHadronEnergyFraction();
	   jet.chargedEmEnergyFraction =pat_jet.chargedEmEnergyFraction();
	   jet.chargedHadronEnergyFraction =pat_jet.chargedHadronEnergyFraction();
	   jet.muonEnergyFraction =pat_jet.muonEnergyFraction();
	   jet.photonEnergyFraction =pat_jet.photonEnergyFraction();
	   jet.chargedMultiplicity =pat_jet.chargedMultiplicity();
	   jet.neutralMultiplicity =pat_jet.neutralMultiplicity();
	   jet.muonMultiplicity =pat_jet.muonMultiplicity();
	   jet.electronMultiplicity =pat_jet.electronMultiplicity();
	   jet.photonMultiplicity =pat_jet.photonMultiplicity();
	 }

	 jecUnc->setJetEta(pat_jet.eta());
	 jecUnc->setJetPt(pat_jet.pt());
	 jet.JEC_uncertainty = jecUnc->getUncertainty(true);

	 jet.btag_simpleSecondaryVertexHighEff=pat_jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	 jet.btag_simpleSecondaryVertexHighPur=pat_jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	 jet.btag_combinedSecondaryVertex=pat_jet.bDiscriminator("combinedSecondaryVertexBJetTags");
	 jet.btag_combinedSecondaryVertexMVA=pat_jet.bDiscriminator("combinedSecondaryVertexMVABJetTags");
	 jet.btag_jetBProbability=pat_jet.bDiscriminator("jetBProbabilityBJetTags");
	 jet.btag_jetProbability=pat_jet.bDiscriminator("jetProbabilityBJetTags");

	 jets[j].push_back(jet);
       }
     }
   }

   // ------------- top jets -------------
   if(doTopJets){
     for(size_t j=0; j< topjet_sources.size(); ++j){
       
       topjets[j].clear();

       edm::Handle<pat::JetCollection> pat_topjets;
       //edm::Handle<std::vector<reco::Jet> > pat_topjets;
       iEvent.getByLabel(topjet_sources[j], pat_topjets);

       for (unsigned int i = 0; i < pat_topjets->size(); i++) {
	 const pat::Jet  pat_topjet =  * dynamic_cast<pat::Jet const *>(&pat_topjets->at(i));
	 if(pat_topjet.pt() < topjet_ptmin) continue;
	 if(fabs(pat_topjet.eta()) > topjet_etamax) continue;

	 TopJet topjet;
	 topjet.charge = pat_topjet.charge();
	 topjet.pt = pat_topjet.pt();
	 topjet.eta = pat_topjet.eta();
	 topjet.phi = pat_topjet.phi();
	 topjet.energy = pat_topjet.energy();
	 topjet.numberOfDaughters =pat_topjet.numberOfDaughters();
	 const reco::TrackRefVector&  topjettracks = pat_topjet.associatedTracks();
	 topjet.nTracks = topjettracks.size();
	 topjet.jetArea = pat_topjet.jetArea();
	 topjet.pileup = pat_topjet.pileup();
//  	 topjet.neutralEmEnergyFraction =pat_topjet.neutralEmEnergyFraction();
//  	 topjet.neutralHadronEnergyFraction =pat_topjet.neutralHadronEnergyFraction();
//  	 topjet.chargedEmEnergyFraction =pat_topjet.chargedEmEnergyFraction();
//  	 topjet.chargedHadronEnergyFraction =pat_topjet.chargedHadronEnergyFraction();
//  	 topjet.muonEnergyFraction =pat_topjet.muonEnergyFraction();
//  	 topjet.photonEnergyFraction =pat_topjet.photonEnergyFraction();
// 	 topjet.chargedMultiplicity =pat_topjet.chargedMultiplicity();
// 	 topjet.neutralMultiplicity =pat_topjet.neutralMultiplicity();
// 	 topjet.muonMultiplicity =pat_topjet.muonMultiplicity();
// 	 topjet.electronMultiplicity =pat_topjet.electronMultiplicity();
// 	 topjet.photonMultiplicity =pat_topjet.photonMultiplicity();

	 jecUnc->setJetEta(pat_topjet.eta());
	 jecUnc->setJetPt(pat_topjet.pt());
	 topjet.JEC_uncertainty = jecUnc->getUncertainty(true);

	 topjet.btag_simpleSecondaryVertexHighEff=pat_topjet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	 topjet.btag_simpleSecondaryVertexHighPur=pat_topjet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	 topjet.btag_combinedSecondaryVertex=pat_topjet.bDiscriminator("combinedSecondaryVertexBJetTags");
	 topjet.btag_combinedSecondaryVertexMVA=pat_topjet.bDiscriminator("combinedSecondaryVertexMVABJetTags");
	 topjet.btag_jetBProbability=pat_topjet.bDiscriminator("jetBProbabilityBJetTags");
	 topjet.btag_jetProbability=pat_topjet.bDiscriminator("jetProbabilityBJetTags");

	 for (unsigned int k = 0; k < pat_topjet.numberOfDaughters(); k++) {
	   Particle subjet_v4;
	   subjet_v4.pt = pat_topjet.daughter(k)->p4().pt();
	   subjet_v4.eta = pat_topjet.daughter(k)->p4().eta();
	   subjet_v4.phi = pat_topjet.daughter(k)->p4().phi(); 
	   subjet_v4.energy = pat_topjet.daughter(k)->p4().E(); 
	   topjet.subjets.push_back(subjet_v4);
	 }
	 topjets[j].push_back(topjet);
       }
     }
   }

   // ------------- photons ------------- 
   if(doPhotons){
     for(size_t j=0; j< photon_sources.size(); ++j){
       phs[j].clear();
       
       edm::Handle< std::vector<pat::Photon> > photon_handle;
       iEvent.getByLabel(photon_sources[j], photon_handle);
       const std::vector<pat::Photon>& pat_photons = *(photon_handle.product());
       
       for (unsigned int i = 0; i < pat_photons.size(); ++i) {
	 pat::Photon pat_photon = pat_photons[i];
	 Photon ph;
	 ph.charge = 0;
	 ph.pt =  pat_photon.pt();
	 ph.eta =  pat_photon.eta();
	 ph.phi =  pat_photon.phi();
	 ph.energy =  pat_photon.energy();
	 ph.vertex_x = pat_photon.vertex().x();
	 ph.vertex_y = pat_photon.vertex().y();
	 ph.vertex_z = pat_photon.vertex().z();
	 ph.supercluster_eta = pat_photon.superCluster()->eta();
	 ph.supercluster_phi = pat_photon.superCluster()->phi();
// 	 ph.neutralHadronIso = pat_photon.neutralHadronIso();
// 	 ph.chargedHadronIso = pat_photon.chargedHadronIso();
	 ph.trackIso = pat_photon.trackIso();
	 phs[j].push_back(ph);
       }
     }
   }

   // ------------- MET -------------
   if(doMET){
     for(size_t j=0; j< met_sources.size(); ++j){
       
       edm::Handle< std::vector<pat::MET> > met_handle;
       iEvent.getByLabel(met_sources[j], met_handle);
       const std::vector<pat::MET>& pat_mets = *(met_handle.product());   
       
       if(pat_mets.size()!=1){
	 std::cout<< "WARNING: number of METs = " << pat_mets.size() <<", should be 1" << std::endl;
       }
       else{
	 pat::MET pat_met = pat_mets[0];
	 
	 met[j].pt=pat_met.pt();
	 met[j].phi=pat_met.phi();
	 met[j].mEtSig= pat_met.mEtSig();
       }
       
     }
   }

   // ------------- trigger -------------
   if(doTrigger){
     edm::InputTag triggerEvent = edm::InputTag("hltTriggerSummaryAOD");
     edm::Handle< trigger::TriggerEvent > dummy_TriggerEvent;
     iEvent.getByLabel( edm::InputTag(triggerEvent.label(), triggerEvent.instance()), dummy_TriggerEvent );
     
     const edm::Provenance *meta = dummy_TriggerEvent.provenance();
     std::string nameProcess = meta->processName();
     edm::InputTag triggerResultTag = edm::InputTag("TriggerResults");
     triggerResultTag = edm::InputTag( triggerResultTag.label(), triggerResultTag.instance(), nameProcess );
     
     edm::Handle<edm::TriggerResults> trigger;
     iEvent.getByLabel(triggerResultTag, trigger);
     const edm::TriggerResults& trig = *(trigger.product());
     
     triggerResults.clear();
     triggerNames.clear();
     L1_prescale.clear();
     HLT_prescale.clear();
     
     edm::Service<edm::service::TriggerNamesService> tns;
     std::vector<std::string> triggerNames_all;
     tns->getTrigPaths(trig,triggerNames_all);
     
     if (trig.size()!=triggerNames_all.size()) std::cout <<"ERROR: length of names and paths not the same: "<<triggerNames_all.size()<<","<<trig.size()<< std::endl;
     for(unsigned int i=0; i<trig.size(); ++i){
       std::vector<std::string>::const_iterator it = trigger_prefixes.begin();
       for(; it!=trigger_prefixes.end(); ++it){
	 if(triggerNames_all[i].substr(0, it->size()) == *it)break;
       }
       if(it==trigger_prefixes.end()) continue;
       
       //triggerResults.insert(std::pair<std::string, bool>(triggerNames[i],trig.accept(i)));
       triggerResults.push_back(trig.accept(i));
       if(newrun) triggerNames.push_back(triggerNames_all[i]);
       if(isRealData){
	 std::pair<int, int> pre=hlt_cfg.prescaleValues(iEvent, iSetup, triggerNames_all[i]);
	 L1_prescale.push_back(pre.first);
	 HLT_prescale.push_back(pre.second);
       }
     }
     //    for(std::map<std::string, bool>::const_iterator iter = triggerResults.begin(); iter!=triggerResults.end(); iter++){
     //      std::cout << (*iter).first << "   " << (*iter).second << std::endl;
     //    }
     newrun=false;
   }

   // ------------- generator info -------------
   
   if(doGenInfo){
     genInfo.weights.clear();
     genInfo.binningValues.clear();
     genps.clear();

     edm::Handle<GenEventInfoProduct> genEventInfoProduct;
     iEvent.getByLabel("generator", genEventInfoProduct);
     const GenEventInfoProduct& genEventInfo = *(genEventInfoProduct.product());
  
     genInfo.binningValues = genEventInfo.binningValues();
     genInfo.weights = genEventInfo.weights();
     genInfo.alphaQCD = genEventInfo.alphaQCD();
     genInfo.alphaQED = genEventInfo.alphaQED();
     genInfo.qScale = genEventInfo.qScale();
     
     const gen::PdfInfo* pdf = genEventInfo.pdf();
     if(pdf){
       genInfo.pdf_id1=pdf->id.first;
       genInfo.pdf_id2=pdf->id.second; 
       genInfo.pdf_x1=pdf->x.first;
       genInfo.pdf_x2=pdf->x.second;
       genInfo.pdf_scalePDF=pdf->scalePDF;
       genInfo.pdf_xPDF1=pdf->xPDF.first;
       genInfo.pdf_xPDF2=pdf->xPDF.second; 
     }
     else{
       genInfo.pdf_id1=-999;
       genInfo.pdf_id2=-999; 
       genInfo.pdf_x1=-999;
       genInfo.pdf_x2=-999;
       genInfo.pdf_scalePDF=-999;
       genInfo.pdf_xPDF1=-999;
       genInfo.pdf_xPDF2=-999;
     }

     edm::Handle<std::vector<PileupSummaryInfo> > pus;
     iEvent.getByLabel(edm::InputTag("addPileupInfo"), pus);
     genInfo.pileup_NumInteractions_intime=0;
     genInfo.pileup_NumInteractions_ootbefore=0;
     genInfo.pileup_NumInteractions_ootafter=0;
     if(pus.isValid()){
       genInfo.pileup_TrueNumInteractions = (float) pus->at(0).getTrueNumInteractions();
       for(size_t i=0; i<pus->size(); ++i){
	 if(pus->at(i).getBunchCrossing() == 0) // intime pileup
	    genInfo.pileup_NumInteractions_intime += pus->at(i).getPU_NumInteractions();
	 else if(pus->at(i).getBunchCrossing() == -1){ // oot pileup before
	    genInfo.pileup_NumInteractions_ootbefore += pus->at(i).getPU_NumInteractions();
	 }
	 else if(pus->at(i).getBunchCrossing() == +1){ // oot pileup before
	   genInfo.pileup_NumInteractions_ootafter += pus->at(i).getPU_NumInteractions();
	 }
       }
     }

     edm::Handle<reco::GenParticleCollection> genPartColl;
     iEvent.getByLabel(edm::InputTag("genParticles"), genPartColl);
     int index=-1;
     for(reco::GenParticleCollection::const_iterator iter = genPartColl->begin(); iter != genPartColl->end(); ++ iter){
       index++;
       
       //write out only top quarks and status 3 particles (works fine only for MadGraph)
       if(abs(iter->pdgId())==6 || iter->status()==3){
	 GenParticle genp;
	 genp.charge = iter->charge();;
	 genp.pt = iter->p4().pt();
	 genp.eta = iter->p4().eta();
	 genp.phi = iter->p4().phi();
	 genp.energy = iter->p4().E();
	 genp.index =index;
	 genp.status = iter->status();
	 genp.pdgId = iter->pdgId();

	 genp.mother1=-1;
	 genp.mother2=-1;
	 genp.daughter1=-1;
	 genp.daughter2=-1;
	 
	 int nm=iter->numberOfMothers();
	 int nd=iter->numberOfDaughters();
	 
	 if (nm>0) genp.mother1 = iter->motherRef(0).key();
	 if (nm>1) genp.mother2 = iter->motherRef(nm-1).key();
	 if (nd>0) genp.daughter1 = iter->daughterRef(0).key();
	 if (nd>1) genp.daughter2 = iter->daughterRef(nd-1).key();
	 
	 genps.push_back(genp);
       }
     }

   }


   tr->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
NtupleWriter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtupleWriter::endJob() 
{
  outfile->cd();
  tr->Write();
  outfile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
NtupleWriter::beginRun(edm::Run const& iRun, edm::EventSetup const&  iSetup)
{
  if(doTrigger){
    bool setup_changed = false;
    hlt_cfg.init(iRun, iSetup, "HLT", setup_changed);
    newrun=true;
  }

  if(doJets || doTopJets){
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    jecUnc = new JetCorrectionUncertainty(JetCorPar);
  }
}

// ------------ method called when ending the processing of a run  ------------
void 
NtupleWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
NtupleWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
NtupleWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtupleWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleWriter);
