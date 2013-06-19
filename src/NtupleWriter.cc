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
// $Id: NtupleWriter.cc,v 1.27 2013/06/19 13:22:06 rkogler Exp $
//
//

#include "UHHAnalysis/NtupleWriter/interface/NtupleWriter.h"
#include "UHHAnalysis/NtupleWriter/interface/JetProps.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"
#include "DataFormats/TrackReco/interface/Track.h"

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
  doGenJets = iConfig.getParameter<bool>("doGenJets");
  doGenTopJets = iConfig.getParameter<bool>("doGenTopJets");  
  doPhotons = iConfig.getParameter<bool>("doPhotons");
  doMET = iConfig.getParameter<bool>("doMET");
  doGenInfo = iConfig.getParameter<bool>("doGenInfo");
  doAllGenParticles = iConfig.getParameter<bool>("doAllGenParticles");
  doLumiInfo = iConfig.getParameter<bool>("doLumiInfo");
  doPV = iConfig.getParameter<bool>("doPV");
  doTopJets = iConfig.getParameter<bool>("doTopJets");
  doTopJetsConstituents = iConfig.getParameter<bool>("doTopJetsConstituents");
  doTrigger = iConfig.getParameter<bool>("doTrigger");
  SVComputer_  = iConfig.getUntrackedParameter<edm::InputTag>("svComputer",edm::InputTag("combinedSecondaryVertex"));
  doTagInfos = iConfig.getUntrackedParameter<bool>("doTagInfos",false);

  // initialization of tree variables

  tr->Branch("run",&run);
  tr->Branch("event",&event);
  tr->Branch("luminosityBlock",&luminosityBlock);
  tr->Branch("isRealData",&isRealData);
  tr->Branch("rho",&rho);
  rho_source = iConfig.getParameter<edm::InputTag>("rho_source");

  //tr->Branch("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  if(doLumiInfo){
    tr->Branch("intgRecLumi",&intgRecLumi);
    tr->Branch("intgDelLumi",&intgDelLumi);
  }
  if(doPV){
    tr->Branch("beamspot_x0",&beamspot_x0);
    tr->Branch("beamspot_y0",&beamspot_y0);
    tr->Branch("beamspot_z0",&beamspot_z0);
  }
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
  if(doGenJets){
    genjet_sources = iConfig.getParameter<std::vector<std::string> >("genjet_sources");
    genjet_ptmin = iConfig.getParameter<double> ("genjet_ptmin");
    genjet_etamax = iConfig.getParameter<double> ("genjet_etamax");
    for(size_t j=0; j< genjet_sources.size(); ++j){  
      tr->Branch( genjet_sources[j].c_str(), "std::vector<Particle>", &genjets[j]);
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
  if(doTopJetsConstituents){
    topjet_constituents_sources = iConfig.getParameter<std::vector<std::string> >("topjet_constituents_sources");
    tr->Branch( "PFParticles", "std::vector<PFParticle>", &pfparticles);
  }
  if(doGenTopJets){
    gentopjet_sources = iConfig.getParameter<std::vector<std::string> >("gentopjet_sources");
    gentopjet_ptmin = iConfig.getParameter<double> ("gentopjet_ptmin");
    gentopjet_etamax = iConfig.getParameter<double> ("gentopjet_etamax");
    for(size_t j=0; j< gentopjet_sources.size(); ++j){  
      tr->Branch( gentopjet_sources[j].c_str(), "std::vector<GenTopJet>", &gentopjets[j]);
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
    genparticle_source= iConfig.getParameter<edm::InputTag>("genparticle_source");
    tr->Branch("genInfo","GenInfo",&genInfo);
    tr->Branch("GenParticles","std::vector<GenParticle>", &genps);
  }
  if(doTrigger){
    trigger_prefixes = iConfig.getParameter<std::vector<std::string> >("trigger_prefixes");
    //tr->Branch("triggerResults","std::map<std::string, bool>",&triggerResults);
    tr->Branch("triggerNames", "std::vector<std::string>", &triggerNames);  
    tr->Branch("triggerResults", "std::vector<bool>", &triggerResults);
    //tr->Branch("L1_prescale", "std::vector<int>", &L1_prescale);
    //tr->Branch("HLT_prescale", "std::vector<int>", &HLT_prescale);
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

   edm::Handle<double> m_rho;
   iEvent.getByLabel(rho_source,m_rho);
   rho=*m_rho;

//    if(isRealData){
//      edm::Handle<bool> bool_handle;
//      iEvent.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),bool_handle);
//      HBHENoiseFilterResult = *(bool_handle.product());
//    }
//    else HBHENoiseFilterResult = false;

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
	 pv.set_x( reco_pv.x());
	 pv.set_y( reco_pv.y());
	 pv.set_z( reco_pv.z());
	 pv.set_nTracks( reco_pv.nTracks());
	 //pv.set_isValid( reco_pv.isValid());
	 pv.set_chi2( reco_pv.chi2());
	 pv.set_ndof( reco_pv.ndof());	 

	 pvs[j].push_back(pv);
       }
     }
      
     edm::Handle<reco::BeamSpot> beamSpot;
     iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpot);
     const reco::BeamSpot & bsp = *beamSpot;
     
     beamspot_x0 = bsp.x0();
     beamspot_y0 = bsp.y0();
     beamspot_z0 = bsp.z0();
   }

 // ------------- generator info -------------
   
   if(doGenInfo){
     genInfo.clear_weights();
     genInfo.clear_binningValues();
     genps.clear();

     edm::Handle<GenEventInfoProduct> genEventInfoProduct;
     iEvent.getByLabel("generator", genEventInfoProduct);
     const GenEventInfoProduct& genEventInfo = *(genEventInfoProduct.product());
  
     for(unsigned int k=0; k<genEventInfo.binningValues().size();++k){
       genInfo.add_binningValue(genEventInfo.binningValues().at(k));
     }
     for(unsigned int k=0; k<genEventInfo.weights().size();++k){
       genInfo.add_weight(genEventInfo.weights().at(k));
     }
     genInfo.set_alphaQCD(genEventInfo.alphaQCD());
     genInfo.set_alphaQED(genEventInfo.alphaQED());
     genInfo.set_qScale(genEventInfo.qScale());
     
     const gen::PdfInfo* pdf = genEventInfo.pdf();
     if(pdf){
       genInfo.set_pdf_id1(pdf->id.first);
       genInfo.set_pdf_id2(pdf->id.second); 
       genInfo.set_pdf_x1(pdf->x.first);
       genInfo.set_pdf_x2(pdf->x.second);
       genInfo.set_pdf_scalePDF(pdf->scalePDF);
       genInfo.set_pdf_xPDF1(pdf->xPDF.first);
       genInfo.set_pdf_xPDF2(pdf->xPDF.second); 
     }
     else{
       genInfo.set_pdf_id1(-999);
       genInfo.set_pdf_id2(-999); 
       genInfo.set_pdf_x1(-999);
       genInfo.set_pdf_x2(-999);
       genInfo.set_pdf_scalePDF(-999);
       genInfo.set_pdf_xPDF1(-999);
       genInfo.set_pdf_xPDF2(-999);
     }

     edm::Handle<std::vector<PileupSummaryInfo> > pus;
     iEvent.getByLabel(edm::InputTag("addPileupInfo"), pus);
     genInfo.set_pileup_NumInteractions_intime(0);
     genInfo.set_pileup_NumInteractions_ootbefore(0);
     genInfo.set_pileup_NumInteractions_ootafter(0);
     if(pus.isValid()){
       genInfo.set_pileup_TrueNumInteractions ( (float) pus->at(0).getTrueNumInteractions());
       for(size_t i=0; i<pus->size(); ++i){
	 if(pus->at(i).getBunchCrossing() == 0) // intime pileup
	   genInfo.set_pileup_NumInteractions_intime( genInfo.pileup_NumInteractions_intime() + pus->at(i).getPU_NumInteractions());
	 else if(pus->at(i).getBunchCrossing() == -1){ // oot pileup before
	   genInfo.set_pileup_NumInteractions_ootbefore( genInfo.pileup_NumInteractions_ootbefore() + pus->at(i).getPU_NumInteractions());
	 }
	 else if(pus->at(i).getBunchCrossing() == +1){ // oot pileup before
	   genInfo.set_pileup_NumInteractions_ootafter( genInfo.pileup_NumInteractions_ootafter() + pus->at(i).getPU_NumInteractions());
	 }
       }
     }

     edm::Handle<reco::GenParticleCollection> genPartColl;
     iEvent.getByLabel(genparticle_source, genPartColl);
     int index=-1;
     for(reco::GenParticleCollection::const_iterator iter = genPartColl->begin(); iter != genPartColl->end(); ++ iter){
       index++;
       
       //write out only top quarks,final state leptons and status 3 particles (works fine only for MadGraph)
       bool islepton = iter->status()==1 && abs(iter->pdgId())>=11 && abs(iter->pdgId())<=16 ;
       if(abs(iter->pdgId())==6 || iter->status()==3 || islepton ||  doAllGenParticles){
	 GenParticle genp;
	 genp.set_charge(iter->charge());
	 genp.set_pt(iter->p4().pt());
	 genp.set_eta(iter->p4().eta());
	 genp.set_phi(iter->p4().phi());
	 genp.set_energy(iter->p4().E());
	 genp.set_index(index);
	 genp.set_status( iter->status());
	 genp.set_pdgId( iter->pdgId());

	 genp.set_mother1(-1);
	 genp.set_mother2(-1);
	 genp.set_daughter1(-1);
	 genp.set_daughter2(-1);
	 
	 int nm=iter->numberOfMothers();
	 int nd=iter->numberOfDaughters();

	 
	 if (nm>0) genp.set_mother1( iter->motherRef(0).key());
	 if (nm>1) genp.set_mother2( iter->motherRef(1).key());
	 if (nd>0) genp.set_daughter1( iter->daughterRef(0).key());
	 if (nd>1) genp.set_daughter2( iter->daughterRef(1).key());

	 genps.push_back(genp);
       }
     }
   }

   // ------------- electrons -------------   
   if(doElectrons){

//      edm::Handle<reco::ConversionCollection> hConversions;
//      iEvent.getByLabel("allConversions", hConversions);

//      edm::Handle<reco::BeamSpot> beamSpot;
//      iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpot);
//      const reco::BeamSpot & bsp = *beamSpot;

     for(size_t j=0; j< electron_sources.size(); ++j){
       eles[j].clear();
       edm::Handle< std::vector<pat::Electron> > ele_handle;
       iEvent.getByLabel(electron_sources[j], ele_handle);
       const std::vector<pat::Electron>& pat_electrons = *(ele_handle.product());
       
       for (unsigned int i = 0; i < pat_electrons.size(); ++i) {
	 pat::Electron pat_ele = pat_electrons[i];
	 Electron ele;
	 
	 ele.set_charge( pat_ele.charge());
	 ele.set_pt( pat_ele.pt());
	 ele.set_eta( pat_ele.eta());
	 ele.set_phi( pat_ele.phi());
	 ele.set_energy( pat_ele.energy());
	 ele.set_vertex_x(pat_ele.vertex().x());
	 ele.set_vertex_y(pat_ele.vertex().y());
	 ele.set_vertex_z(pat_ele.vertex().z());
	 ele.set_supercluster_eta(pat_ele.superCluster()->eta());
	 ele.set_supercluster_phi(pat_ele.superCluster()->phi());
	 ele.set_dB(pat_ele.dB());
	 //ele.set_particleIso(pat_ele.particleIso());
	 ele.set_neutralHadronIso(pat_ele.neutralHadronIso());
	 ele.set_chargedHadronIso(pat_ele.chargedHadronIso());
	 ele.set_trackIso(pat_ele.trackIso());
	 ele.set_photonIso(pat_ele.photonIso());
	 ele.set_puChargedHadronIso(pat_ele.puChargedHadronIso());
	 ele.set_gsfTrack_trackerExpectedHitsInner_numberOfLostHits(pat_ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits());
	 ele.set_gsfTrack_px( pat_ele.gsfTrack()->px());
	 ele.set_gsfTrack_py( pat_ele.gsfTrack()->py());
	 ele.set_gsfTrack_pz( pat_ele.gsfTrack()->pz());
	 ele.set_gsfTrack_vx( pat_ele.gsfTrack()->vx());
	 ele.set_gsfTrack_vy( pat_ele.gsfTrack()->vy());
	 ele.set_gsfTrack_vz( pat_ele.gsfTrack()->vz());
	 //ele.set_passconversionveto(!ConversionTools::hasMatchedConversion(pat_ele,hConversions,bsp.position()));
	 ele.set_passconversionveto(pat_ele.passConversionVeto());
	 ele.set_dEtaIn(pat_ele.deltaEtaSuperClusterTrackAtVtx());
	 ele.set_dPhiIn(pat_ele.deltaPhiSuperClusterTrackAtVtx());
	 ele.set_sigmaIEtaIEta(pat_ele.sigmaIetaIeta());
	 ele.set_HoverE(pat_ele.hadronicOverEm());
	 ele.set_fbrem(pat_ele.fbrem());
	 ele.set_EoverPIn(pat_ele.eSuperClusterOverP());
	 ele.set_EcalEnergy(pat_ele.ecalEnergy());
	 ele.set_mvaTrigV0(pat_ele.electronID("mvaTrigV0"));
	 ele.set_mvaNonTrigV0(pat_ele.electronID("mvaNonTrigV0"));
	 float AEff03 = 0.00;
	 if(isRealData){
	   AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, pat_ele.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
	 }else{
	   AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, pat_ele.superCluster()->eta(), ElectronEffectiveArea::kEleEAFall11MC);
	 }
	 ele.set_AEff(AEff03);

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
	 mu.set_charge( pat_mu.charge());
	 mu.set_pt( pat_mu.pt());
	 mu.set_eta( pat_mu.eta());
	 mu.set_phi( pat_mu.phi());
	 mu.set_energy( pat_mu.energy());
	 mu.set_vertex_x ( pat_mu.vertex().x());
	 mu.set_vertex_y ( pat_mu.vertex().y());
	 mu.set_vertex_z ( pat_mu.vertex().z());
	 mu.set_dB ( pat_mu.dB());
	 //mu.particleIso ( pat_mu.particleIso());
	 mu.set_neutralHadronIso ( pat_mu.neutralHadronIso());
	 mu.set_chargedHadronIso ( pat_mu.chargedHadronIso());
	 mu.set_trackIso ( pat_mu.trackIso());
	 mu.set_photonIso ( pat_mu.photonIso());
	 mu.set_puChargedHadronIso ( pat_mu.puChargedHadronIso());
	 mu.set_isGlobalMuon ( pat_mu.isGlobalMuon());
	 mu.set_isPFMuon ( pat_mu.isPFMuon());
	 mu.set_isStandAloneMuon ( pat_mu.isStandAloneMuon());
	 mu.set_isTrackerMuon ( pat_mu.isTrackerMuon());
	 mu.set_numberOfMatchedStations ( pat_mu.numberOfMatchedStations());
	 reco::TrackRef globalTrack = pat_mu.globalTrack();
	 if(!globalTrack.isNull()){
	   mu.set_globalTrack_chi2 ( globalTrack->chi2());
	   mu.set_globalTrack_ndof ( globalTrack->ndof());
	   mu.set_globalTrack_d0 ( globalTrack->d0());	 
	   mu.set_globalTrack_d0Error ( globalTrack->d0Error());
	   mu.set_globalTrack_numberOfValidHits ( globalTrack->numberOfValidHits());
	   mu.set_globalTrack_numberOfLostHits ( globalTrack->numberOfLostHits());
	   mu.set_globalTrack_numberOfValidMuonHits(globalTrack->hitPattern().numberOfValidMuonHits() );
	 }
	 else{
	   mu.set_globalTrack_chi2 ( 0);
	   mu.set_globalTrack_ndof ( 0);
	   mu.set_globalTrack_d0 ( 0);
	   mu.set_globalTrack_d0Error ( 0);
	   mu.set_globalTrack_numberOfValidHits ( 0);
	   mu.set_globalTrack_numberOfLostHits ( 0);
	 }
	 reco::TrackRef innerTrack = pat_mu.innerTrack();
	 if(!innerTrack.isNull()){
	   mu.set_innerTrack_chi2 ( innerTrack->chi2());
	   mu.set_innerTrack_ndof ( innerTrack->ndof());
	   mu.set_innerTrack_d0 ( innerTrack->d0());	 
	   mu.set_innerTrack_d0Error ( innerTrack->d0Error());
	   mu.set_innerTrack_numberOfValidHits ( innerTrack->numberOfValidHits());
	   mu.set_innerTrack_numberOfLostHits ( innerTrack->numberOfLostHits());
	   mu.set_innerTrack_trackerLayersWithMeasurement ( innerTrack->hitPattern().trackerLayersWithMeasurement());
	   mu.set_innerTrack_numberOfValidPixelHits ( innerTrack->hitPattern().numberOfValidPixelHits());
	 }
	 else{
	   mu.set_innerTrack_chi2 ( 0);
	   mu.set_innerTrack_ndof ( 0);
	   mu.set_innerTrack_d0 ( 0);
	   mu.set_innerTrack_d0Error ( 0);
	   mu.set_innerTrack_numberOfValidHits ( 0);
	   mu.set_innerTrack_numberOfLostHits ( 0);
	   mu.set_innerTrack_trackerLayersWithMeasurement ( 0);
	   mu.set_innerTrack_numberOfValidPixelHits ( 0);
	 }
	 reco::TrackRef outerTrack = pat_mu.outerTrack();
	 if(!outerTrack.isNull()){
	   mu.set_outerTrack_chi2 ( outerTrack->chi2());
	   mu.set_outerTrack_ndof ( outerTrack->ndof());
	   mu.set_outerTrack_d0 ( outerTrack->d0());	 
	   mu.set_outerTrack_d0Error ( outerTrack->d0Error());
	   mu.set_outerTrack_numberOfValidHits ( outerTrack->numberOfValidHits());
	   mu.set_outerTrack_numberOfLostHits ( outerTrack->numberOfLostHits());
	 } 
	 else{
	   mu.set_outerTrack_chi2 ( 0);
	   mu.set_outerTrack_ndof ( 0);
	   mu.set_outerTrack_d0 ( 0);
	   mu.set_outerTrack_d0Error ( 0);
	   mu.set_outerTrack_numberOfValidHits ( 0);
	   mu.set_outerTrack_numberOfLostHits ( 0);
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
	 tau.set_charge( pat_tau.charge());
	 tau.set_pt( pat_tau.pt());
	 tau.set_eta( pat_tau.eta());
	 tau.set_phi( pat_tau.phi());
	 tau.set_energy( pat_tau.energy());
	 tau.set_decayModeFinding ( pat_tau.tauID("decayModeFinding")>0.5); 
	 //tau.set_byVLooseCombinedIsolationDeltaBetaCorr  ( pat_tau.tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byLooseCombinedIsolationDeltaBetaCorr ( pat_tau.tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byMediumCombinedIsolationDeltaBetaCorr ( pat_tau.tauID("byMediumCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byTightCombinedIsolationDeltaBetaCorr ( pat_tau.tauID("byTightCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byLooseIsolationMVA( pat_tau.tauID("byLooseIsolationMVA")>0.5); 
	 tau.set_byMediumIsolationMVA( pat_tau.tauID("byMediumIsolationMVA")>0.5); 
	 tau.set_byTightIsolationMVA( pat_tau.tauID("byTightIsolationMVA")>0.5); 
	 tau.set_byLooseIsolationMVA2( pat_tau.tauID("byLooseIsolationMVA2")>0.5); 
	 tau.set_byMediumIsolationMVA2( pat_tau.tauID("byMediumIsolationMVA2")>0.5); 
	 tau.set_byTightIsolationMVA2( pat_tau.tauID("byTightIsolationMVA2")>0.5);
	 tau.set_byLooseCombinedIsolationDeltaBetaCorr3Hits(  pat_tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>0.5); 
	 tau.set_byMediumCombinedIsolationDeltaBetaCorr3Hits ( pat_tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")>0.5); 
	 tau.set_byTightCombinedIsolationDeltaBetaCorr3Hits ( pat_tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")>0.5); 
	 tau.set_againstElectronLooseMVA3  ( pat_tau.tauID("againstElectronLooseMVA3")>0.5); 
	 tau.set_againstElectronMediumMVA3 ( pat_tau.tauID("againstElectronMediumMVA3")>0.5); 
	 tau.set_againstElectronTightMVA3 ( pat_tau.tauID("againstElectronTightMVA3")>0.5); 
	 tau.set_againstElectronVTightMVA3 ( pat_tau.tauID("againstElectronVTightMVA3")>0.5); 
	 tau.set_againstMuonLoose2 ( pat_tau.tauID("againstMuonLoose2")>0.5); 
	 tau.set_againstMuonMedium2 ( pat_tau.tauID("againstMuonMedium2")>0.5); 
	 tau.set_againstMuonTight2 ( pat_tau.tauID("againstMuonTight2")>0.5); 
	 tau.set_byIsolationMVAraw(  pat_tau.tauID("byIsolationMVAraw"));
	 tau.set_byIsolationMVA2raw(  pat_tau.tauID("byIsolationMVA2raw"));
	 tau.set_decayMode( pat_tau.decayMode() );
	 tau.set_byCombinedIsolationDeltaBetaCorrRaw( pat_tau.tauID("byCombinedIsolationDeltaBetaCorrRaw"));
	 tau.set_byCombinedIsolationDeltaBetaCorrRaw3Hits( pat_tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));

//       std::cout << pat_tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << std::endl;
	 
// 	 reco::PFCandidateRef leadPFCand = pat_tau.leadPFCand();
// 	 if(!leadPFCand.isNull()){
// 	   tau.set_leadPFCand_px ( leadPFCand->px());
// 	   tau.set_leadPFCand_py ( leadPFCand->py());
// 	   tau.set_leadPFCand_pz ( leadPFCand->pz());
// 	 }
// 	 else{
// 	   tau.set_leadPFCand_px ( 0);
// 	   tau.set_leadPFCand_py ( 0);
// 	   tau.set_leadPFCand_pz ( 0);
// 	 }
	 taus[j].push_back(tau);
       }
     }
   }

   //-------------- gen jets -------------

   if(doGenJets){
     for(size_t j=0; j< genjet_sources.size(); ++j){
       
       genjets[j].clear();

       edm::Handle< std::vector<reco::GenJet> > genjet_handle;
       iEvent.getByLabel(genjet_sources[j], genjet_handle);
       const std::vector<reco::GenJet>& gen_jets = *(genjet_handle.product());
  
       for (unsigned int i = 0; i < gen_jets.size(); ++i) {

	 const reco::GenJet* gen_jet = &gen_jets[i];
	 if(gen_jet->pt() < genjet_ptmin) continue;
	 if(fabs(gen_jet->eta()) > genjet_etamax) continue;

	 Particle jet;
	 jet.set_charge(gen_jet->charge());
	 jet.set_pt(gen_jet->pt());
	 jet.set_eta(gen_jet->eta());
	 jet.set_phi(gen_jet->phi());
	 jet.set_energy(gen_jet->energy());

	 // recalculate the jet charge
	 int jet_charge = 0;
	 std::vector<const reco::GenParticle * > jetgenps = gen_jet->getGenConstituents();
	 for(unsigned int l = 0; l<jetgenps.size(); ++l){
	   jet_charge +=  jetgenps[l]->charge();
	 }
	 jet.set_charge(jet_charge);

	 genjets[j].push_back(jet);
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
	 jet.set_charge(pat_jet.charge());
	 jet.set_pt(pat_jet.pt());
	 jet.set_eta(pat_jet.eta());
	 jet.set_phi(pat_jet.phi());
	 jet.set_energy(pat_jet.energy());
	 jet.set_flavor(pat_jet.partonFlavour());
	 jet.set_numberOfDaughters (pat_jet.numberOfDaughters());
	 const reco::TrackRefVector&  jettracks = pat_jet.associatedTracks();
	 jet.set_nTracks ( jettracks.size());
	 jet.set_jetArea(pat_jet.jetArea());
	 jet.set_pileup(pat_jet.pileup());
	 if(pat_jet.isPFJet()){
	   jet.set_neutralEmEnergyFraction (pat_jet.neutralEmEnergyFraction());
	   jet.set_neutralHadronEnergyFraction (pat_jet.neutralHadronEnergyFraction());
	   jet.set_chargedEmEnergyFraction (pat_jet.chargedEmEnergyFraction());
	   jet.set_chargedHadronEnergyFraction (pat_jet.chargedHadronEnergyFraction());
	   jet.set_muonEnergyFraction (pat_jet.muonEnergyFraction());
	   jet.set_photonEnergyFraction (pat_jet.photonEnergyFraction());
	   jet.set_chargedMultiplicity (pat_jet.chargedMultiplicity());
	   jet.set_neutralMultiplicity (pat_jet.neutralMultiplicity());
	   jet.set_muonMultiplicity (pat_jet.muonMultiplicity());
	   jet.set_electronMultiplicity (pat_jet.electronMultiplicity());
	   jet.set_photonMultiplicity (pat_jet.photonMultiplicity());
	 }

	 jet.set_JEC_factor_raw(pat_jet.jecFactor("Uncorrected"));
	 
	 jet.set_btag_simpleSecondaryVertexHighEff(pat_jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	 jet.set_btag_simpleSecondaryVertexHighPur(pat_jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	 jet.set_btag_combinedSecondaryVertex(pat_jet.bDiscriminator("combinedSecondaryVertexBJetTags"));
	 jet.set_btag_combinedSecondaryVertexMVA(pat_jet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	 jet.set_btag_jetBProbability(pat_jet.bDiscriminator("jetBProbabilityBJetTags"));
	 jet.set_btag_jetProbability(pat_jet.bDiscriminator("jetProbabilityBJetTags"));

	 
	 const reco::GenJet *genj = pat_jet.genJet();
	 if(genj){

	   for(unsigned int k=0; k<genjets->size(); ++k){
	     if(genj->pt()==genjets->at(k).pt() && genj->eta()==genjets->at(k).eta()){
	       jet.set_genjet_index(k);
	     }
	   }
// 	   if( jet.genjet_index()<0){
// 	     std::cout<< "genjet not found for " << genj->pt() << "  " << genj->eta() << std::endl;
// 	   }

	   if(doAllGenParticles){
	     std::vector<const reco::GenParticle * > jetgenps = genj->getGenConstituents();
	     for(unsigned int l = 0; l<jetgenps.size(); ++l){
	       for(unsigned int k=0; k< genps.size(); ++k){
		 if(jetgenps[l]->pt() == genps[k].pt() && jetgenps[l]->pdgId() == genps[k].pdgId()){
		   jet.add_genparticles_index(genps[k].index());
		   break;
		 }
	       }
	     }
	     if(jet.genparticles_indices().size()!= jetgenps.size())
	       std::cout << "WARNING: Found only " << jet.genparticles_indices().size() << " from " << jetgenps.size() << " gen particles of this jet"<<std::endl;
	   }
	   
	 }
	 
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
	 topjet.set_charge(pat_topjet.charge());
	 topjet.set_pt(pat_topjet.pt());
	 topjet.set_eta(pat_topjet.eta());
	 topjet.set_phi(pat_topjet.phi());
	 topjet.set_energy(pat_topjet.energy());
	 topjet.set_flavor(pat_topjet.partonFlavour());
	 topjet.set_numberOfDaughters(pat_topjet.numberOfDaughters());
	 const reco::TrackRefVector&  topjettracks = pat_topjet.associatedTracks();
	 topjet.set_nTracks( topjettracks.size());
	 topjet.set_jetArea( pat_topjet.jetArea());
	 topjet.set_pileup( pat_topjet.pileup());

	 topjet.set_JEC_factor_raw( pat_topjet.jecFactor("Uncorrected"));

	 topjet.set_btag_simpleSecondaryVertexHighEff(pat_topjet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	 topjet.set_btag_simpleSecondaryVertexHighPur(pat_topjet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	 topjet.set_btag_combinedSecondaryVertex(pat_topjet.bDiscriminator("combinedSecondaryVertexBJetTags"));
	 topjet.set_btag_combinedSecondaryVertexMVA(pat_topjet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	 topjet.set_btag_jetBProbability(pat_topjet.bDiscriminator("jetBProbabilityBJetTags"));
	 topjet.set_btag_jetProbability(pat_topjet.bDiscriminator("jetProbabilityBJetTags"));

	 /*
	 const reco::GenJet *genj = pat_topjet.genJet();
	 if(genj){
	   topjet.set_genjet_pt ( genj->pt());
	   topjet.set_genjet_eta ( genj->eta());
	   topjet.set_genjet_phi ( genj->phi());
	   topjet.set_genjet_energy ( genj->energy());
	   if(doAllGenParticles){
	     std::vector<const reco::GenParticle * > jetgenps = genj->getGenConstituents();
	     for(unsigned int l = 0; l<jetgenps.size(); ++l){
	       for(unsigned int k=0; k< genps.size(); ++k){
		 if(jetgenps[l]->pt() == genps[k].pt() && jetgenps[l]->pdgId() == genps[k].pdgId()){
		   topjet.add_genparticles_index(genps[k].index());
		   break;
		 }
	       }
	     }
	     if(topjet.genparticles_indices().size()!= jetgenps.size())
	       std::cout << "WARNING: Found only " << topjet.genparticles_indices().size() << " from " << jetgenps.size() << " gen particles of this topjet"<<std::endl;
	   }
	 }
	 */

	 // add constituents to the jet, if requested
	 if (doTopJetsConstituents){

	   if (topjet_constituents_sources.size()>j){ //only add constituents if they are defined

	     edm::Handle<pat::JetCollection> pat_topjets_with_cands;
	     iEvent.getByLabel(topjet_constituents_sources[j], pat_topjets_with_cands);
	     pat::Jet* pat_topjet_wc = NULL;

	     for (unsigned int it = 0; it < pat_topjets_with_cands->size(); it++) {
	       const pat::Jet* cand =  dynamic_cast<pat::Jet const *>(&pat_topjets_with_cands->at(it));
	       double dphi = cand->phi() - pat_topjet.phi();
	       if (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
	       if (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();	       
	       if (fabs(dphi)<0.5 && fabs(cand->eta()-pat_topjet.eta())<0.5){ // be generous: filtering, pruning... can change jet axis
		 pat_topjet_wc = const_cast<pat::Jet*>(cand);
		 break;
	       }
	     }

	     if (pat_topjet_wc){
	       StoreJetConstituents(pat_topjet_wc, &topjet);

	       // now run substructure information
	       JetProps substructure(&topjet);
	       substructure.set_pf_cands(&pfparticles);
	       double tau1 = substructure.GetNsubjettiness(1, Njettiness::onepass_kt_axes, 1., 2.0);
	       double tau2 = substructure.GetNsubjettiness(2, Njettiness::onepass_kt_axes, 1., 2.0);
	       double tau3 = substructure.GetNsubjettiness(3, Njettiness::onepass_kt_axes, 1., 2.0);
	       double qjets = substructure.GetQjetVolatility(iEvent.id().event(), 2.0);
	       topjet.set_tau1(tau1);
	       topjet.set_tau2(tau2);
	       topjet.set_tau3(tau3);
	       topjet.set_qjets_volatility(qjets);

	     }
	   }
	 }



	 for (unsigned int k = 0; k < pat_topjet.numberOfDaughters(); k++) {
	   Particle subjet_v4;

	   reco::Candidate const * subjetd =  pat_topjet.daughter(k);
           pat::Jet const * patsubjetd = dynamic_cast<pat::Jet const *>(subjetd);
	   if(patsubjetd)
	   {
	      subjet_v4.set_pt(patsubjetd->correctedP4(0).pt());
              subjet_v4.set_eta(patsubjetd->correctedP4(0).eta());
              subjet_v4.set_phi(patsubjetd->correctedP4(0).phi()); 
              subjet_v4.set_energy(patsubjetd->correctedP4(0).E());
              topjet.add_subjet(subjet_v4);
	      //btag info
              topjet.add_subFlavour(patsubjetd->partonFlavour());
              topjet.add_subCSV(patsubjetd->bDiscriminator("combinedSecondaryVertexBJetTags"));
	      if (doTagInfos)
		{
		  //ip tag info
		  reco::TaggingVariableList tvlIP=patsubjetd->tagInfoTrackIP("impactParameter")->taggingVariables();
		  topjet.add_subTrackMomentum(tvlIP.getList(reco::btau::trackMomentum,false));
		  topjet.add_subTrackEta(tvlIP.getList(reco::btau::trackEta,false));
		  topjet.add_subTrackEtaRel(tvlIP.getList(reco::btau::trackEtaRel,false));
		  topjet.add_subTrackDeltaR(tvlIP.getList(reco::btau::trackDeltaR,false));
		  topjet.add_subTrackSip3dVal(tvlIP.getList(reco::btau::trackSip3dVal,false));
		  topjet.add_subTrackSip3dSig(tvlIP.getList(reco::btau::trackSip3dSig,false));
		  topjet.add_subTrackSip2dVal(tvlIP.getList(reco::btau::trackSip2dVal,false));
		  topjet.add_subTrackSip2dSig(tvlIP.getList(reco::btau::trackSip2dSig,false));
		  topjet.add_subTrackDecayLenVal(tvlIP.getList(reco::btau::trackDecayLenVal,false));
		  topjet.add_subTrackChi2(tvlIP.getList(reco::btau::trackChi2,false));
		  topjet.add_subTrackNTotalHits(tvlIP.getList(reco::btau::trackNTotalHits,false));
		  topjet.add_subTrackNPixelHits(tvlIP.getList(reco::btau::trackNPixelHits,false));     
		  topjet.add_subTrackPtRel(tvlIP.getList(reco::btau::trackPtRel,false));
		  topjet.add_subTrackPPar(tvlIP.getList(reco::btau::trackPPar,false));
		  topjet.add_subTrackPtRatio(tvlIP.getList(reco::btau::trackPtRatio,false));
		  topjet.add_subTrackPParRatio(tvlIP.getList(reco::btau::trackPParRatio,false));
		  topjet.add_subTrackJetDistVal(tvlIP.getList(reco::btau::trackJetDistVal,false));
		  topjet.add_subTrackJetDistSig(tvlIP.getList(reco::btau::trackJetDistSig,false));
		  topjet.add_subTrackGhostTrackDistVal(tvlIP.getList(reco::btau::trackGhostTrackDistVal,false));
		  topjet.add_subTrackGhostTrackDistSig(tvlIP.getList(reco::btau::trackGhostTrackDistSig,false));
		  topjet.add_subTrackGhostTrackWeight(tvlIP.getList(reco::btau::trackGhostTrackWeight,false));
		  //sv tag info
		  reco::TaggingVariableList tvlSV=patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->taggingVariables();
		  topjet.add_subFlightDistance2dVal(tvlSV.getList(reco::btau::flightDistance2dVal,false));
		  topjet.add_subFlightDistance2dSig(tvlSV.getList(reco::btau::flightDistance2dSig,false));
		  topjet.add_subFlightDistance3dVal(tvlSV.getList(reco::btau::flightDistance3dVal,false));
		  topjet.add_subFlightDistance3dSig(tvlSV.getList(reco::btau::flightDistance3dSig,false));
		  topjet.add_subVertexJetDeltaR(tvlSV.getList(reco::btau::vertexJetDeltaR,false));
		  topjet.add_subJetNSecondaryVertices(tvlSV.get(reco::btau::jetNSecondaryVertices,-9999));
		  topjet.add_subVertexNTracks(tvlSV.get(reco::btau::vertexNTracks,-9999));
		  std::vector<TLorentzVector> vp4; vp4.clear();
		  std::vector<float> vchi2; vchi2.clear();
		  std::vector<float> vndof; vndof.clear();
		  std::vector<float> vchi2ndof; vchi2ndof.clear();
		  std::vector<float> vtrsize; vtrsize.clear();
		  for(unsigned int i=0; i<patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->nVertices(); i++)
		    {
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4 = patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(i).p4();
		      vp4.push_back(TLorentzVector(p4.px(),p4.py(),p4.pz(),p4.e()));
		      vchi2.push_back(patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(i).chi2());  
		      vndof.push_back(patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(i).ndof());  
		      vchi2ndof.push_back(patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(i).normalizedChi2());  
		      vtrsize.push_back(patsubjetd->tagInfoSecondaryVertex("secondaryVertex")->secondaryVertex(i).tracksSize());  
		    }
		  topjet.add_subSecondaryVertex(vp4);
		  topjet.add_subVertexChi2(vchi2);
		  topjet.add_subVertexNdof(vndof);
		  topjet.add_subVertexNormalizedChi2(vchi2ndof);
		  topjet.add_subVertexTracksSize(vtrsize);
		  //try computer
		  edm::ESHandle<JetTagComputer> computerHandle;
		  iSetup.get<JetTagComputerRecord>().get( SVComputer_.label(), computerHandle );
		  const GenericMVAJetTagComputer *computer = dynamic_cast<const GenericMVAJetTagComputer*>( computerHandle.product() );
		  if(computer)
		    {
		      computer->passEventSetup(iSetup);
		      std::vector<const reco::BaseTagInfo*>  baseTagInfos;
		      baseTagInfos.push_back(patsubjetd->tagInfoTrackIP("impactParameter") );
		      baseTagInfos.push_back(patsubjetd->tagInfoSecondaryVertex("secondaryVertex") );      
		      JetTagComputer::TagInfoHelper helper(baseTagInfos);
		      reco::TaggingVariableList vars = computer->taggingVariables(helper);
		      topjet.add_subVertexMassJTC(vars.get(reco::btau::vertexMass,-9999));
		      topjet.add_subVertexCategoryJTC(vars.get(reco::btau::vertexCategory,-9999));
		      topjet.add_subVertexEnergyRatioJTC(vars.get(reco::btau::vertexEnergyRatio,-9999));
		      topjet.add_subTrackSip3dSigAboveCharmJTC(vars.get(reco::btau::trackSip3dSigAboveCharm,-9999));
		    }
		}
	   }
	   else
	     {
	       //filling only standard information in case the subjet has not been pat-tified during the pattuples production
	       subjet_v4.set_pt(pat_topjet.daughter(k)->p4().pt());
	       subjet_v4.set_eta(pat_topjet.daughter(k)->p4().eta());
	       subjet_v4.set_phi(pat_topjet.daughter(k)->p4().phi());
	       subjet_v4.set_energy(pat_topjet.daughter(k)->p4().E());
	       topjet.add_subjet(subjet_v4);
	     }
	   
	   
	 }
	 topjets[j].push_back(topjet);
       }
     }
   }

   
   // ------------- generator top jets -------------
   if(doGenTopJets){
     for(size_t j=0; j< gentopjet_sources.size(); ++j){
       
       gentopjets[j].clear();
       
       edm::Handle<reco::BasicJetCollection> reco_gentopjets;
       iEvent.getByLabel(gentopjet_sources[j], reco_gentopjets);

       for (unsigned int i = 0; i < reco_gentopjets->size(); i++) {
	 
	 const reco::BasicJet  reco_gentopjet =  reco_gentopjets->at(i);
	 if(reco_gentopjet.pt() < gentopjet_ptmin) continue;
	 if(fabs(reco_gentopjet.eta()) > gentopjet_etamax) continue;

	 GenTopJet gentopjet;
	 gentopjet.set_charge(reco_gentopjet.charge());
	 gentopjet.set_pt(reco_gentopjet.pt());
	 gentopjet.set_eta(reco_gentopjet.eta());
	 gentopjet.set_phi(reco_gentopjet.phi());
	 gentopjet.set_energy(reco_gentopjet.energy());

	 for (unsigned int k = 0; k < reco_gentopjet.numberOfDaughters(); k++) {
	   Particle subjet_v4;
	   subjet_v4.set_pt(reco_gentopjet.daughter(k)->p4().pt());
	   subjet_v4.set_eta(reco_gentopjet.daughter(k)->p4().eta());
	   subjet_v4.set_phi(reco_gentopjet.daughter(k)->p4().phi()); 
	   subjet_v4.set_energy(reco_gentopjet.daughter(k)->p4().E()); 
	   gentopjet.add_subjet(subjet_v4);
	 }
	 gentopjets[j].push_back(gentopjet);
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
	 ph.set_charge(0);
	 ph.set_pt( pat_photon.pt());
	 ph.set_eta( pat_photon.eta());
	 ph.set_phi( pat_photon.phi());
	 ph.set_energy( pat_photon.energy());
	 ph.set_vertex_x(pat_photon.vertex().x());
	 ph.set_vertex_y(pat_photon.vertex().y());
	 ph.set_vertex_z(pat_photon.vertex().z());
	 ph.set_supercluster_eta(pat_photon.superCluster()->eta());
	 ph.set_supercluster_phi(pat_photon.superCluster()->phi());
// 	 ph.set_neutralHadronIso(pat_photon.neutralHadronIso());
// 	 ph.set_chargedHadronIso(pat_photon.chargedHadronIso());
	 ph.set_trackIso(pat_photon.trackIso());
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
	 
	 met[j].set_pt(pat_met.pt());
	 met[j].set_phi(pat_met.phi());
	 met[j].set_mEtSig(pat_met.mEtSig());
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
//      L1_prescale.clear();
//      HLT_prescale.clear();
     
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
//        if(isRealData){
// 	 std::pair<int, int> pre=hlt_cfg.prescaleValues(iEvent, iSetup, triggerNames_all[i]);
// 	 L1_prescale.push_back(pre.first);
// 	 HLT_prescale.push_back(pre.second);
// 	 //std::cout <<  triggerNames_all[i] << " " << pre.first << " " <<pre.second << "   " << hlt_cfg.prescaleValue(iEvent, iSetup, triggerNames_all[i]) << std::endl;
//        }
     }
     //    for(std::map<std::string, bool>::const_iterator iter = triggerResults.begin(); iter!=triggerResults.end(); iter++){
     //      std::cout << (*iter).first << "   " << (*iter).second << std::endl;
     //    }
     newrun=false;
   }


   tr->Fill();
   if(doLumiInfo)
     previouslumiblockwasfilled=true;

   // clean up
   m_stored_pfs.clear();
   if(doTopJetsConstituents){
     pfparticles.clear();
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
NtupleWriter::beginJob()
{
  if(doLumiInfo){
    totalRecLumi=0;
    totalDelLumi=0;
    previouslumiblockwasfilled=false;
  }
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
    //bool setup_changed = false;
    //hlt_cfg.init(iRun, iSetup, "HLT", setup_changed);
    newrun=true;
  }

}

// ------------ method called when ending the processing of a run  ------------
void 
NtupleWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
  if(doLumiInfo)
    std::cout << "total integ. luminosity: " << totalDelLumi <<"(del) " << totalRecLumi << "(rec)" << std::endl;
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
NtupleWriter::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const&)
{
  if(doLumiInfo){
    edm::Handle<LumiSummary> l;
    lumi.getByLabel("lumiProducer", l); 
    
    //add lumi of lumi blocks without any event to next lumiblock
    if(previouslumiblockwasfilled){
      intgRecLumi=0;
      intgDelLumi=0; 
    }
    previouslumiblockwasfilled=false;
    
    if (l.isValid()){;
      intgRecLumi+=l->intgRecLumi()*6.37;
      intgDelLumi+=l->intgDelLumi()*6.37;
      totalRecLumi+=l->intgRecLumi()*6.37;
      totalDelLumi+=l->intgDelLumi()*6.37;
    }
    //std::cout << "this lb: " <<l->intgRecLumi()*6.37 <<"   " << l->intgDelLumi()*6.37<<std::endl;
    //std::cout << "summed:  "<< intgRecLumi << "   " << intgDelLumi << std::endl;
  }
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

// ------------ method fills constituents of the pat_jet into the Ntuple and stores a reference 
// ------------ to those in the topjet
// ------------ it is checked if the constituents have been stored already
void
NtupleWriter::StoreJetConstituents(pat::Jet* pat_jet, Jet* jet)
{
  // checks if the pf cand has already been stored, only stores so far missing
  // PF candidates, then stores a reference to the pf candidate in the jet
  // also calculates the jet charge and sets it


  const std::vector<reco::PFCandidatePtr> jconstits = pat_jet->getPFConstituents();

  // loop over all jet constituents and store them 
  for (unsigned int i=0; i<jconstits.size(); ++i){

    PFParticle part;
    const reco::PFCandidate* pf = jconstits[i].get();

    // check if it has already been stored, omit if true
    int is_already_in_list = -1;
    for (unsigned int j=0; j<m_stored_pfs.size(); ++j){
      
      const reco::PFCandidate* spf = m_stored_pfs[j];
      double r2 = pow(pf->eta()-spf->eta(),2) + pow(pf->phi()-spf->phi(),2);
      double dpt = fabs( pf->pt() - spf->pt() );
      if (r2<1e-10 && dpt<1e-10){
	is_already_in_list = j;
	break;
      }            
    }
    
    if (is_already_in_list>-1){
      jet->add_pfconstituents_index(is_already_in_list);
      continue;
    }

    part.set_pt(pf->pt());
    part.set_eta(pf->eta());
    part.set_phi(pf->phi());
    part.set_energy(pf->energy());
    part.set_charge(pf->charge());

    part.set_ecal_en(pf->ecalEnergy());
    part.set_hcal_en(pf->hcalEnergy());
    reco::TrackRef trackref = pf->trackRef();
    if (!trackref.isNull()){
      part.set_track_mom(trackref->p());
    }

    PFParticle::EParticleID id = PFParticle::eX;
    switch ( pf->translatePdgIdToType(pf->pdgId()) ){
    case reco::PFCandidate::X : id = PFParticle::eX; break;
    case reco::PFCandidate::h : id = PFParticle::eH; break;
    case reco::PFCandidate::e : id = PFParticle::eE; break;
    case reco::PFCandidate::mu : id = PFParticle::eMu; break;
    case reco::PFCandidate::gamma : id = PFParticle::eGamma; break;
    case reco::PFCandidate::h0 : id = PFParticle::eH0; break;
    case reco::PFCandidate::h_HF : id = PFParticle::eH_HF; break;
    case reco::PFCandidate::egamma_HF : id = PFParticle::eEgamma_HF; break;
    }
    part.set_particleID(id);

    pfparticles.push_back(part);
    m_stored_pfs.push_back(pf);

    // add a reference to the particle 
    jet->add_pfconstituents_index(pfparticles.size()-1);
    
  }
  
  if(pat_jet->isPFJet()){
    jet->set_charge(pat_jet->charge());
    jet->set_neutralEmEnergyFraction (pat_jet->neutralEmEnergyFraction());
    jet->set_neutralHadronEnergyFraction (pat_jet->neutralHadronEnergyFraction());
    jet->set_chargedEmEnergyFraction (pat_jet->chargedEmEnergyFraction());
    jet->set_chargedHadronEnergyFraction (pat_jet->chargedHadronEnergyFraction());
    jet->set_muonEnergyFraction (pat_jet->muonEnergyFraction());
    jet->set_photonEnergyFraction (pat_jet->photonEnergyFraction());
    jet->set_chargedMultiplicity (pat_jet->chargedMultiplicity());
    jet->set_neutralMultiplicity (pat_jet->neutralMultiplicity());
    jet->set_muonMultiplicity (pat_jet->muonMultiplicity());
    jet->set_electronMultiplicity (pat_jet->electronMultiplicity());
    jet->set_photonMultiplicity (pat_jet->photonMultiplicity());
  }
  
  return;

}

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleWriter);
