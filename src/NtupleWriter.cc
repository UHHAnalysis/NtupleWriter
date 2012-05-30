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
// $Id: NtupleWriter.cc,v 1.17 2012/05/22 09:32:31 peiffer Exp $
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
  doGenTopJets = iConfig.getParameter<bool>("doGenTopJets");  
  doPhotons = iConfig.getParameter<bool>("doPhotons");
  doMET = iConfig.getParameter<bool>("doMET");
  doGenInfo = iConfig.getParameter<bool>("doGenInfo");
  doAllGenParticles = iConfig.getParameter<bool>("doAllGenParticles");
  doLumiInfo = iConfig.getParameter<bool>("doLumiInfo");
  doPV = iConfig.getParameter<bool>("doPV");
  doTopJets = iConfig.getParameter<bool>("doTopJets");
  doTrigger = iConfig.getParameter<bool>("doTrigger");

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
  if(doTopJets){
    topjet_sources = iConfig.getParameter<std::vector<std::string> >("topjet_sources");
    topjet_ptmin = iConfig.getParameter<double> ("topjet_ptmin");
    topjet_etamax = iConfig.getParameter<double> ("topjet_etamax");
    for(size_t j=0; j< topjet_sources.size(); ++j){  
      tr->Branch( topjet_sources[j].c_str(), "std::vector<TopJet>", &topjets[j]);
    }
  }
  if(doGenTopJets){
    gentopjet_sources = iConfig.getParameter<std::vector<std::string> >("gentopjet_sources");
    gentopjet_ptmin = iConfig.getParameter<double> ("gentopjet_ptmin");
    gentopjet_etamax = iConfig.getParameter<double> ("gentopjet_etamax");
    for(size_t j=0; j< gentopjet_sources.size(); ++j){  
      tr->Branch( gentopjet_sources[j].c_str(), "std::vector<TopJet>", &gentopjets[j]);
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
     iEvent.getByLabel(edm::InputTag("genParticles"), genPartColl);
     int index=-1;
     for(reco::GenParticleCollection::const_iterator iter = genPartColl->begin(); iter != genPartColl->end(); ++ iter){
       index++;
       
       //write out only top quarks and status 3 particles (works fine only for MadGraph)
       if(abs(iter->pdgId())==6 || iter->status()==3 || doAllGenParticles){
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

     edm::Handle<reco::ConversionCollection> hConversions;
     iEvent.getByLabel("allConversions", hConversions);

     edm::Handle<reco::BeamSpot> beamSpot;
     iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpot);
     const reco::BeamSpot & bsp = *beamSpot;

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
	 ele.set_passconversionveto(!ConversionTools::hasMatchedConversion(pat_ele,hConversions,bsp.position()));
	 ele.set_dEtaIn(pat_ele.deltaEtaSuperClusterTrackAtVtx());
	 ele.set_dPhiIn(pat_ele.deltaPhiSuperClusterTrackAtVtx());
	 ele.set_sigmaIEtaIEta(pat_ele.sigmaIetaIeta());
	 ele.set_HoverE(pat_ele.hadronicOverEm());
	 ele.set_fbrem(pat_ele.fbrem());
	 ele.set_EoverPIn(pat_ele.eSuperClusterOverP());
	 ele.set_EcalEnergy(pat_ele.ecalEnergy());
	 //ele.set_mvaTrigV0(pat_ele.electronID("mvaTrigV0"));
	 //ele.set_mvaNonTrigV0(pat_ele.electronID("mvaNonTrigV0"));

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
	 tau.set_byVLooseCombinedIsolationDeltaBetaCorr  ( pat_tau.tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byLooseCombinedIsolationDeltaBetaCorr ( pat_tau.tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byMediumCombinedIsolationDeltaBetaCorr ( pat_tau.tauID("byMediumCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_byTightCombinedIsolationDeltaBetaCorr ( pat_tau.tauID("byTightCombinedIsolationDeltaBetaCorr")>0.5); 
	 tau.set_againstElectronLoose  ( pat_tau.tauID("againstElectronLoose")>0.5); 
	 tau.set_againstElectronMedium ( pat_tau.tauID("againstElectronMedium")>0.5); 
	 tau.set_againstElectronTight ( pat_tau.tauID("againstElectronTight")>0.5); 
	 tau.set_againstElectronMVA  ( pat_tau.tauID("againstElectronMVA")>0.5); 
	 tau.set_againstMuonLoose ( pat_tau.tauID("againstMuonLoose")>0.5); 
	 tau.set_againstMuonMedium ( pat_tau.tauID("againstMuonMedium")>0.5); 
	 tau.set_againstMuonTight ( pat_tau.tauID("againstMuonTight")>0.5); 

	 reco::PFCandidateRef leadPFCand = pat_tau.leadPFCand();
	 if(!leadPFCand.isNull()){
	   tau.set_leadPFCand_px ( leadPFCand->px());
	   tau.set_leadPFCand_py ( leadPFCand->py());
	   tau.set_leadPFCand_pz ( leadPFCand->pz());
	 }
	 else{
	   tau.set_leadPFCand_px ( 0);
	   tau.set_leadPFCand_py ( 0);
	   tau.set_leadPFCand_pz ( 0);
	 }
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
	 jet.set_charge(pat_jet.charge());
	 jet.set_pt(pat_jet.pt());
	 jet.set_eta(pat_jet.eta());
	 jet.set_phi(pat_jet.phi());
	 jet.set_energy(pat_jet.energy());
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

	 jecUnc->setJetEta(pat_jet.eta());
	 jecUnc->setJetPt(pat_jet.pt());
	 jet.set_JEC_uncertainty(jecUnc->getUncertainty(true));
	 jet.set_JEC_factor_raw(pat_jet.jecFactor("Uncorrected"));

	 jet.set_btag_simpleSecondaryVertexHighEff(pat_jet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	 jet.set_btag_simpleSecondaryVertexHighPur(pat_jet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	 jet.set_btag_combinedSecondaryVertex(pat_jet.bDiscriminator("combinedSecondaryVertexBJetTags"));
	 jet.set_btag_combinedSecondaryVertexMVA(pat_jet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	 jet.set_btag_jetBProbability(pat_jet.bDiscriminator("jetBProbabilityBJetTags"));
	 jet.set_btag_jetProbability(pat_jet.bDiscriminator("jetProbabilityBJetTags"));

	 const reco::GenJet *genj = pat_jet.genJet();
	 if(genj){
	   jet.set_genjet_pt(genj->pt());
	   jet.set_genjet_eta(genj->eta());
	   jet.set_genjet_phi(genj->phi());
	   jet.set_genjet_energy(genj->energy());
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
	 topjet.set_numberOfDaughters(pat_topjet.numberOfDaughters());
	 const reco::TrackRefVector&  topjettracks = pat_topjet.associatedTracks();
	 topjet.set_nTracks( topjettracks.size());
	 topjet.set_jetArea( pat_topjet.jetArea());
	 topjet.set_pileup( pat_topjet.pileup());
//  	 topjet.set_neutralEmEnergyFraction(pat_topjet.neutralEmEnergyFraction());
//  	 topjet.set_neutralHadronEnergyFraction(pat_topjet.neutralHadronEnergyFraction());
//  	 topjet.set_chargedEmEnergyFraction(pat_topjet.chargedEmEnergyFraction());
//  	 topjet.set_chargedHadronEnergyFraction(pat_topjet.chargedHadronEnergyFraction());
//  	 topjet.set_muonEnergyFraction(pat_topjet.muonEnergyFraction());
//  	 topjet.set_photonEnergyFraction(pat_topjet.photonEnergyFraction());
// 	 topjet.set_chargedMultiplicity(pat_topjet.chargedMultiplicity());
// 	 topjet.set_neutralMultiplicity(pat_topjet.neutralMultiplicity());
// 	 topjet.set_muonMultiplicity(pat_topjet.muonMultiplicity());
// 	 topjet.set_electronMultiplicity(pat_topjet.electronMultiplicity());
// 	 topjet.set_photonMultiplicity(pat_topjet.photonMultiplicity());

	 jecUnc->setJetEta(pat_topjet.eta());
	 jecUnc->setJetPt(pat_topjet.pt());
	 topjet.set_JEC_uncertainty( jecUnc->getUncertainty(true));
	 topjet.set_JEC_factor_raw( pat_topjet.jecFactor("Uncorrected"));

	 topjet.set_btag_simpleSecondaryVertexHighEff(pat_topjet.bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	 topjet.set_btag_simpleSecondaryVertexHighPur(pat_topjet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	 topjet.set_btag_combinedSecondaryVertex(pat_topjet.bDiscriminator("combinedSecondaryVertexBJetTags"));
	 topjet.set_btag_combinedSecondaryVertexMVA(pat_topjet.bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	 topjet.set_btag_jetBProbability(pat_topjet.bDiscriminator("jetBProbabilityBJetTags"));
	 topjet.set_btag_jetProbability(pat_topjet.bDiscriminator("jetProbabilityBJetTags"));

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

	 for (unsigned int k = 0; k < pat_topjet.numberOfDaughters(); k++) {
	   Particle subjet_v4;
	   subjet_v4.set_pt(pat_topjet.daughter(k)->p4().pt());
	   subjet_v4.set_eta(pat_topjet.daughter(k)->p4().eta());
	   subjet_v4.set_phi(pat_topjet.daughter(k)->p4().phi()); 
	   subjet_v4.set_energy(pat_topjet.daughter(k)->p4().E()); 
	   topjet.add_subjet(subjet_v4);
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
       //edm::Handle<std::vector<reco::Jet> > reco_gentopjets;
       iEvent.getByLabel(gentopjet_sources[j], reco_gentopjets);

       for (unsigned int i = 0; i < reco_gentopjets->size(); i++) {
	 
	 const reco::BasicJet  reco_gentopjet =  reco_gentopjets->at(i);
	 if(reco_gentopjet.pt() < gentopjet_ptmin) continue;
	 if(fabs(reco_gentopjet.eta()) > gentopjet_etamax) continue;

	 TopJet gentopjet;
	 gentopjet.set_charge(reco_gentopjet.charge());
	 gentopjet.set_pt(reco_gentopjet.pt());
	 gentopjet.set_eta(reco_gentopjet.eta());
	 gentopjet.set_phi(reco_gentopjet.phi());
	 gentopjet.set_energy(reco_gentopjet.energy());
	 gentopjet.set_numberOfDaughters(reco_gentopjet.numberOfDaughters());

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

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleWriter);
