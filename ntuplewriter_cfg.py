# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("MYNTUPLE")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.options   = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
 fileNames  = cms.untracked.vstring('file:/scratch/hh/lustre/cms/user/peiffer/ttbsm_52x_mc.root'),
 skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

#NtupleWriter
useData = False
writeAllGenParticles=True

process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string('/scratch/hh/lustre/cms/user/peiffer/SFrame_Ntuples/TTbarTest.root'), 
                                  doElectrons = cms.bool(True),
                                  doMuons = cms.bool(True),
                                  doTaus = cms.bool(True),
                                  doJets = cms.bool(True),
                                  doTopJets = cms.bool(True),
                                  doJECUncertainty = cms.bool(False),
                                  doPhotons = cms.bool(False),
                                  doMET = cms.bool(True),
                                  doPV = cms.bool(True),
                                  doGenInfo = cms.bool(not useData),
				  doAllGenParticles = cms.bool(writeAllGenParticles),
				  doLumiInfo = cms.bool(useData),
                                  doTrigger = cms.bool(True),
				  rho_source =  cms.InputTag("kt6PFJetsPFlow", "rho", "PAT"),
                                  electron_sources = cms.vstring("selectedPatElectronsPFlow","selectedPatElectronsLoosePFlow"),
                                  muon_sources = cms.vstring("selectedPatMuonsPFlow","selectedPatMuonsLoosePFlow"),
                                  tau_sources = cms.vstring("selectedPatTausPFlow","selectedPatTaus"),
                                  tau_ptmin = cms.double(0.0),
                                  tau_etamax = cms.double(999.0),
                                  jet_sources = cms.vstring("goodPatJetsPFlow"),
                                  jet_ptmin = cms.double(10.0),
                                  jet_etamax = cms.double(5.0),
				  #photon_sources = cms.vstring("selectedPatPhotons"),
                                  topjet_sources = cms.vstring("goodPatJetsCATopTagPF","goodPatJetsCA8PrunedPF"),
                                  topjet_ptmin = cms.double(150.0), 
                                  topjet_etamax = cms.double(5.0),
				  doGenTopJets = cms.bool(not useData),
                                  gentopjet_sources = cms.vstring("caTopTagGen" ),
                                  gentopjet_ptmin = cms.double(150.0), 
                                  gentopjet_etamax = cms.double(5.0),
                                  met_sources =  cms.vstring("patMETs","patMETsPFlow"),
                                  pv_sources = cms.vstring("goodOfflinePrimaryVertices"),
                                  trigger_prefixes = cms.vstring(#"HLT_IsoMu", "HLT_Mu",
                                                                 #"HLT_L1SingleMu", "HLT_L2Mu",
                                                                 #"HLT_Ele",
                                                                 "HLT_",
                                                                 #"HLT_DoubleMu", "HLT_DoubleEle"
	                                                         ),
                                  
)


process.p = cms.Path(process.MyNtuple)
