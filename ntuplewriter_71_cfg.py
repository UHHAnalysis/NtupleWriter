
# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("MYNTUPLE")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.options   = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
 fileNames  = cms.untracked.vstring('file:/nfs/dust/cms/user/peiffer/TTJets_miniaod/004C6DA7-FB03-E411-96BD-0025905A497A.root'),

 skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

#NtupleWriter
useData = False
writeAllGenParticles=True


process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string("Ntuple.root"), 
                                  doElectrons = cms.bool(True),
                                  doMuons = cms.bool(True),
                                  doTaus = cms.bool(True),
                                  doJets = cms.bool(True),
                                  doJetsConstituents = cms.bool(True),
                                  doTopJets = cms.bool(False),
                                  doTopJetsConstituents = cms.bool(False),
                                  doGenJets = cms.bool(not useData),
                                  doPhotons = cms.bool(False),
                                  doMET = cms.bool(True),
                                  doPV = cms.bool(True),
                                  doRho = cms.untracked.bool(True),
                                  doAllPFConstituents = cms.bool(True),
                                  pf_constituents_sources = cms.vstring("packedPFCandidates"),
				  storePFsAroundLeptons = cms.untracked.bool(True),
                                  doGenInfo = cms.bool(not useData),
				  doAllGenParticles = cms.bool(writeAllGenParticles), #set to true if you want to store all gen particles, otherwise, only tops and status 3 particles are stored
				  doLumiInfo = cms.bool(useData),
                                  doTrigger = cms.bool(True),
                                  doTagInfos = cms.untracked.bool(True), # when set to true crashes for the 'packed' jet collections
                                  svComputer = cms.untracked.InputTag("combinedSecondaryVertex"),
				  rho_source = cms.InputTag("fixedGridRhoFastjetAll"),
                                  genparticle_source = cms.InputTag("prunedGenParticles" ),
                                  stablegenparticle_source = cms.InputTag("packedGenParticles" ),
                                  electron_sources = cms.vstring("slimmedElectrons"),
                                  muon_sources = cms.vstring("slimmedMuons"),
                                  tau_sources = cms.vstring("slimmedTaus" ),
                                  tau_ptmin = cms.double(0.0),
                                  tau_etamax = cms.double(999.0),
                                  jet_sources = cms.vstring("slimmedJets","slimmedJetsAK8"),
                                  jet_ptmin = cms.double(10.0),
                                  jet_etamax = cms.double(5.0),
                                  genjet_sources = cms.vstring("slimmedGenJets"),
                                  genjet_ptmin = cms.double(10.0),
                                  genjet_etamax = cms.double(5.0),
				  #photon_sources = cms.vstring("selectedPatPhotons"),
                                  #missing in miniaod
                                  #topjet_sources = cms.vstring("goodPatJetsCATopTagPFPacked","goodPatJetsCA12TopTagPFPacked", "goodPatJetsCA15FilteredPFPacked", "goodPatJetsCAHEPTopTagPFPacked"),
                                  #topjet_constituents_sources = cms.vstring("goodPatJetsCA8PF","goodPatJetsCA12PF", "goodPatJetsCA15PF"),
                                  #topjet_constituents_sources = cms.vstring(),
                                  #topjet_ptmin = cms.double(30.0), 
                                  #topjet_etamax = cms.double(5.0),
                                  pf_around_leptons_sources = cms.vstring("packedPFCandidates"),
                                  #missing in miniaod
				  doGenTopJets = cms.bool(False),			      
                                  #gentopjet_sources = cms.vstring("caTopTagGen","ca12TopTagGen", "caFilteredGenJetsNoNu", "caHEPTopTagGen"),
                                  #gentopjet_ptmin = cms.double(150.0), 
                                  #gentopjet_etamax = cms.double(5.0),
                                  #missing in miniaod
  				  doGenJetsWithParts = cms.bool(False),
                                  #genjetwithparts_sources = cms.vstring("slimmedGenJets"),
                                  #genjetwithparts_ptmin = cms.double(10.0), 
                                  #genjetwithparts_etamax = cms.double(5.0),
                                  met_sources =  cms.vstring("slimmedMETs"),
                                  pv_sources = cms.vstring("offlineSlimmedPrimaryVertices"),
                                  trigger_bits = cms.InputTag("TriggerResults","","HLT"),
                                  trigger_prefixes = cms.vstring(#"HLT_IsoMu", "HLT_Mu",
                                                                 #"HLT_L1SingleMu", "HLT_L2Mu",
                                                                 #"HLT_Ele",
                                                                 "HLT_",
                                                                 #"HLT_DoubleMu", "HLT_DoubleEle"
	                                                         ),
                                  
)

process.p = cms.Path(process.MyNtuple)




