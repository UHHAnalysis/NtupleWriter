# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("MYGENTEST")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.options   = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",

   fileNames = cms.untracked.vstring('file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ1.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ2.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ3.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ4.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ5.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ6.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ7.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ8.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ9.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ10.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ11.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ12.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ13.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ14.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ15.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ16.root',
                                     'file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ17.root'
                                     ),

   skipEvents = cms.untracked.uint32(0),
   duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *

process.caTopTagGen = cms.EDProducer(
    "CATopJetProducer",
    GenJetParameters.clone(src = cms.InputTag("genParticlesForJetsNoNu"),
                           doAreaFastjet = cms.bool(False),
                           doRhoFastjet = cms.bool(False)),
    AnomalousCellParameters,
    CATopJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )

#NtupleWriter

process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string('/scratch/hh/lustre/cms/user/peiffer/genLQ_test.root'), 
                                  doElectrons = cms.bool(False),
                                  doMuons = cms.bool(False),
                                  doTaus = cms.bool(False),
                                  doJets = cms.bool(False),
                                  doTopJets = cms.bool(False),
                                  doGenTopJets = cms.bool(True),
                                  gentopjet_sources = cms.vstring("caTopTagGen" ),
                                  gentopjet_ptmin = cms.double(150.0), 
                                  gentopjet_etamax = cms.double(5.0),
                                  doPhotons = cms.bool(False),
                                  doMET = cms.bool(False),
                                  doPV = cms.bool(False),
                                  doGenInfo = cms.bool(True),
                                  doAllGenParticles = cms.bool(True),
                                  doTrigger = cms.bool(False),
                                  doLumiInfo = cms.bool(False)
                                  
                                  
)

process.p = cms.Path( process.genParticlesForJetsNoNu*
                      process.caTopTagGen*                        
                      process.MyNtuple)

