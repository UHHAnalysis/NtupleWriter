# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("MYGENTEST")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING'
#process.options   = cms.untracked.PSet(
#  wantSummary = cms.untracked.bool(True)
#)

process.source = cms.Source("PoolSource",

   fileNames = cms.untracked.vstring('file:/scratch/hh/lustre/cms/user/peiffer/LQGen/LQ1.root'),

   skipEvents = cms.untracked.uint32(0),
   duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


#NtupleWriter

process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string('genLQ.root'), 
                                  doElectrons = cms.bool(False),
                                  doMuons = cms.bool(False),
                                  doTaus = cms.bool(False),
                                  doJets = cms.bool(False),
                                  doTopJets = cms.bool(False),
                                  doPhotons = cms.bool(False),
                                  doMET = cms.bool(False),
                                  doPV = cms.bool(False),
                                  doGenInfo = cms.bool(True),
                                  doTrigger = cms.bool(False)
                                  
                                  
)

process.p = cms.Path(process.MyNtuple)

