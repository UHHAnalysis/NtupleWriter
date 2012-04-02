# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.CaloJetParameters_cfi import CaloJetParameters
from RecoJets.JetProducers.GenJetParameters_cfi import GenJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import PFJetParameters
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.CATopJetParameters_cfi import CATopJetParameters
from RecoJets.JetProducers.AnomalousCellParameters_cfi import AnomalousCellParameters

virtualjet_parameters = cms.PSet(jetAlgorithm=cms.string("CambridgeAachen"), rParam=cms.double(0.8), jetCollInstanceName=cms.string("subjets"))

SubJetParameters.nSubjets = cms.int32(2)

caPrunedCaloJets = cms.EDProducer("SubJetProducer",
    SubJetParameters,
    virtualjet_parameters,
    AnomalousCellParameters,
    CaloJetParameters
)

caPrunedPFJets = cms.EDProducer("SubJetProducer",
    SubJetParameters,
    virtualjet_parameters,
    AnomalousCellParameters,
    PFJetParameters
)

caPrunedGenJets = cms.EDProducer("SubJetProducer",
    SubJetParameters,
    virtualjet_parameters,
    AnomalousCellParameters,
    GenJetParameters
)

#2.b. toptag
caTopCaloJets = cms.EDProducer("CATopJetProducer",
       CATopJetParameters,
       CaloJetParameters,
       AnomalousCellParameters,
       jetAlgorithm=cms.string("CambridgeAachen"), rParam=cms.double(0.8)
)

caTopPFJets = cms.EDProducer("CATopJetProducer",
       CATopJetParameters,
       PFJetParameters,
       AnomalousCellParameters,
       jetAlgorithm=cms.string("CambridgeAachen"), rParam=cms.double(0.8)
)

caTopGenJets = cms.EDProducer("CATopJetProducer",
       CATopJetParameters,
       GenJetParameters,
       AnomalousCellParameters,
       jetAlgorithm=cms.string("CambridgeAachen"), rParam=cms.double(0.8)
)

topjet_seq = cms.Sequence(caPrunedCaloJets * caPrunedPFJets * caTopPFJets * caTopCaloJets)
topjet_genseq = cms.Sequence(caPrunedGenJets * caTopGenJets)

