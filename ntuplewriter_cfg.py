# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("MYPATTEST")

isttbar = False
issingletop = False
isvjets = False
isdata = True
nickname = "MC_TTbar"
#nickname= "DATA"
if nickname.startswith('MC_'): isdata = False
if nickname.find("_TTbar") > -1 or nickname.find('_ZP') > -1 : isttbar = True
if nickname.find("_ST_tchan") > -1 : issingletop = True
if nickname.find("_Wjets") > -1 or nickname.find('_Vqq') > -1 or nickname.find('_Wc') > -1 or nickname.find('_Zjets') > -1 : isvjets = True
tmp = os.getenv('CMSSW_BASE')
if 'CMSSW_4_1_' in tmp: cmssw_version = '41'
elif 'CMSSW_4_2_' in tmp: cmssw_version = '42'
elif 'CMSSW_4_4_' in tmp: cmssw_version = '44'
else: raise RuntimeError, 'could not determine cmssw version (cmssw correctly set up?)'

print "detected CMSSW version %s" % cmssw_version
if isdata: print "assuming to run over data"
else: print "assuming to run over MC; (single top, ttbar, vjets) = (", issingletop, isttbar, isvjets, ")"

# for global tag, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
def globaltag():
   if isdata and cmssw_version=='41': return 'GR_R_41_V0::All'
   if not isdata and cmssw_version=='41': return 'START41_V0::All'
   if isdata and cmssw_version=='42': return 'GR_R_42_V25::All'
   if not isdata and cmssw_version=='42': return 'START42_V17::All'
   if isdata and cmssw_version=='44': return 'GR_R_44_V15::All'
   if not isdata and cmssw_version=='44': return 'START44_V13::All'

   raise RuntimeError, "could not determine global tag to use!"


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.options   = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/data/Run2011B/Jet/AOD/19Nov2011-v1/0000/985E38E8-3B15-E111-87A8-0015178C4970.root',
                            #                                  'dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/data/Run2011B/Jet/AOD/19Nov2011-v1/0000/4449A169-1A15-E111-B25F-0015178C4CE8.root'),
   fileNames = cms.untracked.vstring('dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6-START44_V5-v1/0000/68423CA2-BC04-E111-B90B-0018F3D0968C.root'),
   #fileNames = cms.untracked.vstring('dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6-START44_V5-v1/0000/7C7C20EE-BB04-E111-A833-001A92971BDA.root'),
   #fileNames = cms.untracked.vstring('dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/mc/Fall11/W2Jets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START44_V9B-v1/0000/9A4D5BEE-F536-E111-9F71-003048C692CA.root'),                         
   #fileNames = cms.untracked.vstring('dcap://dcache-cms-dcap.desy.de//pnfs/desy.de/cms/tier2/store/mc/Fall11/QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6/AODSIM/PU_S6_START42_V14B-v2/0000/50404FDF-6214-E111-8248-78E7D164BFC8.root'),            
   skipEvents = cms.untracked.uint32(0),
   duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *' , *patEventContent)
                               )

if isdata and cmssw_version == '41': # note: there are no residual corrections in 42(?!);
                                     # needs to be done with recorrector in SKITA
    JECToProcess = cms.vstring(['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual'])
else:
    JECToProcess = cms.vstring(['L1FastJet','L2Relative', 'L3Absolute'])

# for some jet types, no residual correction is available; use l2l3 only:
#l2l3jec = cms.vstring(['L2Relative', 'L3Absolute'])

inputJetCorrLabel = ('AK5PFchs', JECToProcess)

# 1. general setup
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globaltag()
print "using globaltag %s" % str(process.GlobalTag.globaltag)
process.load("Configuration.StandardSequences.MagneticField_cff")



# 3. configure PF2PAT
#process.load("PhysicsTools.PFCandProducer.PF2PAT_cff")

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJets = kt4PFJets.clone()
process.kt6PFJets.rParam = 0.6
process.kt6PFJets.doRhoFastjet = True

process.kt6PFJetsRestrictedRho = process.kt6PFJets.clone() 
process.kt6PFJetsRestrictedRho.Rho_EtaMax = cms.double(1.9)
process.kt6PFJetsRestrictedRho.Ghost_EtaMax = cms.double(2.5)



# Get a list of good primary vertices, in 42x, these are DAF vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )


#top jets

process.load("AnalysisCode.NtupleWriter.Topjets_cfi")

process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8) )



process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
from PhysicsTools.PatAlgos.tools.pfTools import *
usePF2PAT(process, runPF2PAT=True, jetAlgo="AK5", runOnMC = not isdata)
process.pfPileUp.Enable = True
process.pfPileUp.Vertices = 'goodOfflinePrimaryVertices'
process.pfJets.doAreaFastjet = True
process.pfJets.doRhoFastjet = False
process.kt6PFJets.src = process.pfJets.src
process.kt6PFJetsRestrictedRho.src = process.pfJets.src
process.patJetCorrFactors.payload = 'AK5PFchs'
process.patJetCorrFactors.levels = JECToProcess
process.patJetCorrFactors.rho = cms.InputTag("kt6PFJets", "rho")
process.pfNoTau.enable = cms.bool(True) ##change this????
    

from PhysicsTools.PatAlgos.producersLayer1.photonProducer_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
if isdata:
   process.patDefaultSequence *= patPhotons
else:  process.patDefaultSequence *= makePatPhotons
process.patDefaultSequence *= selectedPatPhotons
process.patDefaultSequence *= countPatPhotons

process.pfIsolatedMuons.isolationCut = cms.double(0.2) 


process.patMuons.usePV = cms.bool(False)
process.patElectrons.usePV = cms.bool(False)
process.patJets.addTagInfos = cms.bool(True)

process.pfIsolatedMuonsLoose = process.pfIsolatedMuons.clone(isolationCut = cms.double(float("inf")))
process.patMuonsLoose = process.patMuons.clone(pfMuonSource = cms.InputTag("pfIsolatedMuonsLoose"), genParticleMatch = cms.InputTag("muonMatchLoose"))
# use pf isolation, but do not change matching (hell this "PAT" is really f*** up).
tmp = process.muonMatch.src
adaptPFMuons(process, process.patMuonsLoose, "")
process.muonMatch.src = tmp


process.pfSelectedElectrons.cut = cms.string('et > 15. & abs(eta) < 2.5') 
process.pfSelectedMuons.cut = cms.string('pt > 10. & abs(eta) < 2.5') 
process.pfSelectedPhotons.cut = cms.string('pt > 10.')

#jet and tau cuts are not propagated to ntuplewriter
process.selectedPatJets.cut = cms.string('pt > 10. & abs(eta) < 5.0')
process.selectedPatTaus.cut = cms.string('pt > 10. & abs(eta) < 5.0')

process.muonMatchLoose = process.muonMatch.clone(src = cms.InputTag("pfIsolatedMuonsLoose"))

process.pfIsolatedElectronsLoose = process.pfIsolatedElectrons.clone(isolationCut = cms.double(float("inf")))
process.patElectronsLoose = process.patElectrons.clone(pfElectronSource = cms.InputTag("pfIsolatedElectronsLoose"))
adaptPFElectrons(process, process.patElectronsLoose, "")

process.looseLeptonSequence = cms.Sequence(
    process.pfIsolatedMuonsLoose +
    process.muonMatchLoose +
    process.patMuonsLoose +
    process.pfIsolatedElectronsLoose +
    process.patElectronsLoose
    )
    
if isdata: process.looseLeptonSequence.remove(process.muonMatchLoose)

# note: we cannot use this for the moment: some modules require a primary vertex (even if is is fake!!)
#for module in (process.patMuons, process.patMuonsLoose, process.patElectrons, process.patElectronsLoose):
#    module.pvSrc = "goodOfflinePrimaryVertices"
#for module in (process.impactParameterTagInfos, process.impactParameterTagInfosAOD, process.softMuonTagInfosAOD, process.softMuonTagInfos):
#    module.primaryVertex = "goodOfflinePrimaryVertices"
#for module in (process.patJetCorrFactors,):
#    module.primaryVertices = "goodOfflinePrimaryVertices"
#for module in (process.pfElectronsFromVertex,process.pfMuonsFromVertex):
#    module.vertices = "goodOfflinePrimaryVertices"



process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

process.HBHENoiseFilterResultProducer.minRatio = cms.double(-999)
process.HBHENoiseFilterResultProducer.maxRatio = cms.double(999)
process.HBHENoiseFilterResultProducer.minHPDHits = cms.int32(17)
process.HBHENoiseFilterResultProducer.minRBXHits = cms.int32(999)
process.HBHENoiseFilterResultProducer.minHPDNoOtherHits = cms.int32(10)
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(10)
process.HBHENoiseFilterResultProducer.minHighEHitTime = cms.double(-9999.0)
process.HBHENoiseFilterResultProducer.maxHighEHitTime = cms.double(9999.0)
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(-999.0)
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(999999.)
process.HBHENoiseFilterResultProducer.useTS4TS5 = cms.bool(True)


# add CA PF jets
process.caTopPFJets.src = process.pfJets.src
addJetCollection(process, 
                 cms.InputTag('caTopPFJets'),
                 algoLabel = 'TopTag',  typeLabel = 'PF', doJTA=False,
                 doBTagging=False, jetCorrLabel = inputJetCorrLabel, doType1MET=False,#use ak5 PFjet corrections as described in AN2011_194
                 genJetCollection = cms.InputTag("ca8GenJets"), doJetID = False)

# *** jet pruning ala Jott:
process.caPrunedPFJets.src = process.pfJets.src
addJetCollection(process,
                 cms.InputTag('caPrunedPFJets'),
                 algoLabel = 'Pruned', typeLabel = 'PF', doJTA=True,
                 doBTagging=True, jetCorrLabel = inputJetCorrLabel, doType1MET=False,#use ak5 PFjet corrections as described in AN2011_194
                 genJetCollection = cms.InputTag("ca8GenJets"), doJetID = False)

#configure jet parton flavour determination
for jetcollname in ('TopTagPF', 'PrunedPF'):
   jetcoll = getattr(process, 'patJets'+jetcollname)
   jetcoll.embedGenJetMatch = cms.bool(False)
   jetcoll.addDiscriminators = cms.bool(False)
   # jet flavor stuff: match also the top quark.
   getattr(process, 'patJetPartonAssociation'+jetcollname).coneSizeToAssociate = cms.double(0.8)
   getattr(process, 'patJetPartonAssociation'+jetcollname).doPriority = cms.bool(True)
   getattr(process, 'patJetPartonAssociation'+jetcollname).priorityList = cms.vint32(6)
   #4 = match heaviest flavor
   getattr(process, 'patJetFlavourAssociation'+jetcollname).definition = cms.int32(4)
   getattr(process, 'patJetFlavourAssociation'+jetcollname).physicsDefinition = cms.bool(False)
   getattr(process, 'patJetPartonMatch'+jetcollname).mcPdgId = cms.vint32(1,2,3,4,5,6,21)
   getattr(process, 'patJetPartonMatch'+jetcollname).maxDeltaR = cms.double(0.8)

for icorr in [process.patJetCorrFactorsTopTagPF,
              process.patJetCorrFactorsPrunedPF
              ] :
    icorr.rho = cms.InputTag("kt6PFJets", "rho")   

# 5. compile main sequence

#process.patPF2PATSequence.replace(process.pfJets, process.kt6PFJets * process.kt6PFJetsRestrictedRho * process.pfJets)
process.patPF2PATSequence.replace(process.pfJets, process.kt6PFJets * process.kt6PFJetsRestrictedRho * process.pfJets * process.topjet_seq)

process.main_seq = cms.Sequence(process.goodOfflinePrimaryVertices)

process.main_seq *= process.patPF2PATSequence* process.looseLeptonSequence

process.main_genseq = cms.Sequence(process.genParticlesForJets * process.ca8GenJets * process.topjet_genseq )
#process.main_genseq = cms.Sequence(process.genParticlesForJets * process.ca8GenJets * process.subjet_genseq * process.topjet_genseq)


process.p = cms.Path()

if isdata: process.p *= process.HBHENoiseFilterResultProducer
if not isdata:  process.p *= process.main_genseq
process.p *= process.main_seq

# Apply jet ID to all of the jets upstream. We aren't going to screw around
# with this, most likely. So, we don't really to waste time with it
# at the analysis level. 
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("patJets")
                                   )

process.goodPatJetsPrunedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                              filterParams = pfJetIDSelector.clone(),
                                              src = cms.InputTag("patJetsPrunedPF")
                                              )


process.goodPatJetsTopTagPF  = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                             filterParams = pfJetIDSelector.clone(),
                                             src = cms.InputTag("patJetsTopTagPF")
                                             )



#NtupleWriter

process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string('ttbar.root'), 
                                  doElectrons = cms.bool(True),
                                  doMuons = cms.bool(True),
                                  doTaus = cms.bool(True),
                                  doJets = cms.bool(True),
                                  doTopJets = cms.bool(True),
                                  doPhotons = cms.bool(True),
                                  doMET = cms.bool(True),
                                  doPV = cms.bool(True),
                                  doGenInfo = cms.bool(not isdata),
                                  electron_sources = cms.vstring("patElectrons","patElectronsLoose"),
                                  muon_sources = cms.vstring("patMuons","patMuonsLoose"),
                                  tau_sources = cms.vstring("patTaus"),
                                  tau_ptmin = cms.double(10.0),
                                  tau_etamax = cms.double(5.0),
                                  jet_sources = cms.vstring("goodPatJets"),
                                  jet_ptmin = cms.double(10.0),
                                  jet_etamax = cms.double(5.0),
                                  topjet_sources = cms.vstring("goodPatJetsTopTagPF","goodPatJetsPrunedPF"),
                                  topjet_ptmin = cms.double(150.0), 
                                  topjet_etamax = cms.double(5.0),
                                  photon_sources = cms.vstring("patPhotons"),
                                  met_sources =  cms.vstring("patMETs"),
                                  pv_sources = cms.vstring("goodOfflinePrimaryVertices", "offlinePrimaryVertices"),
                                  trigger_prefixes = cms.vstring(#"HLT_IsoMu", "HLT_Mu",
                                                                 #"HLT_L1SingleMu", "HLT_L2Mu",
                                                                 #"HLT_Ele",
                                                                 "HLT_Jet", "HLT_HT",
                                                                 #"HLT_DoubleMu", "HLT_DoubleEle"
	                                                         ),
                                  
)

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('tree3.root')
#)
process.p *=process.goodPatJets
process.p *=process.goodPatJetsPrunedPF
process.p *=process.goodPatJetsTopTagPF
process.p *= process.MyNtuple

#comment this line in to produce pattuples:
#process.outpath = cms.EndPath(process.out)

process.options.wantSummary = True

#open('dump.py', 'w').write(process.dumpPython())
