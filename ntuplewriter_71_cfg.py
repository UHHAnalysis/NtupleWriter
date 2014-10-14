import FWCore.ParameterSet.Config as cms
###############################################
useMiniAOD = True

# AOD
pfcandidates          = 'particleFlow'
chsstring             = 'pfNoPileUpJME'
genjetparticles       = 'genParticles'
importantgenparticles = 'genParticles'
tracks                = 'generalTracks'
vertices              = 'offlinePrimaryVertices'
mergedvertices        = 'inclusiveMergedVertices' 
mergedvertices2       = '' 
primaryvertices       = 'offlinePrimaryVertices'

#miniAOD
if useMiniAOD:
  pfcandidates          = 'packedPFCandidates'
  genjetparticles       = 'packedGenParticles'
  importantgenparticles = 'prunedGenParticles'
  tracks                = 'unpackedTracksAndVertices'
  vertices              = 'unpackedTracksAndVertices'
  mergedvertices        = 'unpackedTracksAndVertices'
  mergedvertices2       = 'secondary'
  primaryvertices       = 'offlineSlimmedPrimaryVertices'


print 'useMiniAOD = '+str(useMiniAOD)
print ' pfcandidates          = '+pfcandidates         
print ' genjetparticles       = '+genjetparticles      
print ' importantgenparticles = '+importantgenparticles
print ' tracks                = '+tracks               
print ' vertices              = '+vertices             
print ' mergedvertices        = '+mergedvertices       
print ' mergedvertices2       = '+mergedvertices2      
print ' primaryvertices       = '+primaryvertices 

###############################################
# SETUP
process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) , allowUnscheduled = cms.untracked.bool(True) )


###############################################
# SOURCE

#process.source = cms.Source("PoolSource",
#                            fileNames  = cms.untracked.vstring(__FILE_NAMES__),
#                            skipEvents = cms.untracked.uint32(__SKIP_EVENTS__)
#)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(__MAX_EVENTS__))

process.source = cms.Source("PoolSource",
                            fileNames  = cms.untracked.vstring("/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU40bx25_POSTLS170_V7-v2/00000/00800BE3-E826-E411-AD01-20CF3019DEE9.root"),
                            skipEvents = cms.untracked.uint32(0)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))



###############################################
# OUT
outfile = 'TestJetsAOD_706p1.root'
if useMiniAOD:
  outfile = 'TestJetsMiniAOD_706p1.root'

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(outfile),
                               outputCommands = cms.untracked.vstring([
                                #'keep recoCATopJetTagInfos_*_*_*',  #recoCATopJetTagInfos_CATopTagInfos__USER etc.
                                #'keep patJets_patJetsAK15PF_*_*',
                                #'keep patJets_patJetsAK15PFfiltered*_*_*',
                                #'keep patJets_patJetsAK15CHS_*_*',
                                #'keep patJets_patJetsAK15CHSfiltered*_*_*',
                                #'keep patJets_patJetsCA8PF_*_*',
                                #'keep patJets_patJetsCA8PFpruned*_*_*',
                                #'keep patJets_patJetsCA8CHS_*_*',
                                #'keep patJets_patJetsCA8CHSpruned*_*_*',
                                #'keep patJets_patJetsCMSTopTagCHS*_*_*',
                                #'keep patJets_patJetsCMSTopTagFJCHS*_*_*',
                                #'keep patJets_patJetsHEPTopTagCHS*_*_*',
                                'drop *',
                                                                       ])
)

###############################################
# RECO AND GEN SETUP
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START70_V6::All'

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.fixedGridRhoFastjetAll.pfCandidatesTag = pfcandidates


from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.caTopTaggers_cff import *

###############################################
# GEN JETS


from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

process.ak3GenJets = ak5GenJets.clone( rParam = cms.double(0.3),
                                           src = cms.InputTag(genjetparticles))

process.ca8GenJets = ca4GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag(genjetparticles))

process.ak15GenJets = ak5GenJets.clone( rParam = cms.double(1.5),
                                           src = cms.InputTag(genjetparticles))

#ADD NO NU

process.ak15GenJetsFiltered = ak5GenJets.clone(
  rParam = cms.double(1.5),
  src = cms.InputTag(genjetparticles),
  useFiltering = cms.bool(True),
  nFilt = cms.int32(3),
  rFilt = cms.double(0.3),
  writeCompound = cms.bool(True),
  jetCollInstanceName=cms.string("SubJets")
  )

process.ak15GenJetsFilteredNfilt2 = ak5GenJets.clone(
  rParam = cms.double(1.5),
  src = cms.InputTag(genjetparticles),
  useFiltering = cms.bool(True),
  nFilt = cms.int32(2),
  rFilt = cms.double(0.3),
  writeCompound = cms.bool(True),
  jetCollInstanceName=cms.string("SubJets")
  )

process.ca8GenJetsTrimmed = ca4GenJets.clone(
  rParam = cms.double(0.8),
  src = cms.InputTag(genjetparticles),
  useTrimming = cms.bool(True),
  trimPtFracMin = cms.double(0.05),
  rFilt = cms.double(0.2),
  writeCompound = cms.bool(True),
  jetCollInstanceName=cms.string("SubJets")
  )

process.ca8GenJetsPruned = ca4GenJets.clone(
  rParam = cms.double(0.8),
  src = cms.InputTag(genjetparticles),
  usePruning = cms.bool(True),
  zcut = cms.double(0.1),
  rcut_factor = cms.double(0.5),
  nFilt = cms.int32(2),
  writeCompound = cms.bool(True),
  jetCollInstanceName=cms.string("SubJets")
  )

process.cmsTopTagGEN = cms.EDProducer(
    "CATopJetProducer",
    GenJetParameters.clone(src = cms.InputTag(genjetparticles),
                           doAreaFastjet = cms.bool(False),
                           doRhoFastjet = cms.bool(False)),
    AnomalousCellParameters,
    CATopJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )

process.hepTopTagGEN = process.cmsTopTagGEN.clone(
    rParam = cms.double(1.5),
    tagAlgo = cms.int32(2), #2=fastjet heptt
    muCut = cms.double(0.8),
    maxSubjetMass = cms.double(30.0),
    useSubjetMass = cms.bool(False),
)

###############################################
# RECO JETS


if useMiniAOD:
  process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
  chsstring = 'chs'
process.ak4PFJets.src = pfcandidates
process.ca4PFJets.src = pfcandidates

process.ca8PFJets  = process.ca4PFJets.clone(rParam = 0.8,  doAreaFastjet = True)
process.ak15PFJets = process.ak4PFJets.clone(rParam = 1.5,  doAreaFastjet = True)

process.ca8CHSJets  = process.ca8PFJets.clone (src = chsstring)

process.ak15CHSJets = process.ak15PFJets.clone(src = chsstring)

process.load('RecoJets.JetProducers.ak4PFJetsFiltered_cfi')
process.ak15PFJetsFiltered  = process.ak4PFJetsFiltered.clone(rParam=1.5,  doAreaFastjet = True, src = pfcandidates)
process.ak15CHSJetsFiltered = process.ak4PFJetsFiltered.clone(rParam=1.5,  doAreaFastjet = True, src = chsstring)

process.ak15PFJetsFilteredNfilt2  = process.ak4PFJetsFiltered.clone(rParam=1.5, nFilt=2, rFilt = 0.3,  doAreaFastjet = True, src = pfcandidates)
process.ak15CHSJetsFilteredNfilt2 = process.ak4PFJetsFiltered.clone(rParam=1.5, nFilt=2, rFilt = 0.3,  doAreaFastjet = True, src = chsstring)

process.ak15PFJetsMassDropFiltered  = process.ak4PFJetsMassDropFiltered.clone(rParam=1.5,  doAreaFastjet = True, src = pfcandidates)
process.ak15CHSJetsMassDropFiltered = process.ak4PFJetsMassDropFiltered.clone(rParam=1.5,  doAreaFastjet = True, src = chsstring)

process.load('RecoJets.JetProducers.ak4PFJetsTrimmed_cfi')
process.ca8PFJetsTrimmed  = process.ak4PFJetsTrimmed.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), rFilt=0.2, trimPtFracMin=0.05,  doAreaFastjet = True, src = pfcandidates)
process.ca8CHSJetsTrimmed = process.ak4PFJetsTrimmed.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), rFilt=0.2, trimPtFracMin=0.05,  doAreaFastjet = True, src = chsstring)

process.load('RecoJets.JetProducers.ak4PFJetsPruned_cfi')
process.ca8PFJetsPruned  = process.ak4PFJetsPruned.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), doAreaFastjet = True, src = pfcandidates)
process.ca8CHSJetsPruned = process.ak4PFJetsPruned.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), doAreaFastjet = True, src = chsstring)

process.ca8PFJetsPrunedNsub3  = process.ak4PFJetsPruned.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), nFilt=3, doAreaFastjet = True, src = pfcandidates)
process.ca8CHSJetsPrunedNsub3 = process.ak4PFJetsPruned.clone(rParam=0.8, jetAlgorithm = cms.string("CambridgeAachen"), nFilt=3, doAreaFastjet = True, src = chsstring)

# CATopJet PF Jets with adjacency 
#process.cmsTopTagCHS = cmsTopTagPFJetsCHS.clone()
process.cmsTopTagCHS = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag(chsstring),
                           doAreaFastjet = cms.bool(True),
                           doRhoFastjet = cms.bool(False),
                           jetPtMin = cms.double(100.0)
                           ),
    AnomalousCellParameters,
    CATopJetParameters.clone( jetCollInstanceName = cms.string("SubJets"), 
                              verbose = cms.bool(False),
                              algorithm = cms.int32(1),                 # 0 = KT, 1 = CA, 2 = anti-KT
                              tagAlgo = cms.int32(0),                   #0=legacy top
                              useAdjacency = cms.int32(2),              # modified adjacency
                              centralEtaCut = cms.double(2.5),          # eta for defining "central" jets
                              sumEtBins = cms.vdouble(0,1600,2600),     # sumEt bins over which cuts vary. vector={bin 0 lower bound, bin 1 lower bound, ...}
                              rBins = cms.vdouble(0.8,0.8,0.8),         # Jet distance paramter R. R values depend on sumEt bins.
                              ptFracBins = cms.vdouble(0.05,0.05,0.05), # minimum fraction of central jet pt for subjets (deltap)
                              deltarBins = cms.vdouble(0.19,0.19,0.19), # Applicable only if useAdjacency=1. deltar adjacency values for each sumEtBin
                              nCellBins = cms.vdouble(1.9,1.9,1.9), 
                            ),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )

process.cmsTopTagFJCHS = process.cmsTopTagCHS.clone(
    tagAlgo = cms.int32(1), # 1=fastjet cmstagger
    ptFrac = cms.double(0.05),
    rFrac = cms.double(0.4),
    adjacencyParam = cms.double(0.0004),
)

process.hepTopTagCHS = process.cmsTopTagCHS.clone(
    rParam = cms.double(1.5),
    tagAlgo = cms.int32(2), #2=fastjet heptt
    muCut = cms.double(0.8),
    maxSubjetMass = cms.double(30.0),
    useSubjetMass = cms.bool(False),
)

process.ca15CHSJets = process.hepTopTagCHS.clone ( writeCompound = cms.bool(False), jetCollInstanceName = cms.string('') )


process.CATopTagInfos = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("cmsTopTagCHS"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )

process.cmsTopTagFJCHSTagInfos = process.CATopTagInfos.clone(
  src = cms.InputTag("cmsTopTagFJCHS")
)

process.hepTopTagCHSTagInfos = process.CATopTagInfos.clone(
  src = cms.InputTag("hepTopTagCHS")
)

###############################################
# PAT JETS
from PhysicsTools.PatAlgos.tools.jetTools import *

process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

# fix JTA (see https://github.com/cms-sw/cmssw/tree/CMSSW_7_0_X/RecoJets/JetAssociationProducers/python)
process.ak5JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak4PFJets")
process.ak5JetTracksAssociatorAtVertexPF.tracks = cms.InputTag(tracks)
process.impactParameterTagInfos.primaryVertex = cms.InputTag(vertices)
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag(mergedvertices,mergedvertices2,"")
process.combinedSecondaryVertex.trackMultiplicityMin = 1

#AK15PF
addJetCollection(
    process,
    labelName = 'AK15PF',
    jetSource = cms.InputTag('ak15PFJets'),
    algo = 'ak15',
    rParam = 1.5,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],

    )
process.patJetPartonMatchAK15PF.matched= importantgenparticles
process.patJetCorrFactorsAK15PF.primaryVertices = primaryvertices
process.patJetGenJetMatchAK15PF.matched = 'ak15GenJets'#'slimmedGenJets'
process.patJetPartons.particles = importantgenparticles
process.jetTracksAssociatorAtVertexAK15PF=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak15PFJets'), coneSize = 1.5)
process.secondaryVertexTagInfosAK15PF.trackSelection.jetDeltaRMax = cms.double(1.5) # default is 0.3
process.secondaryVertexTagInfosAK15PF.vertexCuts.maxDeltaRToJetAxis = cms.double(1.5) # default is 0.5
process.combinedSecondaryVertexAK15PF= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexAK15PF.trackSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexAK15PF.trackPseudoSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexBJetTagsAK15PF.jetTagComputer = cms.string('combinedSecondaryVertexAK15PF')

#AK15CHS
addJetCollection(
    process,
    labelName = 'AK15CHS',
    jetSource = cms.InputTag('ak15CHSJets'),
    algo = 'ak15',
    rParam = 1.5,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchAK15CHS.matched=importantgenparticles
process.patJetCorrFactorsAK15CHS.primaryVertices = primaryvertices
process.patJetGenJetMatchAK15CHS.matched = 'ak15GenJets'#'slimmedGenJets'
process.jetTracksAssociatorAtVertexAK15CHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak15CHSJets'), coneSize = 1.5)
process.secondaryVertexTagInfosAK15CHS.trackSelection.jetDeltaRMax = cms.double(1.5) # default is 0.3
process.secondaryVertexTagInfosAK15CHS.vertexCuts.maxDeltaRToJetAxis = cms.double(1.5) # default is 0.5
process.combinedSecondaryVertexAK15CHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexAK15CHS.trackSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexAK15CHS.trackPseudoSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexBJetTagsAK15CHS.jetTagComputer = cms.string('combinedSecondaryVertexAK15CHS')


#AK15PF Filt
addJetCollection(
    process,
    labelName = 'AK15PFfiltered',
    jetSource = cms.InputTag('ak15PFJetsFiltered'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchAK15PFfiltered.matched=importantgenparticles
process.patJetCorrFactorsAK15PFfiltered.primaryVertices = primaryvertices
process.patJetGenJetMatchAK15PFfiltered.matched = 'ak15GenJetsFiltered'#'slimmedGenJets'
process.patJetPartonMatchAK15PFfiltered.matched = importantgenparticles
process.jetTracksAssociatorAtVertexAK15PFfiltered=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak15PFJetsFiltered'), coneSize = 1.5)
process.secondaryVertexTagInfosAK15PFfiltered.trackSelection.jetDeltaRMax = cms.double(1.5) # default is 0.3
process.secondaryVertexTagInfosAK15PFfiltered.vertexCuts.maxDeltaRToJetAxis = cms.double(1.5) # default is 0.5
process.combinedSecondaryVertexAK15PFfiltered= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexAK15PFfiltered.trackSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexAK15PFfiltered.trackPseudoSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexBJetTagsAK15PFfiltered.jetTagComputer = cms.string('combinedSecondaryVertexAK15PFfiltered')

addJetCollection(
    process,
    labelName = 'AK15PFfilteredSubjets',
    jetSource = cms.InputTag('ak15PFJetsFiltered','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchAK15PFfilteredSubjets.matched=importantgenparticles
process.patJetCorrFactorsAK15PFfilteredSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchAK15PFfilteredSubjets.matched = 'ak3GenJets'#slimmedGenJets'
process.patJetPartonMatchAK15PFfilteredSubjets.matched = importantgenparticles

process.patJetsAK15PFfilteredPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsAK15PFfiltered" ),
    subjetSrc=cms.InputTag("patJetsAK15PFfilteredSubjets")
      )

#AK15CHS Filt
addJetCollection(
    process,
    labelName = 'AK15CHSfiltered',
    jetSource = cms.InputTag('ak15CHSJetsFiltered'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )

process.patJetPartonMatchAK15CHSfiltered.matched=importantgenparticles
process.patJetCorrFactorsAK15CHSfiltered.primaryVertices = primaryvertices
process.patJetGenJetMatchAK15CHSfiltered.matched = 'ak15GenJetsFiltered'#'slimmedGenJets'
process.patJetPartonMatchAK15CHSfiltered.matched = importantgenparticles
process.jetTracksAssociatorAtVertexAK15CHSfiltered=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak15CHSJetsFiltered'), coneSize = 1.5)
process.secondaryVertexTagInfosAK15CHSfiltered.trackSelection.jetDeltaRMax = cms.double(1.5) # default is 0.3
process.secondaryVertexTagInfosAK15CHSfiltered.vertexCuts.maxDeltaRToJetAxis = cms.double(1.5) # default is 0.5
process.combinedSecondaryVertexAK15CHSfiltered= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexAK15CHSfiltered.trackSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexAK15CHSfiltered.trackPseudoSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexBJetTagsAK15CHSfiltered.jetTagComputer = cms.string('combinedSecondaryVertexAK15CHSfiltered')

addJetCollection(
    process,
    labelName = 'AK15CHSfilteredSubjets',
    jetSource = cms.InputTag('ak15CHSJetsFiltered','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchAK15CHSfilteredSubjets.matched=importantgenparticles
process.patJetCorrFactorsAK15CHSfilteredSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchAK15CHSfilteredSubjets.matched = 'ak3GenJets'#slimmedGenJets'
process.patJetPartonMatchAK15CHSfilteredSubjets.matched = importantgenparticles

process.patJetsAK15CHSfilteredPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsAK15CHSfiltered" ),
    subjetSrc=cms.InputTag("patJetsAK15CHSfilteredSubjets")
      )



# patJetsCA8PF
addJetCollection(
    process,
    labelName = 'CA8PF',
    jetSource = cms.InputTag('ca8PFJets'),
    algo = 'ca8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA8PF.matched=importantgenparticles
process.patJetCorrFactorsCA8PF.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8PF.matched = 'ca8GenJets'#'slimmedGenJets'
process.patJetPartons.particles = importantgenparticles
process.jetTracksAssociatorAtVertexCA8PF=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8PFJets'), coneSize = 0.8)
process.secondaryVertexTagInfosCA8PF.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCA8PF.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCA8PF= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA8PF.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCA8PF.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCA8PF.jetTagComputer = cms.string('combinedSecondaryVertexCA8PF')

#patJetsCA8CHS
addJetCollection(
    process,
    labelName = 'CA8CHS',
    jetSource = cms.InputTag('ca8CHSJets'),
    algo = 'ca8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA8CHS.matched=importantgenparticles
process.patJetCorrFactorsCA8CHS.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8CHS.matched = 'ca8GenJets'#'slimmedGenJets'
process.jetTracksAssociatorAtVertexCA8CHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8CHSJets'), coneSize = 0.8)
process.secondaryVertexTagInfosCA8CHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCA8CHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCA8CHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA8CHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCA8CHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCA8CHS.jetTagComputer = cms.string('combinedSecondaryVertexCA8CHS')


# patJetsCA8PFpruned
addJetCollection(
    process,
    labelName = 'CA8PFpruned',
    jetSource = cms.InputTag('ca8PFJetsPruned'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )

process.patJetPartonMatchCA8PFpruned.matched=importantgenparticles
process.patJetCorrFactorsCA8PFpruned.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8PFpruned.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchCA8PFpruned.matched = importantgenparticles

addJetCollection(
    process,
    labelName = 'CA8PFprunedSubjets',
    jetSource = cms.InputTag('ca8PFJetsPruned','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )

process.patJetPartonMatchCA8PFprunedSubjets.matched=importantgenparticles
process.patJetCorrFactorsCA8PFprunedSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8PFprunedSubjets.matched = 'ak3GenJets'#slimmedGenJets'
process.patJetPartonMatchCA8PFprunedSubjets.matched = importantgenparticles

process.patJetsCA8PFprunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCA8PFpruned" ),
    subjetSrc=cms.InputTag("patJetsCA8PFprunedSubjets")
      )


# patJetsCA8CHSpruned
addJetCollection(
    process,
    labelName = 'CA8CHSpruned',
    jetSource = cms.InputTag('ca8CHSJetsPruned'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA8CHSpruned.matched=importantgenparticles
process.patJetCorrFactorsCA8CHSpruned.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8CHSpruned.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchCA8CHSpruned.matched = importantgenparticles
process.jetTracksAssociatorAtVertexCA8CHSpruned=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8CHSJetsPruned'), coneSize = 0.8)
process.secondaryVertexTagInfosCA8CHSpruned.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCA8CHSpruned.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCA8CHSpruned= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA8CHSpruned.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCA8CHSpruned.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCA8CHSpruned.jetTagComputer = cms.string('combinedSecondaryVertexCA8CHSpruned')

addJetCollection(
    process,
    labelName = 'CA8CHSprunedSubjets',
    jetSource = cms.InputTag('ca8CHSJetsPruned','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )

process.patJetPartonMatchCA8CHSprunedSubjets.matched=importantgenparticles
process.patJetCorrFactorsCA8CHSprunedSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchCA8CHSprunedSubjets.matched = 'ak3GenJets'#slimmedGenJets'
process.patJetPartonMatchCA8CHSprunedSubjets.matched = importantgenparticles

process.patJetsCA8CHSprunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCA8CHSpruned" ),
    subjetSrc=cms.InputTag("patJetsCA8CHSprunedSubjets")
      )

# patJetsCMSTopTagCHS
addJetCollection(
    process,
    labelName = 'CMSTopTagCHS',
    jetSource = cms.InputTag('cmsTopTagCHS'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
        # btagInfos = [
        # 'CATopTagInfos'
        #  ]   
    )
process.patJetPartonMatchCMSTopTagCHS.matched=importantgenparticles
process.patJetCorrFactorsCMSTopTagCHS.primaryVertices = primaryvertices
process.patJetGenJetMatchCMSTopTagCHS.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHS.matched = importantgenparticles
process.jetTracksAssociatorAtVertexCMSTopTagCHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('cmsTopTagCHS'), coneSize = 0.8)
process.secondaryVertexTagInfosCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCMSTopTagCHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCMSTopTagCHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCMSTopTagCHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCMSTopTagCHS.jetTagComputer = cms.string('combinedSecondaryVertexCMSTopTagCHS')
process.patJetsCMSTopTagCHS.addTagInfos = True
process.patJetsCMSTopTagCHS.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopTagInfos')
    )

addJetCollection(
    process,
    labelName = 'CMSTopTagCHSSubjets',
    jetSource = cms.InputTag('cmsTopTagCHS','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchCMSTopTagCHSSubjets.matched=importantgenparticles
process.patJetCorrFactorsCMSTopTagCHSSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchCMSTopTagCHSSubjets.matched = 'cmsTopTagGENSubJets'#slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHSSubjets.matched = importantgenparticles

process.patJetsCMSTopTagCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCMSTopTagCHS" ),
    subjetSrc=cms.InputTag("patJetsCMSTopTagCHSSubjets")
      )


# patJetsCMSTopTagFJCHS
addJetCollection(
    process,
    labelName = 'CMSTopTagFJCHS',
    jetSource = cms.InputTag('cmsTopTagFJCHS'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCMSTopTagFJCHS.matched=importantgenparticles
process.patJetCorrFactorsCMSTopTagFJCHS.primaryVertices = primaryvertices
process.patJetGenJetMatchCMSTopTagFJCHS.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchCMSTopTagFJCHS.matched = importantgenparticles
process.jetTracksAssociatorAtVertexCMSTopTagFJCHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('cmsTopTagFJCHS'), coneSize = 0.8)
process.secondaryVertexTagInfosCMSTopTagFJCHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCMSTopTagFJCHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCMSTopTagFJCHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCMSTopTagFJCHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCMSTopTagFJCHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCMSTopTagFJCHS.jetTagComputer = cms.string('combinedSecondaryVertexCMSTopTagFJCHS')
process.patJetsCMSTopTagFJCHS.addTagInfos = True
process.patJetsCMSTopTagFJCHS.tagInfoSources = cms.VInputTag(
    cms.InputTag('cmsTopTagFJCHSTagInfos')
    )

addJetCollection(
    process,
    labelName = 'CMSTopTagFJCHSSubjets',
    jetSource = cms.InputTag('cmsTopTagFJCHS','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchCMSTopTagFJCHSSubjets.matched=importantgenparticles
process.patJetCorrFactorsCMSTopTagFJCHSSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchCMSTopTagFJCHSSubjets.matched = 'ak3GenJets'#slimmedGenJets'
process.patJetPartonMatchCMSTopTagFJCHSSubjets.matched = importantgenparticles

process.patJetsCMSTopTagFJCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCMSTopTagFJCHS" ),
    subjetSrc=cms.InputTag("patJetsCMSTopTagFJCHSSubjets")
      )

#patJetsCA15CHS
addJetCollection(
    process,
    labelName = 'CA15CHS',
    jetSource = cms.InputTag('ca15CHSJets'),
#    algo = 'ca15',
#    rParam = 1.5,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchCA15CHS.matched=importantgenparticles
process.patJetCorrFactorsCA15CHS.primaryVertices = primaryvertices
process.patJetGenJetMatchCA15CHS.matched = 'ca15GenJets'#'slimmedGenJets'
process.jetTracksAssociatorAtVertexCA15CHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca15CHSJets'), coneSize = 1.5)
process.secondaryVertexTagInfosCA15CHS.trackSelection.jetDeltaRMax = cms.double(1.5) # default is 0.3
process.secondaryVertexTagInfosCA15CHS.vertexCuts.maxDeltaRToJetAxis = cms.double(1.5) # default is 0.5
process.combinedSecondaryVertexCA15CHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCA15CHS.trackSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexCA15CHS.trackPseudoSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexBJetTagsCA15CHS.jetTagComputer = cms.string('combinedSecondaryVertexCA15CHS')


# patJetsHEPTopTagCHS
addJetCollection(
    process,
    labelName = 'HEPTopTagCHS',
    jetSource = cms.InputTag('hepTopTagCHS'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    )
process.patJetPartonMatchHEPTopTagCHS.matched=importantgenparticles
process.patJetCorrFactorsHEPTopTagCHS.primaryVertices = primaryvertices
process.patJetGenJetMatchHEPTopTagCHS.matched = 'ca8GenJetsPruned'#'slimmedGenJets'
process.patJetPartonMatchHEPTopTagCHS.matched = importantgenparticles
process.jetTracksAssociatorAtVertexHEPTopTagCHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('hepTopTagCHS'), coneSize = 1.5)
process.secondaryVertexTagInfosHEPTopTagCHS.trackSelection.jetDeltaRMax = cms.double(1.5) # default is 0.3
process.secondaryVertexTagInfosHEPTopTagCHS.vertexCuts.maxDeltaRToJetAxis = cms.double(1.5) # default is 0.5
process.combinedSecondaryVertexHEPTopTagCHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexHEPTopTagCHS.trackSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexHEPTopTagCHS.trackPseudoSelection.jetDeltaRMax = cms.double(1.5)
process.combinedSecondaryVertexBJetTagsHEPTopTagCHS.jetTagComputer = cms.string('combinedSecondaryVertexHEPTopTagCHS')
process.patJetsHEPTopTagCHS.addTagInfos = True
process.patJetsHEPTopTagCHS.tagInfoSources = cms.VInputTag(
    cms.InputTag('hepTopTagCHSTagInfos')
    )

addJetCollection(
    process,
    labelName = 'HEPTopTagCHSSubjets',
    jetSource = cms.InputTag('hepTopTagCHS','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag(tracks),
    pvSource = cms.InputTag(vertices),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchHEPTopTagCHSSubjets.matched=importantgenparticles
process.patJetCorrFactorsHEPTopTagCHSSubjets.primaryVertices = primaryvertices
process.patJetGenJetMatchHEPTopTagCHSSubjets.matched = 'ak3GenJets'#slimmedGenJets'
process.patJetPartonMatchHEPTopTagCHSSubjets.matched = importantgenparticles

process.patJetsHEPTopTagCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsHEPTopTagCHS" ),
    subjetSrc=cms.InputTag("patJetsHEPTopTagCHSSubjets")
      )

###############################################
# Modifify all

for jetcoll in (process.patJetsAK15PF,
                process.patJetsAK15PFfiltered,
                process.patJetsAK15PFfilteredSubjets,
                process.patJetsAK15CHS,
                process.patJetsAK15CHSfiltered,
                process.patJetsAK15CHSfilteredSubjets,
                process.patJetsCA8PF,
                process.patJetsCA8PFpruned,
                process.patJetsCA8PFprunedSubjets,
                process.patJetsCA8CHS,
                process.patJetsCA15CHS,
                process.patJetsCA8CHSpruned,
                process.patJetsCA8CHSprunedSubjets,
                process.patJetsCMSTopTagCHS,
                process.patJetsCMSTopTagCHSSubjets,
                process.patJetsCMSTopTagFJCHS,
                process.patJetsCMSTopTagFJCHSSubjets,
                process.patJetsHEPTopTagCHS,
                process.patJetsHEPTopTagCHSSubjets,
                ) :
        jetcoll.embedGenJetMatch = False
        jetcoll.getJetMCFlavour = False
        jetcoll.addGenPartonMatch = False


###############################################

#  DOES THIS WORK YET? combinedSecondaryVertexV2BJetTags

###############################################
# N-subjettiness

# #-------------------------------------
# ## N-subjettiness
# from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

# process.Njettiness = Njettiness.clone(
#     src = cms.InputTag("PFJetsCHS"),
#     cone = cms.double(options.jetRadius)
# )

# process.patJets.userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']



###############################################
# HISTOGRAM MAKER
# process.ana = cms.EDAnalyzer('MiniAnalyzer',
#     jets = cms.InputTag("slimmedJets"),
#     fatjets = cms.InputTag("slimmedJetsAK8"),
#     pfCands = cms.InputTag("packedPFCandidates"),
# 	packedGen = cms.InputTag("packedGenParticles"),   #packedGenParticles: all status == 1 particles (for gen jets)
# 	prunedGen = cms.InputTag("prunedGenParticles")    #prunedGenParticles: the interesting particles (for matching)
# )

# process.TFileService = cms.Service("TFileService",
#       fileName = cms.string("CheckMA.root"),
#       closeFileFast = cms.untracked.bool(True)
#   )

#NtupleWriter
useData = False
writeAllGenParticles=True

process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string("Ntuple.root"), 
                                  runOnMiniAOD = cms.bool(useMiniAOD),
                                  doElectrons = cms.bool(True),
                                  doMuons = cms.bool(True),
                                  doTaus = cms.bool(True),
                                  doJets = cms.bool(True),
                                  doJetsConstituents = cms.bool(True),
                                  doTopJets = cms.bool(True),
                                  doTopJetsConstituents = cms.bool(True),
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
                                  doTagInfos = cms.untracked.bool(False), # when set to true crashes for the 'packed' jet collections
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
                                  #topjet_sources = cms.vstring("patJetsHEPTopTagCHSPacked","patJetsCMSTopTagCHSPacked","patJetsCA8CHSprunedPacked"),
                                  #topjet_constituents_sources = cms.vstring("patJetsHEPTopTagCHS","patJetsCMSTopTagCHS", "patJetsCA8CHSpruned"),
                                  topjet_sources = cms.vstring("patJetsCMSTopTagCHSPacked","patJetsHEPTopTagCHSPacked"),
                                  topjet_constituents_sources = cms.vstring("patJetsCA8CHS","patJetsCA15CHS"),
                                  topjet_ptmin = cms.double(100.0), 
                                  topjet_etamax = cms.double(5.0),
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


###############################################
# PATHS

#process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
#  process.ak3GenJets
#  *process.ak15GenJets
#  *process.ak15GenJetsFiltered
#  *process.ca8GenJetsTrimmed
#  *process.ca8GenJetsPruned
#  *process.ak4PFJets
#  *process.ca8PFJets
#  *process.ak15PFJets
#  *process.ca8CHSJets
#  *process.ca15CHSJets
#  *process.ak15CHSJets
#  *process.ak15PFJetsFiltered
#  *process.ak15CHSJetsFiltered
#  *process.ak15PFJetsFilteredNfilt2
#  *process.ak15CHSJetsFilteredNfilt2
#  *process.ak15PFJetsMassDropFiltered
#  *process.ak15CHSJetsMassDropFiltered
#  *process.ca8PFJetsTrimmed
#  *process.ca8CHSJetsTrimmed
#  *process.ca8PFJetsPruned
#  *process.ca8CHSJetsPruned
#  *process.ca8PFJetsPrunedNsub3
#  *process.ca8CHSJetsPrunedNsub3
  process.cmsTopTagCHS
#  *process.cmsTopTagFJCHS
  *process.hepTopTagCHS
  *process.CATopTagInfos
#  *process.cmsTopTagFJCHSTagInfos
  *process.hepTopTagCHSTagInfos
#  *process.patJetsAK15PF
#  *process.patJetsAK15PFfiltered
#  *process.patJetsAK15PFfilteredSubjets
#  *process.patJetsAK15PFfilteredPacked
#  *process.patJetsAK15CHS
#  *process.patJetsAK15CHSfiltered
#  *process.patJetsAK15CHSfilteredSubjets
#  *process.patJetsAK15CHSfilteredPacked
#  *process.patJetsCA8PF
#  *process.patJetsCA8PFpruned
#  *process.patJetsCA8PFprunedSubjets
#  *process.patJetsCA8PFprunedPacked
  *process.patJetsCA8CHS
  *process.patJetsCA15CHS
  *process.patJetsCA8CHSpruned
  *process.patJetsCA8CHSprunedSubjets
  *process.patJetsCA8CHSprunedPacked
  *process.patJetsCMSTopTagCHS
  *process.patJetsCMSTopTagCHSSubjets
  *process.patJetsCMSTopTagCHSPacked
#  *process.patJetsCMSTopTagFJCHS
#  *process.patJetsCMSTopTagFJCHSSubjets
#  *process.patJetsCMSTopTagFJCHSPacked
  *process.patJetsHEPTopTagCHS
  *process.patJetsHEPTopTagCHSSubjets
  *process.patJetsHEPTopTagCHSPacked
  *process.MyNtuple
  #*process.content
  #*process.ana
  )
#process.end = cms.EndPath(process.out)
