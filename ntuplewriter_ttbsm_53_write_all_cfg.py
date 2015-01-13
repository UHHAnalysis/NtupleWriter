# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

###############################
####### Parameters ############
###############################

# workaround: if executed from edmConfigHash, sys.argv is missing, so set it here manually. 
import sys
if not hasattr(sys, 'argv'): sys.argv = ['abc.py'] # VarParsing expects to find some '*.py' in the argument list, so just provide it ...
 

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

options.register ('tlbsmTag',
                  'tlbsm_53x_v3',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'TLBSM tag use in production')

options.register ('useData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  'Run this on real data')

options.register ('globalTag',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  'Overwrite default globalTag')

options.register ('hltProcess',
                  'HLT',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "HLT process name to use.")

options.register ('writePFCands',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Output PF candidates")


options.register ('writeFat',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Output tracks and PF candidates (and GenParticles for MC)")

options.register ('writeSimpleInputs',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Write four-vector and ID of PF candidates")

options.register ('writeGenParticles',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Output GenParticles collection")

options.register ('writeAllGenParticles',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Output full GenParticles collection")

options.register ('forceCheckClosestZVertex',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Force the check of the closest z vertex")


options.register ('useSusyFilter',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Use the SUSY event filter")


options.register ('useExtraJetColls',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Write extra jet collections for substructure studies")


options.register ('usePythia8',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Use status codes from Pythia8 rather than Pythia6")


options.register ('usePythia6andPythia8',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Use status codes from Pythia8 and Pythia6")


options.register ('runOnFastSim',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Option needed to run on fastsim.")


options.register('doJetTauCrossCleaning',
                 False,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Enable cleaning the jet collections based on taus")


options.register ('useExplicitJTA',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  'Run the explicit Jet-track association')

options.parseArguments()


if not options.useData :
    inputJetCorrLabelAK5PFchs = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    inputJetCorrLabelAK7PFchs = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
else :
    inputJetCorrLabelAK5PFchs = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    inputJetCorrLabelAK7PFchs = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])


process.source = cms.Source("PoolSource",
                            fileNames  =
                            cms.untracked.vstring(__FILE_NAMES__),
                            skipEvents =
                            cms.untracked.uint32(__SKIP_EVENTS__)
                            )

#process.source.eventsToProcess = cms.untracked.VEventRange( ['1:86747'] )

#process.source.skipEvents = cms.untracked.uint32(17268) 

print options

print 'Running AK5 jet corrections: '
print inputJetCorrLabelAK5PFchs

print 'Running AK7 jet corrections: '
print inputJetCorrLabelAK7PFchs

import sys


###############################
####### Global Setup ##########
###############################

if options.useData :
    if options.globalTag is '':
        process.GlobalTag.globaltag = cms.string( 'FT_53_V21_AN4::All' )
    else:
        process.GlobalTag.globaltag = cms.string( options.globalTag )
else :
    if options.globalTag is '':
        process.GlobalTag.globaltag = cms.string( 'START53_V24::All' )
    else:
        process.GlobalTag.globaltag = cms.string( options.globalTag )


from PhysicsTools.PatAlgos.patTemplate_cfg import *


## The beam scraping filter __________________________________________________||
process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

## The Good vertices collection needed by the tracking failure filter ________||
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && position.rho < 2")
)

if not options.runOnFastSim:
    ## The tracking failure filter _______________________________________________||
    process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
    process.load('RecoMET.METFilters.trackingPOGFilters_cfi')

    # Tracking coherent noise filter , see https://twiki.cern.ch/twiki/bin/viewauth/CMS/TrackingPOGFilters
    process.manystripclus53X = cms.EDFilter('ByClusterSummaryMultiplicityPairEventFilter',
                                            multiplicityConfig = cms.PSet(
        firstMultiplicityConfig = cms.PSet(
        clusterSummaryCollection = cms.InputTag("clusterSummaryProducer"),
                subDetEnum = cms.int32(5),
                subDetVariable = cms.string("pHits")
                ),
            secondMultiplicityConfig = cms.PSet(
                clusterSummaryCollection = cms.InputTag("clusterSummaryProducer"),
                subDetEnum = cms.int32(0),
                subDetVariable = cms.string("cHits")
                ),
            ),
        cut = cms.string("( mult2 > 20000+7*mult1)")
     )

    process.toomanystripclus53X = cms.EDFilter('ByClusterSummaryMultiplicityPairEventFilter',
        multiplicityConfig = cms.PSet(
            firstMultiplicityConfig = cms.PSet(
                clusterSummaryCollection = cms.InputTag("clusterSummaryProducer"),
                subDetEnum = cms.int32(5),
                subDetVariable = cms.string("pHits")
                ),
            secondMultiplicityConfig = cms.PSet(
                clusterSummaryCollection = cms.InputTag("clusterSummaryProducer"),
                subDetEnum = cms.int32(0),
                subDetVariable = cms.string("cHits")
                ),
            ),
        cut = cms.string("(mult2>50000) && ( mult2 > 20000+7*mult1)")
        )
        




## Add the latest Tau discriminators _________________________________________||
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# load the b-tag configuration ##
process.load('RecoBTag.Configuration.RecoBTag_cff')
svComputer = cms.InputTag("combinedSecondaryVertex"),

# switch on PAT trigger
#from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
#switchOnTrigger( process, hltProcess=options.hltProcess )


###############################
####### DAF PV's     ##########
###############################

pvSrc = 'offlinePrimaryVertices'

## The good primary vertex filter ____________________________________________||
process.primaryVertexFilter = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake & ndof > 4 & abs(z) <= 24 & position.Rho <= 2"),
    filter = cms.bool(True)
    )


from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( maxZ = cms.double(24.0),
                                     minNdof = cms.double(4.0) # this is >= 4
                                     ),
    src=cms.InputTag(pvSrc)
    )


###############################
########## Gen Setup ##########
###############################


process.load("RecoJets.Configuration.GenJetParticles_cff")
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ca8GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.ca8TopGenJetsNoNu = ca4GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.ca15GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(1.5),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.ca15TopGenJetsNoNu = ca4GenJets.clone( rParam = cms.double(1.5),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.ca12GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(1.2),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.ca12TopGenJetsNoNu = ca4GenJets.clone( rParam = cms.double(1.2),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))

process.ak8GenJetsNoNu = ak5GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))


process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

# add the flavor history
process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")


# prune gen particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
                                            src = cms.InputTag("genParticles"),
                                            select = cms.vstring(
                                                "drop  *"
                                                ,"keep status = 3" #keeps  particles from the hard matrix element
                                                ,"keep (abs(pdgId) >= 11 & abs(pdgId) <= 16) & status = 1" #keeps e/mu and nus with status 1
                                                ,"keep (abs(pdgId)  = 15) & status = 3" #keeps taus
                                                )
                                            )

if options.usePythia8 :
    process.prunedGenParticles.select = cms.vstring(
                                                "drop  *"
                                                ,"keep status = 21" #keeps  particles from the hard matrix element
                                                ,"keep status = 22" #keeps  particles from the hard matrix element
                                                ,"keep status = 23" #keeps  particles from the hard matrix element
                                                ,"keep (abs(pdgId) >= 11 & abs(pdgId) <= 16) & status = 1" #keeps e/mu and nus with status 1
                                                ,"keep (abs(pdgId)  = 15) & (status = 21 || status = 22 || status = 23) " #keeps taus
                                                )
if options.usePythia6andPythia8 :
    process.prunedGenParticles.select = cms.vstring(
                                                "drop  *"
                                                ,"keep status = 3" #keeps  particles from the hard matrix element
                                                ,"keep status = 21" #keeps  particles from the hard matrix element
                                                ,"keep status = 22" #keeps  particles from the hard matrix element
                                                ,"keep status = 23" #keeps  particles from the hard matrix element
                                                ,"keep (abs(pdgId) >= 11 & abs(pdgId) <= 16) & status = 1" #keeps e/mu and nus with status 1
                                                ,"keep (abs(pdgId)  = 15) & (status = 3 || status = 21 || status = 22 || status = 23)" #keeps taus
                                                )                                      


## process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
##                                             src = cms.InputTag("genParticles"),
##                                             select = cms.vstring(
##                                                 "drop  *"
##                                                 ,"keep++ (abs(pdgId) =6) "
##                                                 )
##                                             )

###############################
#### Jet RECO includes ########
###############################

from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.CaloJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *


###############################
########## PF Setup ###########
###############################


# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff. 
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=not options.useData, postfix=postfix,
	  jetCorrections=inputJetCorrLabelAK5PFchs, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True)
if not options.forceCheckClosestZVertex :
    process.pfPileUpPFlow.checkClosestZVertex = False

# change the cone size of electron isolation to 0.3 as default.
process.pfIsolatedElectronsPFlow.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFlow"))
process.pfIsolatedElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFlow")
process.pfIsolatedElectronsPFlow.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPFlow"), cms.InputTag("elPFIsoValueGamma03PFIdPFlow"))

process.pfElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFlow"))
process.pfElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFlow" )
process.pfElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag( "elPFIsoValueNeutral03PFIdPFlow"), cms.InputTag("elPFIsoValueGamma03PFIdPFlow"))

process.patElectronsPFlow.isolationValues = cms.PSet(
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFlow"),
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFlow"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFlow"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFlow"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFlow")
        )

postfixLoose = "PFlowLoose"
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=not options.useData, postfix=postfixLoose,
	  jetCorrections=inputJetCorrLabelAK5PFchs, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True)
if not options.forceCheckClosestZVertex :
    process.pfPileUpPFlowLoose.checkClosestZVertex = False
    
    
# change the cone size of electron isolation to 0.3 as default.
process.pfIsolatedElectronsPFlowLoose.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFlowLoose"))
process.pfIsolatedElectronsPFlowLoose.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFlowLoose")
process.pfIsolatedElectronsPFlowLoose.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPFlowLoose"), cms.InputTag("elPFIsoValueGamma03PFIdPFlowLoose"))

process.pfElectronsPFlowLoose.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFlowLoose"))
process.pfElectronsPFlowLoose.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFlowLoose" )
process.pfElectronsPFlowLoose.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag( "elPFIsoValueNeutral03PFIdPFlowLoose"), cms.InputTag("elPFIsoValueGamma03PFIdPFlowLoose"))

process.patElectronsPFlowLoose.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFlowLoose"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFlowLoose"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFlowLoose"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFlowLoose"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFlowLoose")
    )

# enable/disable tau cleaning
if not options.doJetTauCrossCleaning:
    # if jetCrossCleaning is false, we want to disable
    # the cross cleaning (which is on by default)
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoTau"+postfixLoose).enable = False
else:
    getattr(process,"pfNoTau"+postfix).enable = False
    getattr(process,"pfNoTau"+postfixLoose).enable = False

# Set up "loose" leptons. 

process.pfIsolatedMuonsPFlowLoose.isolationCut = cms.double(999.0) 
process.pfIsolatedElectronsPFlowLoose.isolationCut = cms.double(999.0)
process.patMuonsPFlowLoose.pfMuonSource = "pfMuonsPFlowLoose"
process.patElectronsPFlowLoose.pfElectronSource = "pfElectronsPFlowLoose"

# Keep additional PF information for taus
# embed in AOD externally stored leading PFChargedHadron candidate
process.patTausPFlow.embedLeadPFChargedHadrCand = cms.bool(True)  
# embed in AOD externally stored signal PFChargedHadronCandidates
process.patTausPFlow.embedSignalPFChargedHadrCands = cms.bool(True)  
# embed in AOD externally stored signal PFGammaCandidates
process.patTausPFlow.embedSignalPFGammaCands = cms.bool(True) 
# embed in AOD externally stored isolation PFChargedHadronCandidates
process.patTausPFlow.embedIsolationPFChargedHadrCands = cms.bool(True) 
# embed in AOD externally stored isolation PFGammaCandidates
process.patTausPFlow.embedIsolationPFGammaCands = cms.bool(True)
# embed in AOD externally stored leading PFChargedHadron candidate
process.patTaus.embedLeadPFChargedHadrCand = cms.bool(True)  
# embed in AOD externally stored signal PFChargedHadronCandidates 
process.patTaus.embedSignalPFChargedHadrCands = cms.bool(True)  
# embed in AOD externally stored signal PFGammaCandidates
process.patTaus.embedSignalPFGammaCands = cms.bool(True) 
# embed in AOD externally stored isolation PFChargedHadronCandidates 
process.patTaus.embedIsolationPFChargedHadrCands = cms.bool(True) 
# embed in AOD externally stored isolation PFGammaCandidates
process.patTaus.embedIsolationPFGammaCands = cms.bool(True)

# turn to false when running on data
if options.useData :
    removeMCMatching( process, ['All'] )

###############################
###### Electron ID ############
###############################


process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi') 
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
#Electron ID
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )

process.patElectronsPFlowLoose.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlowLoose.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePFlowLoose.replace( process.patElectronsPFlowLoose, process.eidMVASequence * process.patElectronsPFlowLoose )

#Convesion Rejection
# this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
process.patConversionsPFlow = cms.EDProducer("PATConversionProducer",
                                             electronSource = cms.InputTag("selectedPatElectronsPFlow")      
                                             )
process.patPF2PATSequencePFlow += process.patConversionsPFlow
process.patConversionsPFlowLoose = cms.EDProducer("PATConversionProducer",
                                                  electronSource = cms.InputTag("selectedPatElectronsPFlowLoose")  
                                                  )
process.patPF2PATSequencePFlowLoose += process.patConversionsPFlowLoose


###############################
###### Bare KT 0.6 jets #######store/mc/Summer12/W1JetsToLNu_matchingdown_TuneZ2Star_8TeV-madgraph/AODSIM/START53_V7C_FSIM-v1/10000/0013D85A-7D9E-E211-B3E2-80000048FE80.root
###############################

from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
process.kt6PFJetsForIsolation =  kt6PFJets.clone(
    rParam = 0.6,
    doRhoFastjet = True,
    Rho_EtaMax = cms.double(2.5),
    )

###############################
###### Bare CA 0.8 jets #######
###############################
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca8PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(6.0),
    Ghost_EtaMax = cms.double(8.0)
    )

###############################
###### Bare CA 1.2 jets #######
###############################
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca12PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(1.2),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(6.0),
    Ghost_EtaMax = cms.double(8.0)
    )

###############################
###### Bare CA 1.5 jets #######
###############################
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.ca15PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(1.5),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(100.0)
    )



###############################
###### AK 0.7 jets ############
###############################
process.ak7PFlow = process.pfJetsPFlow.clone(
	rParam = cms.double(0.7)
    )


###############################
###### AK 0.8 jets ############
###############################
process.ak8PFlow = process.pfJetsPFlow.clone(
	rParam = cms.double(0.8)
    )


###############################
###### AK 0.5 jets groomed ####
###############################

from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
process.ak5TrimmedPFlow = ak5PFJetsTrimmed.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True)
    )

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.ak5FilteredPFlow = ak5PFJetsFiltered.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True)
    )

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ak5PrunedPFlow = ak5PFJetsPruned.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True)
    )



###############################
###### AK 0.7 jets groomed ####
###############################

process.ak7TrimmedPFlow = process.ak5TrimmedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.7)
    )

process.ak7FilteredPFlow = process.ak5FilteredPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.7)
	)

process.ak7PrunedPFlow = process.ak5PrunedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.7)
    )


process.ak7TrimmedGenJetsNoNu = ak5GenJets.clone(
	rParam = cms.double(0.7),
	src = cms.InputTag("genParticlesForJetsNoNu"),
	useTrimming = cms.bool(True),
	rFilt = cms.double(0.2),
	trimPtFracMin = cms.double(0.03),
	)

process.ak7FilteredGenJetsNoNu = ak5GenJets.clone(
	rParam = cms.double(0.7),
	src = cms.InputTag("genParticlesForJetsNoNu"),
	useFiltering = cms.bool(True),
	nFilt = cms.int32(3),
	rFilt = cms.double(0.3),
	writeCompound = cms.bool(True),
	jetCollInstanceName=cms.string("SubJets")
	)



process.ak7PrunedGenJetsNoNu = ak5GenJets.clone(
	SubJetParameters,
	rParam = cms.double(0.7),
	src = cms.InputTag("genParticlesForJetsNoNu"),
	usePruning = cms.bool(True),
	writeCompound = cms.bool(True),
	jetCollInstanceName=cms.string("SubJets")
	)



###############################
###### AK 0.8 jets groomed ####
###############################

process.ak8TrimmedPFlow = process.ak5TrimmedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.8)
    )

process.ak8FilteredPFlow = process.ak5FilteredPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.8)
	)

process.ak8PrunedPFlow = process.ak5PrunedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.8)
    )

###############################
###### CA8 Pruning Setup ######
###############################


# Pruned PF Jets
process.caPrunedPFlow = process.ak5PrunedPFlow.clone(
	jetAlgorithm = cms.string("CambridgeAachen"),
	rParam       = cms.double(0.8)
)


process.caPrunedGen = process.ca8GenJetsNoNu.clone(
	SubJetParameters,
	usePruning = cms.bool(True),
	useExplicitGhosts = cms.bool(True),
	writeCompound = cms.bool(True),
	jetCollInstanceName=cms.string("SubJets")
)

###############################
###### CA8 Filtered Setup #####
###############################


# Filtered PF Jets
process.caFilteredPFlow = ak5PFJetsFiltered.clone(
	src = cms.InputTag('pfNoElectron'+postfix),
	jetAlgorithm = cms.string("CambridgeAachen"),
	rParam       = cms.double(1.5),
	writeCompound = cms.bool(True),
	doAreaFastjet = cms.bool(True),
	jetPtMin = cms.double(100.0)
)

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsMassDropFiltered
process.caMassDropFilteredPFlow = ak5PFJetsMassDropFiltered.clone(
	src = cms.InputTag('pfNoElectron'+postfix),
	jetAlgorithm = cms.string("CambridgeAachen"),
	rParam       = cms.double(1.5),
	writeCompound = cms.bool(True),
	doAreaFastjet = cms.bool(True),
	jetPtMin = cms.double(100.0)
)

process.caFilteredGenJetsNoNu = process.ca8GenJetsNoNu.clone(
	nFilt = cms.int32(3),
	rFilt = cms.double(0.3),
	useFiltering = cms.bool(True),
	useExplicitGhosts = cms.bool(True),
	writeCompound = cms.bool(True),
	rParam       = cms.double(1.5),
	jetCollInstanceName=cms.string("SubJets"),
	jetPtMin = cms.double(150.0)
)

process.caMassDropFilteredGenJetsNoNu = process.caFilteredGenJetsNoNu.clone(
        src = cms.InputTag('genParticlesForJetsNoNu'),
	useMassDropTagger = cms.bool(True),
	muCut = cms.double(0.667),
	yCut = cms.double(0.08)
)

###############################
#### CATopTag Setup ###########
###############################

# CATopJet PF Jets
# with adjacency 
process.caTopTagPFlow = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag('pfNoElectron'+postfix),
                           doAreaFastjet = cms.bool(True),
                           doRhoFastjet = cms.bool(False),
                           Rho_EtaMax = cms.double(6.0),
                           Ghost_EtaMax = cms.double(7.0),
			   jetPtMin = cms.double(100.0)
                           ),
    AnomalousCellParameters,
    CATopJetParameters,
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )

process.ca12TopTagPFlow = process.caTopTagPFlow.clone(
	rParam = cms.double(1.2)
)


process.caHEPTopTagPFlow = process.caTopTagPFlow.clone(
	rParam = cms.double(1.5),
	tagAlgo = cms.int32(2)
)


process.CATopTagInfosPFlow = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("caTopTagPFlow"),
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


process.CA12TopTagInfosPFlow = process.CATopTagInfosPFlow.clone(
	src = cms.InputTag("ca12TopTagPFlow")
)


process.CATopTagInfosHEPTopTagPFlow = process.CATopTagInfosPFlow.clone(
	src = cms.InputTag("caHEPTopTagPFlow")
)

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

process.ca12TopTagGen = process.caTopTagGen.clone(
	rParam = cms.double(1.2)
)

process.caHEPTopTagGen = process.caTopTagGen.clone(
	rParam = cms.double(1.5),
        jetPtMin = cms.double(150.0)
)

process.CATopTagInfosGen = cms.EDProducer("CATopJetTagger",
                                          src = cms.InputTag("caTopTagGen"),
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


process.CA12TopTagInfosGen = process.CATopTagInfosGen.clone(
    src = cms.InputTag("ca12TopTagGen")
    )


# CATopJet PF Jets

for ipostfix in [postfix] :
    for module in (
        getattr(process,"ca8PFJets" + ipostfix),
        getattr(process,"ca12PFJets" + ipostfix),
        getattr(process,"CATopTagInfos" + ipostfix),
        getattr(process,"CATopTagInfosHEPTopTag" + ipostfix),
        getattr(process,"CA12TopTagInfos" + ipostfix),
        getattr(process,"caTopTag" + ipostfix),
        getattr(process,"ca12TopTag" + ipostfix),
        getattr(process,"caHEPTopTag" + ipostfix),
        getattr(process,"caPruned" + ipostfix),
        getattr(process,"ca15PFJets" + ipostfix),
        getattr(process,"caFiltered" + ipostfix),
        getattr(process,"caMassDropFiltered" + ipostfix)
        ) :
        getattr(process,"patPF2PATSequence"+ipostfix).replace( getattr(process,"pfNoElectron"+ipostfix), getattr(process,"pfNoElectron"+ipostfix)*module )


    if options.useExtraJetColls : 
	    for module in (
		getattr(process,"ak5Trimmed" + ipostfix),
		getattr(process,"ak5Filtered" + ipostfix),
		getattr(process,"ak5Pruned" + ipostfix),
		getattr(process,"ak7Trimmed" + ipostfix),
		getattr(process,"ak7Filtered" + ipostfix),
		getattr(process,"ak7Pruned" + ipostfix),
		getattr(process,"ak7" + ipostfix),
		getattr(process,"ak8Trimmed" + ipostfix),
		getattr(process,"ak8Filtered" + ipostfix),
		getattr(process,"ak8Pruned" + ipostfix),
		getattr(process,"ak8" + ipostfix)
		) :
		    getattr(process,"patPF2PATSequence"+ipostfix).replace( getattr(process,"pfNoElectron"+ipostfix), getattr(process,"pfNoElectron"+ipostfix)*module )



# Use the good primary vertices everywhere. 
for imod in [process.patMuonsPFlow,
             process.patMuonsPFlowLoose,
             process.patElectronsPFlow,
             process.patElectronsPFlowLoose,
             process.patMuons,
             process.patElectrons] :
    imod.pvSrc = "goodOfflinePrimaryVertices"
    imod.embedTrack = True
    

addJetCollection(process, 
                 cms.InputTag('ca8PFJetsPFlow'),
                 'CA8', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
                 doJetID = False
                 )

addJetCollection(process, 
                 cms.InputTag('ca12PFJetsPFlow'),
                 'CA12', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca12GenJetsNoNu"),
                 doJetID = False
                 )


addJetCollection(process, 
                 cms.InputTag('caPrunedPFlow'),
                 'CA8Pruned', 'PF',
                 doJTA=False,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
                 doJetID = False
                 )


addJetCollection(process,
                 cms.InputTag('caPrunedPFlow','SubJets'),
                 'CA8PrunedSubjets', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK5PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection=cms.InputTag('caPrunedGen','SubJets'),
                 doJetID=False
                 )

addJetCollection(process, 
                 cms.InputTag('caTopTagPFlow'),
                 'CATopTag', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
                 doJetID = False
                 )

addJetCollection(process, 
                 cms.InputTag('ca12TopTagPFlow'),
                 'CA12TopTag', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = cms.InputTag("ca12GenJetsNoNu"),
                 doJetID = False
                 )

addJetCollection(process,
                 cms.InputTag('caTopTagPFlow', 'caTopSubJets'),
                 'CATopTagSubjets', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK5PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )

addJetCollection(process,
                 cms.InputTag('ca12TopTagPFlow', 'caTopSubJets'),
                 'CA12TopTagSubjets', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK5PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )


addJetCollection(process, 
                 cms.InputTag('caHEPTopTagPFlow'),
                 'CAHEPTopTag', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )

addJetCollection(process, 
                 cms.InputTag('caHEPTopTagPFlow', 'caTopSubJets'),
                 'CAHEPTopTagSubjets', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK5PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )

addJetCollection(process, 
                 cms.InputTag('ca15PFJetsPFlow'),
                 'CA15', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=True,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )

addJetCollection(process, 
                 cms.InputTag('caFilteredPFlow'),
                 'CA15Filtered', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=True,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )


addJetCollection(process, 
                 cms.InputTag('caMassDropFilteredPFlow'),
                 'CA15MassDropFiltered', 'PF',
                 doJTA=True,
                 doBTagging=True,
                 jetCorrLabel=inputJetCorrLabelAK7PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                 )


addJetCollection(process, 
                 cms.InputTag('caMassDropFilteredPFlow', 'SubJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
                 'CA15MassDropFilteredSubjets', 'PF',
                 doJTA=True,            # Run Jet-Track association & JetCharge
                 doBTagging=True,       # Run b-tagging
                 jetCorrLabel=inputJetCorrLabelAK5PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                     )

addJetCollection(process, 
                 cms.InputTag('caFilteredPFlow', 'SubJets'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
                 'CA15FilteredSubjets', 'PF',
                 doJTA=True,            # Run Jet-Track association & JetCharge
                 doBTagging=True,       # Run b-tagging
                 jetCorrLabel=inputJetCorrLabelAK5PFchs,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 genJetCollection = None,
                 doJetID = False
                     )


##############################################################
### For subjet b tagging with explicit jet-track association
###
### This requires the following additional packages
###
### addpkg RecoJets/JetAssociationAlgorithms V03-01-01-00
### addpkg RecoJets/JetAssociationProducers  V03-02-01

### Substitute the standard jet-track association with the explicit jet-track association
### (this will keep the original module names unchanged so might be a bit misleading at first glance)

if options.useExplicitJTA : 
    for xtrplabel in ['CA8PrunedSubjets', 'CATopTagSubjets', 'CA12TopTagSubjets','CAHEPTopTagSubjets' , 'CA15MassDropFilteredSubjets' , 'CA15FilteredSubjets'] :
        if hasattr( process, 'jetTracksAssociatorAtVertex' + xtrplabel + 'PF' ):
            from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorExplicit
            m = 'jetTracksAssociatorAtVertex' + xtrplabel + 'PF'
            print 'Switching ' + m + ' to explicit jet-track association'
            setattr( process, m, ak5JetTracksAssociatorExplicit.clone(jets = getattr(getattr(process,m),'jets')) )

###
##############################################################



if options.useExtraJetColls: 


	addJetCollection(process, 
			 cms.InputTag('ak5PrunedPFlow'),
			 'AK5Pruned', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK5PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
			 doJetID = False
			 )


	addJetCollection(process, 
			 cms.InputTag('ak5FilteredPFlow'),
			 'AK5Filtered', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK5PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
			 doJetID = False
			 )

	addJetCollection(process, 
			 cms.InputTag('ak5TrimmedPFlow'),
			 'AK5Trimmed', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK5PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
			 doJetID = False
			 )


	addJetCollection(process, 
			 cms.InputTag('ak7PFlow'),
			 'AK7', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )

	addJetCollection(process, 
			 cms.InputTag('ak7PrunedPFlow'),
			 'AK7Pruned', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )


	addJetCollection(process, 
			 cms.InputTag('ak7FilteredPFlow'),
			 'AK7Filtered', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )

	addJetCollection(process, 
			 cms.InputTag('ak7TrimmedPFlow'),
			 'AK7Trimmed', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )





	addJetCollection(process, 
			 cms.InputTag('ak8PFlow'),
			 'AK8', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )

	addJetCollection(process, 
			 cms.InputTag('ak8PrunedPFlow'),
			 'AK8Pruned', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )


	addJetCollection(process, 
			 cms.InputTag('ak8FilteredPFlow'),
			 'AK8Filtered', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )

	addJetCollection(process, 
			 cms.InputTag('ak8TrimmedPFlow'),
			 'AK8Trimmed', 'PF',
			 doJTA=False,
			 doBTagging=False,
			 jetCorrLabel=inputJetCorrLabelAK7PFchs,
			 doType1MET=True,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )

switchJetCollection(process,cms.InputTag('ak5PFJets'),
		    doJTA        = False,
		    doBTagging   = False,
		    jetCorrLabel = inputJetCorrLabelAK5PFchs,
		    doType1MET   = True,
		    genJetCollection=cms.InputTag("ak5GenJetsNoNu"),
		    doJetID      = False
		    )

for icorr in [process.patJetCorrFactors,
	      process.patJetCorrFactorsCATopTagPF,
	      process.patJetCorrFactorsCA12TopTagPF,
	      process.patJetCorrFactorsCAHEPTopTagPF,
              process.patJetCorrFactorsCA8PrunedPF,
              process.patJetCorrFactorsCA12PF,
              process.patJetCorrFactorsCA15PF,
              process.patJetCorrFactorsCA15FilteredPF,
              process.patJetCorrFactorsCA15MassDropFilteredPF,
	      process.patJetCorrFactorsCATopTagSubjetsPF,
	      process.patJetCorrFactorsCA12TopTagSubjetsPF,
	      process.patJetCorrFactorsCAHEPTopTagSubjetsPF,
              process.patJetCorrFactorsCA8PrunedSubjetsPF,
              process.patJetCorrFactorsCA15MassDropFilteredSubjetsPF,
              process.patJetCorrFactorsCA15FilteredSubjetsPF,
              process.patJetCorrFactorsCA8PF ] :
    icorr.rho = cms.InputTag("kt6PFJets", "rho")


if options.useExtraJetColls: 
	for icorr in [process.patJetCorrFactorsAK5PrunedPF,
		      process.patJetCorrFactorsAK5FilteredPF,
		      process.patJetCorrFactorsAK5TrimmedPF,
		      process.patJetCorrFactorsAK7PF,
		      process.patJetCorrFactorsAK7PrunedPF,
		      process.patJetCorrFactorsAK7FilteredPF,
		      process.patJetCorrFactorsAK7TrimmedPF,
		      process.patJetCorrFactorsAK8PF,
		      process.patJetCorrFactorsAK8PrunedPF,
		      process.patJetCorrFactorsAK8FilteredPF,
		      process.patJetCorrFactorsAK8TrimmedPF] :
	    icorr.rho = cms.InputTag("kt6PFJets", "rho")



###############################
### TagInfo and Matching Setup#
###############################

# Do some configuration of the jet substructure things
for jetcoll in (process.patJetsPFlow,
		process.patJets,
                process.patJetsCA8PF,
                process.patJetsCA12PF,
                process.patJetsCA8PrunedPF,
                process.patJetsCATopTagPF,
                process.patJetsCA12TopTagPF,
                process.patJetsCAHEPTopTagPF,
                process.patJetsCA15PF,
                process.patJetsCA15FilteredPF,
                process.patJetsCA15MassDropFilteredPF,
                process.patJetsCA8PrunedSubjetsPF,
                process.patJetsCATopTagSubjetsPF,
                process.patJetsCA12TopTagSubjetsPF,
                process.patJetsCAHEPTopTagSubjetsPF,
                process.patJetsCA15MassDropFilteredSubjetsPF,
                process.patJetsCA15FilteredSubjetsPF
                ) :
    if options.useData == False :
        jetcoll.embedGenJetMatch = False
        jetcoll.getJetMCFlavour = True
        jetcoll.addGenPartonMatch = True
    # Add the calo towers and PFCandidates.
    # I'm being a little tricksy here, because I only
    # actually keep the products if the "writeFat" switch
    # is on. However, this allows for overlap checking
    # with the Refs so satisfies most use cases without
    # having to add to the object size
    jetcoll.addBTagInfo = False
    jetcoll.embedCaloTowers = True
    if not options.writeFat and not options.writePFCands : 
        jetcoll.embedPFCandidates = True

        

# Add CATopTag and b-tag info... piggy-backing on b-tag functionality
process.patJetsPFlow.addBTagInfo = True
process.patJetsCATopTagPF.addBTagInfo = True
process.patJetsCA12TopTagPF.addBTagInfo = True
process.patJetsCAHEPTopTagPF.addBTagInfo = True
process.patJetsCA8PrunedPF.addBTagInfo = True
process.patJetsCA8PrunedSubjetsPF.addBTagInfo = True
process.patJetsCA15MassDropFilteredSubjetsPF.addBTagInfo = True
process.patJetsCA15FilteredSubjetsPF.addBTagInfo = True
process.patJetsCA15FilteredPF.addBTagInfo = True
process.patJetsCA15PF.addBTagInfo = True
process.patJetsCATopTagSubjetsPF.addBTagInfo = True
process.patJetsCA12TopTagSubjetsPF.addBTagInfo = True
process.patJetsCAHEPTopTagSubjetsPF.addBTagInfo = True


process.patJetsCA8PrunedSubjetsPF.embedPFCandidates = False

#Add b-tagging taginfos on subjets
process.patJetsCAHEPTopTagSubjetsPF.addTagInfos = True
process.patJetsCATopTagSubjetsPF.addTagInfos = True
process.patJetsCA12TopTagSubjetsPF.addTagInfos = True
process.patJetsCA15FilteredSubjetsPF.addTagInfos = True

# Do some configuration of the jet substructure things
if options.useExtraJetColls: 
	for jetcoll in (process.patJetsAK5TrimmedPF,
			process.patJetsAK5PrunedPF,
			process.patJetsAK5FilteredPF,
			process.patJetsAK7PF,
			process.patJetsAK7TrimmedPF,
			process.patJetsAK7PrunedPF,
			process.patJetsAK7FilteredPF,
			process.patJetsAK8PF,
			process.patJetsAK8TrimmedPF,
			process.patJetsAK8PrunedPF,
			process.patJetsAK8FilteredPF
			) :
	    if options.useData == False :
		jetcoll.embedGenJetMatch = False
		jetcoll.getJetMCFlavour = True
		jetcoll.addGenPartonMatch = True
	    # Add the calo towers and PFCandidates.
	    # I'm being a little tricksy here, because I only
	    # actually keep the products if the "writeFat" switch
	    # is on. However, this allows for overlap checking
	    # with the Refs so satisfies most use cases without
	    # having to add to the object size
	    jetcoll.addBTagInfo = False
	    jetcoll.embedCaloTowers = True
            if not options.writeFat and not options.writePFCands : 
                jetcoll.embedPFCandidates = True




#################################################
#### Fix the PV collections for the future ######
#################################################
for module in [process.patJetCorrFactors,
               process.patJetCorrFactorsPFlow,
               process.patJetCorrFactorsCATopTagPF,
               process.patJetCorrFactorsCA12TopTagPF,
               process.patJetCorrFactorsCAHEPTopTagPF,
               process.patJetCorrFactorsCA8PrunedPF,
               process.patJetCorrFactorsCA15PF,
               process.patJetCorrFactorsCA15FilteredPF,
               process.patJetCorrFactorsCA15MassDropFilteredPF,
               process.patJetCorrFactorsCATopTagSubjetsPF,
               process.patJetCorrFactorsCA12TopTagSubjetsPF,
               process.patJetCorrFactorsCAHEPTopTagSubjetsPF,
               process.patJetCorrFactorsCA8PrunedSubjetsPF,
               process.patJetCorrFactorsCA15MassDropFilteredSubjetsPF,
               process.patJetCorrFactorsCA15FilteredSubjetsPF,
               process.patJetCorrFactorsCA8PF,
               process.patJetCorrFactorsCA12PF
               ]:
    module.primaryVertices = "goodOfflinePrimaryVertices"

    
if options.useExtraJetColls: 
	for module in [process.patJetCorrFactorsAK5TrimmedPF,
		       process.patJetCorrFactorsAK5PrunedPF,
		       process.patJetCorrFactorsAK5FilteredPF,
		       process.patJetCorrFactorsAK7PF,
		       process.patJetCorrFactorsAK7TrimmedPF,
		       process.patJetCorrFactorsAK7PrunedPF,
		       process.patJetCorrFactorsAK7FilteredPF,
		       process.patJetCorrFactorsAK8PF,
		       process.patJetCorrFactorsAK8TrimmedPF,
		       process.patJetCorrFactorsAK8PrunedPF,
		       process.patJetCorrFactorsAK8FilteredPF
		       ]:
	    module.primaryVertices = "goodOfflinePrimaryVertices"


###############################
#### Selections Setup #########
###############################

# AK5 Jets
process.selectedPatJetsPFlow.cut = cms.string("pt > 5")
process.selectedPatJetsPFlowLoose.cut = cms.string("pt > 20")
process.patJetsPFlow.addTagInfos = True
process.patJetsPFlow.tagInfoSources = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODPFlow")
    )
process.patJetsPFlow.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
                                                      "tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0")
process.patJetsPFlow.userData.userFunctionLabels = cms.vstring('secvtxMass')

# CA8 jets
process.selectedPatJetsCA8PF.cut = cms.string("pt > 30 & abs(rapidity) < 2.5")
#process.selectedPatJetsCA8PF.cut = cms.string("pt > 20")

# CA8 Pruned jets
process.selectedPatJetsCA8PrunedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
process.patJetsCA8PrunedSubjetsPF.addTagInfos = False
#process.selectedPatJetsCA8PrunedSubjetsPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")

                                                    
# CA8 TopJets
process.selectedPatJetsCATopTagPF.cut = cms.string("pt > 30 & abs(rapidity) < 2.5")
#process.selectedPatJetsCATopTagSubjetsPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
process.patJetsCATopTagPF.addTagInfos = True
process.patJetsCATopTagPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopTagInfosPFlow')
    )

# CA1.5 HEPTopTagTopJets
process.selectedPatJetsCAHEPTopTagPF.cut = cms.string("pt > 150 & abs(rapidity) < 2.5")
#process.selectedPatJetsCAHEPTopTagSubjetsPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
process.patJetsCAHEPTopTagPF.addTagInfos = True
process.patJetsCAHEPTopTagPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopTagInfosHEPTopTagPFlow')
    )

# CA15 jets
process.selectedPatJetsCA15PF.cut = cms.string("pt > 120 & abs(rapidity) < 2.5")

# CA15 Filtered jets
process.selectedPatJetsCA15FilteredPF.cut = cms.string("pt > 150 & abs(rapidity) < 2.5")
process.selectedPatJetsCA15MassDropFilteredPF.cut = cms.string("pt > 150 & abs(rapidity) < 2.5")
#process.selectedPatJetsCA15MassDropFilteredSubjetsPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")

# CA12 jets
process.selectedPatJetsCA12PF.cut = cms.string("pt > 30 & abs(rapidity) < 2.5")

# CA12 TopJets
process.selectedPatJetsCA12TopTagPF.cut = cms.string("pt > 30 & abs(rapidity) < 2.5")
#process.selectedPatJetsCATopTagSubjetsPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
process.patJetsCA12TopTagPF.addTagInfos = True
process.patJetsCA12TopTagPF.tagInfoSources = cms.VInputTag(
    cms.InputTag('CA12TopTagInfosPFlow')
    )
if options.useExtraJetColls: 

	# AK5 groomed jets
	process.selectedPatJetsAK5PrunedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK5TrimmedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK5FilteredPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")


	# AK7 groomed jets
	process.selectedPatJetsAK7PF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK7PrunedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK7TrimmedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK7FilteredPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")


	# AK8 groomed jets
	process.selectedPatJetsAK8PF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK8PrunedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK8TrimmedPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	process.selectedPatJetsAK8FilteredPF.cut = cms.string("pt > 20 & abs(rapidity) < 2.5")
	


# electrons
process.selectedPatElectrons.cut = cms.string('pt > 5.0 & abs(eta) < 2.5')
process.patElectrons.embedTrack = cms.bool(True)
process.selectedPatElectronsPFlow.cut = cms.string('pt > 5.0 & abs(eta) < 2.5')
process.patElectronsPFlow.embedTrack = cms.bool(True)
process.selectedPatElectronsPFlowLoose.cut = cms.string('pt > 5.0 & abs(eta) < 2.5')
process.patElectronsPFlowLoose.embedTrack = cms.bool(True)
# muons
process.selectedPatMuons.cut = cms.string('pt > 5.0 & abs(eta) < 2.5')
process.patMuons.embedTrack = cms.bool(True)
process.selectedPatMuonsPFlow.cut = cms.string("pt > 5.0 & abs(eta) < 2.5")
process.patMuonsPFlow.embedTrack = cms.bool(True)
process.selectedPatMuonsPFlowLoose.cut = cms.string("pt > 5.0 & abs(eta) < 2.5")
process.patMuonsPFlowLoose.embedTrack = cms.bool(True)
# taus
process.selectedPatTausPFlow.cut = cms.string("pt > 10.0 & abs(eta) < 3")
process.selectedPatTaus.cut = cms.string("pt > 10.0 & abs(eta) < 3")
process.patTausPFlow.isoDeposits = cms.PSet()
process.patTaus.isoDeposits = cms.PSet()
# photons
process.patPhotonsPFlow.isoDeposits = cms.PSet()
process.patPhotons.isoDeposits = cms.PSet()

# Apply jet ID to all of the jets upstream. We aren't going to screw around
# with this, most likely. So, we don't really to waste time with it
# at the analysis level. 
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
process.goodPatJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJetsPFlow")
                                        )
process.goodPatJetsCA8PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJetsCA8PF")
                                        )
process.goodPatJetsCA12PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                        filterParams = pfJetIDSelector.clone(),
                                        src = cms.InputTag("selectedPatJetsCA12PF")
                                        )
process.goodPatJetsCA8PrunedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                              filterParams = pfJetIDSelector.clone(),
                                              src = cms.InputTag("selectedPatJetsCA8PrunedPF")
                                              )

process.goodPatJetsCATopTagPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                             filterParams = pfJetIDSelector.clone(),
                                             src = cms.InputTag("selectedPatJetsCATopTagPF")
                                             )

process.goodPatJetsCA12TopTagPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                             filterParams = pfJetIDSelector.clone(),
                                             src = cms.InputTag("selectedPatJetsCA12TopTagPF")
                                             )                                             


process.goodPatJetsCAHEPTopTagPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                             filterParams = pfJetIDSelector.clone(),
                                             src = cms.InputTag("selectedPatJetsCAHEPTopTagPF")
                                             )

process.goodPatJetsCA15PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                         filterParams = pfJetIDSelector.clone(),
                                         src = cms.InputTag("selectedPatJetsCA15PF")
                                         )

process.goodPatJetsCA15FilteredPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                                 filterParams = pfJetIDSelector.clone(),
                                                 src = cms.InputTag("selectedPatJetsCA15FilteredPF")
                                                 )

process.goodPatJetsCA15MassDropFilteredPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                                         filterParams = pfJetIDSelector.clone(),
                                                         src = cms.InputTag("selectedPatJetsCA15MassDropFilteredPF")
                                                         )

if options.useExtraJetColls:


	process.goodPatJetsAK5PrunedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK5PrunedPF")
						      )
	process.goodPatJetsAK5FilteredPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK5FilteredPF")
						      )
	process.goodPatJetsAK5TrimmedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK5TrimmedPF")
						      )

	process.goodPatJetsAK7PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK7PF")
						      )
	process.goodPatJetsAK7PrunedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK7PrunedPF")
						      )
	process.goodPatJetsAK7FilteredPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK7FilteredPF")
						      )
	process.goodPatJetsAK7TrimmedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK7TrimmedPF")
						      )



	process.goodPatJetsAK8PF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK8PF")
						      )
	process.goodPatJetsAK8PrunedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK8PrunedPF")
						      )
	process.goodPatJetsAK8FilteredPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK8FilteredPF")
						      )
	process.goodPatJetsAK8TrimmedPF = cms.EDFilter("PFJetIDSelectionFunctorFilter",
						      filterParams = pfJetIDSelector.clone(),
						      src = cms.InputTag("selectedPatJetsAK8TrimmedPF")
						      )



process.goodPatJetsCA8PrunedPFPacked = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("goodPatJetsCA8PrunedPF"),
                                                      subjetSrc=cms.InputTag("selectedPatJetsCA8PrunedSubjetsPF")
    )

process.goodPatJetsCATopTagPFPacked = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("goodPatJetsCATopTagPF"),
                                                      subjetSrc=cms.InputTag("selectedPatJetsCATopTagSubjetsPF")
    )

process.goodPatJetsCA12TopTagPFPacked = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("goodPatJetsCA12TopTagPF"),
                                                      subjetSrc=cms.InputTag("selectedPatJetsCA12TopTagSubjetsPF")
    )

process.goodPatJetsCAHEPTopTagPFPacked = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("goodPatJetsCAHEPTopTagPF"),
                                                      subjetSrc=cms.InputTag("selectedPatJetsCAHEPTopTagSubjetsPF")
    )

process.goodPatJetsCA15MassDropFilteredPFPacked = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("goodPatJetsCA15MassDropFilteredPF"),
                                                      subjetSrc=cms.InputTag("selectedPatJetsCA15MassDropFilteredSubjetsPF")
    )

process.goodPatJetsCA15FilteredPFPacked = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("goodPatJetsCA15FilteredPF"),
                                                      subjetSrc=cms.InputTag("selectedPatJetsCA15FilteredSubjetsPF")
    )

if options.writeSimpleInputs :
	process.pfInputs = cms.EDProducer(
	    "CandViewNtpProducer", 
	    src = cms.InputTag('selectedPatJetsCA8PF', 'pfCandidates'),
	    lazyParser = cms.untracked.bool(True),
	    eventInfo = cms.untracked.bool(False),
	    variables = cms.VPSet(
		cms.PSet(
		    tag = cms.untracked.string("px"),
		    quantity = cms.untracked.string("px")
		    ),
		cms.PSet(
		    tag = cms.untracked.string("py"),
		    quantity = cms.untracked.string("py")
		    ),
		cms.PSet(
		    tag = cms.untracked.string("pz"),
		    quantity = cms.untracked.string("pz")
		    ),
		cms.PSet(
		    tag = cms.untracked.string("energy"),
		    quantity = cms.untracked.string("energy")
		    ),
		cms.PSet(
		    tag = cms.untracked.string("pdgId"),
		    quantity = cms.untracked.string("pdgId")
		    )
		)
	)


if options.useExtraJetColls:
	process.ak5Lite = cms.EDProducer(
	    "CandViewNtpProducer", 
	    src = cms.InputTag('goodPatJetsPFlow'),
	    lazyParser = cms.untracked.bool(True),
	    eventInfo = cms.untracked.bool(False),
	    variables = cms.VPSet(
			cms.PSet(
				tag = cms.untracked.string("px"),
				quantity = cms.untracked.string("px")
				),
			cms.PSet(
				tag = cms.untracked.string("py"),
				quantity = cms.untracked.string("py")
				),
			cms.PSet(
				tag = cms.untracked.string("pz"),
				quantity = cms.untracked.string("pz")
				),
			cms.PSet(
				tag = cms.untracked.string("energy"),
				quantity = cms.untracked.string("energy")
				),
			cms.PSet(
				tag = cms.untracked.string("jetArea"),
				quantity = cms.untracked.string("jetArea")
				),
			cms.PSet(
				tag = cms.untracked.string("jecFactor"),
				quantity = cms.untracked.string("jecFactor(0)")
				)
				)
	)


	process.ak5TrimmedLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK5TrimmedPF')
		)

	process.ak5PrunedLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK5PrunedPF')
		)

	process.ak5FilteredLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK5FilteredPF')
		)

	process.ak7Lite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK7PF')
		)

	process.ak7TrimmedLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK7TrimmedPF')
		)

	process.ak7PrunedLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK7PrunedPF')
		)

	process.ak7FilteredLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK7FilteredPF')
		)




	process.ak7TrimmedGenLite = cms.EDProducer(
	    "CandViewNtpProducer", 
	    src = cms.InputTag('ak7TrimmedGenJetsNoNu'),
	    lazyParser = cms.untracked.bool(True),
	    eventInfo = cms.untracked.bool(False),
	    variables = cms.VPSet(
			cms.PSet(
				tag = cms.untracked.string("px"),
				quantity = cms.untracked.string("px")
				),
			cms.PSet(
				tag = cms.untracked.string("py"),
				quantity = cms.untracked.string("py")
				),
			cms.PSet(
				tag = cms.untracked.string("pz"),
				quantity = cms.untracked.string("pz")
				),
			cms.PSet(
				tag = cms.untracked.string("energy"),
				quantity = cms.untracked.string("energy")
				)
				)
	)


	process.ak7PrunedGenLite = process.ak7TrimmedGenLite.clone(
		src = cms.InputTag('ak7PrunedGenJetsNoNu')
		)

	process.ak7FilteredGenLite = process.ak7TrimmedGenLite.clone(
		src = cms.InputTag('ak7FilteredGenJetsNoNu')
		)

        process.ca8PrunedGenLite = process.ak7TrimmedGenLite.clone(
                src = cms.InputTag('caPrunedGen')
                )

        process.ca12FilteredGenLite = process.ak7TrimmedGenLite.clone(
                src = cms.InputTag('caFilteredGenJetsNoNu')
                )

        process.ca12MassDropFilteredGenLite = process.ak7TrimmedGenLite.clone(
                src = cms.InputTag('caMassDropFilteredGenJetsNoNu')
                )



	process.ak8Lite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK8PF')
		)

	process.ak8TrimmedLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK8TrimmedPF')
		)

	process.ak8PrunedLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK8PrunedPF')
		)

	process.ak8FilteredLite = process.ak5Lite.clone(
		src = cms.InputTag('goodPatJetsAK8FilteredPF')
		)


## IVF and BCandidate producer for Vbb cross check analysis
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')

## don't cluster e, mu, tau into gen jets
process.genParticlesForJetsNoNu.excludeFromResonancePids = cms.vuint32(12, 14, 16)

# let it run
if options.runOnFastSim:
    process.filterSeq = cms.Sequence(process.primaryVertexFilter *process.goodVertices)

else:
 process.filtersSeq = cms.Sequence(
   process.primaryVertexFilter *
   process.noscraping *
   process.HBHENoiseFilter *
   process.CSCTightHaloFilter *
   process.hcalLaserEventFilter *
   process.EcalDeadCellTriggerPrimitiveFilter *
   process.goodVertices * process.trackingFailureFilter *
   ~process.manystripclus53X *
   ~process.toomanystripclus53X *
   ~process.logErrorTooManyClusters *
   ~process.logErrorTooManyTripletsPairs *
   ~process.logErrorTooManySeeds *
   process.eeBadScFilter
 )


# remove not needed collections from pat sequence
process.patDefaultSequence.remove( process.selectedPatElectrons )
process.patDefaultSequence.remove( process.selectedPatMuons )
process.patDefaultSequence.remove( process.selectedPatTaus )
process.patDefaultSequence.remove( process.selectedPatPhotons )
process.patDefaultSequence.remove( process.cleanPatElectrons )
process.patDefaultSequence.remove( process.cleanPatMuons )
process.patDefaultSequence.remove( process.cleanPatTaus )
process.patDefaultSequence.remove( process.cleanPatPhotons )
process.patDefaultSequence.remove( process.countPatElectrons )
process.patDefaultSequence.remove( process.countPatMuons )
process.patDefaultSequence.remove( process.countPatTaus )
process.patDefaultSequence.remove( process.countPatPhotons )
process.patDefaultSequence.remove( process.countPatLeptons )
process.patDefaultSequence.remove( process.makePatElectrons )
process.patDefaultSequence.remove( process.makePatMuons )
process.patDefaultSequence.remove( process.makePatTaus )
process.patDefaultSequence.remove( process.makePatMETs )
process.patDefaultSequence.remove( process.makePatPhotons )
process.patDefaultSequence.remove( process.cleanPatJets )
process.patDefaultSequence.remove( process.cleanPatCandidateSummary )
process.patPF2PATSequencePFlowLoose.remove (process.makePatTausPFlowLoose)
process.patPF2PATSequencePFlowLoose.remove (process.selectedPatTausPFlowLoose)
process.patPF2PATSequencePFlowLoose.remove (process.countPatTausPFlowLoose)
process.patPF2PATSequencePFlowLoose.remove (process.countPatLeptonsPFlowLoose)
process.PFBRECOPFlowLoose.remove( process.pfMETPFlowLoose )
process.patCandidatesPFlowLoose.remove (process.makePatMETsPFlowLoose)
process.patPF2PATSequencePFlowLoose.remove (process.producePatPFMETCorrectionsPFlowLoose) 
process.patPF2PATSequencePFlowLoose.remove (process.patMETsPFlowLoose)
process.PFBRECOPFlowLoose.remove( process.pfJetSequencePFlowLoose )
process.PFBRECOPFlowLoose.remove( process.pfTauSequencePFlowLoose)
process.PFBRECOPFlowLoose.remove( process.pfNoTauPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patJetCorrectionsPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.jetTracksAssociatorAtVertexPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.btaggingAODPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patJetChargePFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patJetPartonMatchPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.genForPF2PATSequence )
process.patPF2PATSequencePFlowLoose.remove ( process.patJetGenJetMatchPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patJetFlavourIdPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patJetsPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.selectedPatJetsPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.countPatJetsPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patHPSPFTauDiscriminationUpdatePFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patPFTauIsolationPFlowLoose )
process.patPF2PATSequencePFlowLoose.remove ( process.patTausPFlowLoose )


if options.runOnFastSim:
    process.patseq = cms.Sequence(
        process.goodOfflinePrimaryVertices*
        process.inclusiveVertexing*
        process.genParticlesForJetsNoNu*
        process.ca8GenJetsNoNu*
        process.ca12GenJetsNoNu*
        process.ca15GenJetsNoNu*
        process.ca8TopGenJetsNoNu*
        process.ca12TopGenJetsNoNu*
        process.ca15TopGenJetsNoNu*
        process.ak8GenJetsNoNu*
        process.caFilteredGenJetsNoNu*
        process.caMassDropFilteredGenJetsNoNu*
        process.caPrunedGen*
        process.caTopTagGen*
        process.ca12TopTagGen*
        process.CATopTagInfosGen*
        process.CA12TopTagInfosGen*
        process.caHEPTopTagGen*
        getattr(process,"patPF2PATSequence"+postfix)*
        process.patDefaultSequence*
        process.goodPatJetsPFlow*
        process.goodPatJetsCA8PF*
        process.goodPatJetsCA12PF*
        process.goodPatJetsCA8PrunedPF*
        process.goodPatJetsCATopTagPF*
        process.goodPatJetsCA12TopTagPF*
        process.goodPatJetsCAHEPTopTagPF*
        process.goodPatJetsCA15PF*
        process.goodPatJetsCA15MassDropFilteredPF*
        process.goodPatJetsCA15FilteredPF*
        process.goodPatJetsCA8PrunedPFPacked*
        process.goodPatJetsCATopTagPFPacked*
        process.goodPatJetsCA12TopTagPFPacked*
        process.goodPatJetsCAHEPTopTagPFPacked*
        process.goodPatJetsCA15MassDropFilteredPFPacked*
        process.goodPatJetsCA15FilteredPFPacked*
        process.flavorHistorySeq*
        process.prunedGenParticles*
        process.kt6PFJetsForIsolation*
        getattr(process,"patPF2PATSequence"+postfixLoose)#*
        #    process.miniPFLeptonSequence
        )
else:    
    process.patseq = cms.Sequence(
        process.filtersSeq*
        process.goodOfflinePrimaryVertices*
        process.inclusiveVertexing*
        process.genParticlesForJetsNoNu*
        process.ca8GenJetsNoNu*
        process.ca12GenJetsNoNu*
        process.ca15GenJetsNoNu*
        process.ca8TopGenJetsNoNu*
        process.ca12TopGenJetsNoNu*
        process.ca15TopGenJetsNoNu*
        process.ak8GenJetsNoNu*
        process.caFilteredGenJetsNoNu*
        process.caMassDropFilteredGenJetsNoNu*
        process.caPrunedGen*
        process.caTopTagGen*
        process.ca12TopTagGen*
        process.CATopTagInfosGen*
        process.CA12TopTagInfosGen*
        process.caHEPTopTagGen*
        getattr(process,"patPF2PATSequence"+postfix)*
        process.patDefaultSequence*
        process.goodPatJetsPFlow*
        process.goodPatJetsCA8PF*
        process.goodPatJetsCA12PF*
        process.goodPatJetsCA8PrunedPF*
        process.goodPatJetsCATopTagPF*
        process.goodPatJetsCA12TopTagPF*
        process.goodPatJetsCAHEPTopTagPF*
        process.goodPatJetsCA15PF*
        process.goodPatJetsCA15MassDropFilteredPF*
        process.goodPatJetsCA15FilteredPF*
        process.goodPatJetsCA8PrunedPFPacked*
        process.goodPatJetsCATopTagPFPacked*
        process.goodPatJetsCA12TopTagPFPacked*
        process.goodPatJetsCAHEPTopTagPFPacked*
        process.goodPatJetsCA15MassDropFilteredPFPacked*
        process.goodPatJetsCA15FilteredPFPacked*
        process.flavorHistorySeq*
        process.prunedGenParticles*
        process.kt6PFJetsForIsolation*
        getattr(process,"patPF2PATSequence"+postfixLoose)#*
        #    process.miniPFLeptonSequence
        )

if options.useExtraJetColls:
	process.extraJetSeq = cms.Sequence(
		process.ak7TrimmedGenJetsNoNu*
		process.ak7FilteredGenJetsNoNu*
		process.ak7PrunedGenJetsNoNu*
                process.goodPatJetsCA15PF*
		process.goodPatJetsCA15FilteredPF*
		process.goodPatJetsCA15MassDropFilteredPF*
		process.goodPatJetsAK5TrimmedPF*
		process.goodPatJetsAK5FilteredPF*
		process.goodPatJetsAK5PrunedPF*
		process.goodPatJetsAK7PF*
		process.goodPatJetsAK7TrimmedPF*
		process.goodPatJetsAK7FilteredPF*
		process.goodPatJetsAK7PrunedPF*
		process.goodPatJetsAK8PF*
		process.goodPatJetsAK8TrimmedPF*
		process.goodPatJetsAK8FilteredPF*
		process.goodPatJetsAK8PrunedPF*
		process.ak5Lite*
		process.ak5TrimmedLite*
		process.ak5FilteredLite*
		process.ak5PrunedLite*
		process.ak7Lite*
		process.ak7TrimmedLite*
		process.ak7FilteredLite*
		process.ak7PrunedLite*
		process.ak7TrimmedGenLite*
		process.ak7FilteredGenLite*
		process.ak7PrunedGenLite*
		process.ak8Lite*
		process.ak8TrimmedLite*
		process.ak8FilteredLite*
		process.ak8PrunedLite*
                process.ca8PrunedGenLite*
                process.ca12FilteredGenLite*
                process.ca12MassDropFilteredGenLite
	)
	process.patseq *= process.extraJetSeq


if options.useData :
    process.patseq.remove( process.genParticlesForJetsNoNu )
    process.patseq.remove( process.genJetParticles )
    process.patseq.remove( process.ak8GenJetsNoNu )
    process.patseq.remove( process.ca8GenJetsNoNu )
    process.patseq.remove( process.ca12GenJetsNoNu )
    process.patseq.remove( process.ca15GenJetsNoNu )
    process.patseq.remove( process.ca8TopGenJetsNoNu )
    process.patseq.remove( process.ca12TopGenJetsNoNu )
    process.patseq.remove( process.ca15TopGenJetsNoNu )
    process.patseq.remove( process.caFilteredGenJetsNoNu )
    process.patseq.remove( process.flavorHistorySeq )
    process.patseq.remove( process.caPrunedGen )
    process.patseq.remove( process.caHEPTopTagGen)
    process.patseq.remove( process.caTopTagGen )
    process.patseq.remove( process.ca12TopTagGen )
    process.patseq.remove( process.CATopTagInfosGen )
    process.patseq.remove( process.CA12TopTagInfosGen )
    process.patseq.remove( process.prunedGenParticles )
    process.patseq.remove( process.caMassDropFilteredGenJetsNoNu )
    process.patseq.remove( process.caFilteredGenJetsNoNu )

    if options.useExtraJetColls:
	    process.patseq.remove( process.ak8GenJetsNoNu )
	    process.patseq.remove( process.caFilteredGenJetsNoNu )
	    process.patseq.remove( process.ak7TrimmedGenJetsNoNu )
	    process.patseq.remove( process.ak7FilteredGenJetsNoNu )
	    process.patseq.remove( process.ak7PrunedGenJetsNoNu )
	    process.patseq.remove( process.ak7TrimmedGenLite )
	    process.patseq.remove( process.ak7FilteredGenLite )
	    process.patseq.remove( process.ak7PrunedGenLite )
            process.patseq.remove( process.ca8PrunedGenLite )
            process.patseq.remove( process.ca12FilteredGenLite )
            process.patseq.remove( process.ca12MassDropFilteredGenLite )

if options.runOnFastSim:
    process.patseq.remove( process.HBHENoiseFilter )
    process.patseq.remove( process.CSCTightHaloFilter ) 

if options.writeSimpleInputs :
	process.patseq *= cms.Sequence(process.pfInputs)

if options.useSusyFilter :
	process.patseq.remove( process.HBHENoiseFilter )
	process.load( 'PhysicsTools.HepMCCandAlgos.modelfilter_cfi' )
	process.modelSelector.parameterMins = [500.,    0.] # mstop, mLSP
	process.modelSelector.parameterMaxs  = [7000., 200.] # mstop, mLSP
	process.p0 = cms.Path(
		process.modelSelector *
		process.patseq
	)



else :
	process.p0 = cms.Path(
		process.patseq
	)

process.MyNtuple = cms.EDAnalyzer('NtupleWriter',
                                  fileName = cms.string("Ntuple.root"), 
                                  doElectrons = cms.bool(True),
                                  doMuons = cms.bool(True),
                                  doTaus = cms.bool(True),
                                  doJets = cms.bool(True),
                                  doTopJets = cms.bool(True),
                                  doTopJetsConstituents = cms.bool(True),
                                  doGenJets = cms.bool(not options.useData),
                                  doPhotons = cms.bool(False),
                                  doMET = cms.bool(True),
                                  doPV = cms.bool(True),
                                  doAllPFParticles = cms.bool(options.writePFCands),
				  storePFsAroundLeptons = cms.untracked.bool(False),
                                  doGenInfo = cms.bool(not options.useData),
				  doAllGenParticles = cms.bool(options.writeAllGenParticles), #set to true if you want to store all gen particles, otherwise, only tops and status 3 particles are stored
				  doLumiInfo = cms.bool(options.useData),
                                  doTrigger = cms.bool(True),
                                  doTagInfos = cms.untracked.bool(True), # when set to true crashes for the 'packed' jet collections
                                  svComputer = cms.untracked.InputTag("combinedSecondaryVertex"),
				  rho_source = cms.InputTag("kt6PFJets", "rho"),
                                  genparticle_source = cms.InputTag("prunedGenParticles" ),
                                  electron_sources = cms.vstring("selectedPatElectronsPFlow","selectedPatElectronsPFlowLoose"),
                                  muon_sources = cms.vstring("selectedPatMuonsPFlow","selectedPatMuonsPFlowLoose"),
                                  tau_sources = cms.vstring("selectedPatTausPFlow" ),
                                  tau_ptmin = cms.double(0.0),
                                  tau_etamax = cms.double(999.0),
                                  jet_sources = cms.vstring("goodPatJetsPFlow","goodPatJetsCA15PF"),
                                  jet_ptmin = cms.double(10.0),
                                  jet_etamax = cms.double(5.0),
                                  genjet_sources = cms.vstring("ak5GenJetsNoNu", "ca8GenJetsNoNu", "ca12GenJetsNoNu", "ca15GenJetsNoNu"),
                                  genjet_ptmin = cms.double(10.0),
                                  genjet_etamax = cms.double(5.0),
				  #photon_sources = cms.vstring("selectedPatPhotons"),
                                  topjet_sources = cms.vstring("goodPatJetsCATopTagPFPacked","goodPatJetsCA12TopTagPFPacked", "goodPatJetsCA15FilteredPFPacked", "goodPatJetsCAHEPTopTagPFPacked"),
                                  topjet_constituents_sources = cms.vstring("goodPatJetsCA8PF","goodPatJetsCA12PF", "goodPatJetsCA15PF"),
                                  #topjet_constituents_sources = cms.vstring(),
                                  topjet_ptmin = cms.double(30.0), 
                                  topjet_etamax = cms.double(5.0),
                                  pf_around_leptons_sources = cms.vstring("pfNoPileUpIsoPFlow", "pfPileUpIsoPFlow"),
                                  pf_collection_source = cms.string("pfNoElectronPFlow"),
				  doGenTopJets = cms.bool(not options.useData),			      
                                  gentopjet_sources = cms.vstring("caTopTagGen","ca12TopTagGen", "caFilteredGenJetsNoNu", "caHEPTopTagGen"),
                                  gentopjet_ptmin = cms.double(150.0), 
                                  gentopjet_etamax = cms.double(5.0),
  				  doGenJetsWithParts = cms.bool(not options.useData),
                                  genjetwithparts_sources = cms.vstring("ca8TopGenJetsNoNu","ca12TopGenJetsNoNu", "ca15TopGenJetsNoNu"),
                                  genjetwithparts_ptmin = cms.double(150.0), 
                                  genjetwithparts_etamax = cms.double(5.0),
                                  met_sources =  cms.vstring("patMETsPFlow"),
                                  pv_sources = cms.vstring("goodOfflinePrimaryVertices"),
                                  trigger_prefixes = cms.vstring(#"HLT_IsoMu", "HLT_Mu",
                                                                 #"HLT_L1SingleMu", "HLT_L2Mu",
                                                                 #"HLT_Ele",
                                                                 "HLT_",
                                                                 #"HLT_DoubleMu", "HLT_DoubleEle"
	                                                         ),
                                  
)

process.p0 *= process.MyNtuple

process.outpath = cms.EndPath()


# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(50)


# process all the events
process.maxEvents = cms.untracked.PSet( input =cms.untracked.int32(__MAX_EVENTS__))
process.options.wantSummary = False
process.out.dropMetaData = cms.untracked.string("DROPPED")

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

#open('junk.py','w').write(process.dumpPython())


