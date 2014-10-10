#Set up the pat environment
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

#Import the pat configurations
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#Setting up various environmental stuff that makes all of this jazz actually work.

###############################
####### Global Setup ##########
###############################


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
process.load("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")

process.load("FWCore.Framework.test.cmsExceptionsFatal_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.destinations = ['cerr']
process.MessageLogger.statistics = []
process.MessageLogger.fwkJobReports = []
process.MessageLogger.categories=cms.untracked.vstring('FwkJob'
                                                       ,'FwkReport'
                                                       ,'FwkSummary'
                                                       )

process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(
                     wantSummary = cms.untracked.bool(True)
                     )

process.GlobalTag.globaltag = cms.string('FT_53_V21_AN6::All')

#There's a bit in here about some btau tags that the code looks for. I don't know if this is significant, however. I'm going to ignore it for now.


#Import jet reco things. Apparently this makes cmsRun crash.
process.load('RecoJets.Configuration.RecoPFJets_cff')

#Now do cool fast jet correction things!

process.kt6PFJets.doRhoFastjet = True

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )



from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = cms.PSet( minNdof = cms.double( 4. ) ,
                             maxZ = cms.double( 24. ) ,
                             maxRho = cms.double( 2. ) ) ,
    filter = cms.bool( True) ,
    src = cms.InputTag( 'offlinePrimaryVertices' ) )

## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")

## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

## The ECAL laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

## The tracking POG filters __________________________________________________||
process.load('RecoMET.METFilters.trackingPOGFilters_cff')

process.goodVertices = cms.EDFilter(
      "VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
      )

process.filtersSeq = cms.Sequence(
    process.HBHENoiseFilter *
    process.CSCTightHaloFilter *
    process.hcalLaserEventFilter *
    process.EcalDeadCellTriggerPrimitiveFilter *
    process.goodVertices *
    process.trackingFailureFilter *
    process.eeBadScFilter *
    process.ecalLaserCorrFilter *
    process.trkPOGFilters
    )


#Gen Setup - I'm unsure what this does, and I can't actually do it anyway as I don't think TopQuarkAnalysis actually exists? Or I may have to import it. I don't really know what this does. I assume this is related to that bit of global tags that I didn't do, so I might just ignore this for now.
#process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")





###############################
####### PF2PAT Setup ##########
###############################

# Default PF2PAT with AK5 jets. Make sure to turn ON the L1fastjet stuff.
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PF2PAT"
usePF2PAT(process,runPF2PAT=True, jetAlgo="AK5", runOnMC=False, postfix=postfix, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True)


getattr(process,"pfNoPileUp"  +postfix).enable = True
getattr(process,"pfNoMuon"    +postfix).enable = False
getattr(process,"pfNoElectron"+postfix).enable = False
getattr(process,"pfNoTau"     +postfix).enable = False
getattr(process,"pfNoJet"     +postfix).enable = False

#process.patJetsPF2PAT.discriminatorSources = cms.VInputTag(
#            cms.InputTag("combinedSecondaryVertexBJetTagsAODPF2PAT"),
#            cms.InputTag("combinedSecondaryVertexMVABJetTagsAODPF2PAT"),
#            )



#Gsf Electrons
#useGsfElectrons(process, postfix, "03")

#process.pfIsolatedMuonsPF2PAT.doDeltaBetaCorrection = True
process.pfIsolatedMuonsPF2PAT.isolationCut = cms.double(9999.)
process.pfIsolatedMuonsPF2PAT.combinedIsolationCut = cms.double(9999.)
#process.pfSelectedMuonsPF2PAT.cut = cms.string('pt > 10. && abs(eta) < 2.5')
process.pfIsolatedMuonsPF2PAT.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("muPFIsoValueCharged04PF2PAT"))
process.pfIsolatedMuonsPF2PAT.deltaBetaIsolationValueMap = cms.InputTag("muPFIsoValuePU04PF2PAT")
process.pfIsolatedMuonsPF2PAT.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("muPFIsoValueNeutral04PF2PAT"), cms.InputTag("muPFIsoValueGamma04PF2PAT"))

process.pfIsolatedElectronsPF2PAT.isolationCut = cms.double(9999.)
process.pfIsolatedElectronsPF2PAT.combinedIsolationCut = cms.double(9999.)
#process.pfIsolatedElectronsPF2PAT.doDeltaBetaCorrection = True
#process.pfSelectedElectronsPF2PAT.cut = cms.string('pt > 15. && abs(eta) < 2.5')
process.pfIsolatedElectronsPF2PAT.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"))
process.pfIsolatedElectronsPF2PAT.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT")
process.pfIsolatedElectronsPF2PAT.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPF2PAT"), cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT"))


#I don't think I need these right now, but it can't hurt to include them. Unless it breaks my code. Which is entirely possible. It did! I knew it!
process.patElectronsPF2PAT.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPF2PAT"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPF2PAT"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPF2PAT"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPF2PAT"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPF2PAT")
    )

#process.patElectronsPF2PAT.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
#process.patElectronsPF2PAT.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
#process.patPF2PATSequencePF2PAT.replace( process.patElectronsPF2PAT, process.eidMVASequence * process.patElectronsPF2PAT )

#Now do a bit of JEC
process.patJetCorrFactorsPF2PAT.payload = 'AK5PFchs'
process.patJetCorrFactorsPF2PAT.levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])
process.pfPileUpPF2PAT.checkClosestZVertex = False



###############################
###### Electron ID ############
###############################


#process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
#process.eidMVASequence = cms.Sequence( process.mvaTrigV0 + process.mvaNonTrigV0 )



process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
#Electron ID
process.patElectronsPF2PAT.electronIDSources.mvaTrigV0	 = cms.InputTag("mvaTrigV0")
process.patElectronsPF2PAT.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0") 
process.patPF2PATSequencePF2PAT.replace( process.patElectronsPF2PAT, process.eidMVASequence * process.patElectronsPF2PAT )


#Convesion Rejection
# this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
process.patConversionsPF2PAT = cms.EDProducer("PATConversionProducer",
                                             electronSource = cms.InputTag("selectedPatElectronsPF2PAT")      
                                             )
					     
process.patPF2PATSequencePF2PAT += process.patConversionsPF2PAT




###############################
###### Bare KT 0.6 jets ####### Because why not?
###############################

#from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
# For electron (effective area) isolation
#process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#mprocess.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

###############################
#### Selections Setup #########
###############################

# AK5 Jets
#   PF
process.selectedPatJetsPF2PAT.cut = cms.string("pt > 5.0")

# Flavor history stuff - don't really know what this is, but it was in the other one too so I guess I need to include it.
#process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")
#process.flavorHistoryFilter.pathToSelect = cms.int32(-1)



#Letting pat run
process.patseq = cms.Sequence(
  #  process.kt6PFJetsForIsolation*
    process.goodOfflinePrimaryVertices*
    process.primaryVertexFilter * #removes events with no good pv (but if cuts to determine good pv change...)
    process.filtersSeq *
    getattr(process,"patPF2PATSequence"+postfix) # main PF2PAT
#   * process.flavorHistorySeq
    )

####
# The N-tupliser/cutFlow
####

triggerStringName = 'HLT'

process.load("NTupliser.SingleTop.MakeTopologyNtuple_cfi")
process.makeTopologyNtuple.flavorHistoryTag=cms.bool(False) # change to false at your convenience
process.makeTopologyNtuple.runMCInfo=cms.bool(False) # prevent checking gen info
process.makeTopologyNtuple.doJERSmear=cms.bool(False)
#process.makeTopologyNtuple.doCuts=cms.bool(True) # if set to false will skip ALL cuts. Z veto still applies electron cuts.
process.makeTopologyNtuple.triggerTag = cms.InputTag("TriggerResults","",triggerStringName) # or HLT, depends on file   

#settings to apply tight selection:
process.makeTopologyNtuple.minJetPt=cms.double(30)
process.makeTopologyNtuple.maxJetEta=cms.double(2.5)
process.makeTopologyNtuple.bDiscCut=cms.double(0.679)
process.makeTopologyNtuple.minEleEt=cms.double(20)
process.makeTopologyNtuple.maxEleEta=cms.double(2.5)
process.makeTopologyNtuple.ignoreElectronID=cms.bool(False)
process.makeTopologyNtuple.eleCombRelIso=cms.double(0.15)
process.makeTopologyNtuple.maxEled0=cms.double(0.04)
process.makeTopologyNtuple.eleInterECALEtaLow=cms.double(1.4442)
process.makeTopologyNtuple.eleInterECALEtaHigh=cms.double(1.5660)
process.makeTopologyNtuple.minEleEtLooseZVeto=cms.double(15)
process.makeTopologyNtuple.minEleEtaLooseZVeto=cms.double(2.5)
process.makeTopologyNtuple.eleCombRelIsoLooseZVeto=cms.double(0.15)
process.makeTopologyNtuple.dREleJetCrossClean=cms.double(0.3)
process.makeTopologyNtuple.maxMuonEta=cms.double(2.4)
process.makeTopologyNtuple.minMuonPt=cms.double(20)
process.makeTopologyNtuple.maxMuonD0=cms.double(0.02)
process.makeTopologyNtuple.muonRelIsoTight=cms.double(0.2)
process.makeTopologyNtuple.muoNormalizedChi2=cms.double(10)
process.makeTopologyNtuple.muoNTrkHits=cms.double(11)
process.makeTopologyNtuple.muonECalIso=cms.double(4)
process.makeTopologyNtuple.muonHCalIso=cms.double(6)
process.makeTopologyNtuple.dREleGeneralTrackMatchForPhotonRej=cms.double(0.3)
process.makeTopologyNtuple.maxDistForPhotonRej=cms.double(0.04)
process.makeTopologyNtuple.maxDcotForPhotonRej=cms.double(0.03)
process.makeTopologyNtuple.fillAll=cms.bool(True)
process.makeTopologyNtuple.processingLoose=cms.bool(False)
#process.makeTopologyNtuple.btagParameterizationList = cms.vstring()
#process.makeTopologyNtuple.btagParameterizationMode = cms.vstring()
process.makeTopologyNtuple.runSwissCross = cms.bool(False)
#Don't actually do cuts
process.makeTopologyNtuple.doCuts = cms.bool(False)
process.makeTopologyNtuple.runCutFlow=cms.double(0)

#Make the inputs for the n-tupliser right.
process.makeTopologyNtuple.electronPFTag = cms.InputTag("selectedPatElectronsPF2PAT")
process.makeTopologyNtuple.tauPFTag = cms.InputTag("selectedPatTausPF2PAT")
process.makeTopologyNtuple.muonPFTag = cms.InputTag("selectedPatMuonsPF2PAT")
process.makeTopologyNtuple.jetPFTag = cms.InputTag("selectedPatJetsPF2PAT")
process.makeTopologyNtuple.metPFTag = cms.InputTag("patType1CorrectedPFMetPF2PAT")                                                                                  
#For now this is just the patseq, but soon this will also involve the ntupliser. And then minor corrections for the data version which will include more filters and such.
process.p = cms.Path(
    process.patseq
*    process.makeTopologyNtuple
    )

process.source.fileNames = [
    #Danny's selection
#    "file:/afs/cern.ch/work/l/leggat/FA8766A2-38EA-E111-9624-001A92811738.root",
#    "file:eeMetDump/res/pickevents_1_1_NM8.root",
    "file:synchFiles/data1.root",
    "file:synchFiles/data2.root",
    "file:synchFiles/data3.root"
#    "file:synchFiles/synch1.root",
#    "file:synchFiles/synch2.root",
#    "file:synchFiles/synch3.root",
#    "file:synchFiles/synch4.root",
#    "file:synchFiles/synch5.root",
#    "file:synchFiles/synch6.root"
#        "file:/afs/cern.ch/work/j/jandrea/public/147DB408-446A-E311-8E63-00259073E32A.root"
    ]

process.maxEvents.input = cms.untracked.int32(-1)

from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands += patEventContent
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += cms.untracked.vstring('keep *_flavorHistoryFilter_*_*','keep *_TriggerResults_*_*','keep *_selectedPat*_*_*', 'keep *_*goodOfflinePrimaryVertices*_*_*','keep double_*_rho_*', 'keep patMuons_*_*_*', 'keep *MET*_*_*_*', 'keep *_*MET*_*_*')

#PAT output and various other outpath stuff which is a bit dumb coz I'm probably not even gonna use the outpath. Nevermind.
process.out.fileName = cms.untracked.string('Data_out.root')

#NTuple output
process.TFileService = cms.Service("TFileService", fileName = cms.string('Data_Ntuple_out.root') )
process.options.wantSummary = False
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))

#To run with pat output:
#process.Fin = cms.EndPath(process.out)
#process.schedule = cms.Schedule(process.p, process.Fin)

#Removing pat output (coz we really don't need it now)
del process.out
del process.outpath

process.schedule = cms.Schedule(process.p)
