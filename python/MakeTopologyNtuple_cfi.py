import FWCore.ParameterSet.Config as cms

makeTopologyNtuple = cms.EDAnalyzer('MakeTopologyNtuple',
                                    # "Calo"
                                    electronTag = cms.InputTag("cleanPatElectrons"),
                                    tauTag      = cms.InputTag("cleanPatTaus"),
                                    muonTag     = cms.InputTag("cleanPatMuons"),
                                    jetTag      = cms.InputTag("cleanPatJets"),
                                    genJetTag   = cms.InputTag("ak5GenJetsNoNu"),
                                    photonTag   = cms.InputTag("cleanPatPhotons"),
                                    metTag      = cms.InputTag("patMETs"),
                                    # PF
                                    electronPFTag = cms.InputTag("selectedPatElectronsPFlow"),
                                    tauPFTag      = cms.InputTag("selectedPatTausPFlow"),
                                    muonPFTag     = cms.InputTag("selectedPatMuonsPFlow"),
                                    jetPFTag      = cms.InputTag("selectedPatJetsPFlow"),
                                    jetPFRecoTag  = cms.InputTag("selectedPatJetsAK5PF"),
#                                    photonPFTag   = cms.InputTag("cleanPatPhotons"), #
                                    metPFTag      = cms.InputTag("patType1CorrectedPFMETPF2PAT"),
                                    # JPT
                                    jetJPTTag         = cms.InputTag("selectedPatJetsAK5JPT"),
                                    metJPTTag      = cms.InputTag("patMETsTC"),
                                    primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
                                    rho              = cms.InputTag("kt6PFJets", "rho"),
                                    triggerTag  = cms.InputTag("TriggerResults","","HLT8E29"),
                                    fakeTriggerList = cms.vstring(), # empty. You can add fake triggers that are run on the fly to this list. No check on the process name is made so when duplicates are available only the latest one is added.

                                    triggerList = cms.vstring(                                                              #Updated Triggers
                                                              #Menu 5E33
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8',
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9',
                                                              'HLT_Mu17_Mu8_v13',
                                                              'HLT_Mu17_Mu8_v14',
                                                              'HLT_Mu17_Mu8_v15',
                                                              'HLT_Mu17_Mu8_v16',
                                                              'HLT_Mu17_Mu8_v17',
                                                              'HLT_Mu17_Mu8_v18',
                                                              'HLT_Mu17_Mu8_v19',
                                                              'HLT_Mu17_Mu8_v15',
                                                              'HLT_Mu17_Mu8_v14',
                                                              'HLT_Mu17_Mu8_v20',
                                                              'HLT_Mu17_Mu8_v21',
                                                              'HLT_Mu17_Mu8_v22',
                                                              'HLT_Mu17_TkMu8_v6',
                                                              'HLT_Mu17_TkMu8_v7',
                                                              'HLT_Mu17_TkMu8_v8',
                                                              'HLT_Mu17_TkMu8_v9',
                                                              'HLT_Mu17_TkMu8_v10',
                                                              'HLT_Mu17_TkMu8_v11',
                                                              'HLT_Mu17_TkMu8_v12',
                                                              'HLT_Mu17_TkMu8_v13',
                                                              'HLT_Mu17_TkMu8_v14',
                                                              #MC Menu
                                                              'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v9',
                                                              'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v9',
                                                              'HLT_Mu17_Mu8_v12',
                                                              'HLT_Mu17_TkMu8_v5',
                                                              'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11'
                                                              #Menu 7E33
                                                               ),
                                    l1TriggerTag = cms.InputTag("gtDigis"),                                    
                                    checkTriggers = cms.bool(True),
                                    genParticles = cms.InputTag("genParticles"),
                                    runMCInfo = cms.bool(True), # if set to true will skip MCInfo section
                                    doJERSmear = cms.bool(True), # as run MC is true, may as well be true too.
                                    runPUReWeight = cms.bool(False), #Run pile-up reweighting. Don't do if this is data I guess.
                                    doCuts = cms.bool(True), # if set to true will skip ALL cuts. Z veto still applies electron cuts.
                                    # default preselection settings! see https://twiki.cern.ch/twiki/bin/view/CMS/VplusJets for inspiration


                                    #Some jet cuts.
                                    minJetPt = cms.double(0), #min jet pT in GeV/c
                                    maxJetEta = cms.double(5), # jet |eta|
                                    isPF = cms.bool(False), #Particle flow or not
                                    jetMinConstituents = cms.double(2),
                                    jetNHEF = cms.double(0.99), 
                                    jetNEEF = cms.double(0.99),
                                    ecalEndRejectAngle = cms.double(2.4), #Angle before which we care about the charged fraction
                                    jetCEF = cms.double(0.99),
                                    jetCHF = cms.double(0.),
                                    jetNCH = cms.double(0.),
                                    jetPtCutLoose = cms.double(20.),
                                    
                                    #electron ID (tight, for analysis)
                                    runSwissCross = cms.bool(True),
                                    runReweightTest = cms.bool(False), #This is just a little test to see what happens when using the default reweight class, i.e. whether it breaks like my one does.
                                    runPDFUncertainties = cms.bool(False),
                                    useResidualJEC = cms.bool(False),
                                    electronID = cms.string('eidRobustTight'),
                                    ignoreElectronID = cms.bool(True), # if set to true will save all electrons, also those not passing electronID.
                                    minEleEt = cms.double(0), #  electron ET in GeV
                                    eleMvaCut=cms.double(0.0), #mva minimum. Maximum mva value is hard-coded as 1, as I think that's the highest it can be.
                                    maxEleEta = cms.double(5), #  electron |eta|
                                    eleCombRelIso = cms.double(1.0), # V+jets adviced cut: 0.1. 
                                    maxEled0 = cms.double(500), # D0 (beam spot corrected) in cm, V+jets adviced cut: 0.02 for prompt electrons
                                    eleInterECALEtaLow = cms.double(1.4442),
                                    eleInterECALEtaHigh = cms.double(1.5660),
                                    # electron ID (loose, for Z Veto)
                                    electronIDLooseZVeto = cms.string('eidRobustLoose'),
                                    minEleEtLooseZVeto = cms.double(0),
                                    maxEleEtaLooseZVeto = cms.double(5),
                                    eleCombRelIsoLooseZVeto = cms.double(1.0), # always accept the electron
                                    maxEled0LooseZVeto = cms.double(100.), # always accept the electron
                                    # cross-cleaning parameter (jet is rejected if inside electron cone0
                                    dREleJetCrossClean = cms.double(-1), # cone distance that an electron-jet match needs to have to reject the jet (the electron is always kept)
                                    # muon identification
                                    maxMuonEta = cms.double(5),
                                    minMuonPt = cms.double(0.),
                                    maxMuonD0 = cms.double(50), # D0 in cm, already corrected for Beam position. V+jets def: 0.02 (for prompt muons)
                                    muoCombRelIso = cms.double(1.), #combined track isolation,
                                    muoNormalizedChi2 = cms.double(5000), #normalized chi2 (Chi2/NDOF)
                                    muoNTrkHits = cms.double(0), # minimal number of track hits
                                    muonECalIso = cms.double(7000),
                                    muonHCalIso = cms.double(7000),
                                    flavorHistoryTag = cms.bool(True),
                                    muoValidHits = cms.double(1), #at least one valid muon hit
                                    muonMatchedStations = cms.double(2),
                                    muonDZCut = cms.double(0.5),
                                    muonDBCut = cms.double(0.2),
                                    muonPixelHits = cms.double(1),#minimum of one
                                    muonTrackLayersWithHits = cms.double(6), # 5 or less is skipped.
                                    muonRelIsoTight = cms.double(0.12),
                                    muonPtLoose = cms.double(10.),
                                    muonEtaLoose = cms.double(2.5),
                                    muoRelIsoLoose = cms.double(0.2),
                                    metCut = cms.double(30.0),
                                    fillAll = cms.bool(False),
                                    # photon rejection:
                                    dREleGeneralTrackMatchForPhotonRej=cms.double(0.3),
                                    magneticFieldForPhotonRej=cms.double(3.8),
                                    correctFactorForPhotonRej=cms.double(-0.003),
                                    maxDistForPhotonRej=cms.double(0),
                                    maxDcotForPhotonRej=cms.double(0),
                                    ebRecHits=cms.InputTag('reducedEcalRecHitsEB'),
                                    eeRecHits=cms.InputTag('reducedEcalRecHitsEE'),
                                    isMCatNLO=cms.bool(False),
                                    #New B-tagging info
                                    bDiscCut=cms.double(0),
                                    bDiscName=cms.string('combinedSecondaryVertexBJetTags'),

                                    # Btagging algorithms to look at (default discriminant is used). The pat::Jet::bDiscriminator(string) function is used.


                                    
                                    btagAlgorithmsToNtuple = cms.vstring("trackCountingHighPurBJetTags","trackCountingHighEffBJetTags",
                                                                         "jetProbabilityBJetTags", "jetBProbabilityBJetTags", "softElectronBJetTags","softMuonBJetTags","softMuonNoIPBJetTags",
                                                                         "simpleSecondaryVertexHighEffBJetTags","simpleSecondaryVertexNegativeBJetTags",
                                                                         "combinedSecondaryVertexMVABJetTags", "simpleSecondaryVertexHighPurBJetTags"
    
                                                                         ),
                                    # Btagging parameterizations to look at (the vectors btagParameterizationList and btagParameterizationMode are coupled!). Documentation on algo names (go in btagParamerizationList) and parameterizations (go in btagParameterizationMode) are available on this twiki:
                                    # https://twiki.cern.ch/twiki/bin/view/CMS/BtagOctober09ExerciseUsePayload
                                    btagParameterizationList = cms.vstring(
                                        #MC Measurements
                                        "MCCaloSSVHPTb","MCCaloSSVHPTc","MCCaloSSVHPTl",
                                        "MCCaloSSVHEMb","MCCaloSSVHEMc","MCCaloSSVHEMl",
                                        "MCCaloSSVHETb","MCCaloSSVHETc","MCCaloSSVHETl",
                                        "MCCaloTCHELb","MCCaloTCHELc","MCCaloTCHELl",
                                        "MCCaloTCHEMb","MCCaloTCHEMc","MCCaloTCHEMl",
                                        "MCCaloTCHETb","MCCaloTCHETc","MCCaloTCHETl",
                                        #MC Measurements Errors
                                        "MCCaloSSVHPTb","MCCaloSSVHPTc","MCCaloSSVHPTl",
                                        "MCCaloSSVHEMb","MCCaloSSVHEMc","MCCaloSSVHEMl",
                                        "MCCaloSSVHETb","MCCaloSSVHETc","MCCaloSSVHETl",
                                        "MCCaloTCHELb","MCCaloTCHELc","MCCaloTCHELl",
                                        "MCCaloTCHEMb","MCCaloTCHEMc","MCCaloTCHEMl",
                                        "MCCaloTCHETb","MCCaloTCHETc","MCCaloTCHETl",
                                        #Mistag fall10
                                        "MISTAGSSVHEM", "MISTAGSSVHEM", "MISTAGSSVHEM", "MISTAGSSVHEM",
                                        "MISTAGSSVHPT", "MISTAGSSVHPT", "MISTAGSSVHPT", "MISTAGSSVHPT",
                                        "MISTAGTCHEL", "MISTAGTCHEL", "MISTAGTCHEL", "MISTAGTCHEL",
                                        "MISTAGTCHEM", "MISTAGTCHEM", "MISTAGTCHEM", "MISTAGTCHEM"
                                        ),
                                    btagParameterizationMode = cms.vstring(
                                        "BTAGBEFF", "BTAGCEFF", "BTAGLEFF",
                                        "BTAGBEFF", "BTAGCEFF", "BTAGLEFF",
                                        "BTAGBEFF", "BTAGCEFF", "BTAGLEFF",
                                        "BTAGBEFF", "BTAGCEFF", "BTAGLEFF",
                                        "BTAGBEFF", "BTAGCEFF", "BTAGLEFF",
                                        "BTAGBEFF", "BTAGCEFF", "BTAGLEFF",
                                        "BTAGBERR", "BTAGCERR", "BTAGLERR",
                                        "BTAGBERR", "BTAGCERR", "BTAGLERR",
                                        "BTAGBERR", "BTAGCERR", "BTAGLERR",
                                        "BTAGBERR", "BTAGCERR", "BTAGLERR",
                                        "BTAGBERR", "BTAGCERR", "BTAGLERR",
                                        "BTAGBERR", "BTAGCERR", "BTAGLERR",
                                        #Mistag fall10
                                        "BTAGLEFF", "BTAGLERR", "BTAGLEFFCORR", "BTAGLERRCORR",
                                        "BTAGLEFF", "BTAGLERR", "BTAGLEFFCORR", "BTAGLERRCORR",
                                        "BTAGLEFF", "BTAGLERR", "BTAGLEFFCORR", "BTAGLERRCORR",
                                        "BTAGLEFF", "BTAGLERR", "BTAGLEFFCORR", "BTAGLERRCORR"
                                        ),
                                    runCutFlow = cms.double(0), #0 is no cut flow, 1 is ee, 2 is emu, 3 mumu.
                                    isttBar = cms.bool(False),# This affects reweighting things. If set to false, then has a weight of 1.
                                    ttGenEvent = cms.InputTag("null")
                                    )# end of MakeTopologyNtuple
