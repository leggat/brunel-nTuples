/// -*- C++ -*-
//
// Package:    MakeTopologyNtuple
// Class:      MakeTopologyNtuple
// %
/**\class MakeTopologyNtuple MakeTopologyNtuple.cc FreyaAnalysis/MakeTopologyNtuple/src/MakeTopologyNtuplecc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Wed April 22 19:23:10 CET 2009
// $Id: MakeTopologyNtuple.h,v 1.68 2010/11/05 15:32:16 chadwick Exp $
//
//

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

using namespace cu_ejetmet;

class MakeTopologyNtuple : public edm::EDAnalyzer {
public:
  explicit MakeTopologyNtuple(const edm::ParameterSet&);
  ~MakeTopologyNtuple();


private:
//  virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;
  std::map<std::string,TH1D*> histocontainer_; // simple map to contain all histograms. Histograms are booked in the beginJob() method
  std::map<std::string,TH2D*> histocontainer2D_; // simple map to contain all histograms. Histograms are booked in the beginJob() method (2D)
  
  edm::InputTag eleLabel_;
  edm::InputTag muoLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag tauLabel_;
  edm::InputTag metLabel_;
  edm::InputTag phoLabel_;
  edm::InputTag electronPFTag_;	
  edm::InputTag tauPFTag_;	
  edm::InputTag muonPFTag_;	
  edm::InputTag jetPFTag_;	
  edm::InputTag genJetTag_;
  edm::InputTag jetPFRecoTag_;	
  edm::InputTag metPFTag_;	
  edm::InputTag jetJPTTag_;	
  edm::InputTag metJPTTag_;      
  edm::InputTag trigLabel_;
  edm::InputTag ttGenEvent_;
  std::vector<std::string> fakeTrigLabelList_;
  std::vector<std::string> triggerList_;
  edm::InputTag l1TrigLabel_;
  edm::InputTag genParticles_;
  edm::InputTag pvLabel_;
  edm::InputTag rho_;

  std::map<std::string,int> hltpasses_;
  std::vector<std::string> hltnames_;
  std::vector<std::string> btaggingparamnames_;
  std::vector<std::string> btaggingparaminputtypes_;
  std::vector<std::string> btaggingtontuplenames_;// these are the discriminant outputs (so the data) that are saved into the ntuple...
  std::map<std::string,PerformanceResult::ResultType> btaggingparamtype_;
  std::vector<std::string> eleIDsToNtuple_;

  bool filledBIDInfo_;
  bool runMCInfo_;
  bool runPUReWeight_;
  bool doCuts_;

  //jet cuts
  double jetPtCut_;
  double jetEtaCut_;
  double jetMinConstituents_;
  double jetNHEF_;
  double jetNEEF_;
  double ecalEndRejectAngle_;
  double jetCEF_;
  double jetCHF_;
  double jetNCH_;
  std::string bDiscName_; //The name of the bdsicriminator to be used for the analysis.
  double bDiscCut_; //The cut applied for the bDiscriminator;
  double jetPtCutLoose_;

  //bool runSwissCross_;
  bool runReweightingTests_;
  bool runPDFUncertainties_;
  bool useResidualJEC_;
  std::string eleIDquality_; 
  std::string eleIDqualityLoose_;// only used for Z rejection
  bool ignore_emIDtight_; // possibility to completely ignore EM id

  double eleEtCut_;
  double eleEtaCut_;
  double eleIsoCut_;
  double eled0Cut_;
  double eleECALbadLo_;
  double eleECALbadHi_;
  double eleEtCutLoose_;
  double eleEtaCutLoose_;
  double eleIsoCutLoose_;
  double eled0CutLoose_;
  double eleMvaCut_;
  double dREleJetCrossClean_;
  double muoEtaCut_;
  double muoPtCut_;
  double muoD0Cut_;
  double muoNTkHitsCut_;
  double muoIsoCut_;
  double muoNormChi2_;
  double muoHCalIso_;
  double muoECalIso_;
  double muoVldHits_;
  double muoMtchdStns_;
  double muoDB_;
  double muoDZCut_;
  double muoPxlHits_;
  double muoTkLyrsWthHts_;
  double muoRelIsoTight_;
  double muoPtLoose_;
  double muoEtaLoose_;
  double muoRelIsoLoose_;
  double metCut_;
  double rhoIso;
  
  bool ran_jetloop_;
  bool ran_eleloop_;
  bool ran_muonloop_;
  bool ran_mcloop_;
  bool ran_postloop_;
  bool ran_PV_;
  bool ran_tracks_;
  bool ran_photonTau_;
  bool check_triggers_;
  std::string muoIDquality_;
  bool flavorHistoryTag_;
  bool bookedBDBDisc_;
  double  dREleGeneralTrackMatch_;
  double  magneticField_;
  double  correctFactor_;
  double  maxDist_;
  double  maxDcot_;
  edm::InputTag ebRecHits_;
  edm::InputTag eeRecHits_;
  bool isMCatNLO_;


  // and an ntuple (filling in the methods)
  void fillBeamSpot(const edm::Event&, const edm::EventSetup&);
  void fillJets(const edm::Event&, const edm::EventSetup&, edm::InputTag, std::string);
  void fillBTagInfo(const pat::Jet &jet, const size_t jetindex, std::string ID);
  void fillBTagInfoNew(const pat::Jet &jet, const size_t jetindex, std::string ID);
  void fillBIDParameters(const edm::EventSetup&, std::string);
  void makeBIDInfoHistos(const edm::EventSetup&);
  void fillBIDInfoHistos( const edm::EventSetup&, TH2D *, PerformanceResult::ResultType &, BinningVariables::BinningVariablesType, BinningVariables::BinningVariablesType, const BtagPerformance &);
  void fillOtherJetInfo(const pat::Jet &jet, const size_t jetindex, std::string ID, const edm::Event& iEvent);
  void fillMCJetInfo(const reco::GenJet &jet, const size_t jetindex, std::string ID, bool fillMC);
  void fillMCJetInfo(int empty, const size_t jetindex, std::string ID, bool fillMC);
  void fillLooseJetInfo(const pat::Jet &jet, const size_t jetindex, float jetPt, std::string ID);
  void fillZVeto(const edm::Event &, const edm::EventSetup&, edm::InputTag, std::string);
  void fillMuons(const edm::Event&, const edm::EventSetup&, edm::InputTag, std::string);
  void fillPhotons(const edm::Event&, const edm::EventSetup&, edm::InputTag, std::string);
  void fillTaus(const edm::Event&, const edm::EventSetup&, edm::InputTag, std::string);
  void fillElectrons(const edm::Event&, const edm::EventSetup&, edm::InputTag, std::string);
  void fillMissingET(const edm::Event&, const edm::EventSetup&,edm::InputTag, std::string);
  void fillEventInfo(const edm::Event&, const edm::EventSetup&);
  void fillMCInfo(const edm::Event&, const edm::EventSetup&);
  void fillTriggerData(const edm::Event&);
  void fillSummaryVariables(void);// should only be called after all other functions.
  void fillGeneralTracks(const edm::Event&, const edm::EventSetup&);
  void fillFlavorHistory(const edm::Event&, const edm::EventSetup&);
  // ID functions
  bool looseElectronID(const pat::Electron &, bool = false);
  bool tightElectronID(const pat::Electron &, bool = false);
  bool jetID(const pat::Jet &,const size_t jetindex, std::string, float);
  bool jetIDLoose(const pat::Jet &, float);
  bool muonID(const pat::Muon &);
  bool muonIDLoose(const pat::Muon &);
  bool photonConversionVeto(const pat::Electron &, float &, float &);

  void bookBranches(void);// does all the branching.
  void bookJetBranches(std::string ID, std::string name);// called by bookBranches, makes jet branches.
  void bookBIDInfoBranches(std::string, std::string);// called by bookJetBranches, makes branches for B-ID.
  void bookCaloJetBranches(std::string ID, std::string name);// called by bookBranches, makes jet branches.
  void bookPFJetBranches(std::string ID, std::string name);// called by bookBranches, makes jet branches.
  void bookJPTJetBranches(std::string ID, std::string name);// called by bookBranches, makes jet branches.
  void bookTauBranches(std::string ID, std::string name);
  void bookPhotonBranches(std::string ID, std::string name);
  void bookElectronBranches(std::string ID, std::string name);// called by bookBranches, makes electron branches.
  void bookMuonBranches(std::string ID, std::string name);// called by bookBranches, makes muon branches.
  void bookMETBranches(std::string ID, std::string name); //called by bookBranches , makes MET branches
  void bookCaloMETBranches(std::string ID, std::string name); //called by bookBranches , makes MET branches
  void bookMCBranches(void); // called by bookBranches, makes MC branches.
  void bookGeneralTracksBranches(void); // called by bookBranches, makes generalTracks branches.
  
  float getAEff03(float); //Used to get the effective area of the electron. Because ElectronEffectiveArea.h is outdated. FIXME in future releases.
  
  TTree *mytree_;

  int processId_; int genMyProcId;
  float processPtHat_;
  
  double weight_;
  double topPtReweight;

  int numGeneralTracks;
    std::map< std::string, int > numJet;
  std::map< std::string, int > numLooseBJets;
    std::map< std::string, int > numEle;
    std::map< std::string, int > numLooseEle;
    std::map< std::string, int > numMuo;
  std::map< std::string, int > numLooseMuo;
 
  math::XYZPoint beamSpotPoint_; 
  math::XYZPoint vertexPoint_; 
  
  unsigned int flavorHistory;

  template <class C>
  struct IndexSorter {
    IndexSorter (const C& values, bool decreasing = true) : 
      values_(values), decrease_(decreasing) {}
    std::vector<int> operator() () const {
      std::vector<int> result;
      result.reserve(values_.size());
      for ( size_t i=0; i<values_.size(); ++i )  result.push_back(i);
      sort(result.begin(),result.end(),*this);
      return result;
    }
    bool operator() (int a, int b) {
      if ( decrease_ )
	return values_[a]>values_[b];
      else
	return values_[a]<values_[b];
    }
    const C& values_;
    bool decrease_;
  };

  void cleararrays(void);// used to set everything in the following arrays to zero or unphysical numbers
  void clearjetarrays(std::string);// clearing jet info, used by cleararrays]
  void clearLooseJetarrays(std::string);
  void clearTauArrays(std::string);
  void clearPhotonArrays(std::string);
  void clearelectronarrays(std::string);//clearing electron info, used by cleararrays
  void clearmuonarrays(std::string);//clearing muon info, used by cleararrays
  void clearMetArrays(std::string); //clearing met info
  void clearMCarrays(void); //clearing MC info
  void clearGeneralTracksarrays(void);//clearing generalTracks info, used by cleararrays

  std::vector<float> electronEts; //just used for sorting

  float beamSpotX;
  float beamSpotY;
  float beamSpotZ;

  int numPv;
  float pvX;
  float pvY;
  float pvZ;
  float pvDX;
  float pvDY;
  float pvDZ;
  float pvRho;
  int pvIsFake;
  float pvChi2;
  float pvNdof;

  std::map< std::string, int >nzcandidates; 
  std::map< std::string, std::vector<float> > zcandidatesvector;// stores the Z candidates

  size_t NELECTRONSMAX;// set to array size in constructor
  std::map< std::string, std::vector<float> > electronSortedE;
  std::map< std::string, std::vector<float> > electronSortedEt;
  std::map< std::string, std::vector<float> > electronSortedEta;
  std::map< std::string, std::vector<float> > electronSortedPt; 
  std::map< std::string, std::vector<float> > electronSortedTheta;
  std::map< std::string, std::vector<float> > electronSortedPhi;
  std::map< std::string, std::vector<float> > electronSortedPx;
  std::map< std::string, std::vector<float> > electronSortedPy;
  std::map< std::string, std::vector<float> > electronSortedPz;
  std::map< std::string, std::vector<int> > electronSortedCharge;
  std::map< std::string, std::vector<float> > electronSortedMVA;
  //  std::map< std::string, std::vector<int> > electronSortedIDQuality;
  //std::map< std::string, std::vector<int> > electronSortedIDQualityLoose;
  std::map< std::string, std::vector<float> > electronSortedChargedHadronIso;
  std::map< std::string, std::vector<float> > electronSortedNeutralHadronIso;
  std::map< std::string, std::vector<float> > electronSortedPhotonIso;
  std::map< std::string, std::vector<float> > electronSortedTrackPt;
  std::map< std::string, std::vector<float> > electronSortedTrackEta;
  std::map< std::string, std::vector<float> > electronSortedTrackPhi;
  std::map< std::string, std::vector<float> > electronSortedTrackChi2;
  std::map< std::string, std::vector<float> > electronSortedTrackNDOF;
  std::map< std::string, std::vector<float> > electronSortedTrackD0;
  std::map< std::string, std::vector<float> > electronSortedBeamSpotCorrectedTrackD0;
  std::map< std::string, std::vector<float> > electronSortedDBBeamSpotCorrectedTrackD0;

  //std::map< std::string, std::vector<float> > electronSortedDBInnerTrackD0;

  std::map< std::string, std::vector<float> > electronSortedTrackDz;
  std::map< std::string, std::vector<float> > electronSortedTrackD0PV;
  std::map< std::string, std::vector<float> > electronSortedTrackDZPV;
  std::map< std::string, std::vector<float> > electronSortedVtxZ;
  std::map< std::string, std::vector<float> > electronSortedBeamSpotCorrectedTrackDz;
  std::map< std::string, std::vector<int> > electronSortedIsGsf;
  std::map< std::string, std::vector<float> > electronSortedGsfPx;
  std::map< std::string, std::vector<float> > electronSortedGsfPy;
  std::map< std::string, std::vector<float> > electronSortedGsfPz;
  std::map< std::string, std::vector<float> > electronSortedGsfE;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterEta;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterE;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterPhi;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterSigmaEtaEta;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterE1x5;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterE2x5max;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterE5x5;
  std::map< std::string, std::vector<float> > electronSortedSuperClusterSigmaIEtaIEta;
  std::map< std::string, std::vector<float> > electronSortedTrackIso04;
  std::map< std::string, std::vector<float> > electronSortedECalIso04;
  std::map< std::string, std::vector<float> > electronSortedHCalIso04;
  std::map< std::string, std::vector<float> > electronSortedTrackIso03;
  std::map< std::string, std::vector<float> > electronSortedECalIso03;
  std::map< std::string, std::vector<float> > electronSortedHCalIso03;
  std::map< std::string, std::vector<float> > electronSorteddr04EcalRecHitSumEt;
  std::map< std::string, std::vector<float> > electronSorteddr03EcalRecHitSumEt;
  std::map< std::string, std::vector<float> > electronSortedECalIsoDeposit;
  std::map< std::string, std::vector<float> > electronSortedHCalIsoDeposit;
  std::map< std::string, std::vector<float> > electronSortedCaloIso;
  std::map< std::string, std::vector<float> > electronSortedTriggerMatch;
  std::map< std::string, std::vector<float> > electronSortedJetOverlap;
  std::map< std::string, std::vector<float> > electronSortedComRelIso;
  std::map< std::string, std::vector<float> > electronSortedComRelIsodBeta;
  std::map< std::string, std::vector<float> > electronSortedComRelIsoRho;
  std::map< std::string, std::vector<float> > electronSortedChHadIso;
  std::map< std::string, std::vector<float> > electronSortedNtHadIso;
  std::map< std::string, std::vector<float> > electronSortedGammaIso;
  std::map< std::string, std::vector<float> > electronSortedRhoIso;
  std::map< std::string, std::vector<float> > electronSortedAEff03;
  std::map< std::string, std::vector<int> > electronSortedMissingInnerLayers;
  std::map< std::string, std::vector<float> > electronSortedHoverE;
  std::map< std::string, std::vector<float> > electronSortedDeltaPhiSC;
  std::map< std::string, std::vector<float> > electronSortedDeltaEtaSC;
  std::map< std::string, std::vector<int> > electronSortedIsBarrel;
  std::map< std::string, std::vector<int> > electronSortedPhotonConversionTag;
  std::map< std::string, std::vector<int> > electronSortedPhotonConversionTagCustom;
  std::map< std::string, std::vector<float> > electronSortedPhotonConversionDcot;
  std::map< std::string, std::vector<float> > electronSortedPhotonConversionDist;
  std::map< std::string, std::vector<int> > electronSortedPhotonConversionVeto;
  std::map< std::string, std::vector<float> > electronSortedPhotonConversionDcotCustom;
  std::map< std::string, std::vector<float> > electronSortedPhotonConversionDistCustom;
  //std::map< std::string, std::vector<float> > electronSortedSwissCross;

  std::map< std::string, std::vector<float> > electronSortedImpactTransDist;
  std::map< std::string, std::vector<float> > electronSortedImpactTransError;
  std::map< std::string, std::vector<float> > electronSortedImpactTransSignificance;
  std::map< std::string, std::vector<float> > electronSortedImpact3DDist;
  std::map< std::string, std::vector<float> > electronSortedImpact3DError;
  std::map< std::string, std::vector<float> > electronSortedImpact3DSignificance;

  //  std::map< std::string, std::vector<float> > electronSortedIDResults_;

  std::map< std::string, std::vector<float> > genElectronSortedEt;
  std::map< std::string, std::vector<float> > genElectronSortedEta;
  std::map< std::string, std::vector<float> > genElectronSortedTheta;
  std::map< std::string, std::vector<float> > genElectronSortedPhi;
  std::map< std::string, std::vector<float> > genElectronSortedPx;
  std::map< std::string, std::vector<float> > genElectronSortedPy;
  std::map< std::string, std::vector<float> > genElectronSortedPz;
  std::map< std::string, std::vector<int> > genElectronSortedCharge;
  
  //Information for loose electrons
  std::map< std::string, std::vector<float> > looseElectronSortedEt;
  std::map< std::string, std::vector<float> > looseElectronSortedPt;
  std::map< std::string, std::vector<float> > looseElectronSortedEta;
  std::map< std::string, std::vector<float> > looseElectronSortedMVA;
  std::map< std::string, std::vector<float> > looseElectronSortedRelIso;
  



  // MC Truth
  int nT;
  int nThadronic, nb, nWhadronic;
  int nTleptonic, nWleptonic;
  int VQQBosonAbsId; 

  int NTOPMCINFOSMAX;
  float T_hadronicMCTruthE[20];
  float T_hadronicMCTruthEt[20];
  float T_hadronicMCTruthPx[20];
  float T_hadronicMCTruthPy[20];
  float T_hadronicMCTruthPz[20];
  int T_hadronicMotherIndex[20];

  float T_leptonicMCTruthE[20];
  float T_leptonicMCTruthEt[20];
  float T_leptonicMCTruthPx[20];
  float T_leptonicMCTruthPy[20];
  float T_leptonicMCTruthPz[20];
  int T_leptonicMotherIndex[20];

  float bMCTruthE[20];
  float bMCTruthEt[20];
  float bMCTruthPx[20];
  float bMCTruthPy[20];
  float bMCTruthPz[20];
  int bMCTruthMother[20];

  float W_hadronicMCTruthE[20];
  float W_hadronicMCTruthEt[20];
  float W_hadronicMCTruthPx[20];
  float W_hadronicMCTruthPy[20];
  float W_hadronicMCTruthPz[20];
  int W_hadronicMCTruthPID[20];
  int W_hadronicMCTruthMother[20];

  float W_leptonicMCTruthE[20];
  float W_leptonicMCTruthEt[20];
  float W_leptonicMCTruthPx[20];
  float W_leptonicMCTruthPy[20];
  float W_leptonicMCTruthPz[20];
  int W_leptonicMCTruthPID[20];
  int W_leptonicMCTruthMother[20];

  int isElePlusJets;
  
  //  float remainingEnergy[20];

  std::map< std::string, double > metEt;
  std::map< std::string, double > metEtRaw;  
  std::map< std::string, double > metPhi;
  std::map< std::string, double > metPt;
  std::map< std::string, double > metPx; 
  std::map< std::string, double > metPy;
  std::map< std::string, float > metSignificance;
  std::map< std::string, float > metScalarEt;
  std::map< std::string, float > metEtUncorrected;
  std::map< std::string, float > metPhiUncorrected;
  std::map< std::string, float > metMaxEtEM;
  std::map< std::string, float > metMaxEtHad;
  std::map< std::string, float > metEtFracHad;
  std::map< std::string, float > metEtFracEM;
  std::map< std::string, float > metHadEtHB;
  std::map< std::string, float > metHadEtHO;
  std::map< std::string, float > metHadEtHE;
  std::map< std::string, float > metEmEtEE;
  std::map< std::string, float > metEmEtEB;
  std::map< std::string, float > metEmEtHF;
  std::map< std::string, float > metHadEtHF;
  std::map< std::string, float > genMetEt; 
  std::map< std::string, float > genMetPhi;
  std::map< std::string, float > genMetPt; 
  std::map< std::string, float > genMetPx; 
  std::map< std::string, float > genMetPy; 

  edm::LumiReWeighting LumiWeightsA;
  edm::LumiReWeighting LumiWeightsB;
  edm::LumiReWeighting LumiWeightsC;

  int numVert;
  double pileUpWeightA;
  double pileUpWeightB;
  double pileUpWeightC;

  float mhtPx;
  float mhtPy;
  float mhtPt;
  float mhtPhi;
  float mhtSumEt;
  float mhtSignif;
  
  size_t NMUONSMAX;// max array size, set in constructor.
  std::vector<float> muonEts;
  std::map< std::string, std::vector<float> > muonSortedE;
  std::map< std::string, std::vector<float> > muonSortedEt;
  std::map< std::string, std::vector<float> > muonSortedPt;
  std::map< std::string, std::vector<float> > muonSortedEta;
  std::map< std::string, std::vector<float> > muonSortedTheta;
  std::map< std::string, std::vector<float> > muonSortedPhi;
  std::map< std::string, std::vector<float> > muonSortedPx;
  std::map< std::string, std::vector<float> > muonSortedPy;
  std::map< std::string, std::vector<float> > muonSortedPz;
  std::map< std::string, std::vector<int> > muonSortedCharge;
  std::map< std::string, std::vector<float> > muonSortedGlobalID;
  std::map< std::string, std::vector<float> > muonSortedTrackID;
  std::map< std::string, std::vector<float> > muonSortedChi2;
  std::map< std::string, std::vector<float> > muonSortedD0;
  std::map< std::string, std::vector<float> > muonSortedDBBeamSpotCorrectedTrackD0;

  std::map< std::string, std::vector<float> > muonSortedDBInnerTrackD0;

  std::map< std::string, std::vector<float> > muonSortedBeamSpotCorrectedD0;
  std::map< std::string, std::vector<int> > muonSortedTrackNHits;
  std::map< std::string, std::vector<int> > muonSortedValidHitsGlobal;
  std::map< std::string, std::vector<float> > muonSortedNDOF; //n_d.o.f
  //Extra muon variables used for ID and stuff
  std::map< std::string, std::vector<int> > muonSortedTkLysWithMeasurements;
  std::map< std::string, std::vector<float> > muonSortedGlbTkNormChi2;
  std::map< std::string, std::vector<float> > muonSortedDBPV;
  std::map< std::string, std::vector<float> > muonSortedDZPV;
  std::map< std::string, std::vector<int> > muonSortedVldPixHits;
  std::map< std::string, std::vector<int> > muonSortedMatchedStations;

  //Vertex location information. For dZ cuts.
  std::map< std::string, std::vector<float> > muonSortedVertX;
  std::map< std::string, std::vector<float> > muonSortedVertY;
  std::map< std::string, std::vector<float> > muonSortedVertZ;
  

  std::map< std::string, std::vector<float> > muonSortedChargedHadronIso;
  std::map< std::string, std::vector<float> > muonSortedNeutralHadronIso;
  std::map< std::string, std::vector<float> > muonSortedPhotonIso;

  std::map< std::string, std::vector<float> > muonSortedTrackIso;
  std::map< std::string, std::vector<float> > muonSortedECalIso;
  std::map< std::string, std::vector<float> > muonSortedHCalIso;
  std::map< std::string, std::vector<float> > muonSortedComRelIso;
  std::map< std::string, std::vector<float> > muonSortedComRelIsodBeta;
  std::map< std::string, std::vector<int> > muonSortedIsPFMuon;

  std::map< std::string, std::vector<int> > muonSortedNumChambers;
  std::map< std::string, std::vector<int> > muonSortedNumMatches;

  std::map< std::string, std::vector<float> > genMuonSortedEt;
  std::map< std::string, std::vector<float> > genMuonSortedEta;
  std::map< std::string, std::vector<float> > genMuonSortedTheta;
  std::map< std::string, std::vector<float> > genMuonSortedPhi;
  std::map< std::string, std::vector<float> > genMuonSortedPx;
  std::map< std::string, std::vector<float> > genMuonSortedPy;
  std::map< std::string, std::vector<float> > genMuonSortedPz;
  std::map< std::string, std::vector<int> > genMuonSortedCharge;


  //Loose muon information

  std::map< std::string, std::vector<float> > looseMuonSortedEt;
  std::map< std::string, std::vector<float> > looseMuonSortedPt;
  std::map< std::string, std::vector<float> > looseMuonSortedEta;
  std::map< std::string, std::vector<float> > looseMuonSortedRelIso;
  std::map< std::string, std::vector<float> > looseMuonSortedisGlb;
  std::map< std::string, std::vector<float> > looseMuonSortedisTrk;

  size_t NJETSMAX; // max number of jets, set in constructor;

//JEC to be initialised once per collection.
    FactorizedJetCorrector *jecCalo;
    FactorizedJetCorrector *jecPF;
    FactorizedJetCorrector *jecJPT;
    JetCorrectionUncertainty *jecCaloUncertainty;
    JetCorrectionUncertainty *jecPFUncertainty;
    JetCorrectionUncertainty *jecJPTUncertainty;

    std::vector<float> correctedJetEts;
    std::map< std::string, std::vector<double> > jetSortedE;
    std::map< std::string, std::vector<double> > jetSortedEt;
    std::map< std::string, std::vector<double> > jetSortedPt;
    std::map< std::string, std::vector<double> > jetSortedPtRaw;
    std::map< std::string, std::vector<double> > jetSortedUnCorEt;
    std::map< std::string, std::vector<double> > jetSortedUnCorPt;
    std::map< std::string, std::vector<double> > jetSortedEta;
    std::map< std::string, std::vector<double> > jetSortedTheta;
    std::map< std::string, std::vector<double> > jetSortedPhi;
    std::map< std::string, std::vector<double> > jetSortedPx;
    std::map< std::string, std::vector<double> > jetSortedPy;
    std::map< std::string, std::vector<double> > jetSortedPz;
    std::map< std::string, std::vector<double> > jetSortedClosestLepton;
    std::map< std::string, std::vector<int> >   jetSortedNtracksInJet;
    std::map< std::string, std::vector<float> > jetSortedJetCharge;
    std::map< std::string, std::vector<float> > jetSortedfHPD;
    std::map< std::string, std::vector<float> > jetSortedCorrFactor;// only include full JES corrections for now.
    std::map< std::string, std::vector<float> > jetSortedCorrResidual;
    std::map< std::string, std::vector<float> > jetSortedL2L3ResErr;
    std::map< std::string, std::vector<float> > jetSortedCorrErrLow;// JES uncertainty
    std::map< std::string, std::vector<float> > jetSortedCorrErrHi;
    std::map< std::string, std::vector<float> > jetSortedN90Hits;
    std::map< std::string, std::vector<float> > jetSortedTriggered;
    std::map< std::string, std::vector<float> > jetSortedSVPT;
    std::map< std::string, std::vector<float> > jetSortedSVL2D;
    std::map< std::string, std::vector<float> > jetSortedSVL2Dxy;
    std::map< std::string, std::vector<float> > jetSortedSVL2DxyErr;
    std::map< std::string, std::vector<float> > jetSortedSVL2DxySig;
    std::map< std::string, std::vector<float> > jetSortedSVL3D;
    std::map< std::string, std::vector<float> > jetSortedSVL3DErr;
    std::map< std::string, std::vector<float> > jetSortedSVL3DSig;
    std::map< std::string, std::vector<float> > jetSortedSVMass;
    std::map< std::string, std::vector<int> >   jetSortedSVNtracks;
    std::map< std::string, std::vector<float> >  jetSortedSVX;
    std::map< std::string, std::vector<float> >  jetSortedSVY;
    std::map< std::string, std::vector<float> >  jetSortedSVZ;
    std::map< std::string, std::vector<float> >  jetSortedSVDX;
    std::map< std::string, std::vector<float> >  jetSortedSVDY;
    std::map< std::string, std::vector<float> > jetSortedSVDZ;
    std::map< std::string, std::vector<int> > jetSortedNConstituents;
    std::map< std::string, std::vector<float> > jetSortedBDiscriminator;

//Calo Jet
    std::map< std::string, std::vector<float> > jetSortedEMEnergyInEB;
    std::map< std::string, std::vector<float> > jetSortedEMEnergyInEE;
    std::map< std::string, std::vector<float> > jetSortedEMEnergyFraction;
    std::map< std::string, std::vector<float> > jetSortedEMEnergyInHF;
    std::map< std::string, std::vector<float> > jetSortedHadEnergyInHB;
    std::map< std::string, std::vector<float> > jetSortedHadEnergyInHE;
    std::map< std::string, std::vector<float> > jetSortedHadEnergyInHF;
    std::map< std::string, std::vector<float> > jetSortedHadEnergyInHO;
    std::map< std::string, std::vector<float> > jetSortedN60;
    std::map< std::string, std::vector<float> > jetSortedN90;
//PF Specific
    std::map< std::string, std::vector<float> > jetSortedNeutralEmEnergy;
    std::map< std::string, std::vector<float> > jetSortedMuEnergy;
    std::map< std::string, std::vector<float> > jetSortedMuEnergyFraction;
    std::map< std::string, std::vector<int> > jetSortedChargedMultiplicity;
    std::map< std::string, std::vector<float> > jetSortedNeutralHadEnergy;

    std::map< std::string, std::vector<int> > jetSortedNeutralMultiplicity;
    std::map< std::string, std::vector<float> > jetSortedChargedHadronEnergyFraction;
    std::map< std::string, std::vector<float> > jetSortedNeutralHadronEnergyFraction;
    std::map< std::string, std::vector<float> > jetSortedChargedEmEnergyFraction;
    std::map< std::string, std::vector<float> > jetSortedNeutralEmEnergyFraction;
    std::map< std::string, std::vector<float> > jetSortedChargedHadronEnergyFractionCorr;
    std::map< std::string, std::vector<float> > jetSortedNeutralHadronEnergyFractionCorr;
    std::map< std::string, std::vector<float> > jetSortedChargedEmEnergyFractionCorr;
    std::map< std::string, std::vector<float> > jetSortedNeutralEmEnergyFractionCorr;

    // more detailed BID info for a few algorithms.
    std::map< std::string, std::vector<float> > jetSortedBtagSoftMuonPtRel;
    std::map< std::string, std::vector<float> > jetSortedBtagSoftMuonQuality;
    std::map< std::string,std::vector<float> > jetSortedBtagDiscriminants_; // stores the data output
    std::map< std::string,std::vector<float> > jetSortedBIDParams_; // stores the parameter (db) output
    std::map< std::string,float > bidParamsDiscCut_;

    std::map< std::string,std::vector<float> > genJetSortedEt;
    std::map< std::string,std::vector<float> > genJetSortedPt;
    std::map< std::string,std::vector<float> > genJetSortedEta;
    std::map< std::string,std::vector<float> > genJetSortedTheta;
    std::map< std::string,std::vector<float> > genJetSortedPhi;
    std::map< std::string,std::vector<float> > genJetSortedPx;
    std::map< std::string,std::vector<float> > genJetSortedPy;
    std::map< std::string,std::vector<float> > genJetSortedPz;
    std::map< std::string, std::vector<int> > jetSortedPID;
    std::map< std::string, std::vector<int> > genJetSortedPID;
    std::map< std::string,std::vector<float> > genJetSortedClosestB;
    std::map< std::string,std::vector<float> > genJetSortedClosestC;
    std::map< std::string,std::vector<float> > genJetSortedBtag;

  //Loose jet info
  std::map< std::string, std::vector<float> > jetLooseSortedPt;
  std::map< std::string, std::vector<float> > jetLooseSortedEt;
  std::map< std::string, std::vector<float> > jetLooseSortedEta;                    
  std::map< std::string, std::vector<float> > jetLooseSortedBDisc;                    




  //generalTracks are used to subtract photon conversion background
  size_t NTRACKSMAX;
  float generalTracksPt[1000];
  float generalTracksEta[1000];
  float generalTracksTheta[1000];
  float generalTracksPhi[1000];
  float generalTracksBeamSpotCorrectedD0[1000];
  int   generalTracksCharge[1000];
  

  // gen particle vars
  int NGENPARMAX;
  int nGenPar;
  float genParEta[50];
  float genParPhi[50];
  float genParE[50];
  float genParPt[50];
  int genParId[50];
  int genParCharge[50];
//PDF info
    float genPDFScale;
    float genPDFx1;
    float genPDFx2;
    int genPDFf1;
    int genPDFf2;
//CTEQ_6.6 general purpose
    float genCTEQ66_Weight[44];
//MRST98 NLO
    float genMRST2006nnlo_Weight[31];

  // basic 4-vectors for photons,taus as we're not interested in them.
  int NTAUSMAX; // set to array size in constructor
  int NPHOTONSMAX; // set to array size in constructor
  std::map< std::string, int > ntaus;
  std::map< std::string, int > nphotons;
  std::map< std::string, std::vector<float> > photon_e;
  std::map< std::string, std::vector<float> > photon_phi;
  std::map< std::string, std::vector<float> > photon_eta;
  std::map< std::string, std::vector<float> > photon_pt;
  std::map< std::string, std::vector<float> > tau_e;
  std::map< std::string, std::vector<float> > tau_phi;
  std::map< std::string, std::vector<float> > tau_eta;
  std::map< std::string, std::vector<float> > tau_pt;
  
  std::vector<int> triggerRes;
  std::vector<int> HLT_fakeTriggerValues;
  int nTriggerBits;
  int TriggerBits[700];
  
  float topo_sphericity;
  float topo_aplanarity;
  float topo_sphericitye;
  float topo_aplanaritye;
  float topo_oblateness;
  float topo_sqrts;
  float topo_sqrtse;
  float topo_ht3;
  float topo_hte;
  float topo_ht;
  
  int evtRun;
  int evtnum;
  float evtlumiblock;
  int eventCount;
  
  int bTags;
  int softTags;

  //Debugging variables used to see how many electrons are found at each cut.
  int eleDebugNumberTotal;
  int eleDebugNumbermvaID;
  int eleDebugNumberEt;
  int eleDebugNumberEta;
  int eleDebugNumberCrack;
  int eleDebugNumberIso;
  int eleDebugNumberD0;
  int eleDebugNumberConV;


  //Debugging variables for muons
  int vanillaMuons;
  int globalPFMuons;
  int ptMuons;
  int validHitsMuons;
  int chi2Muons;
  int tkHitsMuons;
  int dbMuons;
  int dzMuons;
  int pixelHitsMuons;
  int trackerLayersMuons;
  int mvaTrig;
  int mvaAsFunc;

  //cut flow related variables.
  int runCutFlow_;
  bool doCutFlow;
  bool doJERSmear_;
  bool fillAll_;
  bool processingLoose_;
  
  //Sets whether the sample is ttbar or not. Default is false. This affects top pt reweighting of the sample.
  bool isttbar_;

};

namespace LHAPDF {
    enum SetType {
	EVOLVE = 0, LHPDF = 0, 
	INTERPOLATE = 1, LHGRID = 1
    };
    enum Verbosity { SILENT=0, LOWKEY=1, DEFAULT=2 };
    void setVerbosity(Verbosity noiselevel);
    void initPDFSet(int nset, const std::string& filename, int member=0);
    void initPDFSet(const std::string& name, LHAPDF::SetType type, int member=0);
    void initPDFSet(int nset, const std::string& name, LHAPDF::SetType type, int member=0);
    int numberPDF(int nset);
    void usePDFMember(int nset, int member);
    void usePDFMember(int member);
    double xfx(int nset, double x, double Q, int fl);
    double xfx(double x, double Q, int fl);
    double getXmin(int nset, int member);
    double getXmax(int nset, int member);
    double getQ2min(int nset, int member);
    double getQ2max(int nset, int member);
    void extrapolate(bool extrapolate=true);
}
