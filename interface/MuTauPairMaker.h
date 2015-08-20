#ifndef MuTauPairMaker_h
#define MuTauPairMaker_h


#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "TRandom2.h"


namespace mtp {

  /// a simple struct to hold muon and tau candidate pairs
  struct MuTauPair{
    // reco::CandidateBaseRef cand1, cand2, pair;
    pat::CompositeCandidateRef cand1, cand2, pair;
    float mass;
    MuTauPair() {}
    // MuTauPair(const reco::CandidateBaseRef &c1, const reco::CandidateBaseRef &c2, const reco::CandidateBaseRef &p, float m) : cand1(c1), cand2(c2), pair(p), mass(m) {} 
    MuTauPair(const pat::CompositeCandidateRef &c1, const pat::CompositeCandidateRef &c2, const pat::CompositeCandidateRef &p, float m) : cand1(c1), cand2(c2), pair(p), mass(m) {} 

  };
  typedef std::vector<MuTauPair> MuTauPairs;

  class MuTauPairMaker{
      public: 
          MuTauPairMaker( const edm::ParameterSet &iConfig, edm::ConsumesCollector && iC );
          ~MuTauPairMaker() {}
          /// fill in the MuTau pairs for this event
          MuTauPairs run(const edm::Event &iEvent) const ;

      private: 
	  // edm::EDGetTokenT<reco::CandidateView> srcToken_;
	  edm::EDGetTokenT<pat::CompositeCandidate> srcToken_;  
          enum Arbitration {None};
          Arbitration arbitration_;
          void arbitrate(MuTauPairs &pairs) const;
 
  };



}

#endif
