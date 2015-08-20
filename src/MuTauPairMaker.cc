#include "MYANALYZER/MYANALYZER/interface/MuTauPairMaker.h"

mtp::MuTauPairMaker::MuTauPairMaker(const edm::ParameterSet &iConfig, edm::ConsumesCollector && iC):
  // srcToken_(iC.consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("muTauPairs")))
  srcToken_(iC.consumes<pat::CompositeCandidate>(iConfig.getParameter<edm::InputTag>("muTauPairs")))
{
  std::string arbitration = iConfig.getParameter<std::string>("arbitration");
  if(arbitration == "None") {
    arbitration_ = None;
  } else throw cms::Exception("Configuration") << "TagProbePairMakerOnTheFly: the only currently "
	                                       << "allowed values for 'arbitration' are "
                                               << "'None'\n";
}

// THIS FUNCTION RETURN A VECTOR OF TYPE mtp::MuTauPairs
mtp::MuTauPairs
mtp::MuTauPairMaker::run(const edm::Event &iEvent) const
{
  // declare output
  mtp::MuTauPairs pairs;

  // read from event
  // edm::Handle<reco::CandidateView> src;
  // iEvent.getByToken(srcToken_, src);

  // read from event
  edm::Handle<pat::CompositeCandidateCollection> src;
  iEvent.getByToken(srcToken_, src);



  /*
  // convert
  for(reco::CandidateView::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) { 
    const reco::Candidate & mother = *it;
    if (mother.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << "Tag&Probe pair with " << mother.numberOfDaughters() << " daughters\n";
    pairs.push_back( mtp::MuTauPair(mother.daughter(0)->masterClone(), mother.daughter(1)->masterClone(),src->refAt(it - src->begin()), mother.mass() ));  
 
  }
  */

  // convert
  for(pat::CompositeCandidateCollection::const_iterator it = src->begin(), ed = src->end(); it != ed; ++it) { 
      const pat::CompositeCandidate & mother = *it;
      // const pat::CompositeCandidate* mother = &(it);
      // if (mother->numberOfDaughters() != 2) throw cms::Exception("CorruptData") << "Tag&Probe pair with " << mother->numberOfDaughters() << " daughters\n";
      // pairs.push_back( mtp::MuTauPair(mother->daughter(0)->masterClone(), mother->daughter(1)->masterClone(),src->refAt(it - src->begin()), mother->mass() ));  
 
      if (mother.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << "Tag&Probe pair with " << mother.numberOfDaughters() << " daughters\n";                     
      pairs.push_back( mtp::MuTauPair(mother.daughter(0)->masterClone(), mother.daughter(1)->masterClone(),src->refAt(it - src->begin()), mother.mass() ));    

  }


  if ((arbitration_ != None) && (pairs.size() > 1)) {
        // might need to clean up
        arbitrate(pairs);
     }

  // return
  return pairs;

}

void 
mtp::MuTauPairMaker::arbitrate(MuTauPairs &pairs) const
{

  std::cout<< "DO NOTHING AT THE MOMENT"<< std::endl;


}


