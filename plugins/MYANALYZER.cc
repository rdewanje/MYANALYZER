// -*- C++ -*-
//
// Package:    MYANALYZER/MYANALYZER
// Class:      MYANALYZER
// 
/**\class MYANALYZER MYANALYZER.cc MYANALYZER/MYANALYZER/plugins/MYANALYZER.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  ram dewanjee
//         Created:  Mon, 15 Jun 2015 19:25:36 GMT
//
//


// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <TNtuple.h>
#include <TTree.h>
#include <iterator>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// MY HEADERS
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/TauReco/interface/PFTau.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

// FOR MUTAU BASELINE
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"    // NEEDED IF RUNNING ON MINIAOD (since it has pat::Muon in place of reco::Muon )
#include "DataFormats/PatCandidates/interface/Tau.h"     // NEEDED IF RUNNING ON MINIAOD (since it has pat::Tau in place of reco::Tau )
#include "DataFormats/PatCandidates/interface/Electron.h" 
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// #include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimatorRun2Phys14NonTrig.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class declaration
//

class MYANALYZER : public edm::EDAnalyzer {
   public:
      explicit MYANALYZER(const edm::ParameterSet&);
      ~MYANALYZER();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      std::string TauLabel;    
      std::string MuonLabel;
      // std::string VertexLabel; // MUTAU BASELINE SELECTIONS
      TString theFileName;
      TTree *myTree;      
      edm::InputTag srcTauCandidates_;
      typedef std::vector<edm::InputTag> vInputTag;
      vInputTag srcDiscriminators_;

      // BRANCH VARIBLES
     

     // EVENT BRANCHES
     std::vector<Int_t> _evt;      
     std::vector<Int_t> _run;      
     std::vector<Int_t> _lumi;      


      // TAU BRANCHES
      std::vector<Float_t> _tau_pt;
      std::vector<Float_t> _tau_px;
      std::vector<Float_t> _tau_py;
      std::vector<Float_t> _tau_pz;
      std::vector<Float_t> _tau_e;
      std::vector<Float_t> _tau_eta;
      std::vector<Float_t> _tau_phi;
      std::vector<Bool_t> _tau_passed;
      std::vector<Float_t> _tau_chg_iso;
      std::vector<Float_t> _tau_neutral_iso;
      std::vector<Float_t> _tau_pucorr_iso;

      // MUON BRANCHES
      std::vector<Float_t> _mu_pt;
      std::vector<Float_t> _mu_eta;
      std::vector<Float_t> _mu_phi;
      std::vector<Float_t> _mu_px;
      std::vector<Float_t> _mu_py;
      std::vector<Float_t> _mu_pz;
      std::vector<Float_t> _mu_e;
      std::vector<Bool_t> _mu_isglobal;
      std::vector<Bool_t> _mu_istracker;
      std::vector<Bool_t> _mu_iso;
      std::vector<Float_t> _mvis;

     // MUTAU BASELINE QUANTITIES
     std::vector<Bool_t> _mu_isPFMuon ;
     std::vector<Bool_t> _mu_MediumId ;
     std::vector<Float_t> _mu_dxy;
     std::vector<Float_t> _mu_dz;
     std::vector<Float_t> _mu_chg;      
     std::vector<Float_t> _tau_chg;      
     std::vector<Float_t> _tau_dz;      
     std::vector<Float_t> _tau_ID_NewDMs;
     std::vector<Float_t> _tau_ID_Elec_VLooseMVA5;
     std::vector<Float_t> _tau_ID_MuonTight3;
     std::vector<Float_t> _tau_ID_DBCorrRaw3Hits;



};

/*
void MYANALYZER::Intialize(){
  _tau_pt.clear();
  _mu_pt.clear();

}
*/


//
// constants, enums and typedefs
//

//
// static data member definitions
//


bool isMediumMuon(const reco::Muon & recoMu) 
{
  bool goodGlob = recoMu.isGlobalMuon() && 
    recoMu.globalTrack()->normalizedChi2() < 3 && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.8 && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}





//
// constructors and destructor
//
MYANALYZER::MYANALYZER(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  TauLabel           = iConfig.getUntrackedParameter<std::string>("TauCollection");
  MuonLabel          = iConfig.getUntrackedParameter<std::string>("MuonCollection");
  // VertexLabel        = iConfig.getUntrackedParameter<std::string>("VertexCollection");
  theFileName        = iConfig.getUntrackedParameter<std::string>("fileName");
  srcTauCandidates_  = iConfig.getParameter<edm::InputTag>("srcTauCandidates");
  srcDiscriminators_ = iConfig.getParameter<vInputTag>("srcDiscriminators");

 
  //  Initialize(); 
  _evt.clear();
  _run.clear();
  _lumi.clear();
  _tau_pt.clear();
  _tau_px.clear();
  _tau_py.clear();
  _tau_pz.clear();
  _tau_e.clear();
  _tau_eta.clear();
  _tau_phi.clear();
  _tau_passed.clear();
  _tau_chg_iso.clear();
  _tau_neutral_iso.clear();
  _tau_pucorr_iso.clear();
  _mu_pt.clear();
  _mu_px.clear();
  _mu_py.clear();
  _mu_pz.clear();
  _mu_e.clear();
  _mu_eta.clear();
  _mu_phi.clear();
  _mu_isglobal.clear();
  _mu_istracker.clear();
  _mu_iso.clear();
  _mvis.clear(); // NOT USED RIGHT NOW !

  // MUTAU BASELINE QUANTITIES
  _mu_MediumId.clear();  
  _mu_dxy.clear();
  _mu_dz.clear();
  _mu_chg.clear();
  _tau_chg.clear();
  _tau_dz.clear();
  _mu_isPFMuon.clear(); 
  _tau_ID_NewDMs.clear();
  _tau_ID_Elec_VLooseMVA5.clear();
  _tau_ID_MuonTight3.clear();
  _tau_ID_DBCorrRaw3Hits.clear();
}


MYANALYZER::~MYANALYZER()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}





// ------------ method called once each job just before starting event loop  ------------
void 
MYANALYZER::beginJob()
{
  edm::Service<TFileService> fs;
  myTree = fs->make<TTree>("MuTauTree","MuTauTree");

  myTree->Branch("Event",&_evt);
  myTree->Branch("Run",&_run);
  myTree->Branch("Lumi",&_lumi);
  myTree->Branch("TauPt",&_tau_pt);
  myTree->Branch("TauPx",&_tau_px);
  myTree->Branch("TauPy",&_tau_py);
  myTree->Branch("TauPz",&_tau_pz);
  myTree->Branch("TauE",&_tau_e);
  myTree->Branch("TauEta",&_tau_eta);
  myTree->Branch("TauPhi",&_tau_phi);
  myTree->Branch("TauPassed",&_tau_passed);
  myTree->Branch("TauChgIso",&_tau_chg_iso);
  myTree->Branch("TauNeutralIso",&_tau_neutral_iso);
  myTree->Branch("TauPUCorrIso",&_tau_pucorr_iso);
  myTree->Branch("MuPt",&_mu_pt);
  myTree->Branch("MuPx",&_mu_px);
  myTree->Branch("MuPy",&_mu_py);
  myTree->Branch("MuPz",&_mu_pz);
  myTree->Branch("MuE",&_mu_e);
  myTree->Branch("MuEta",&_mu_eta);
  myTree->Branch("MuPhi",&_mu_phi);
  myTree->Branch("isMuGlobal",&_mu_isglobal);
  myTree->Branch("isMuTracker",&_mu_istracker);
  myTree->Branch("MuIso",&_mu_iso);
  myTree->Branch("Mvis",&_mvis);
  myTree->Branch("MuonMediumId",&_mu_MediumId);
  myTree->Branch("Mu_dxy",&_mu_dxy);
  myTree->Branch("Mu_dz",&_mu_dz);
  myTree->Branch("Mu_chg",&_mu_chg);
  myTree->Branch("Tau_chg",&_tau_chg);
  myTree->Branch("Tau_dz",&_tau_dz);
  myTree->Branch("isMuPF",&_mu_isPFMuon);
  myTree->Branch("TauID_NewDMs",&_tau_ID_NewDMs);
  myTree->Branch("TauID_Elec_VLooseMVA5",&_tau_ID_Elec_VLooseMVA5);
  myTree->Branch("TauID_MuonTight3",&_tau_ID_MuonTight3);
  myTree->Branch("TauID_DBCorrRaw3Hits",&_tau_ID_DBCorrRaw3Hits);

}




//
// member functions
//

// ------------ method called for each event  ------------
void
MYANALYZER::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;
   using namespace pat;

   // ******************* EVENT ***********************
      _evt.push_back(iEvent.id().event());
      _run.push_back(iEvent.run());
      _lumi.push_back(iEvent.luminosityBlock());
   // ************************************************



   // *********************** VERTICES ******************************
      // edm::Handle<vector<reco::Vertex> >  vertices;

      edm::Handle<edm::View<reco::Vertex>>  vertices;
      iEvent.getByLabel("offlineSlimmedPrimaryVertices",vertices); // FOR MINIAOD
      edm::View<reco::Vertex>::const_iterator iVertex = vertices->begin(); 
      // for(edm::View<reco::Vertex>::const_iterator iVertex = vertices->begin(); iVertex != vertices->end();iVertex++){
      //     if(std::distance(vertices->begin(), iVertex) == 0  ){ std::cout<< "vertices->size(): " << vertices->size() << std::endl;}
      // }     
   // ***************************************************************


   // ******************** ELECTRONS ******************************
      // Get the electron ID data from the event stream.
      // Note: this implies that the VID ID modules have been run upstream.
      // If you need more info, check with the EGM group.
      edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
      edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
      // iEvent.getByToken("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80",medium_id_decisions);
      // iEvent.getByToken("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90",tight_id_decisions);
      iEvent.getByLabel("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80",medium_id_decisions);
      iEvent.getByLabel("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90",tight_id_decisions);

      // Get MVA values and categories (optional)
      edm::Handle<edm::ValueMap<float> > mvaValues;
      edm::Handle<edm::ValueMap<int> > mvaCategories;
      // iEvent.getByToken("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues",mvaValues);
      // iEvent.getByToken("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories",mvaCategories); 

      iEvent.getByLabel("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues",mvaValues);
      iEvent.getByLabel("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories",mvaCategories); 




      // FOR MINIAODSIM                                                                                            
      edm::Handle<pat::ElectronCollection>ElectronHandle;
      iEvent.getByLabel("slimmedElectrons",ElectronHandle);
      const pat::ElectronCollection*  Electrons = ElectronHandle.product();
      for(pat::ElectronCollection::const_iterator iElectron = Electrons->begin(); iElectron != Electrons->end();iElectron++){
          std::cout<< "ELECTRON PT: "<< iElectron->pt() << std::endl;
          std::cout<< "ELECTRON ETA: "<< iElectron->eta() << std::endl;
          std::cout<< "ELECTRON DXY: "<< iElectron->gsfTrack()->dxy(iVertex->position()) << std::endl;
          std::cout<< "ELECTRON DZ: "<< iElectron->gsfTrack()->dz(iVertex->position()) << std::endl;  
	  std::cout<< "ELECTRON PASS CONVERSION VETO: "<< iElectron->passConversionVeto() << std::endl;
	  std::cout<< "ELECTRON MISSING INNER HITS: "<< iElectron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) << std::endl;
          std::cout<< "ELECTRON ISOLATION: "<<  (iElectron->pfIsolationVariables().sumChargedHadronPt + max(iElectron->pfIsolationVariables().sumNeutralHadronEt + iElectron->pfIsolationVariables().sumPhotonEt - 0.5 * iElectron->pfIsolationVariables().sumPUPt, 0.0)) / iElectron->pt()<< std::endl;
      }


      edm::Handle<edm::View<reco::GsfElectron> > electrons; 
      iEvent.getByLabel("slimmedElectrons",electrons);
      // iEvent.getByToken("slimmedElectrons",electrons);
      

      for(size_t i = 0; i < electrons->size(); ++i){
          const auto el = electrons->ptrAt(i);
          bool isPassMedium = (*medium_id_decisions)[el];
	  bool isPassTight  = (*tight_id_decisions)[el]; 
          std::cout<< "ELECTRON SC ETA: "<< el->superCluster()->eta() << std::endl;
          std::cout<< "ELECTRON PASS MEDIUM: "<< (int)isPassMedium << std::endl;    
          std::cout<< "ELECTRON PASS TIGHT: "<< (int)isPassTight << std::endl;    
	  std::cout<< "ELECTRON MVA VALUE: "<< (*mvaValues)[el] << std::endl;        
	  std::cout<< "ELECTRON MVA CATS: "<< (*mvaCategories)[el]<< std::endl;        
      }



   // *************************************************************






   // ************************** MUONS *******************************
   /*
   // FOR AODSIM
   edm::Handle<edm::View<reco::Muon>>MuonHandle;
   iEvent.getByLabel("muons",MuonHandle);  
   const edm::View<reco::Muon>*  Muons = MuonHandle.product();

   
   for(edm::View<reco::Muon>::const_iterator iMuon = Muons->begin(); iMuon != Muons->end();iMuon++){
     // if( std::fabs(iMuon->eta()) < 2.1 && iMuon->isGlobalMuon() && iMuon->isTrackerMuon() && ( iMuon->pfIsolationR03().sumChargedHadronPt ) < (0.3*iMuon->pt()) ){
     // std::cout<< "Passed Muon pt(): " << iMuon->pt() << std::endl;
         _mu_pt.push_back( iMuon->pt() );
         _mu_px.push_back( iMuon->px() );
         _mu_py.push_back( iMuon->py() );
         _mu_pz.push_back( iMuon->pz() );
         _mu_e.push_back( iMuon->energy() );
         _mu_eta.push_back( iMuon->eta() );
         _mu_phi.push_back( iMuon->phi() );
         _mu_isglobal.push_back( iMuon->isGlobalMuon() );
         _mu_istracker.push_back( iMuon->isTrackerMuon() );
         _mu_iso.push_back(iMuon->pfIsolationR03().sumChargedHadronPt );
         _mu_MediumId.push_back(isMediumMuon( *iMuon  ));    

    //    }
   }
  */




   // FOR MINIAODSIM
   edm::Handle<pat::MuonCollection>MuonHandle;
   iEvent.getByLabel("slimmedMuons",MuonHandle);  
   const pat::MuonCollection*  Muons = MuonHandle.product();
   
   for(pat::MuonCollection::const_iterator iMuon = Muons->begin(); iMuon != Muons->end();iMuon++){
     // if( std::fabs(iMuon->eta()) < 2.1 && iMuon->isGlobalMuon() && iMuon->isTrackerMuon() && ( iMuon->pfIsolationR03().sumChargedHadronPt ) < (0.3*iMuon->pt()) ){
     // std::cout<< "Passed Muon pt(): " << iMuon->pt() << std::endl;
         _mu_pt.push_back( iMuon->pt() );
         _mu_px.push_back( iMuon->px() );
         _mu_py.push_back( iMuon->py() );
         _mu_pz.push_back( iMuon->pz() );
         _mu_e.push_back( iMuon->energy() );
         _mu_eta.push_back( iMuon->eta() );
         _mu_phi.push_back( iMuon->phi() );
         _mu_isglobal.push_back( iMuon->isGlobalMuon() );
         _mu_istracker.push_back( iMuon->isTrackerMuon() );
         _mu_isPFMuon.push_back( iMuon->isPFMuon() );
         _mu_iso.push_back(iMuon->pfIsolationR03().sumChargedHadronPt );
         _mu_MediumId.push_back(isMediumMuon( *iMuon  ));
         _mu_dxy.push_back( iMuon->muonBestTrack()->dxy(iVertex->position()) );
         _mu_dz.push_back(iMuon->muonBestTrack()->dz(iVertex->position()) );
         _mu_chg.push_back(iMuon->charge());
         std::cout<< "MU PX: "<< iMuon->px() << std::endl;
	 std::cout<< "MU PT: "<< iMuon->pt() << std::endl;
	 // std::cout<< "MU dxy: "<< iMuon->muonBestTrack()->dxy(iVertex->position()) << std::endl;
	 // std::cout<< "MU dz: "<< iMuon->muonBestTrack()->dz(iVertex->position()) << std::endl;
	 // std::cout<< "Muon charge: "<< iMuon->charge() << std::endl; 

    //    }
   }
   // **********************************************************************************





   // *********************************** TAUS *****************************************
   edm::Handle<pat::TauCollection>TauHandle;
   iEvent.getByLabel("slimmedTaus",TauHandle);
   const pat::TauCollection*  Taus = TauHandle.product();
   for(pat::TauCollection::const_iterator iTau = Taus->begin(); iTau != Taus->end();iTau++){
     _tau_pt.push_back( iTau->pt() );                                                                                                                                                    
     _tau_px.push_back( iTau->px() );                                                                                                                                                    
     _tau_py.push_back( iTau->py() );                                                                                                                                                    
     _tau_pz.push_back( iTau->pz() );                                                                                                                                                    
     _tau_e.push_back( iTau->energy() );                                                                                                                                                 
     _tau_eta.push_back( iTau->eta() );                                                                                                                                                  
     _tau_phi.push_back( iTau->phi() );
     _tau_chg.push_back(iTau->charge());   
     _tau_ID_NewDMs.push_back( iTau->tauID("decayModeFindingNewDMs")   );
     _tau_ID_Elec_VLooseMVA5.push_back( iTau->tauID("againstElectronVLooseMVA5") );
     _tau_ID_MuonTight3.push_back( iTau->tauID("againstMuonTight3")  );
     _tau_ID_DBCorrRaw3Hits.push_back( iTau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits")  );

     pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(iTau->leadChargedHadrCand().get());
     _tau_dz.push_back(packedLeadTauCand->dz()); 
     std::cout<< "Tau ID: "<< iTau->tauID("decayModeFindingNewDMs")<< std::endl;


     // std::cout<< "TAU DZ: "<< packedLeadTauCand->dz() << std::endl;
     std::cout<< "TAU PX: "<< iTau->px() << std::endl;
     std::cout<< "TAU PT: "<< iTau->pt() << std::endl;
   }
   // ************************************************************************************





   /*
   // OBSOLETE 
   edm::Handle<edm::View<reco::PFTau>>TauHandle;
   iEvent.getByLabel("hpsPFTauProducer",TauHandle);  
   const edm::View<reco::PFTau>* Taus = TauHandle.product();
   math::XYZTLorentzVector pfour_tau, pfour_muon;
   for(edm::View<reco::PFTau>::const_iterator iTau = Taus->begin(); iTau != Taus->end();iTau++){
       _tau_pt.push_back( iTau->pt() ); 
       pfour_tau = iTau->p4();
       std::cout<< "iTau pt(): " << iTau->pt() << " Tau mass: "<< pfour_tau.mass() << std::endl; 
       // if(iTau->tauID("hpsPFTauDiscriminationByDecayModeFindingNewDMs")){std::cout<< "NewDMs passed" << std::endl;}

      for(edm::View<reco::Muon>::const_iterator iMuon = Muons->begin(); iMuon != Muons->end();iMuon++){
          pfour_muon = iMuon->p4();
          pfour_muon += pfour_tau;
          _mvis.push_back( pfour_muon.mass() );
        } // Muon Loop ends
   } // Tau Loop ends
   */


   


/*
   edm::Handle<reco::PFTauCollection> tauCandidates;
   iEvent.getByLabel(srcTauCandidates_, tauCandidates);

   size_t numTauCandidates = tauCandidates->size();
   for ( size_t idxTauCandidate = 0; idxTauCandidate < numTauCandidates; ++idxTauCandidate ) { // NEW TAU LOOP
     reco::PFTauRef tauCandidate(tauCandidates, idxTauCandidate);
 
     bool passesDiscriminators = true;
     double tauDiscriminatorValue_chg_iso = 0.;
     double tauDiscriminatorValue_neutral_iso = 0.;
     double tauDiscriminatorValue_pucorr_iso = 0.;

     for ( vInputTag::const_iterator srcDiscriminator = srcDiscriminators_.begin();
	   srcDiscriminator != srcDiscriminators_.end(); ++srcDiscriminator ) { // LOOP OVER TAU DISCRIM.S
       edm::Handle<reco::PFTauDiscriminator> tauDiscriminators;
       iEvent.getByLabel(*srcDiscriminator, tauDiscriminators);
       double minValue = 0.5;
       double maxValue = 1.e+3;
       if ( srcDiscriminator->label() == "hpsPFTauMVA3IsolationChargedIsoPtSum" ) { // special treatment of "discriminators" that are in fact isolation pT-sums
	 minValue = -1.;
	 maxValue = 2.;
       }

       if ( srcDiscriminator->label() == "hpsPFTauMVA3IsolationChargedIsoPtSum" ||
            srcDiscriminator->label() == "hpsPFTauDiscriminationByDecayModeFindingNewDMs" ||
            srcDiscriminator->label() == "hpsPFTauDiscriminationByLooseMuonRejection3" ) {
           double tauDiscriminatorValue = (*tauDiscriminators)[tauCandidate];
       if( !(tauDiscriminatorValue > minValue && tauDiscriminatorValue < maxValue) ){ passesDiscriminators = false;}
       }


       if( srcDiscriminator->label() == "hpsPFTauMVA3IsolationChargedIsoPtSum" ){
         tauDiscriminatorValue_chg_iso = (*tauDiscriminators)[tauCandidate];   
	 // std::cout << "Chg Iso. value: "<< tauDiscriminatorValue_chg_iso << std::endl;
       }


       if( srcDiscriminator->label() == "hpsPFTauMVA3IsolationNeutralIsoPtSum" ){
         tauDiscriminatorValue_neutral_iso = (*tauDiscriminators)[tauCandidate];   
	 // std::cout << "Neutral Iso. value: "<< tauDiscriminatorValue_neutral_iso << std::endl;
       }


      if( srcDiscriminator->label() == "hpsPFTauMVA3IsolationPUcorrPtSum" ){
         tauDiscriminatorValue_pucorr_iso = (*tauDiscriminators)[tauCandidate];   
	 // std::cout << "PU corr. Iso. value: "<< tauDiscriminatorValue_pucorr_iso << std::endl;
       }
    } // LOOP OVER TAU DISC ENDS
   


     // std::cout << "Chg Iso. value: "<< tauDiscriminatorValue_chg_iso << " Neutral Iso. value: "<< tauDiscriminatorValue_neutral_iso << " PU corr. Iso. value: "<< tauDiscriminatorValue_pucorr_iso << std::endl;


// if( passesDiscriminators && std::fabs(tauCandidate->eta()) < 2.3  ){
     // std::cout<< "Passed Tau pT: "<< tauCandidate->pt() << std::endl;
         _tau_pt.push_back( tauCandidate->pt() ); 
         _tau_px.push_back( tauCandidate->px() );
         _tau_py.push_back( tauCandidate->py() );
         _tau_pz.push_back( tauCandidate->pz() );
         _tau_e.push_back( tauCandidate->energy() );
         _tau_eta.push_back( tauCandidate->eta() );
         _tau_phi.push_back( tauCandidate->phi() );
         _tau_passed.push_back(passesDiscriminators);
         _tau_chg_iso.push_back( tauDiscriminatorValue_chg_iso );  
         _tau_neutral_iso.push_back( tauDiscriminatorValue_neutral_iso );  
         _tau_pucorr_iso.push_back( tauDiscriminatorValue_pucorr_iso );  

//      }



   } // NEW TAU LOOP ENDS

*/





   /*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   */
   myTree->Fill();

   /*
   //  Clear all vectors
   _evt.clear();
   _run.clear();
   _lumi.clear();
   _tau_pt.clear();
   _tau_px.clear();
   _tau_py.clear();
   _tau_pz.clear();
   _tau_e.clear();
   _tau_eta.clear();
   _tau_phi.clear();
   _tau_passed.clear();
   _tau_chg_iso.clear();
   _tau_neutral_iso.clear();
   _tau_pucorr_iso.clear();
   _mu_pt.clear();
   _mu_px.clear();
   _mu_py.clear();
   _mu_pz.clear();
   _mu_e.clear();
   _mu_eta.clear();
   _mu_phi.clear();
   _mu_isglobal.clear();
   _mu_istracker.clear();
   _mu_iso.clear();
   _mvis.clear(); // NOT USED RIGHT NOW !          

   // MUTAU BASELINE QUANTITIES
   _mu_MediumId.clear();
   */




}




// ------------ method called once each job just after ending the event loop  ------------
void 
MYANALYZER::endJob() 
{

  //  Clear all vectors                                                                                                                                                                             
  _evt.clear();
  _run.clear();
  _lumi.clear();
  _tau_pt.clear();
  _tau_px.clear();
  _tau_py.clear();
  _tau_pz.clear();
  _tau_e.clear();
  _tau_eta.clear();
  _tau_phi.clear();
  _tau_passed.clear();
  _tau_chg_iso.clear();
  _tau_neutral_iso.clear();
  _tau_pucorr_iso.clear();
  _mu_pt.clear();
  _mu_px.clear();
  _mu_py.clear();
  _mu_pz.clear();
  _mu_e.clear();
  _mu_eta.clear();
  _mu_phi.clear();
  _mu_isglobal.clear();
  _mu_istracker.clear();
  _mu_iso.clear();
  _mvis.clear(); // NOT USED RIGHT NOW !                                                                                                                                                            

  // MUTAU BASELINE QUANTITIES                                                                                                                                                                      
  _mu_MediumId.clear();
  _mu_dxy.clear();
  _mu_dz.clear();
  _mu_chg.clear();
  _tau_chg.clear();
  _tau_dz.clear();
  _mu_isPFMuon.clear(); 
  _tau_ID_NewDMs.clear();
  _tau_ID_Elec_VLooseMVA5.clear();
  _tau_ID_MuonTight3.clear();
  _tau_ID_DBCorrRaw3Hits.clear();


}

// ------------ method called when starting to processes a run  ------------
/*
void 
MYANALYZER::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MYANALYZER::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MYANALYZER::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MYANALYZER::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MYANALYZER::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MYANALYZER);
