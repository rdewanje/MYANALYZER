#import FWCore.ParameterSet.Config as cms
#process = cms.Process("Demo")
#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
#    )
#)

#process.demo = cms.EDAnalyzer('PAIRANALYZER'
#)
#process.p = cms.Path(process.demo)


import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff") ### NEEDED FOR ELECTRON ID       

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "IDEAL_V9::All"                                                                                   
process.GlobalTag.globaltag = "MCRUN2_74_V9::All"

process.source = cms.Source("PoolSource",
   # replace 'myfile.root' with the source file you want to use                                                                                                                                     
    fileNames = cms.untracked.vstring(
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root'
    )
)

##### ELECTRON MVA ID JARGON #########                                                                                                                                                              
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
## dataFormat = DataFormat.AOD                                                                                                                    
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce                                                                                                                                                               
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']

#add them to the VID producer                                                                                                                                                                       
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
######################################                                                   


#from PAIRANALYZER.PAIRANALYZER.MuTau_pair_cff import *
#pairsequence = MutaupairSequence.clone()
#process.load("PAIRANALYZER.PAIRANALYZER.MuTau_pair_cff")


#process.demo = cms.EDAnalyzer('PAIRANALYZER'                                                                           
#)                                                   
#process.p = cms.Path( pairsequence * process.demo )

#process.p = cms.Path(MutaupairSequence)

################################################



### MUON LEG SELECTIONS                                                                                           
process.selectedPatMuons = cms.EDFilter("PATMuonSelector",                                                              
        src = cms.InputTag("slimmedMuons"),                                                   
       cut = cms.string("pt > 18. && abs(eta) < 2.1 && isGlobalMuon & isTrackerMuon")                                                                           
)           


process.PATMuSkimmedBy1 = cms.EDFilter("CandViewCountFilter",                                                          
                               src = cms.InputTag('selectedPatMuons'),                                                                             
                               minNumber = cms.uint32(1)                                                                                                           
)      



## TAU LEG SELECTIONS                                                                                                             
process.selectedPatTaus = cms.EDFilter("PATTauSelector",
     src = cms.InputTag("slimmedTaus"),
#     cut = cms.string("pt > 18. && abs(eta) < 2.3 && decayModeFindingNewDMs > 0.5")  ### decaymode not working
     cut = cms.string("pt > 18. && abs(eta) < 2.3")  ### decaymode not working

)

process.PATTauSkimmedBy1 = cms.EDFilter("CandViewCountFilter",                                                                        
                                src = cms.InputTag('selectedPatTaus'),                                                                           
                                minNumber = cms.uint32(1)                                                                                                             
)       

### MUTAU PAIR NON OVERLAP CRITERION                                                                      
process.MuTauPairs = cms.EDProducer("CandViewShallowCloneCombiner",                                  
                          decay = cms.string("selectedPatMuons selectedPatTaus"),                                                                                                          
                                    checkCharge = cms.bool(False),
                                    cut         = cms.string("sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+ min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi) ) * min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi) ) )>0.3"),                                            

)                




process.demo = cms.EDAnalyzer('MYANALYZER',
               muTauPairs        = cms.InputTag("MuTauPairs"),
               arbitration       = cms.string("None"),
               eleMediumIdMap    = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80"),
               eleTightIdMap     = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90"),
               mvaValuesMap      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues"),
               mvaCategoriesMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigCategories"),
               TauCollection     = cms.untracked.string("hpsPFTauProducer"),
               MuonCollection    = cms.untracked.string("muons"),
               fileName          = cms.untracked.string("Output.root"),
               srcTauCandidates  = cms.InputTag('hpsPFTauProducer'),
               srcDiscriminators = cms.VInputTag(
                    'hpsPFTauDiscriminationByDecayModeFindingNewDMs',
                    'hpsPFTauMVA3IsolationChargedIsoPtSum',
                    'hpsPFTauDiscriminationByLooseMuonRejection3',
                    'hpsPFTauMVA3IsolationNeutralIsoPtSum',
                    'hpsPFTauMVA3IsolationPUcorrPtSum'
               )

               )




#process.p = cms.Path( process.selectedPatTaus * process.demo)
process.p = cms.Path( ( (process.selectedPatMuons)  +  (process.selectedPatTaus) ) * process.MuTauPairs * ( process.egmGsfElectronIDSequence * process.demo))





'''
process.PATTauSkimmedBy1 = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedPatTaus'),
    minNumber = cms.uint32(1)
)

### MUON LEG SELECTIONS                                                                                                       
process.selectedPatMuons = cms.EDFilter("PATMuonSelector",
     src = cms.InputTag("slimmedMuons"),
     cut = cms.string("pt > 18. && abs(eta) < 2.1 && isGlobalMuon & isTrackerMuon")
)

process.PATMuSkimmedBy1 = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag('selectedPatMuons'),
    minNumber = cms.uint32(1)
)

### MUTAU PAIR NON OVERLAP CRITERION                                                                                          
process.MuTauPairs = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("PATTauSkimmedBy1 PATMuSkimmedBy1"),
                                    checkCharge = cms.bool(False),
                                    cut         = cms.string("sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+ min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi) ) * min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi) ) )>0.3"),

)

process.MuTauPairSkimmedBy1 = cms.EDFilter("CandViewCountFilter",
   src = cms.InputTag('MuTauPairs'),
   minNumber = cms.uint32(1)
)

## Make the tree                                                                                                  
process.pairTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
     # pairs                                                                                                           
     tagProbePairs = cms.InputTag("MuTauPairSkimmedBy1"),
     arbitration   = cms.string("None"),
     # variables to use                                                                                                                         
     variables = cms.PSet(
         ## methods of reco::Candidate                                                                                                
         eta = cms.string("eta"),
         pt  = cms.string("pt"),
         ## a method of the reco::Muon object (thanks to the 3.4.X StringParser)   
         #nsegm = cms.string("numberOfMatches"),
         ## this one is an external variable                    
         #drj = cms.InputTag("drToNearestJet"),                                                                                                                                   
    ),
     # choice of what defines a 'passing' probe                                                                           
     flags = cms.PSet(),                                                        
     ## one defined by an external collection of passing probes                                                                                                                                  
     #    passingCal = cms.InputTag("probesPassingCal"),                                            
     ## two defined by simple string cuts                                                                                                                                                        
     #    passingGlb = cms.string("isGlobalMuon"),      
     #    passingIso = cms.string("(isolationR03.hadEt+isolationR03.emEt+isolationR03.sumPt) < 0.1 * pt"),                                                           
     #),                                                                                                                                              
     # mc-truth info                                                                                                             
     isMC = cms.bool(False),
     #motherPdgId = cms.vint32(22,23),                                                                                                              
     #makeMCUnbiasTree = cms.bool(True),                                                                                      
     #checkMotherInUnbiasEff = cms.bool(True),                      
     #tagMatches = cms.InputTag("muMcMatch"),                                                                        
     #probeMatches  = cms.InputTag("muMcMatch"),                                                               
     #allProbes     = cms.InputTag("probeMuons"),                                                                                                                                   
 )


process.p = cms.Path( process.selectedPatTaus * process.PATTauSkimmedBy1 * process.selectedPatMuons * process.PATMuSkimmedBy1 * process.MuTauPairs * process.MuTauPairSkimmedBy1 * process.pairTree )
'''

################################################


process.TFileService = cms.Service("TFileService", fileName = cms.string("testTagProbeFitTreeProducer_ZMuTau.root"))
