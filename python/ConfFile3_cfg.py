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
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/002F7FDD-BA13-E511-AA63-0026189437F5.root', ## 2
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00610EE7-C213-E511-842C-00304833529A.root',
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00623DCC-A813-E511-A302-0025905B85A2.root',  ## 3
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0071047F-E813-E511-B38C-842B2B2922E2.root',  ## 4
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/008C704F-C313-E511-B6FC-0025905B85F6.root',
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02147476-C813-E511-984D-90B11C27F8B2.root', ## 5
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02A85CDF-BA13-E511-AF45-00259073E38A.root',
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02B3B2D1-BC13-E511-A895-008CFA110B10.root',
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02C3CC1F-1714-E511-8569-0025905A60E0.root',
#    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02F7DC38-BF13-E511-AD70-6C3BE5B594A0.root'


    )
)





###### ELECTRON MVA ID JARGON #########
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# dataFormat = DataFormat.AOD
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
######################################







#process.TFileService=cms.Service('TFileService',fileName=cms.string('/eos/uscms/store/user/rkd123/MuTau_skim_collections_DYJETS_AODSIM_ALL.root'))
process.TFileService=cms.Service('TFileService',fileName=cms.string('MuTau_skim_collections_DYJETS_MINIAODSIM_ALL.root'))


process.demo = cms.EDAnalyzer('MYANALYZER',
               TauCollection  = cms.untracked.string("hpsPFTauProducer"),                    
               MuonCollection = cms.untracked.string("muons"),       
               fileName       = cms.untracked.string("Output.root"),
               srcTauCandidates = cms.InputTag('hpsPFTauProducer'),
               srcDiscriminators = cms.VInputTag(
                    'hpsPFTauDiscriminationByDecayModeFindingNewDMs',
                    'hpsPFTauMVA3IsolationChargedIsoPtSum',
                    'hpsPFTauDiscriminationByLooseMuonRejection3',
                    'hpsPFTauMVA3IsolationNeutralIsoPtSum', 
                    'hpsPFTauMVA3IsolationPUcorrPtSum'
               )       

               )


## DEFAULT LINES  (BEFORE ELECTRON ID)
# process.p = cms.Path(process.demo)

# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmGsfElectronIDSequence * process.demo)
