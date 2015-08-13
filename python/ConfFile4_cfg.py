import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/16FCAB85-D527-E511-BBB5-0025905A6090.root',
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/5CC262BC-BD27-E511-B09E-002618943836.root',
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/A0A80B20-E627-E511-9245-00259059642E.root',
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/A4A4849D-A927-E511-B9FF-00261894396E.root',
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/B620EFF3-9B27-E511-8119-003048FFCB84.root',
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/CA5826B2-7228-E511-816C-002590596486.root',
    '/store/relval/CMSSW_7_4_6_patch6/SingleMu/RAW-RECO/MuTau-GR_H_V57A_skimTest_RelVal_mu2012A-v1/00000/E0E37E47-AC28-E511-BFEA-00261894389A.root'


    )
)

process.TFileService=cms.Service('TFileService',fileName=cms.string('MuTau_skim_SingleMu_RELVAL_TEST.root'))

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


process.p = cms.Path(process.demo)
