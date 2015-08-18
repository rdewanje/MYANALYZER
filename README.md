# MYANALYZER
My sample EDAnalyzer to read MuTau skim related quantities in the form of MuTau pairs per event on the lines of "PhysicsTools/TagAndProbe" (in CMSSW_7_4_8_pre1)
The actual code is implemented inside the plugins/MYANALYZER.cc while the config file to run the analyzer is 
python/ConfFile_MuTaupair_cfg.py

The relevant files for the MuTauPairMaker are: 

interface/MuTauPairMaker.h
src/MuTauPairMaker.cc