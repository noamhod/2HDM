# 2HDM
2HDM reweighting code based on Madgraph and Pythia8
Prequisits:
- MG5_aMC_v2_3_3 or later
- Higgs_Effective_Couplings_FormFact model
- model parameters ntuple from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsBSM2HDMRecommendations
- pythia8215 or later
- Python (with numpy and scipy)
- ROOT v5 (not yet 6), copmpatile with the available Python version

Put all files found here in the pythia8215/examples directory and have fun :)

To run,
- edit the number of events to generate in main34.cc
- make main34
- execute `./main34 SM`
- execute `./main34 H`
- etc...
- edit the t2HDM class in THDM.py with the parameters of the models or load it with other parameters later
- execute `python makeMatrix.py` to make all libraries for this configuration
- execute `python run2.py` to reweight the SM sample into this configuration and store the weights in the tree
- execute `python template.py` to plot the tempaltes
- execute `python run3.py` to validate a single template
