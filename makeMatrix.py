#!/usr/bin/python

import ROOT
from ROOT import std, gROOT, gStyle, gPad, TCanvas, TH1, TH1D, TH2D, TLegend, TLine, TFile, TTree, TLorentzVector, TMath, TVirtualPad, TEventList
import kinematics
import THDM
import sys
import os

##################
### Begin code ###
##################

### make the definitions
type=THDM.model.type # 2
sba=THDM.model.sba #1
mX=THDM.model.mX #500
nameX=THDM.model.nameX #"H"
cuts = THDM.model.cuts #"m"+nameX+"=="+str(mX)+" && sba==1 && TMath::ATan(tanb)>0. && TMath::ATan(tanb)<TMath::Pi()/2. && TMath::Abs(cba)<=1. && type=="+str(type)+" && (status&3)==0"
mgpath = THDM.model.mgpath #"/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/"
### make the libraries
THDM.testTHDM(mgpath,nameX,mX)

# testTHDM("/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/","H",750,"(width_H/mH>0.05 && width_H/mH<0.1)")
# testTHDM("/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/","H",450)
# testTHDM("/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/","A",450)