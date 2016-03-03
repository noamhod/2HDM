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
mX     = THDM.model.mX
nameX  = THDM.model.nameX
mgpath = THDM.model.mgpath
### make the libraries
THDM.testTHDM(mgpath,nameX,mX)