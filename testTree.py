#!/usr/bin/python
import ROOT
from ROOT import std, gROOT, gStyle, gPad, TCanvas, TH1, TH1D, TH2D, TLegend, TLine, TFile, TTree, TLorentzVector, TMath, TVirtualPad, TEventList


fXX = TFile("tops.SMIA.tmp.10k.root","READ")
tXX = fXX.Get("SMIA")
pXX = ROOT.vector(TLorentzVector)()
iXX = std.vector(int)()
tXX.SetBranchAddress("p4",pXX);
tXX.SetBranchAddress("id",iXX);

hmSMXgen = TH1D("hmSMXgen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSMXgen.Sumw2()
hmSMXgen.SetLineWidth(2)
hmSMXgen.SetLineColor(ROOT.kAzure+9)
hmSMXgen.SetMarkerColor(ROOT.kAzure+9)
hmSMXgen.SetMarkerStyle(24)

n=1
for event in tXX:
   if(n%1000==0): print "processed |SM+X|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()
   hmSMXgen.Fill(mtt)
   n+=1

cnv = TCanvas("cnv","",600,600)
cnv.Draw()
cnv.cd()
hmSMXgen.Draw()
cnv.Update()
cnv.RedrawAxis()
cnv.SaveAs("testTree.pdf")
