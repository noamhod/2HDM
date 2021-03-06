#!/usr/bin/python

import ROOT
from ROOT import *
import kinematics as kin
import THDM
import sys
import os
import math
from array import *

def normalizeToBinWidth(h, unit=0): ## normalize to nEvts/40 GeV
   for b in xrange(h.GetNbinsX()+1):
      sf = unit/h.GetBinWidth(b) if(unit>0) else 1./h.GetBinWidth(b)
      h.SetBinContent(b, h.GetBinContent(b)*sf)
      h.SetBinError(b,   h.GetBinError(b)*sf)
   return h


def initdecaychain():
   chain = {
      "t"  : {"t":-1,"g1":-1,"g2":-1,"w+":-1,"fd":-1,"fu":-1,"b":-1},
      "t~" : {"t~":-1,"g1":-1,"g2":-1,"w-":-1,"fd":-1,"fu":-1,"b~":-1},
   }
   return chain

def setdecaychain(chain,pdgid,parents,children):
   for j in xrange(pdgid.size()):
      if(pdgid[j]==21 and chain["t"]["g1"]<0):
         chain["t"]["g1"]  = j
         chain["t~"]["g1"] = j
      if(pdgid[j]==21 and chain["t"]["g2"]<0 and chain["t"]["g1"]!=j):
         chain["t"]["g2"]  = j
         chain["t~"]["g2"] = j
      if(pdgid[j]==6):
         p1 = parents[j][0]
         p2 = parents[j][1]
         c1 = children[j][0]
         c2 = children[j][1]
         c11 = children[c1][0]
         c12 = children[c1][1]
         c21 = children[c2][0]
         c22 = children[c2][1]
         chain["t"]["t"] = j
         if(pdgid[c1]==24):
            chain["t"]["w+"] = c1
            chain["t"]["b"] = c2
            chain["t"]["fd"] = c11 if(pdgid[c11]%2!=0) else c12
            chain["t"]["fu"] = c11 if(pdgid[c11]%2==0) else c12
         if(pdgid[c2]==24):
            chain["t"]["w+"] = c2
            chain["t"]["b"] = c1
            chain["t"]["fd"] = c21 if(pdgid[c21]%2!=0) else c22
            chain["t"]["fu"] = c21 if(pdgid[c21]%2==0) else c22

      if(pdgid[j]==-6):
         p1 = parents[j][0]
         p2 = parents[j][1]
         c1 = children[j][0]
         c2 = children[j][1]
         c11 = children[c1][0]
         c12 = children[c1][1]
         c21 = children[c2][0]
         c22 = children[c2][1]
         chain["t~"]["t~"] = j
         if(pdgid[c1]==-24):
            chain["t~"]["w-"] = c1
            chain["t~"]["b~"]  = c2
            chain["t~"]["fd"] = c11 if(pdgid[c11]%2!=0) else c12
            chain["t~"]["fu"] = c11 if(pdgid[c11]%2==0) else c12
         if(pdgid[c2]==-24):
            chain["t~"]["w1"] = c2
            chain["t~"]["b~"]  = c1
            chain["t~"]["fd"] = c21 if(pdgid[c21]%2!=0) else c22
            chain["t~"]["fu"] = c21 if(pdgid[c21]%2==0) else c22
      if(chain["t"]["b"]>=0 and chain["t~"]["b~"]>=0): break
   return chain


def printdcaychain(chian,pdgid,parents,children):
   for j in xrange(pdgid.size()):
      print '[%g] pdgid=%g' % (j,pdgid[j])
      if(parents[j].size()>1):  print '  parent[%g]=%g, parent[%g]=%g' % (parents[j][0],pdgid[parents[j][0]],parents[j][1],pdgid[parents[j][1]])
      if(children[j].size()>1): print '    child[%g]=%g, child[%g]=%g'   % (children[j][0],pdgid[children[j][0]],children[j][1],pdgid[children[j][1]])
   print 'chain["t"]["t"]   = %g, pdgid=%g'  % (chain["t"]["t"],   pdgid[chain["t"]["t"]])
   print 'chain["t"]["g1"]  = %g, pdgid=%g'  % (chain["t"]["g1"],  pdgid[chain["t"]["g1"]])
   print 'chain["t"]["g2"]  = %g, pdgid=%g'  % (chain["t"]["g2"],  pdgid[chain["t"]["g2"]])
   print 'chain["t"]["b"]  = %g, pdgid=%g'  % (chain["t"]["b"],  pdgid[chain["t"]["b"]])
   print 'chain["t"]["w+"]  = %g, pdgid=%g'  % (chain["t"]["w+"],  pdgid[chain["t"]["w+"]])
   print 'chain["t"]["fd"]  = %g, pdgid=%g'  % (chain["t"]["fd"],  pdgid[chain["t"]["fd"]])
   print 'chain["t"]["fu"]  = %g, pdgid=%g'  % (chain["t"]["fu"],  pdgid[chain["t"]["fu"]])
   print 'chain["t~"]["t~"] = %g, pdgid=%g'  % (chain["t~"]["t~"], pdgid[chain["t~"]["t~"]])
   print 'chain["t~"]["g1"] = %g, pdgid=%g' % (chain["t~"]["g1"],  pdgid[chain["t~"]["g1"]])
   print 'chain["t~"]["g2"] = %g, pdgid=%g' % (chain["t~"]["g2"],  pdgid[chain["t~"]["g1"]])
   print 'chain["t~"]["b~"]  = %g, pdgid=%g'  % (chain["t~"]["b~"],  pdgid[chain["t~"]["b~"]])
   print 'chain["t~"]["w-"] = %g, pdgid=%g' % (chain["t~"]["w-"],  pdgid[chain["t~"]["w-"]])
   print 'chain["t~"]["fd"] = %g, pdgid=%g' % (chain["t~"]["fd"],  pdgid[chain["t~"]["fd"]])
   print 'chain["t~"]["fu"] = %g, pdgid=%g' % (chain["t~"]["fu"],  pdgid[chain["t~"]["fu"]])


def setStyle():
   gROOT.Reset()
   icol=0;  # WHITE
   font=42; # Helvetica
   tsize=0.05;
   gStyle.SetFrameBorderMode(icol);
   gStyle.SetFrameFillColor(icol);
   gStyle.SetCanvasBorderMode(icol);
   gStyle.SetCanvasColor(icol);
   gStyle.SetPadBorderMode(icol);
   gStyle.SetPadColor(icol);
   gStyle.SetStatColor(icol);
   gStyle.SetPaperSize(20,26);
   gStyle.SetPadTopMargin(0.05);
   gStyle.SetPadRightMargin(0.08);
   gStyle.SetPadBottomMargin(0.15);
   gStyle.SetPadLeftMargin(0.12);
   gStyle.SetTitleXOffset(1.05);
   gStyle.SetTitleYOffset(0.95);
   gStyle.SetTextFont(font);
   gStyle.SetTextSize(tsize);
   gStyle.SetLabelFont(font,"x");
   gStyle.SetTitleFont(font,"x");
   gStyle.SetLabelFont(font,"y");
   gStyle.SetTitleFont(font,"y");
   gStyle.SetLabelFont(font,"z");
   gStyle.SetTitleFont(font,"z");
   gStyle.SetLabelSize(tsize*0.85,"x");
   gStyle.SetTitleSize(tsize*1.10,"x");
   gStyle.SetLabelSize(tsize*0.85,"y");
   gStyle.SetTitleSize(tsize*1.10,"y");
   gStyle.SetLabelSize(tsize*0.85,"z");
   gStyle.SetTitleSize(tsize*1.10,"z");
   gStyle.SetMarkerStyle(20);
   gStyle.SetMarkerSize(1.);
   gStyle.SetHistLineWidth(2);
   gStyle.SetLineStyleString(2,"[12 12]"); # postscript dashes
   gStyle.SetEndErrorSize(0.);
   # gStyle.SetOptTitle(0);
   gStyle.SetOptStat(0);
   gStyle.SetOptFit(0);
   gStyle.SetPadTickX(1);
   gStyle.SetPadTickY(1);
   ROOT.TGaxis.SetMaxDigits(4)

def plot(h1, h2, h3, h4, tanb, wX, nkevents, cme, lumi, fname, ymin=-1, ymax=-1):
   cnv = TCanvas("cnv","",600,600);
   cnv.Divide(1,3);
   tvp_hists = cnv.cd(1);
   tvp_subtr = cnv.cd(2);
   tvp_ratio = cnv.cd(3);
   cnv.Draw();
   tvp_hists.SetPad(0.00, 0.40, 1.00, 1.00);
   tvp_subtr.SetPad(0.00, 0.20, 1.00, 0.40);
   tvp_ratio.SetPad(0.00, 0.02, 1.00, 0.20);
   tvp_hists.SetBottomMargin(0.012);
   tvp_subtr.SetBottomMargin(0.20);
   tvp_subtr.SetTopMargin(0.012);
   tvp_ratio.SetBottomMargin(0.40);
   tvp_ratio.SetTopMargin(0.012);
   tvp_hists.SetTicks(1,1);
   tvp_subtr.SetTicks(1,1);	
   tvp_ratio.SetTicks(1,1);	
	
   sXtitle = h1.GetXaxis().GetTitle();
   cloneName_n = h1.GetName();
   cloneName_d = h2.GetName();
   th1n_tmp = h1.Clone(cloneName_n+"_th1n_tmp");
   th1d_tmp = h2.Clone(cloneName_d+"_th1d_tmp");
   th1n_tmp.SetBinErrorOption(TH1.kPoisson);
   th1d_tmp.SetBinErrorOption(TH1.kPoisson);
	
   # hs = h2.Clone("subtr");
   # hs.SetTitle(";"+sXtitle+";|SM+#it{"+nameX+"}|^{2}/|SM|^{2}-1 [%]");
   # hs.Add(h1,-1.)
   # hs.Divide(h1)
   # hs.Scale(100)

   rmin=-13;
   rmax=+13;
   h4.Divide(h1)
   h4.Scale(100)
   h4.SetTitle(";"+sXtitle+";(SM+#it{"+nameX+"})/SM-1 [%]");
   xLabelSize = h4.GetXaxis().GetLabelSize()*2;
   yLabelSize = h4.GetYaxis().GetLabelSize()*2;
   xTitleSize = h4.GetXaxis().GetTitleSize()*2;
   yTitleSize = h4.GetYaxis().GetTitleSize()*2;
   titleSize  = h4.GetTitleSize()           *2;
   h4.GetXaxis().SetLabelSize(xLabelSize);
   h4.GetYaxis().SetLabelSize(yLabelSize);
   h4.GetXaxis().SetTitleSize(xTitleSize);
   h4.GetYaxis().SetTitleSize(yTitleSize);
   h4.SetTitleSize(titleSize);
   h4.GetYaxis().SetTitleOffset(0.42);
   h4.GetXaxis().SetTitleOffset(0.83);
   h4.SetMinimum(rmin);
   h4.SetMaximum(rmax);
   lineS = TLine(h4.GetXaxis().GetXmin(),0.,h4.GetXaxis().GetXmax(),0.);

   rmin=+0.83;
   rmax=+1.17;
   hr = th1d_tmp.Clone("ratio")
   hr.SetTitle(";"+sXtitle+";(SM+#it{"+nameX+"})_{rwt}/(SM+#it{"+nameX+"})_{gen}")
   hr.Divide(h3)
   hr.SetMarkerStyle(h3.GetMarkerStyle());
   hr.SetMarkerSize(h3.GetMarkerSize());
   hr.SetMarkerColor(h3.GetMarkerColor());
   hr.SetLineColor(h3.GetLineColor());
   hr.SetLineStyle(h3.GetLineStyle());
   hr.SetLineWidth(h3.GetLineWidth());
   xLabelSize = hr.GetXaxis().GetLabelSize()*1.9
   yLabelSize = hr.GetYaxis().GetLabelSize()*1.9
   xTitleSize = hr.GetXaxis().GetTitleSize()*1.9
   yTitleSize = hr.GetYaxis().GetTitleSize()*1.9
   titleSize  = hr.GetTitleSize()           *1.9
   hr.GetXaxis().SetLabelSize(xLabelSize);
   hr.GetYaxis().SetLabelSize(yLabelSize);
   hr.GetXaxis().SetTitleSize(xTitleSize);
   hr.GetYaxis().SetTitleSize(yTitleSize);
   hr.SetTitleSize(titleSize);
   hr.GetYaxis().SetTitleOffset(0.4);
   hr.GetXaxis().SetTitleOffset(0.83);
   hr.SetMinimum(rmin);
   hr.SetMaximum(rmax);
   lineR = TLine(hr.GetXaxis().GetXmin(),1.,hr.GetXaxis().GetXmax(),1.);

   tvp_hists.SetBottomMargin(0);
   # tvp_subtr.SetTopMargin(0);
   tvp_subtr.SetBottomMargin(0);
   tvp_ratio.SetTopMargin(0);
   tvp_ratio.SetBottomMargin(0.2);
   
   leg = TLegend(0.5,0.4,0.87,0.9,"","brNDC")
   leg.SetFillStyle(4000); # will be transparent
   leg.SetFillColor(0)
   leg.SetTextFont(42)
   leg.SetBorderSize(0)
   leg.AddEntry(0, "MadGraph+Pythia8", "")
   leg.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} ("+str(nkevents)+"k events)", "")
   leg.AddEntry(0, cme+", "+lumi, "")
   leg.AddEntry(h1,"SM","ple")
   leg.AddEntry(h2,"SM+#it{"+nameX+"} reweighted","ple")
   leg.AddEntry(h3,"SM+#it{"+nameX+"} generated","ple")
   leg.AddEntry(0, "sin(#beta-#alpha)=1", "")
   leg.AddEntry(0, "tan#beta="+str(tanb), "")
   leg.AddEntry(0, "#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV", "")
   leg.AddEntry(0, "#Gamma_{#it{"+nameX+"}}/#it{m}_{#it{"+nameX+"}}="+str(wX)+" [%]", "")
   
   tvp_hists.cd()
   # tvp_hists.SetLogy()
   if(ymin>-1): h2.SetMinimum(ymin)
   else:        h2.SetMinimum(1)
   if(ymax>-1): h2.SetMaximum(ymax)
   h2.Draw();
   h1.Draw("same");
   h3.Draw("same");
   leg.Draw("same");
   tvp_hists.Update();
   tvp_hists.RedrawAxis();
   
   tvp_subtr.cd();
   tvp_subtr.SetGridy();
   h4.Draw("e1p");
   lineS.Draw("same");
   h4.Draw("e1p same");
   tvp_subtr.Update();
   tvp_subtr.RedrawAxis();

   tvp_ratio.cd();
   tvp_ratio.SetGridy();
   hr.Draw("e1p");
   lineR.Draw("same");
   hr.Draw("same");
   tvp_ratio.Update();
   tvp_ratio.RedrawAxis();
   
   cnv.Update();
   cnv.SaveAs(fname.replace("(","").replace(")","").replace(".pdf",".single.pdf"));
   cnv.SaveAs(fname);


##################
### Begin code ###
##################

setStyle()

type   = THDM.model.type
sba    = THDM.model.sba
mX     = THDM.model.mX
nameX  = THDM.model.nameX
cuts   = THDM.model.cuts
mgpath = THDM.model.mgpath
alphaS = THDM.model.alphaS
nhel   = THDM.model.nhel
libmatrix = "/Users/hod/GitHub/2HDM/matrix/"+nameX+"/"+str(mX)+"/"+str(sba)+"/"
THDM.setParameters(nameX,mX,cuts,type,sba)
###############
index=0 
THDM.setModules(libmatrix,nameX,len(THDM.parameters),"All")
print THDM.modules
###############


listbins = [0,80,160,240,320,360,400,440,500,560,600,640,680,720,760,800,860,920,1040,1160,1280]
arrbins = array("d", listbins)
nbins = len(listbins)-1

### mass histos
hmSMgen    = TH1D("hmSMgen",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXXrwt    = TH1D("hmXXrwt",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0rwt    = TH1D("hmX0rwt",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXIrwt    = TH1D("hmXIrwt",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXI0rwt   = TH1D("hmXI0rwt",  ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXXgen    = TH1D("hmXXgen",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0gen    = TH1D("hmX0gen",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXIabs    = TH1D("hmXIabs",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmX0abs    = TH1D("hmX0abs",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmX0absGen = TH1D("hmX0absGen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmXIabsStd = TH1D("hmXIabsStd", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0absStd = TH1D("hmX0absStd", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0absStdGen = TH1D("hmX0absStdGen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)

### cosphi histos
hcSMgen = TH1D("hcSMgen",   ";Truth cos#phi;Normalised",25,-1,+1)
# hcXXrwt = TH1D("hcXXrwt",   ";Truth cos#phi;Normalised",25,-1,+1)
hcX0rwt = TH1D("hcX0rwt",   ";Truth cos#phi;Normalised",25,-1,+1)
# hcXXgen = TH1D("hcXXgen",   ";Truth cos#phi;Normalised",25,-1,+1)
hcX0gen = TH1D("hcX0gen",   ";Truth cos#phi;Normalised",25,-1,+1)



hmSMgen.Sumw2()
hmSMgen.SetLineWidth(2)
hmSMgen.SetLineColor(ROOT.kBlack)
hmSMgen.SetMarkerColor(ROOT.kBlack)
hmSMgen.SetMarkerStyle(24)
hcSMgen.Sumw2()
hcSMgen.SetLineWidth(2)
hcSMgen.SetLineColor(ROOT.kBlack)
hcSMgen.SetMarkerColor(ROOT.kBlack)
hcSMgen.SetMarkerStyle(24)

hmXXrwt.Sumw2()
hmXXrwt.SetLineWidth(2)
hmXXrwt.SetLineColor(ROOT.kRed)
hmXXrwt.SetMarkerColor(ROOT.kRed)
hmXXrwt.SetMarkerStyle(20)
# hcXXrwt.Sumw2()
# hcXXrwt.SetLineWidth(2)
# hcXXrwt.SetLineColor(ROOT.kRed)
# hcXXrwt.SetMarkerColor(ROOT.kRed)
# hcXXrwt.SetMarkerStyle(20)

hmX0rwt.Sumw2()
hmX0rwt.SetLineWidth(2)
hmX0rwt.SetLineColor(ROOT.kRed)
hmX0rwt.SetMarkerColor(ROOT.kRed)
hmX0rwt.SetMarkerStyle(20)
hcX0rwt.Sumw2()
hcX0rwt.SetLineWidth(2)
hcX0rwt.SetLineColor(ROOT.kRed)
hcX0rwt.SetMarkerColor(ROOT.kRed)
hcX0rwt.SetMarkerStyle(20)

hmXIrwt.Sumw2()
hmXIrwt.SetLineWidth(2)
hmXIrwt.SetLineColor(ROOT.kRed)
hmXIrwt.SetMarkerColor(ROOT.kRed)
hmXIrwt.SetMarkerStyle(20)

hmXI0rwt.Sumw2()
hmXI0rwt.SetLineWidth(2)
hmXI0rwt.SetLineColor(ROOT.kRed)
hmXI0rwt.SetMarkerColor(ROOT.kRed)
hmXI0rwt.SetMarkerStyle(20)

hmXXgen.Sumw2()
hmXXgen.SetLineWidth(2)
hmXXgen.SetLineColor(ROOT.kAzure+9)
hmXXgen.SetMarkerColor(ROOT.kAzure+9)
hmXXgen.SetMarkerStyle(24)
# hcXXgen.Sumw2()
# hcXXgen.SetLineWidth(2)
# hcXXgen.SetLineColor(ROOT.kAzure+9)
# hcXXgen.SetMarkerColor(ROOT.kAzure+9)
# hcXXgen.SetMarkerStyle(24)

hmX0gen.Sumw2()
hmX0gen.SetLineWidth(2)
hmX0gen.SetLineColor(ROOT.kAzure+9)
hmX0gen.SetMarkerColor(ROOT.kAzure+9)
hmX0gen.SetMarkerStyle(24)
hcX0gen.Sumw2()
hcX0gen.SetLineWidth(2)
hcX0gen.SetLineColor(ROOT.kAzure+9)
hcX0gen.SetMarkerColor(ROOT.kAzure+9)
hcX0gen.SetMarkerStyle(24)

hmXIabs.Sumw2()
hmXIabs.SetLineWidth(2)
hmXIabs.SetLineColor(ROOT.kRed)

hmX0abs.Sumw2()
hmX0abs.SetLineWidth(2)
hmX0abs.SetLineColor(ROOT.kBlue)
hmX0abs.SetFillColor(ROOT.kBlue)

hmX0absGen.Sumw2()
hmX0absGen.SetLineWidth(2)
hmX0absGen.SetLineColor(ROOT.kGreen)

hmXIabsStd.Sumw2()
hmXIabsStd.SetLineWidth(2)
hmXIabsStd.SetLineColor(ROOT.kRed)

hmX0absStd.Sumw2()
hmX0absStd.SetLineWidth(2)
hmX0absStd.SetLineColor(ROOT.kBlue)
hmX0absStd.SetFillColor(ROOT.kBlue)

hmX0absStdGen.Sumw2()
hmX0absStdGen.SetLineWidth(2)
hmX0absStdGen.SetLineColor(ROOT.kBlue)
hmX0absStdGen.SetFillColor(ROOT.kBlue)


#############################
lumiData = 1 # 20.3 fb-1
sqrts = 13 #or 8 TeV
docuts = False
#############################

gROOT.ProcessLine(".L Loader.C+")

fpath = "/Users/hod/MC/Pythia/pythia8215/examples/"
fSM = TFile(fpath+"tops.SM_nobmass_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.SM_nobmass.root","READ")
tSM = fSM.Get("SM_nobmass_8TeV") if(sqrts==8) else fSM.Get("SM_nobmass")

# fXX = TFile(fpath+"tops.SMIA_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.SMIA.root","READ")
# tXX = fXX.Get("SMIA_8TeV") if(sqrts==8) else fXX.Get("SMIA")
fXX = TFile(fpath+"tops.A_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.A.root","READ")
tXX = fXX.Get("A_8TeV") if(sqrts==8) else fXX.Get("A")

fX0 = TFile(fpath+"tops.A_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.A.root","READ")
tX0 = fX0.Get("A_8TeV") if(sqrts==8) else fX0.Get("A")

nkevents = tXX.GetEntries()/1000
legcme = "#sqrt{#it{s}} = 8 TeV" if("8TeV" in fSM.GetName()) else "#sqrt{#it{s}} = 13 TeV"
scme = "8TeV" if("8TeV" in fSM.GetName()) else "13TeV"
leglumi = "1 fb^{-1}" # "20.3 fb^{-1}"


n=1
sumofweights_X0 = 0
sumofweights_SM = 0
for event in tSM:
   if(n%10000==0): print "processed |SM|^2 and reweighting ", n
   chain = initdecaychain()
   chian = setdecaychain(chain,event.id,event.parents,event.children)
   # printdcaychain(chian,event.id,event.parents,event.children)
   g1 = event.p4[chain["t"]["g1"]]
   g2 = event.p4[chain["t"]["g2"]]
   t1 = event.p4[chain["t"]["t"]]
   t2 = event.p4[chain["t~"]["t~"]]
   b1 = event.p4[chain["t"]["b"]]
   b2 = event.p4[chain["t~"]["b~"]]
   w1 = event.p4[chain["t"]["w+"]]
   w2 = event.p4[chain["t~"]["w-"]]
   fd1 = event.p4[chain["t"]["fd"]]
   fd2 = event.p4[chain["t~"]["fd"]]
   fu1 = event.p4[chain["t"]["fu"]]
   fu2 = event.p4[chain["t~"]["fu"]]

   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   p = [[ g1.E(), g1.Px(), g1.Py(), g1.Pz() ],
        [ g2.E(), g2.Px(), g2.Py(), g2.Pz() ],
        [ t1.E(), t1.Px(), t1.Py(), t1.Pz() ],
        [ t2.E(), t2.Px(), t2.Py(), t2.Pz() ]]
   P=THDM.invert_momenta(p)

   # Q = 0.5*(math.sqrt(t1.M()*t1.M()+t1.Pt()*t1.Pt())+math.sqrt(t2.M()*t2.M()+t2.Pt()*t2.Pt()))
   alphaS = event.aS #THDM.model.AlphaS.alphasQ(Q)

   ## the ME^2 and the weight
   # me2SM = THDM.modules['matrix2SMpy'].get_me(P,alphaS,nhel)                           ### calculate the SM ME^2
   # me2XX = THDM.modules['matrix2'+nameX+str(index)+'py'].get_me(P,alphaS,nhel)         ### calculate the X ME^2
   me2SM = THDM.modules['matrix2SMttxpy'].get_me(P,alphaS,nhel)                          ### calculate the SM ME^2
   me2XX = THDM.modules['matrix2'+nameX+str(index)+'ttxpy'].get_me(P,alphaS,nhel)        ### calculate the X ME^2
   me2X0 = THDM.modules['matrix2'+nameX+"only"+str(index)+'ttxpy'].get_me(P,alphaS,nhel) ### calculate the X ME^2
   weightXX = me2XX/me2SM                                                                 ### calculate the weight
   weightX0 = me2X0/me2SM                                                                 ### calculate the weight

   hcSMgen.Fill(cosphi)
   hmSMgen.Fill(mtt)

   # hcXXrwt.Fill(cosphi,weightXX)
   hmXXrwt.Fill(mtt,weightXX)

   hcX0rwt.Fill(cosphi,weightX0)
   hmX0rwt.Fill(mtt,weightX0)

   hmXIrwt.Fill(mtt,weightXX-1)
   hmXI0rwt.Fill(mtt,weightX0-1)
   hmXIabs.Fill(mtt,weightXX-1)
   hmX0abs.Fill(mtt,weightX0)
   hmXIabsStd.Fill(mtt,weightXX-1)
   hmX0absStd.Fill(mtt,weightX0)

   sumofweights_SM += 1
   sumofweights_X0 += weightX0
 
   n+=1


n=1
for event in tXX:
   if(n%10000==0): print "processed |SM+X|^2 generated ", n
   # t1=event.p4[2]
   # t2=event.p4[3]

   chain = initdecaychain()
   chian = setdecaychain(chain,event.id,event.parents,event.children)
   # printdcaychain(chian,event.id,event.parents,event.children)
   g1 = event.p4[chain["t"]["g1"]]
   g2 = event.p4[chain["t"]["g2"]]
   t1 = event.p4[chain["t"]["t"]]
   t2 = event.p4[chain["t~"]["t~"]]
   b1 = event.p4[chain["t"]["b"]]
   b2 = event.p4[chain["t~"]["b~"]]
   w1 = event.p4[chain["t"]["w+"]]
   w2 = event.p4[chain["t~"]["w-"]]
   fd1 = event.p4[chain["t"]["fd"]]
   fd2 = event.p4[chain["t~"]["fd"]]
   fu1 = event.p4[chain["t"]["fu"]]
   fu2 = event.p4[chain["t~"]["fu"]]

   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   # hcXXgen.Fill(cosphi)
   hmXXgen.Fill(mtt)

   n+=1

n=1
for event in tX0:
   if(n%10000==0): print "processed |X|^2 generated ", n
   # t1=event.p4[2]
   # t2=event.p4[3]

   chain = initdecaychain()
   chian = setdecaychain(chain,event.id,event.parents,event.children)
   # printdcaychain(chian,event.id,event.parents,event.children)
   g1 = event.p4[chain["t"]["g1"]]
   g2 = event.p4[chain["t"]["g2"]]
   t1 = event.p4[chain["t"]["t"]]
   t2 = event.p4[chain["t~"]["t~"]]
   b1 = event.p4[chain["t"]["b"]]
   b2 = event.p4[chain["t~"]["b~"]]
   w1 = event.p4[chain["t"]["w+"]]
   w2 = event.p4[chain["t~"]["w-"]]
   fd1 = event.p4[chain["t"]["fd"]]
   fd2 = event.p4[chain["t~"]["fd"]]
   fu1 = event.p4[chain["t"]["fu"]]
   fu2 = event.p4[chain["t~"]["fu"]]

   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   hcX0gen.Fill(cosphi)
   hmX0gen.Fill(mtt)

   hmX0absGen.Fill(mtt)
   hmX0absStdGen.Fill(mtt)

   n+=1


hmSMgen = normalizeToBinWidth(hmSMgen,40)
hmXXgen = normalizeToBinWidth(hmXXgen,40)
hmX0gen = normalizeToBinWidth(hmX0gen,40)
hmXXrwt = normalizeToBinWidth(hmXXrwt,40)
hmX0rwt = normalizeToBinWidth(hmX0rwt,40)
hmXIrwt = normalizeToBinWidth(hmXIrwt,40)
hmXI0rwt = normalizeToBinWidth(hmXI0rwt,40)
# hmXIabs = normalizeToBinWidth(hmXIabs)
# hmX0abs = normalizeToBinWidth(hmX0abs)
# hmX0absGen = normalizeToBinWidth(hmX0absGen)
hmXIabsStd = normalizeToBinWidth(hmXIabsStd,40)
hmX0absStd = normalizeToBinWidth(hmX0absStd,40)
hmX0absStdGen = normalizeToBinWidth(hmX0absStdGen,40)

hcSMgen.Scale(1./hcSMgen.Integral())
# hcXXgen.Scale(1./hcXXgen.Integral())
hcX0gen.Scale(1./hcX0gen.Integral())
# hcXXrwt.Scale(1./hcXXrwt.Integral())
hcX0rwt.Scale(1./hcX0rwt.Integral())


mb2fb = 1.e12
pb2fb = 1.e3
sigmaSM8TeV = 1.285e-07*mb2fb  # 253*pb2fb
sigmaSM13TeV = 4.458e-07*mb2fb # 831*pb2fb
sigmaXonly8TeV = 3.978e-09*mb2fb
sigmaSM = sigmaSM8TeV if("8TeV" in fSM.GetName()) else sigmaSM13TeV
neventsSM = tXX.GetEntries()
neventsXonly = tX0.GetEntries()
lumiSM = neventsSM/sigmaSM
lumiXonly = neventsXonly/sigmaXonly8TeV

hmSMgen.Scale(lumiData/lumiSM)
hmXXgen.Scale(lumiData/lumiSM)
hmX0gen.Scale(lumiData/lumiXonly)
hmXXrwt.Scale(lumiData/lumiSM)
hmX0rwt.Scale(lumiData/lumiSM)

print "hmX0gen.Integral()/hmX0rwt.Integral() = ",hmX0gen.Integral()/hmX0rwt.Integral()
#########################################################################
#########################################################################
factor = 1/math.sqrt(2) #hmX0gen.Integral()/hmX0rwt.Integral() ##########
#########################################################################
#########################################################################
if("8TeV" in fSM.GetName()):
   xsecReweighted = (sumofweights_X0*factor/sumofweights_SM)*sigmaSM
   print "sum of weights X0 = ",sumofweights_X0*factor
   print "sum of weights SM = ",sumofweights_SM
   print "cross section X0 reweighted = ",xsecReweighted
   print "cross section X0 generated  = ",sigmaXonly8TeV
   print "difference in percent  = ",(sigmaXonly8TeV-xsecReweighted)/sigmaXonly8TeV*100


hmX0rwt.Scale(factor)
hmXIrwt.Scale(lumiData/lumiSM)
hmXI0rwt.Scale(lumiData/lumiSM)
hmXIabs.Scale(lumiData/lumiSM)
hmX0abs.Scale(lumiData/lumiSM * factor)
hmX0absGen.Scale(lumiData/lumiXonly)
hmXIabsStd.Scale(lumiData/lumiSM)
hmX0absStd.Scale(lumiData/lumiSM * factor)
hmX0absStdGen.Scale(lumiData/lumiXonly)



tanb = '%.2f' % THDM.parameters[index].get("tanb")
wXprcnt = '%.1f' % (THDM.parameters[index].get("w"+nameX)/mX*100.)
stanb = tanb.replace(".","")
pdfname = "2HDM.reweighting.sqrts"+scme+"."+nameX+"."+str(mX)+"GeV.tanb"+stanb+".pdf"



plot(hmSMgen,hmXXrwt,hmXXgen,hmXIrwt,tanb,wXprcnt,nkevents,legcme,leglumi,pdfname+"(")
plot(hmSMgen,hmX0rwt,hmX0gen,hmXI0rwt,tanb,wXprcnt,nkevents,legcme,leglumi,pdfname)



legabs = TLegend(0.5,0.5,0.87,0.9,"","brNDC")
legabs.SetFillStyle(4000); # will be transparent
legabs.SetFillColor(0)
legabs.SetTextFont(42)
legabs.SetBorderSize(0)
legabs.AddEntry(0, "MadGraph+Pythia8", "")
legabs.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} ("+str(nkevents)+"k events)", "")
legabs.AddEntry(0, legcme+", "+leglumi, "")
legabs.AddEntry(hmXIabs,"|SM+#it{"+nameX+"}|^{2}-|SM|^{2} reweighted","l")
legabs.AddEntry(0, "sin(#beta-#alpha)=1", "")
legabs.AddEntry(0, "tan#beta="+str(tanb), "")
legabs.AddEntry(0, "#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV", "")
legabs.AddEntry(0, "#Gamma_{#it{"+nameX+"}}/#it{m}_{#it{"+nameX+"}}="+str(wXprcnt)+" [%]", "")


if(sqrts==13):
   hmXIabs.SetMinimum(hmXIabs.GetMinimum()*1.1)
   hmXIabs.SetMaximum(hmXIabs.GetMaximum()*2.3)
   hmX0abs.SetMinimum(hmXIabs.GetMinimum())
   hmX0abs.SetMaximum(hmXIabs.GetMaximum())
   hmX0absGen.SetMinimum(hmXIabs.GetMinimum())
   hmX0absGen.SetMaximum(hmXIabs.GetMaximum())
if(sqrts==8):
   hmXIabs.SetMinimum(-5000)
   hmXIabs.SetMaximum(+21500)
   hmX0abs.SetMinimum(-5000)
   hmX0abs.SetMaximum(+21500)
   hmX0absGen.SetMinimum(-5000)
   hmX0absGen.SetMaximum(+21500)
cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hmXIabs.Draw("hist")
hmX0abs.Draw("hist same")
hmXIabs.Draw("hist same")
hmX0absGen.Draw("hist same")
legabs.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".abs1.pdf"))
cnv.SaveAs(pdfname)

hmXIabsStd.SetMinimum(hmXIabsStd.GetMinimum()*1.1)
hmXIabsStd.SetMaximum(hmXIabsStd.GetMaximum()*2.3)
hmX0absStd.SetMinimum(hmXIabsStd.GetMinimum())
hmX0absStd.SetMaximum(hmXIabsStd.GetMaximum())
cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hmXIabsStd.Draw("hist")
hmX0absStd.Draw("hist same")
hmXIabsStd.Draw("hist same")
legabs.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".abs2.pdf"))
cnv.SaveAs(pdfname)


legcosphi = TLegend(0.15,0.5,0.4,0.9,"","brNDC")
legcosphi.SetFillStyle(4000); # will be transparent
legcosphi.SetFillColor(0)
legcosphi.SetTextFont(42)
legcosphi.SetBorderSize(0)
legcosphi.AddEntry(0, "MadGraph+Pythia8", "")
legcosphi.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} ("+str(nkevents)+"k events)", "")
legcosphi.AddEntry(0, legcme+", "+leglumi, "")
legcosphi.AddEntry(hcSMgen,"SM","ple")
legcosphi.AddEntry(hcX0rwt,"#it{"+nameX+"} reweighted","ple")
legcosphi.AddEntry(hcX0gen,"#it{"+nameX+"} generated","ple")
legcosphi.AddEntry(0, "sin(#beta-#alpha)=1", "")
legcosphi.AddEntry(0, "tan#beta="+str(tanb), "")
legcosphi.AddEntry(0, "#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV", "")
legcosphi.AddEntry(0, "#Gamma_{#it{"+nameX+"}}/#it{m}_{#it{"+nameX+"}}="+str(wXprcnt)+" [%]", "")

cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hcSMgen.SetMinimum(0)
hcSMgen.SetMaximum(hcSMgen.GetMaximum()*2)
hcSMgen.Draw("hist")
hcX0gen.Draw("hist same")
hcX0rwt.Draw("hist same")
legcosphi.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".cosphi.pdf"))
cnv.SaveAs(pdfname+")")
