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
            chain["t"]["b"]  = c2
            chain["t"]["fd"] = c11 if(pdgid[c11]%2!=0) else c12
            chain["t"]["fu"] = c11 if(pdgid[c11]%2==0) else c12
         if(pdgid[c2]==24):
            chain["t"]["w+"] = c2
            chain["t"]["b"]  = c1
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
            chain["t~"]["b~"] = c2
            chain["t~"]["fd"] = c11 if(pdgid[c11]%2!=0) else c12
            chain["t~"]["fu"] = c11 if(pdgid[c11]%2==0) else c12
         if(pdgid[c2]==-24):
            chain["t~"]["w1"] = c2
            chain["t~"]["b~"] = c1
            chain["t~"]["fd"] = c21 if(pdgid[c21]%2!=0) else c22
            chain["t~"]["fu"] = c21 if(pdgid[c21]%2==0) else c22
      if(chain["t"]["b"]>=0 and chain["t~"]["b~"]>=0): break
   return chain


def printdcaychain(chian,pdgid,parents,children):
   for j in xrange(pdgid.size()):
      print '[%g] pdgid=%g' % (j,pdgid[j])
      if(parents[j].size()>1):  print '  parent[%g]=%g, parent[%g]=%g' % (parents[j][0],pdgid[parents[j][0]],parents[j][1],pdgid[parents[j][1]])
      if(children[j].size()>1): print '    child[%g]=%g, child[%g]=%g' % (children[j][0],pdgid[children[j][0]],children[j][1],pdgid[children[j][1]])
   print 'chain["t"]["t"]   = %g, pdgid=%g' % (chain["t"]["t"],   pdgid[chain["t"]["t"]])
   print 'chain["t"]["g1"]  = %g, pdgid=%g' % (chain["t"]["g1"],  pdgid[chain["t"]["g1"]])
   print 'chain["t"]["g2"]  = %g, pdgid=%g' % (chain["t"]["g2"],  pdgid[chain["t"]["g2"]])
   print 'chain["t"]["b"]   = %g, pdgid=%g' % (chain["t"]["b"],   pdgid[chain["t"]["b"]])
   print 'chain["t"]["w+"]  = %g, pdgid=%g' % (chain["t"]["w+"],  pdgid[chain["t"]["w+"]])
   print 'chain["t"]["fd"]  = %g, pdgid=%g' % (chain["t"]["fd"],  pdgid[chain["t"]["fd"]])
   print 'chain["t"]["fu"]  = %g, pdgid=%g' % (chain["t"]["fu"],  pdgid[chain["t"]["fu"]])
   print 'chain["t~"]["t~"] = %g, pdgid=%g' % (chain["t~"]["t~"], pdgid[chain["t~"]["t~"]])
   print 'chain["t~"]["g1"] = %g, pdgid=%g' % (chain["t~"]["g1"], pdgid[chain["t~"]["g1"]])
   print 'chain["t~"]["g2"] = %g, pdgid=%g' % (chain["t~"]["g2"], pdgid[chain["t~"]["g1"]])
   print 'chain["t~"]["b~"] = %g, pdgid=%g' % (chain["t~"]["b~"], pdgid[chain["t~"]["b~"]])
   print 'chain["t~"]["w-"] = %g, pdgid=%g' % (chain["t~"]["w-"], pdgid[chain["t~"]["w-"]])
   print 'chain["t~"]["fd"] = %g, pdgid=%g' % (chain["t~"]["fd"], pdgid[chain["t~"]["fd"]])
   print 'chain["t~"]["fu"] = %g, pdgid=%g' % (chain["t~"]["fu"], pdgid[chain["t~"]["fu"]])


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

def plot(h1, h2, h3, h4, tanb, wX, nkevents, cme, lumi, fname, fnamesingle, ymin=-1, ymax=-1):
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
   cnv.SaveAs(fname.replace("(","").replace(")","").replace(".pdf","."+fnamesingle+".pdf"));
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


############### Temporary stuff ###############
matrixpathSM    = "/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/ggtt-SMtest-partonlevel/SubProcesses/P1_gg_ttx_t_bwp_wp_udx_tx_bxwm_wm_mumvmx"
matrixpathAonly = "/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/ggtt-AonlyTopsDecay-partonlevel/SubProcesses/P1_gg_h1_ttx_no_h_t_bwp_wp_udx_tx_bxwm_wm_mumvmx"
matrixpathSMIA  = "/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/ggtt-Atest-partonlevel/SubProcesses/P1_gg_ttx_no_h_t_bwp_wp_udx_tx_bxwm_wm_mumvmx"
sys.path.append(matrixpathAonly)
sys.path.append(matrixpathSMIA)
sys.path.append(matrixpathSM)
import matrixA2py    as matrixX0
import matrixSMIA2py as matrixXX
import matrixSM2py   as matrixSM
matrixX0.initialise(matrixpathAonly+'/../../Cards/param_card.dat')
matrixXX.initialise(matrixpathSMIA+'/../../Cards/param_card.dat')
matrixSM.initialise(matrixpathSM+'/../../Cards/param_card.dat')
###############################################


listbins = [0,80,160,240,320,360,400,440,500,560,600,640,680,720,760,800,860,920,1040,1160,1280]
arrbins = array("d", listbins)
nbins = len(listbins)-1

### mass histos
hmSMgen    = TH1D("hmSMgen",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXXrwt    = TH1D("hmXXrwt",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0rwt    = TH1D("hmX0rwt",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXIrwt    = TH1D("hmXIrwt",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXI0rwt   = TH1D("hmXI0rwt",         ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXXgen    = TH1D("hmXXgen",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0gen    = TH1D("hmX0gen",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX1gen    = TH1D("hmX1gen",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXIabs    = TH1D("hmXIabs",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmX0abs    = TH1D("hmX0abs",          ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmX0absGen = TH1D("hmX0absGen",       ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmXIabsStd = TH1D("hmXIabsStd",       ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0absStd = TH1D("hmX0absStd",       ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmX0absStdGen = TH1D("hmX0absStdGen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)

### costheta histos
haSMgen = TH1D("haSMgen",   ";Truth cos#it{#theta} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
haXXrwt = TH1D("haXXrwt",   ";Truth cos#it{#theta} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
haX0rwt = TH1D("haX0rwt",   ";Truth cos#it{#theta} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
haXXgen = TH1D("haXXgen",   ";Truth cos#it{#theta} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
haX0gen = TH1D("haX0gen",   ";Truth cos#it{#theta} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
haX1gen = TH1D("haX1gen",   ";Truth cos#it{#theta} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)

### cosphi histos
hcSMgen = TH1D("hcSMgen",   ";Truth cos#it{#phi} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
hcXXrwt = TH1D("hcXXrwt",   ";Truth cos#it{#phi} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
hcX0rwt = TH1D("hcX0rwt",   ";Truth cos#it{#phi} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
hcXXgen = TH1D("hcXXgen",   ";Truth cos#it{#phi} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
hcX0gen = TH1D("hcX0gen",   ";Truth cos#it{#phi} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)
hcX1gen = TH1D("hcX1gen",   ";Truth cos#it{#phi} in #it{m}_{#it{t}#bar{#it{t}}}<600 GeV;Normalised",40,-1,+1)

### pt histos
hpSMgen = TH1D("hpSMgen",   ";Truth lepton #it{p}_{T} [GeV] in |#it{#eta}|<2.5;Normalised",50,0,200)
hpXXrwt = TH1D("hpXXrwt",   ";Truth lepton #it{p}_{T} [GeV] in |#it{#eta}|<2.5;Normalised",50,0,200)
hpX0rwt = TH1D("hpX0rwt",   ";Truth lepton #it{p}_{T} [GeV] in |#it{#eta}|<2.5;Normalised",50,0,200)
hpXXgen = TH1D("hpXXgen",   ";Truth lepton #it{p}_{T} [GeV] in |#it{#eta}|<2.5;Normalised",50,0,200)
hpX0gen = TH1D("hpX0gen",   ";Truth lepton #it{p}_{T} [GeV] in |#it{#eta}|<2.5;Normalised",50,0,200)
hpX1gen = TH1D("hpX1gen",   ";Truth lepton #it{p}_{T} [GeV] in |#it{#eta}|<2.5;Normalised",50,0,200)

### eta histos
htSMgen = TH1D("htSMgen",   ";Truth lepton #it{#eta};Normalised",35,-4.5,+4.5)
htXXrwt = TH1D("htXXrwt",   ";Truth lepton #it{#eta};Normalised",35,-4.5,+4.5)
htX0rwt = TH1D("htX0rwt",   ";Truth lepton #it{#eta};Normalised",35,-4.5,+4.5)
htXXgen = TH1D("htXXgen",   ";Truth lepton #it{#eta};Normalised",35,-4.5,+4.5)
htX0gen = TH1D("htX0gen",   ";Truth lepton #it{#eta};Normalised",35,-4.5,+4.5)
htX1gen = TH1D("htX1gen",   ";Truth lepton #it{#eta};Normalised",35,-4.5,+4.5)



hmSMgen.Sumw2()
hmSMgen.SetLineWidth(2)
hmSMgen.SetLineColor(ROOT.kBlack)
hmSMgen.SetMarkerColor(ROOT.kBlack)
hmSMgen.SetMarkerStyle(24)
haSMgen.Sumw2()
haSMgen.SetLineWidth(2)
haSMgen.SetLineColor(ROOT.kBlack)
haSMgen.SetMarkerColor(ROOT.kBlack)
haSMgen.SetMarkerStyle(24)
hcSMgen.Sumw2()
hcSMgen.SetLineWidth(2)
hcSMgen.SetLineColor(ROOT.kBlack)
hcSMgen.SetMarkerColor(ROOT.kBlack)
hcSMgen.SetMarkerStyle(24)
hpSMgen.Sumw2()
hpSMgen.SetLineWidth(2)
hpSMgen.SetLineColor(ROOT.kBlack)
hpSMgen.SetMarkerColor(ROOT.kBlack)
hpSMgen.SetMarkerStyle(24)
htSMgen.Sumw2()
htSMgen.SetLineWidth(2)
htSMgen.SetLineColor(ROOT.kBlack)
htSMgen.SetMarkerColor(ROOT.kBlack)
htSMgen.SetMarkerStyle(24)

hmXXrwt.Sumw2()
hmXXrwt.SetLineWidth(2)
hmXXrwt.SetLineColor(ROOT.kRed)
hmXXrwt.SetMarkerColor(ROOT.kRed)
hmXXrwt.SetMarkerStyle(20)
haXXrwt.Sumw2()
haXXrwt.SetLineWidth(2)
haXXrwt.SetLineColor(ROOT.kRed)
haXXrwt.SetLineStyle(3)
haXXrwt.SetMarkerColor(ROOT.kRed)
haXXrwt.SetMarkerStyle(20)
hcXXrwt.Sumw2()
hcXXrwt.SetLineWidth(2)
hcXXrwt.SetLineColor(ROOT.kRed)
hcXXrwt.SetLineStyle(3)
hcXXrwt.SetMarkerColor(ROOT.kRed)
hcXXrwt.SetMarkerStyle(20)
hpXXrwt.Sumw2()
hpXXrwt.SetLineWidth(2)
hpXXrwt.SetLineColor(ROOT.kRed)
hpXXrwt.SetLineStyle(2)
hpXXrwt.SetMarkerColor(ROOT.kRed)
hpXXrwt.SetMarkerStyle(20)
htXXrwt.Sumw2()
htXXrwt.SetLineWidth(2)
htXXrwt.SetLineColor(ROOT.kRed)
htXXrwt.SetLineStyle(2)
htXXrwt.SetMarkerColor(ROOT.kRed)
htXXrwt.SetMarkerStyle(20)

hmX0rwt.Sumw2()
hmX0rwt.SetLineWidth(2)
hmX0rwt.SetLineColor(ROOT.kRed)
hmX0rwt.SetMarkerColor(ROOT.kRed)
hmX0rwt.SetMarkerStyle(20)
haX0rwt.Sumw2()
haX0rwt.SetLineWidth(2)
haX0rwt.SetLineColor(ROOT.kRed)
haX0rwt.SetMarkerColor(ROOT.kRed)
haX0rwt.SetMarkerStyle(20)
hcX0rwt.Sumw2()
hcX0rwt.SetLineWidth(2)
hcX0rwt.SetLineColor(ROOT.kRed)
hcX0rwt.SetMarkerColor(ROOT.kRed)
hcX0rwt.SetMarkerStyle(20)
hpX0rwt.Sumw2()
hpX0rwt.SetLineWidth(2)
hpX0rwt.SetLineColor(ROOT.kRed)
hpX0rwt.SetMarkerColor(ROOT.kRed)
hpX0rwt.SetMarkerStyle(20)
htX0rwt.Sumw2()
htX0rwt.SetLineWidth(2)
htX0rwt.SetLineColor(ROOT.kRed)
htX0rwt.SetMarkerColor(ROOT.kRed)
htX0rwt.SetMarkerStyle(20)

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
haXXgen.Sumw2()
haXXgen.SetLineWidth(2)
haXXgen.SetLineColor(ROOT.kAzure+9)
haXXgen.SetLineStyle(2)
haXXgen.SetMarkerColor(ROOT.kAzure+9)
haXXgen.SetMarkerStyle(24)
hcXXgen.Sumw2()
hcXXgen.SetLineWidth(2)
hcXXgen.SetLineColor(ROOT.kAzure+9)
hcXXgen.SetLineStyle(2)
hcXXgen.SetMarkerColor(ROOT.kAzure+9)
hcXXgen.SetMarkerStyle(24)
hpXXgen.Sumw2()
hpXXgen.SetLineWidth(2)
hpXXgen.SetLineColor(ROOT.kAzure+9)
hpXXgen.SetLineStyle(2)
hpXXgen.SetMarkerColor(ROOT.kAzure+9)
hpXXgen.SetMarkerStyle(24)
htXXgen.Sumw2()
htXXgen.SetLineWidth(2)
htXXgen.SetLineColor(ROOT.kAzure+9)
htXXgen.SetLineStyle(2)
htXXgen.SetMarkerColor(ROOT.kAzure+9)
htXXgen.SetMarkerStyle(24)

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
haX0gen.Sumw2()
haX0gen.SetLineWidth(2)
haX0gen.SetLineColor(ROOT.kAzure+9)
haX0gen.SetMarkerColor(ROOT.kAzure+9)
haX0gen.SetMarkerStyle(24)
hpX0gen.Sumw2()
hpX0gen.SetLineWidth(2)
hpX0gen.SetLineColor(ROOT.kAzure+9)
hpX0gen.SetMarkerColor(ROOT.kAzure+9)
hpX0gen.SetMarkerStyle(24)
htX0gen.Sumw2()
htX0gen.SetLineWidth(2)
htX0gen.SetLineColor(ROOT.kAzure+9)
htX0gen.SetMarkerColor(ROOT.kAzure+9)
htX0gen.SetMarkerStyle(24)

hmX1gen.Sumw2()
hmX1gen.SetLineWidth(2)
hmX1gen.SetLineColor(ROOT.kGreen+1)
hmX1gen.SetMarkerColor(ROOT.kGreen+1)
hmX1gen.SetMarkerStyle(24)
haX1gen.Sumw2()
haX1gen.SetLineWidth(2)
haX1gen.SetLineColor(ROOT.kGreen+1)
haX1gen.SetMarkerColor(ROOT.kGreen+1)
haX1gen.SetMarkerStyle(24)
hcX1gen.Sumw2()
hcX1gen.SetLineWidth(2)
hcX1gen.SetLineColor(ROOT.kGreen+1)
hcX1gen.SetMarkerColor(ROOT.kGreen+1)
hcX1gen.SetMarkerStyle(24)
hpX1gen.Sumw2()
hpX1gen.SetLineWidth(2)
hpX1gen.SetLineColor(ROOT.kGreen+1)
hpX1gen.SetMarkerColor(ROOT.kGreen+1)
hpX1gen.SetMarkerStyle(24)
htX1gen.Sumw2()
htX1gen.SetLineWidth(2)
htX1gen.SetLineColor(ROOT.kGreen+1)
htX1gen.SetMarkerColor(ROOT.kGreen+1)
htX1gen.SetMarkerStyle(24)

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

fXX = TFile(fpath+"tops.SMIA_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.SMIA.root","READ")
tXX = fXX.Get("SMIA_8TeV") if(sqrts==8) else fXX.Get("SMIA")

fX0 = TFile(fpath+"tops.A_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.A.root","READ")
tX0 = fX0.Get("A_8TeV") if(sqrts==8) else fX0.Get("A")

fX1 = TFile(fpath+"tops.H.root","READ")
tX1 = fX1.Get("H")

nkevents = tXX.GetEntries()/1000
legcme = "#sqrt{#it{s}} = 8 TeV" if("8TeV" in fSM.GetName()) else "#sqrt{#it{s}} = 13 TeV"
scme = "8TeV" if("8TeV" in fSM.GetName()) else "13TeV"
leglumi = "1 fb^{-1}" # "20.3 fb^{-1}"


def tos(p):
   return "PxPyPzE=(%.2f,%.2f,%.2f,%.2f)" % (p.Px(),p.Py(),p.Pz(),p.E())

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

   costheta = kin.costheta(t1+t2,t1)
   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()
   # cosphi = kin.cosphi(b1+fd1+fu1,fd1,b2+fd2+fu2,fd2)
   # mtt = (b1+fd1+fu1+b2+fd2+fu2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   p = [[ g1.E(), g1.Px(), g1.Py(), g1.Pz() ],
        [ g2.E(), g2.Px(), g2.Py(), g2.Pz() ],
        [ t1.E(), t1.Px(), t1.Py(), t1.Pz() ],
        [ t2.E(), t2.Px(), t2.Py(), t2.Pz() ]]
   P=THDM.invert_momenta(p)

   ############### Temporary stuff ###############
   # gg_h_ttx_t_bwp_wp_udx_tx_bxwm_wm_mumvmx
   # gg_h1_ttx_no_h_t_bwp_wp_udx_tx_bxwm_wm_mumvmx
   p1 = [[ g1.E(), g1.Px(), g1.Py(), g1.Pz()     ],
         [ g2.E(), g2.Px(), g2.Py(), g2.Pz()     ],
         [ b1.E(), b1.Px(), b1.Py(), b1.Pz()     ],
         [ fd1.E(), fd1.Px(), fd1.Py(), fd1.Pz() ],
         [ fu1.E(), fu1.Px(), fu1.Py(), fu1.Pz() ],
         [ b2.E(), b2.Px(), b2.Py(), b2.Pz()     ],
         [ fd2.E(), fd2.Px(), fd2.Py(), fd2.Pz() ],
         [ fu2.E(), fu2.Px(), fu2.Py(), fu2.Pz() ]]
   P1=THDM.invert_momenta(p1)
   ################################################

   # Q = 0.5*(math.sqrt(t1.M()*t1.M()+t1.Pt()*t1.Pt())+math.sqrt(t2.M()*t2.M()+t2.Pt()*t2.Pt()))
   alphaS = event.aS #THDM.model.AlphaS.alphasQ(Q)

   ## the ME^2 and the weight
   me2SM = THDM.modules['matrix2SMttxpy'].get_me(P,alphaS,nhel)                          ### calculate the SM ME^2
   # me2SM = matrixSM.get_me(P1,alphaS,nhel)
   me2XX = THDM.modules['matrix2'+nameX+str(index)+'ttxpy'].get_me(P,alphaS,nhel)        ### calculate the X ME^2
   # me2XX = matrixXX.get_me(P1,alphaS,nhel)
   me2X0 = THDM.modules['matrix2'+nameX+"only"+str(index)+'ttxpy'].get_me(P,alphaS,nhel) ### calculate the X ME^2
   # me2X0 = matrixX0.get_me(P1,alphaS,nhel)

   # print "me2SM=%g, me2SMx=%g | me2XX=%g, me2XXx=%g | me2X0=%g, me2X0x=%g" % (me2SM, me2SMx, me2XX, me2XXx, me2X0, me2X0x)

   weightXX = me2XX/me2SM                                                                 ### calculate the weight
   weightX0 = me2X0/me2SM                                                                 ### calculate the weight

   htSMgen.Fill(fd2.Eta())
   if(abs(fd2.Eta())<2.5): hpSMgen.Fill(fd2.Pt())
   if(mtt<600): haSMgen.Fill(costheta)
   if(mtt<600): hcSMgen.Fill(cosphi)
   hmSMgen.Fill(mtt)

   htXXrwt.Fill(fd2.Eta(),weightXX)
   if(abs(fd2.Eta())<2.5): hpXXrwt.Fill(fd2.Pt(),weightXX)
   if(mtt<600): haXXrwt.Fill(costheta,weightXX)
   if(mtt<600): hcXXrwt.Fill(cosphi,weightXX)
   hmXXrwt.Fill(mtt,weightXX)

   htX0rwt.Fill(fd2.Eta(),weightX0)
   if(abs(fd2.Eta())<2.5): hpX0rwt.Fill(fd2.Pt(),weightX0)
   if(mtt<600): haX0rwt.Fill(costheta,weightX0)
   if(mtt<600): hcX0rwt.Fill(cosphi,weightX0)
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

   costheta = kin.costheta(t1+t2,t1)
   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()
   # cosphi = kin.cosphi(b1+fd1+fu1,fd1,b2+fd2+fu2,fd2)
   # mtt = (b1+fd1+fu1+b2+fd2+fu2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   htXXgen.Fill(fd2.Eta())
   if(abs(fd2.Eta())<2.5): hpXXgen.Fill(fd2.Pt())
   if(mtt<600): haXXgen.Fill(costheta)
   if(mtt<600): hcXXgen.Fill(cosphi)
   hmXXgen.Fill(mtt)

   n+=1

n=1
for event in tX0:
   if(n%10000==0): print "processed |X|^2 generated ", n
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

   costheta = kin.costheta(t1+t2,t1)
   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()
   # cosphi = kin.cosphi(b1+fd1+fu1,fd1,b2+fd2+fu2,fd2)
   # mtt = (b1+fd1+fu1+b2+fd2+fu2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   htX0gen.Fill(fd2.Eta())
   if(abs(fd2.Eta())<2.5): hpX0gen.Fill(fd2.Pt())
   if(mtt<600): haX0gen.Fill(costheta)
   if(mtt<600): hcX0gen.Fill(cosphi)
   hmX0gen.Fill(mtt)

   hmX0absGen.Fill(mtt)
   hmX0absStdGen.Fill(mtt)

   n+=1


n=1
for event in tX1:
   if(n%10000==0): print "processed |X|^2 generated ", n
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

   costheta = kin.costheta(t1+t2,t1)
   cosphi = kin.cosphi(t1,fd1,t2,fd2)
   mtt = (t1+t2).M()
   # cosphi = kin.cosphi(b1+fd1+fu1,fd1,b2+fd2+fu2,fd2)
   # mtt = (b1+fd1+fu1+b2+fd2+fu2).M()

   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue

   htX1gen.Fill(fd2.Eta())
   if(abs(fd2.Eta())<2.5): hpX1gen.Fill(fd2.Pt())
   if(mtt<600): haX1gen.Fill(costheta)
   if(mtt<600): hcX1gen.Fill(cosphi)
   hmX1gen.Fill(mtt)

   n+=1


## clone before normalising
hcrXXgen = hcXXgen.Clone("hcrXXgen")
hcrXXgen.Divide(hcSMgen)
harXXgen = haXXgen.Clone("harXXgen")
harXXgen.Divide(haSMgen)


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

haSMgen.Scale(1./haSMgen.Integral())
haXXgen.Scale(1./haXXgen.Integral())
haX0gen.Scale(1./haX0gen.Integral())
haX1gen.Scale(1./haX1gen.Integral())
haXXrwt.Scale(1./haXXrwt.Integral())
haX0rwt.Scale(1./haX0rwt.Integral())

hcSMgen.Scale(1./hcSMgen.Integral())
hcXXgen.Scale(1./hcXXgen.Integral())
hcX0gen.Scale(1./hcX0gen.Integral())
hcX1gen.Scale(1./hcX1gen.Integral())
hcXXrwt.Scale(1./hcXXrwt.Integral())
hcX0rwt.Scale(1./hcX0rwt.Integral())

hpSMgen.Scale(1./hpSMgen.Integral())
hpXXgen.Scale(1./hpXXgen.Integral())
hpX0gen.Scale(1./hpX0gen.Integral())
hpX1gen.Scale(1./hpX1gen.Integral())
hpXXrwt.Scale(1./hpXXrwt.Integral())
hpX0rwt.Scale(1./hpX0rwt.Integral())

htSMgen.Scale(1./htSMgen.Integral())
htXXgen.Scale(1./htXXgen.Integral())
htX0gen.Scale(1./htX0gen.Integral())
htX1gen.Scale(1./htX1gen.Integral())
htXXrwt.Scale(1./htXXrwt.Integral())
htX0rwt.Scale(1./htX0rwt.Integral())


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



plot(hmSMgen,hmXXrwt,hmXXgen,hmXIrwt,tanb,wXprcnt,nkevents,legcme,leglumi,pdfname+"(", "mttFull")
plot(hmSMgen,hmX0rwt,hmX0gen,hmXI0rwt,tanb,wXprcnt,nkevents,legcme,leglumi,pdfname,"mttRes")



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





legcostheta = TLegend(0.30,0.40,0.75,0.92,"","brNDC")
legcostheta.SetFillStyle(4000); # will be transparent
legcostheta.SetFillColor(0)
legcostheta.SetTextFont(42)
legcostheta.SetBorderSize(0)
legcostheta.AddEntry(0, "MadGraph+Pythia8", "")
legcostheta.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}}#rightarrow#it{b}#it{W}^{+}(#it{u}#bar{#it{d}})#bar{#it{b}}#it{W}^{-}(#it{#mu}^{-}#bar{#it{#nu}}_{#mu})", "")
legcostheta.AddEntry(0, legcme+" ("+str(nkevents)+"k events)", "")
legcostheta.AddEntry(haSMgen,"SM (full range)","l")
legcostheta.AddEntry(haXXgen,"SM+#it{"+nameX+"} generated","l")
legcostheta.AddEntry(haXXrwt,"SM+#it{"+nameX+"} reweighted","l")
legcostheta.AddEntry(haX0gen,"#it{"+nameX+"} (only) generated","l")
legcostheta.AddEntry(haX0rwt,"#it{"+nameX+"} (only) reweighted","l")
legcostheta.AddEntry(haX1gen,"#it{H} (only) generated","l")
legcostheta.AddEntry(0, "sin(#beta-#alpha)=1", "")
legcostheta.AddEntry(0, "tan#beta="+str(tanb), "")
legcostheta.AddEntry(0, "#it{m}_{#it{X}}="+str(mX)+" GeV", "")

cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
haSMgen.SetMinimum(0)
haSMgen.SetMaximum(haSMgen.GetMaximum()*2.5)
haSMgen.Draw("hist")
haXXgen.Draw("hist same")
haXXrwt.Draw("hist same")
haX0gen.Draw("hist same")
haX1gen.Draw("hist same")
haX0rwt.Draw("hist same")
legcostheta.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".costheta.pdf"))
cnv.SaveAs(pdfname)







legcosphi = TLegend(0.30,0.40,0.75,0.92,"","brNDC")
legcosphi.SetFillStyle(4000); # will be transparent
legcosphi.SetFillColor(0)
legcosphi.SetTextFont(42)
legcosphi.SetBorderSize(0)
legcosphi.AddEntry(0, "MadGraph+Pythia8", "")
legcosphi.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}}#rightarrow#it{b}#it{W}^{+}(#it{u}#bar{#it{d}})#bar{#it{b}}#it{W}^{-}(#it{#mu}^{-}#bar{#it{#nu}}_{#mu})", "")
legcosphi.AddEntry(0, legcme+" ("+str(nkevents)+"k events)", "")
legcosphi.AddEntry(hcSMgen,"SM (full range)","l")
legcosphi.AddEntry(hcXXgen,"SM+#it{"+nameX+"} generated","l")
legcosphi.AddEntry(hcXXrwt,"SM+#it{"+nameX+"} reweighted","l")
legcosphi.AddEntry(hcX0gen,"#it{"+nameX+"} (only) generated","l")
legcosphi.AddEntry(hcX0rwt,"#it{"+nameX+"} (only) reweighted","l")
legcosphi.AddEntry(hcX1gen,"#it{H} (only) generated","l")
legcosphi.AddEntry(0, "sin(#beta-#alpha)=1", "")
legcosphi.AddEntry(0, "tan#beta="+str(tanb), "")
legcosphi.AddEntry(0, "#it{m}_{#it{X}}="+str(mX)+" GeV", "")

cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hcSMgen.SetMinimum(0)
hcSMgen.SetMaximum(hcSMgen.GetMaximum()*2.5)
hcSMgen.Draw("hist")
hcXXgen.Draw("hist same")
hcXXrwt.Draw("hist same")
hcX0gen.Draw("hist same")
hcX1gen.Draw("hist same")
hcX0rwt.Draw("hist same")
legcosphi.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".cosphi.pdf"))
cnv.SaveAs(pdfname)


legpt = TLegend(0.45,0.40,0.90,0.92,"","brNDC")
legpt.SetFillStyle(4000); # will be transparent
legpt.SetFillColor(0)
legpt.SetTextFont(42)
legpt.SetBorderSize(0)
legpt.AddEntry(0, "MadGraph+Pythia8", "")
legpt.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}}#rightarrow#it{b}#it{W}^{+}(#it{u}#bar{#it{d}})#bar{#it{b}}#it{W}^{-}(#it{#mu}^{-}#bar{#it{#nu}}_{#mu})", "")
legpt.AddEntry(0, legcme+" ("+str(nkevents)+"k events)", "")
legpt.AddEntry(hpSMgen,"SM (full range)","l")
legpt.AddEntry(hpXXgen,"SM+#it{"+nameX+"} generated","l")
legpt.AddEntry(hpXXrwt,"SM+#it{"+nameX+"} reweighted","l")
legpt.AddEntry(hpX0gen,"#it{"+nameX+"} (only) generated","l")
legpt.AddEntry(hpX0rwt,"#it{"+nameX+"} (only) reweighted","l")
legpt.AddEntry(hpX1gen,"#it{H} (only) generated","l")
legpt.AddEntry(0, "sin(#beta-#alpha)=1", "")
legpt.AddEntry(0, "tan#beta="+str(tanb), "")
legpt.AddEntry(0, "#it{m}_{#it{X}}="+str(mX)+" GeV", "")

cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hpSMgen.SetMinimum(0)
hpSMgen.SetMaximum(hpSMgen.GetMaximum()*2)
hpSMgen.Draw("hist")
hpXXgen.Draw("hist same")
hpXXrwt.Draw("hist same")
hpX0gen.Draw("hist same")
hpX1gen.Draw("hist same")
hpX0rwt.Draw("hist same")
legpt.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".pTlepton.pdf"))
cnv.SaveAs(pdfname)


legeta = TLegend(0.15,0.40,0.50,0.92,"","brNDC")
legeta.SetFillStyle(4000); # will be transparent
legeta.SetFillColor(0)
legeta.SetTextFont(42)
legeta.SetBorderSize(0)
legeta.AddEntry(0, "MadGraph+Pythia8", "")
legeta.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}}#rightarrow#it{b}#it{W}^{+}(#it{u}#bar{#it{d}})#bar{#it{b}}#it{W}^{-}(#it{#mu}^{-}#bar{#it{#nu}}_{#mu})", "")
legeta.AddEntry(0, legcme+" ("+str(nkevents)+"k events)", "")
legeta.AddEntry(htSMgen,"SM (full range)","l")
legeta.AddEntry(htXXgen,"SM+#it{"+nameX+"} generated","l")
legeta.AddEntry(htXXrwt,"SM+#it{"+nameX+"} reweighted","l")
legeta.AddEntry(htX0gen,"#it{"+nameX+"} (only) generated","l")
legeta.AddEntry(htX0rwt,"#it{"+nameX+"} (only) reweighted","l")
legeta.AddEntry(htX1gen,"#it{H} (only) generated","l")
legeta.AddEntry(0, "sin(#beta-#alpha)=1", "")
legeta.AddEntry(0, "tan#beta="+str(tanb), "")
legeta.AddEntry(0, "#it{m}_{#it{X}}="+str(mX)+" GeV", "")

cnv = TCanvas("cnv","",600,445)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
htSMgen.SetMinimum(0)
htSMgen.SetMaximum(htSMgen.GetMaximum()*2)
htSMgen.Draw("hist")
htXXgen.Draw("hist same")
htXXrwt.Draw("hist same")
htX0gen.Draw("hist same")
htX1gen.Draw("hist same")
htX0rwt.Draw("hist same")
legeta.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".etalepton.pdf"))
cnv.SaveAs(pdfname)

### Ratios of angular distributions
cnv = TCanvas("cnv","",1200,445)
cnv.Divide(2,1)
cnv.Draw()
pad1 = cnv.cd(1)
pad1.SetGrid(1,1)
hcrXXgen.SetMinimum(0)
hcrXXgen.SetMaximum(2)
hcrXXgen.Draw()
pad1.RedrawAxis()
pad2 = cnv.cd(2)
pad2.SetGrid(1,1)
harXXgen.SetMinimum(0)
harXXgen.SetMaximum(2)
harXXgen.Draw()
pad2.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname.replace("(","").replace(")","").replace(".pdf",".anglesratio.pdf"))
cnv.SaveAs(pdfname+")")


