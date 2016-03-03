#!/usr/bin/python

import ROOT
from ROOT import std, gROOT, gStyle, gPad, TCanvas, TH1, TH2, TH1D, TH2D, TLegend, TLine, TFile, TTree, TLorentzVector, TMath, TVirtualPad, TEventList, TPaveText, TGraph, TVectorD
from array import *
from sortedcontainers import SortedDict
import numpy as NUMPY
import sys
import os
import THDM


def setStyle():
   gROOT.Reset()
   icol=0; # WHITE
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
   gStyle.SetPaintTextFormat("4.1f")


def makeTemplate(sinba,nameX,mX,h1dict,hname,htitle,ybins): #,nybins,ymin,ymax):
   h2 = TH2D(hname,htitle+";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];tan#beta;Events",25,350,850,len(ybins)-1,ybins) #,nybins,ymin,ymax)
   for stanb, h1 in h1dict.iteritems():
      tanb = float(stanb)
      ################################################
      doSkip = skipPoint(stanb,sinba,nameX,mX)
      if(doSkip): continue
      ################################################
      y = h2.GetYaxis().FindBin(tanb)
      for x in range(1,h1.GetNbinsX()+1):
         val = h1.GetBinContent(x)
         # if(val<0): print "mtt="+str(h1.GetBinCenter(x))+", tan(beta)="+stanb+" --> A+I="+str(val)
         h2.SetBinContent(x,y,val)
   return h2


def plotWidth(dwdictX,fname,nameX,mX,cuts):
   sorted_dwdictX = SortedDict(dwdictX)
   n = len(sorted_dwdictX)-1
   x = array('d',sorted_dwdictX.keys())
   y = array('d',sorted_dwdictX.values())
   gwX = TGraph(n,x,y)
   gwX.SetName("gwX")
   gwX.SetTitle("")
   gwX.GetXaxis().SetTitle("tan#beta")
   gwX.GetYaxis().SetTitle("#Gamma_{#it{"+nameX+"}}/#it{m}_{#it{"+nameX+"}} [%]")
   gwX.SetLineColor(ROOT.kBlack)
   gwX.SetMarkerColor(ROOT.kBlack)
   gwX.SetMarkerStyle(20)
   gwX.SetMarkerSize(0.5)

   ptxt = TPaveText(0.62,0.70,0.87,0.87,"NDC")
   ptxt.SetFillStyle(4000) #will be transparent
   ptxt.SetFillColor(0)
   ptxt.SetTextFont(42)
   ptxt.SetBorderSize(0)
   ptxt.AddText("sin(#beta-#alpha)=1")
   ptxt.AddText("#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV")

   c = TCanvas("c","c",600,600)
   c.cd()
   c.SetLogx()
   c.SetLogy()
   c.SetGridx()
   c.SetGridy()
   c.SetTicks(1,1)
   c.Draw()
   # gwX.Draw("p")
   gwX.Draw()
   ptxt.Draw("same")
   c.Modified()
   c.Update()
   c.SaveAs(fname)




def plotYMX(dydictX,fname,nameX,mX,cuts):
   sorted_dydictX = SortedDict(dydictX)
   n = len(sorted_dydictX)-1
   x = array('d',sorted_dydictX.keys())
   y = array('d',sorted_dydictX.values("ymt")) ## ??

   ymt = TGraph(n,x,y)
   ymt.SetName("ymt")
   ymt.GetXaxis().SetTitle("tan#beta")
   ymt.GetYaxis().SetTitle("#it{ymt-"+nameX+"}")
   ymt.SetLineColor(ROOT.kBlack)
   ymt.SetMarkerColor(ROOT.kBlack)
   ymt.SetMarkerStyle(20)
   ymt.SetMarkerSize(0.8)

   ptxt = TPaveText(0.62,0.70,0.87,0.87,"NDC")
   ptxt.SetFillStyle(4000) #will be transparent
   ptxt.SetFillColor(0)
   ptxt.SetTextFont(42)
   ptxt.SetBorderSize(0)
   ptxt.AddText("sin(#beta-#alpha)=1")
   ptxt.AddText("#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV")

   c = TCanvas("c","c",600,600)
   c.cd()
   c.SetLogx()
   c.SetLogy()
   c.SetGridx()
   c.SetGridy()
   c.SetTicks(1,1)
   c.Draw()
   ymt.Draw()
   ptxt.Draw("same")
   c.Modified()
   c.Update()
   c.SaveAs(fname)


def plot1d(hSM,hXX,hIX,wX,fname,stanb,nameX,mX,ymin=-1,ymax=-1):
   cnv = TCanvas("cnv","",600,600);
   cnv.Divide(1,2);
   tvp_hists = cnv.cd(1);
   tvp_ratio = cnv.cd(2);
   cnv.Draw();
   tvp_hists.SetPad(0.00, 0.35, 1.00, 1.00);
   tvp_ratio.SetPad(0.00, 0.02, 1.00, 0.35);
   tvp_hists.SetBottomMargin(0.012);
   tvp_ratio.SetBottomMargin(0.20);
   tvp_ratio.SetTopMargin(0.012);
   tvp_hists.SetTicks(1,1);
   tvp_ratio.SetTicks(1,1);	

   sXtitle = hXX.GetXaxis().GetTitle();
   cloneName_n = hXX.GetName();
   cloneName_d = hSM.GetName();
   th1n_tmp = hXX.Clone(cloneName_n+"_th1n_tmp");
   th1d_tmp = hSM.Clone(cloneName_d+"_th1d_tmp");
   th1n_tmp.SetBinErrorOption(TH1.kPoisson);
   th1d_tmp.SetBinErrorOption(TH1.kPoisson);

   hr = hIX.Clone("ratio")
   hr.Sumw2()
   hr.SetTitle(";"+sXtitle+";|SM+#it{"+nameX+"}|^{2}/|SM|^{2}-1 [%]")
   hr.Divide(hSM)
   hr.Scale(100.)

   rmin=-11;
   rmax=+11;
   hr.SetMarkerStyle(20);
   hr.SetMarkerSize(0.8);
   hr.SetMarkerColor(ROOT.kBlack);
   hr.SetLineColor(ROOT.kBlack);
   hr.SetLineStyle(1);
   hr.SetLineWidth(2);
   xLabelSize = hr.GetXaxis().GetLabelSize()*1.85;
   yLabelSize = hr.GetYaxis().GetLabelSize()*1.85;
   xTitleSize = hr.GetXaxis().GetTitleSize()*1.85;
   yTitleSize = hr.GetYaxis().GetTitleSize()*1.85;
   titleSize  = hr.GetTitleSize()           *1.85;
   hr.GetXaxis().SetLabelSize(xLabelSize);
   hr.GetYaxis().SetLabelSize(yLabelSize);
   hr.GetXaxis().SetTitleSize(xTitleSize);
   hr.GetYaxis().SetTitleSize(yTitleSize);
   hr.SetTitleSize(titleSize);
   hr.GetYaxis().SetTitleOffset(0.55);
   hr.GetXaxis().SetTitleOffset(0.83);
   hr.SetMinimum(rmin);
   hr.SetMaximum(rmax);
   # lineR = TLine(hr.GetXaxis().GetXmin(),1.,hr.GetXaxis().GetXmax(),1.);

   tvp_hists.SetBottomMargin(0);
   tvp_ratio.SetTopMargin(0);
   tvp_ratio.SetBottomMargin(0.20);

   leg = TLegend(0.5,0.50,0.87,0.9,"","brNDC");
   leg.SetFillStyle(4000); # will be transparent
   leg.SetFillColor(0);
   leg.SetTextFont(42);
   leg.SetBorderSize(0);
   leg.AddEntry(0, "MadGraph+Pythia8", "");
   leg.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} (50k events)", "");
   leg.AddEntry(hSM,"|SM+#it{"+nameX+"}|^{2} reweighted","ple");
   leg.AddEntry(hXX,"|SM|^{2}","ple");
   leg.AddEntry(0, "tan#beta="+stanb, "");
   leg.AddEntry(0, "sin(#beta-#alpha)=1", "");
   leg.AddEntry(0, "#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV", "");
   leg.AddEntry(0, "#Gamma_{#it{"+nameX+"}}="+'%.4f' % wX+"%", "");

   tvp_hists.cd();
   if(ymin>-1): hSM.SetMinimum(ymin);
   if(ymax>-1): hSM.SetMaximum(ymax);
   hXX.Draw();
   hSM.Draw("same");
   leg.Draw("same");
   tvp_hists.Update();
   tvp_hists.RedrawAxis();

   tvp_ratio.cd();
   tvp_ratio.SetGridy();
   hr.Draw("e1p");
   # lineR.Draw("same");
   tvp_ratio.Update();
   tvp_ratio.RedrawAxis();

   cnv.Update();
   cnv.SaveAs(fname)



# def plot2d(h,nsmooth,fname,drawopt,nameX,mX,hR=0,minz=+1.e11,maxz=-1.e11,phi=-110,theta=30,dologx=False,dology=True,dologz=False):
def plot2d(h,nsmooth,fname,drawopt,nameX,mX,hR=0,minz=+1.e11,maxz=-1.e11,phi=-50,theta=25,dologx=False,dology=True,dologz=False):
   for i in range(0,nsmooth):
      h.Smooth()
   c = TCanvas("c"+h.GetName(),"c"+h.GetName(),1000,600)
   c.cd()
   c.Draw()
   c.SetTopMargin(0.1)

   if(dologx): c.SetLogx()
   if(dology): c.SetLogy()
   if(dologz): c.SetLogz()

   if(minz<+1.e10): h.SetMinimum(minz)
   if(maxz>-1.e10): h.SetMaximum(maxz)

   if(hR):
      h.Divide(hR)
      h.Scale(100)
      h.GetZaxis().SetTitle("(SM+#it{"+nameX+"})/SM-1 [%]")
      h.GetZaxis().SetTitleOffset(-0.1)

   # h.GetYaxis().SetLimits(0.,10.)
   # ymin = 0.
   # ymax = 10.
   # nPrimitiveBins = 10
   # nBins = nPrimitiveBins*int(ymax-ymin)
   # h.GetYaxis().Set(nBins,ymin,ymax)

   # h.GetYaxis().SetNdivisions((Int_t)(10.*16),kFALSE)
   # h.RebinY(5)

   h.SetTitleOffset(0.5)
   if(("SURF" in drawopt) or ("LEGO" in drawopt)):
      c.SetPhi(phi)
      c.SetTheta(theta)
      h.GetZaxis().SetTitleOffset(1.10)
      h.GetXaxis().SetTitleOffset(1.60)
      h.GetYaxis().SetTitleOffset(1.30)
   h.Draw(drawopt);
   c.Modified()
   c.Update()
   c.SaveAs(fname)
   c.SaveAs(fname.replace(".pdf",".root"))


def skipPoint(stanb,sinba,nameX,mX):
   if(nameX=="H" and mX==500 and sinba==1 and stanb=="0.22"): return True
   if(nameX=="H" and mX==500 and sinba==1 and stanb=="0.54"): return True
   if(nameX=="H" and mX==500 and sinba==1 and stanb=="0.74"): return True ## not sure if this is needed
   return False


##################
### Begin code ###
##################



#### style
setStyle()


### make the libraries
# type=2
# sba=1
# mX=500
# nameX="H"
# cuts = "m"+nameX+"=="+str(mX)+" && sba==1 && TMath::ATan(tanb)>0. && TMath::ATan(tanb)<TMath::Pi()/2. && TMath::Abs(cba)<=1. && type=="+str(type)+" && (status&3)==0"
type   = THDM.model.type
sba    = THDM.model.sba
mX     = THDM.model.mX
nameX  = THDM.model.nameX
cuts   = THDM.model.cuts
THDM.setParameters(nameX,mX,cuts,type,sba)


### get the tree
# fSM = TFile("2HDM."+nameX+"."+str(mX)+"GeV.tree.50k.root","READ")
fSM = TFile("2HDM."+nameX+"."+str(mX)+"GeV.tree.500k.root","READ")
tSM = fSM.Get("SM")
pSM = ROOT.vector(TLorentzVector)()
iSM = std.vector(int)()
wgt = std.vector(float)()
tnb = std.vector(float)()
wdtA = std.vector(float)()
wdtH = std.vector(float)()
ymt   = std.vector(float)()
ymb   = std.vector(float)()
ymc   = std.vector(float)()
ymtau = std.vector(float)()
ymmu  = std.vector(float)()
tSM.SetBranchAddress("p4",pSM);
tSM.SetBranchAddress("id",iSM);
tSM.SetBranchAddress("wgt",wgt);
tSM.SetBranchAddress("tnb",tnb);
tSM.SetBranchAddress("wdtA",wdtA);
tSM.SetBranchAddress("wdtH",wdtH);
tSM.SetBranchAddress("ymt",ymt);
tSM.SetBranchAddress("ymb",ymb);
tSM.SetBranchAddress("ymc",ymc);
tSM.SetBranchAddress("ymtau",ymtau);
tSM.SetBranchAddress("ymmu",ymmu);


### build the dictionaries
h1dictXX = {}
h1dictSM = {}
h1dictIX = {}
swdictX = {}
dwdictX = {}
sydictX = {}
dydictX = {}
nparameters = len(THDM.parameters)
for i in range(0,nparameters):
   tanb = '%.2f' % THDM.parameters[i].get("tanb")
   wX = THDM.parameters[i].get("w"+nameX)/mX*100

   h1XX = TH1D("h1XX."+tanb,";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
   h1SM = TH1D("h1SM."+tanb,";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
   h1IX = TH1D("h1IX."+tanb,";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
 
   h1SM.Sumw2()
   h1SM.SetLineColor(ROOT.kBlack)
   h1SM.SetMarkerColor(ROOT.kBlack)
   h1SM.SetMarkerStyle(24)

   h1XX.Sumw2()
   h1XX.SetLineColor(ROOT.kRed)
   h1XX.SetMarkerColor(ROOT.kRed)
   h1XX.SetMarkerStyle(20)

   h1IX.Sumw2()
   h1IX.SetLineColor(ROOT.kBlack)
   h1IX.SetMarkerColor(ROOT.kBlack)
   h1IX.SetMarkerStyle(20)

   ### Fill the dictionaries
   h1dictXX.update({tanb:h1XX})
   h1dictSM.update({tanb:h1SM})
   h1dictIX.update({tanb:h1IX})
   swdictX.update({tanb:wX})
   dwdictX.update({float(tanb):wX})

   # ymx = {"ymt":THDM.parameters[i].get("ymt"), "ymb":THDM.parameters[i].get("ymb"), "ymc":THDM.parameters[i].get("ymc"), "ymtau":THDM.parameters[i].get("ymtau"), "ymmu":THDM.parameters[i].get("ymmu")}
   # sydictX.update({tanb:ymx})
   # dydictX.update({float(tanb):ymx})


n=1
for event in tSM:
   if(n%5000==0): print "processed |SM|^2 and reweighting ", n
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()
   for j in range(0,wgt.size()):
      weight = wgt[j]
      tanbet = '%.2f' % tnb[j]
      h1dictXX[tanbet].Fill(mtt,weight)
      h1dictSM[tanbet].Fill(mtt,1)
      h1dictIX[tanbet].Fill(mtt,weight-1.)
   n+=1


### the binning for the tan(beta) axis !!
tanblist = []
for i in range(0,len(THDM.parameters)):
   ################################################
   stanb = '%.2f' % THDM.parameters[i].get("tanb")
   doSkip = skipPoint(stanb,sba,nameX,mX)
   if(doSkip): continue
   ################################################
   tanblist.append(THDM.parameters[i].get("tanb"))
tanblist_sorted = sorted(tanblist)
print "tan(beta) points:"
for tb in tanblist_sorted:
   sys.stdout.write('%.2f' % tb)
   sys.stdout.write(", ")
sys.stdout.write("\n\n")
tanbbins = []
ntb = len(tanblist_sorted)
tanbbins.append(tanblist_sorted[0] - (tanblist_sorted[1]-tanblist_sorted[0])/2.)
for i in range(0,ntb):
   point = -1
   if(i<(ntb-1)): point = tanblist_sorted[i]+(tanblist_sorted[i+1]-tanblist_sorted[i])/2.
   else:          point = tanblist_sorted[i]+(tanblist_sorted[i]-tanblist_sorted[i-1])/2.
   tanbbins.append(point)
print "tan(beta) bins:"
for tb in tanbbins:
   sys.stdout.write('%.2f' % tb)
   sys.stdout.write(", ")
sys.stdout.write("\n\n")
n = len(tanbbins)
bins = array('d',tanbbins)



# tanbbins = 25
# mintanb  = 0
# maxtanb  = 8
hXX = makeTemplate(sba,nameX,mX,h1dictXX,"hXX","SM+#it{"+nameX+"} #it{gg}#rightarrow#it{t}#bar{#it{t}} at #it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV",bins)#tanbbins,mintanb,maxtanb)
hSM = makeTemplate(sba,nameX,mX,h1dictSM,"hSM","SM only #it{gg}#rightarrow#it{t}#bar{#it{t}}",bins)#tanbbins,mintanb,maxtanb)
hIX = makeTemplate(sba,nameX,mX,h1dictIX,"hIX",'"I"+#it{'+nameX+"} #it{gg}#rightarrow#it{t}#bar{#it{t}} at #it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV",bins)#tanbbins,mintanb,maxtanb)

fname = "2HDM.reweighting.2d.histos."+nameX+"."+str(mX)+".GeV.pdf"
nsmooth = 3
c = TCanvas("c","c",600,600)
c.SaveAs(fname+"(")
plot2d(hSM,nsmooth,fname,"SURF3",nameX,mX)
plot2d(hXX,nsmooth,fname,"SURF3",nameX,mX)
plot2d(hIX,nsmooth,fname,"SURF3",nameX,mX,hSM)
plotWidth(dwdictX,fname,nameX,mX,cuts)
# plotYMX(dydictX,fname,nameX,mX,cuts)

for tanb in tanblist_sorted:
   stanb = '%.2f' % tanb
# for stanb, h1 in h1dictXX.iteritems():
   plot1d(h1dictSM[stanb],h1dictXX[stanb],h1dictIX[stanb],swdictX[stanb],fname,stanb,nameX,mX)
c = TCanvas("c","c",600,600)
c.SaveAs(fname+")")