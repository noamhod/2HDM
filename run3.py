#!/usr/bin/python

import ROOT
from ROOT import *
import kinematics
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

def plot(h1, h2, h3, h4, tanb, wX, nkevents, cme, fname, ymin=-1, ymax=-1):
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
   tvp_ratio.SetBottomMargin(0.20);
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
   xLabelSize = h4.GetXaxis().GetLabelSize()*1.85;
   yLabelSize = h4.GetYaxis().GetLabelSize()*1.85;
   xTitleSize = h4.GetXaxis().GetTitleSize()*1.85;
   yTitleSize = h4.GetYaxis().GetTitleSize()*1.85;
   titleSize  = h4.GetTitleSize()           *1.85;
   h4.GetXaxis().SetLabelSize(xLabelSize);
   h4.GetYaxis().SetLabelSize(yLabelSize);
   h4.GetXaxis().SetTitleSize(xTitleSize);
   h4.GetYaxis().SetTitleSize(yTitleSize);
   h4.SetTitleSize(titleSize);
   h4.GetYaxis().SetTitleOffset(0.55);
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
   lineR = TLine(hr.GetXaxis().GetXmin(),1.,hr.GetXaxis().GetXmax(),1.);

   tvp_hists.SetBottomMargin(0);
   tvp_subtr.SetTopMargin(0);
   tvp_subtr.SetBottomMargin(0);
   tvp_ratio.SetTopMargin(0);
   tvp_ratio.SetBottomMargin(0.0);
   
   leg = TLegend(0.5,0.4,0.87,0.9,"","brNDC")
   leg.SetFillStyle(4000); # will be transparent
   leg.SetFillColor(0)
   leg.SetTextFont(42)
   leg.SetBorderSize(0)
   leg.AddEntry(0, "MadGraph+Pythia8", "")
   leg.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} ("+str(nkevents)+"k events)", "")
   leg.AddEntry(0, cme, "")
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
hmSMgen = TH1D("hmSMgen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXXrwt = TH1D("hmXXrwt", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXIrwt = TH1D("hmXIrwt", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXXgen = TH1D("hmXXgen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)
hmXIabs = TH1D("hmXIabs", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 10 GeV",70,300,1000)
hmXIabsStd = TH1D("hmXIabsStd", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events / 40 GeV",nbins,arrbins)

hmSMgen.Sumw2()
hmSMgen.SetLineWidth(2)
hmSMgen.SetLineColor(ROOT.kBlack)
hmSMgen.SetMarkerColor(ROOT.kBlack)
hmSMgen.SetMarkerStyle(24)

hmXXrwt.Sumw2()
hmXXrwt.SetLineWidth(2)
hmXXrwt.SetLineColor(ROOT.kRed)
hmXXrwt.SetMarkerColor(ROOT.kRed)
hmXXrwt.SetMarkerStyle(20)

hmXIrwt.Sumw2()
hmXIrwt.SetLineWidth(2)
hmXIrwt.SetLineColor(ROOT.kRed)
hmXIrwt.SetMarkerColor(ROOT.kRed)
hmXIrwt.SetMarkerStyle(20)

hmXXgen.Sumw2()
hmXXgen.SetLineWidth(2)
hmXXgen.SetLineColor(ROOT.kAzure+9)
hmXXgen.SetMarkerColor(ROOT.kAzure+9)
hmXXgen.SetMarkerStyle(24)

hmXIabs.Sumw2()
hmXIabs.SetLineWidth(2)
hmXIabs.SetLineColor(ROOT.kRed)

hmXIabsStd.Sumw2()
hmXIabsStd.SetLineWidth(2)
hmXIabsStd.SetLineColor(ROOT.kRed)




fSM = TFile("/Users/hod/MC/Pythia/pythia8215/examples/tops.SM_nobmass_8TeV.root","READ")
tSM = fSM.Get("SM_nobmass_8TeV")
# fSM = TFile("/Users/hod/MC/Pythia/pythia8215/examples/tops.SM_8TeV.root","READ")
# tSM = fSM.Get("SM_8TeV")
pSM = ROOT.vector(TLorentzVector)()
iSM = std.vector(int)()
tSM.SetBranchAddress("p4",pSM);
tSM.SetBranchAddress("id",iSM);

fXX = TFile("/Users/hod/MC/Pythia/pythia8215/examples/tops.SMIA_8TeV.root","READ")
tXX = fXX.Get("SMIA_8TeV")
pXX = ROOT.vector(TLorentzVector)()
iXX = std.vector(int)()
tXX.SetBranchAddress("p4",pXX);
tXX.SetBranchAddress("id",iXX);

nkevents = tXX.GetEntries()/1000
cme = "#sqrt{#it{s}} = 8 TeV" if("8TeV" in fSM.GetName()) else "#sqrt{#it{s}} = 13 TeV"

n=1
for event in tSM:
   if(n%10000==0): print "processed |SM|^2 and reweighting ", n
   g1=event.p4[0]
   g2=event.p4[1]
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()

   p = [[ g1.E(), g1.Px(), g1.Py(), g1.Pz() ],
        [ g2.E(), g2.Px(), g2.Py(), g2.Pz() ],
        [ t1.E(), t1.Px(), t1.Py(), t1.Pz() ],
        [ t2.E(), t2.Px(), t2.Py(), t2.Pz() ]]
   P=THDM.invert_momenta(p)

   # Q = 0.5*(math.sqrt(t1.M()*t1.M()+t1.Pt()*t1.Pt())+math.sqrt(t2.M()*t2.M()+t2.Pt()*t2.Pt()))
   alphaS = event.aS #THDM.model.AlphaS.alphasQ(Q)

   ## the ME^2 and the weight
   # me2SM = THDM.modules['matrix2SMpy'].get_me(P,alphaS,nhel)                   ### calculate the SM ME^2
   # me2XX = THDM.modules['matrix2'+nameX+str(index)+'py'].get_me(P,alphaS,nhel) ### calculate the X ME^2
   me2SM = THDM.modules['matrix2SMttxpy'].get_me(P,alphaS,nhel)                   ### calculate the SM ME^2
   me2XX = THDM.modules['matrix2'+nameX+str(index)+'ttxpy'].get_me(P,alphaS,nhel) ### calculate the X ME^2
   weightX = me2XX/me2SM                                                       ### calculate the weight

   hmSMgen.Fill(mtt)
   hmXXrwt.Fill(mtt,weightX)
   hmXIrwt.Fill(mtt,weightX-1)
   hmXIabs.Fill(mtt,weightX-1)
   hmXIabsStd.Fill(mtt,weightX-1)

   n+=1



n=1
for event in tXX:
   if(n%10000==0): print "processed |SM+X|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()
   hmXXgen.Fill(mtt)
   n+=1

hmSMgen = normalizeToBinWidth(hmSMgen,40)
hmXXgen = normalizeToBinWidth(hmXXgen,40)
hmXXrwt = normalizeToBinWidth(hmXXrwt,40)
hmXIrwt = normalizeToBinWidth(hmXIrwt,40)
# hmXIabs = normalizeToBinWidth(hmXIabs)
hmXIabsStd = normalizeToBinWidth(hmXIabsStd,40)


tanb = '%.2f' % THDM.parameters[index].get("tanb")
wXprcnt = '%.1f' % (THDM.parameters[index].get("w"+nameX)/mX*100.)
stanb = tanb.replace(".","")
pdfname = "2HDM.reweighting."+nameX+"."+str(mX)+"GeV.tanb"+stanb+".pdf"
plot(hmSMgen,hmXXrwt,hmXXgen,hmXIrwt,tanb,wXprcnt,nkevents,cme,pdfname+"(")



mb2fb = 1.e12
pb2fb = 1.e3
sigmaSM = 1.285e-07*mb2fb # 252.89*pb2fb
neventsSM = tXX.GetEntries()
lumiSM = neventsSM/sigmaSM
lumiData = 20.3 ## fb-1

legabs = TLegend(0.5,0.5,0.87,0.9,"","brNDC")
legabs.SetFillStyle(4000); # will be transparent
legabs.SetFillColor(0)
legabs.SetTextFont(42)
legabs.SetBorderSize(0)
legabs.AddEntry(0, "MadGraph+Pythia8", "")
legabs.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} ("+str(nkevents)+"k events)", "")
legabs.AddEntry(0, cme+", "+str(lumiData)+" fb^{-1}", "")
legabs.AddEntry(hmXIabs,"|SM+#it{"+nameX+"}|^{2}-|SM|^{2} reweighted","l")
legabs.AddEntry(0, "sin(#beta-#alpha)=1", "")
legabs.AddEntry(0, "tan#beta="+str(tanb), "")
legabs.AddEntry(0, "#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV", "")
legabs.AddEntry(0, "#Gamma_{#it{"+nameX+"}}/#it{m}_{#it{"+nameX+"}}="+str(wXprcnt)+" [%]", "")



hmXIabs.Scale(lumiData/lumiSM)
hmXIabs.SetMinimum(hmXIabs.GetMinimum()*1.1)
hmXIabs.SetMaximum(hmXIabs.GetMaximum()*2.3)
cnv = TCanvas("cnv","",600,500)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hmXIabs.Draw("hist")
legabs.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname)

hmXIabsStd.Scale(lumiData/lumiSM)
hmXIabsStd.SetMinimum(hmXIabsStd.GetMinimum()*1.1)
hmXIabsStd.SetMaximum(hmXIabsStd.GetMaximum()*2.3)
cnv = TCanvas("cnv","",600,500)
cnv.Draw()
cnv.cd()
cnv.SetGrid(1,1)
hmXIabsStd.Draw("hist")
legabs.Draw("same")
cnv.RedrawAxis()
cnv.Update()
cnv.SaveAs(pdfname+")")
