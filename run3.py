#!/usr/bin/python

import ROOT
from ROOT import std, gROOT, gStyle, gPad, TCanvas, TH1, TH1D, TH2D, TLegend, TLine, TFile, TTree, TLorentzVector, TMath, TVirtualPad, TEventList
import kinematics
import THDM
import sys
import os
from array import *
# import subprocess # just to call an arbitrary command e.g. 'ls'
# from contextlib import contextmanager
# from pipes import quote
# from pprint import pprint
# import imp

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

def plot(h1, h2, h3, tanb, wX, fname, ymin=-1, ymax=-1):
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
	
   hs = h2.Clone("subtr");
   hs.SetTitle(";"+sXtitle+";|SM+#it{"+nameX+"}|^{2}/|SM|^{2}-1 [%]");
   hs.Add(h1,-1.)
   hs.Divide(h1)
   hs.Scale(100)

   rmin=-11;
   rmax=+11;
   hs.SetMarkerStyle(20);
   hs.SetMarkerSize(0.8);
   hs.SetMarkerColor(ROOT.kBlack);
   hs.SetLineColor(ROOT.kBlack);
   hs.SetLineStyle(1);
   hs.SetLineWidth(2);
   xLabelSize = hs.GetXaxis().GetLabelSize()*1.85;
   yLabelSize = hs.GetYaxis().GetLabelSize()*1.85;
   xTitleSize = hs.GetXaxis().GetTitleSize()*1.85;
   yTitleSize = hs.GetYaxis().GetTitleSize()*1.85;
   titleSize  = hs.GetTitleSize()           *1.85;
   hs.GetXaxis().SetLabelSize(xLabelSize);
   hs.GetYaxis().SetLabelSize(yLabelSize);
   hs.GetXaxis().SetTitleSize(xTitleSize);
   hs.GetYaxis().SetTitleSize(yTitleSize);
   hs.SetTitleSize(titleSize);
   hs.GetYaxis().SetTitleOffset(0.55);
   hs.GetXaxis().SetTitleOffset(0.83);
   hs.SetMinimum(rmin);
   hs.SetMaximum(rmax);
   lineS = TLine(hs.GetXaxis().GetXmin(),0.,hs.GetXaxis().GetXmax(),0.);

   hr = th1d_tmp.Clone("ratio")
   hr.Sumw2()
   hr.SetTitle(";"+sXtitle+";|SM+#it{"+nameX+"}|^{2}_{rwt}/|SM+#it{"+nameX+"}|^{2}_{gen}")
   hr.Divide(h3)

   rmin=+0.80;
   rmax=+1.20;
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
   lineR = TLine(hr.GetXaxis().GetXmin(),1.,hr.GetXaxis().GetXmax(),1.);

   tvp_hists.SetBottomMargin(0);
   tvp_subtr.SetTopMargin(0);
   tvp_subtr.SetBottomMargin(0);
   tvp_ratio.SetTopMargin(0);
   tvp_ratio.SetBottomMargin(0.20);
   
   leg = TLegend(0.5,0.4,0.87,0.9,"","brNDC");
   leg.SetFillStyle(4000); # will be transparent
   leg.SetFillColor(0);
   leg.SetTextFont(42);
   leg.SetBorderSize(0);
   leg.AddEntry(0, "MadGraph+Pythia8", "");
   leg.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} (500k events)", "");
   leg.AddEntry(0, "Truth level for now", "");
   leg.AddEntry(h1,"|SM|^{2}","ple");
   leg.AddEntry(h2,"|SM+#it{"+nameX+"}|^{2} reweighted","ple");
   leg.AddEntry(h3,"|SM+#it{"+nameX+"}|^{2} generated","ple");
   leg.AddEntry(0, "#it{m}_{#it{"+nameX+"}}="+str(mX)+" GeV", "");
   leg.AddEntry(0, "#Gamma_{#it{"+nameX+"}}/#it{m}_{#it{"+nameX+"}}="+str(wX)+" [%]", "");
   leg.AddEntry(0, "sin(#beta-#alpha)=1", "");
   leg.AddEntry(0, "tan#beta="+str(tanb), "");
   
   tvp_hists.cd();
   # tvp_hists.SetLogy()
   if(ymin>-1): h3.SetMinimum(ymin);
   if(ymax>-1): h3.SetMaximum(ymax);
   h2.Draw();
   h1.Draw("same");
   h3.Draw("same");
   leg.Draw("same");
   tvp_hists.Update();
   tvp_hists.RedrawAxis();
   
   tvp_subtr.cd();
   tvp_subtr.SetGridy();
   hs.Draw("e1p");
   lineS.Draw("same");
   tvp_subtr.Update();
   tvp_subtr.RedrawAxis();

   tvp_ratio.cd();
   tvp_ratio.SetGridy();
   hr.Draw("e1p");
   lineR.Draw("same");
   tvp_ratio.Update();
   tvp_ratio.RedrawAxis();
   
   cnv.Update();
   cnv.SaveAs(fname);


# def getME2(libpath,g1,g2,t1,t2,printerr=False,printout=False):
#    sg1 = "-g1e "+str(g1.E())+" -g1px "+str(g1.Px())+" -g1py "+str(g1.Py())+" -g1pz "+str(g1.Pz())
#    sg2 = "-g2e "+str(g2.E())+" -g2px "+str(g2.Px())+" -g2py "+str(g2.Py())+" -g2pz "+str(g2.Pz())
#    st1 = "-t1e "+str(t1.E())+" -t1px "+str(t1.Px())+" -t1py "+str(t1.Py())+" -t1pz "+str(t1.Pz())
#    st2 = "-t2e "+str(t2.E())+" -t2px "+str(t2.Px())+" -t2py "+str(t2.Py())+" -t2pz "+str(t2.Pz())
#    momenta = sg1+" "+sg2+" "+st1+" "+st2
#    args    = "-path "+libpath+" "+momenta
#    me2 = -1
#    proc = subprocess.Popen("python matrix.py "+args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    out, err = proc.communicate()
#    if(printerr): print "err=",err
#    if(printout): print "out=",out
#    me2 = float(out) 
#    return me2

##################
### Begin code ###
##################

setStyle()


### make the libraries
# type=2
# sba=1
# mX=500
# nameX="H"
# cuts = "m"+nameX+"=="+str(mX)+" && sba==1 && TMath::ATan(tanb)>0. && TMath::ATan(tanb)<TMath::Pi()/2. && TMath::Abs(cba)<=1. && type=="+str(type)+" && (status&3)==0"
# mgpath = "/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/"
# alphaS = 0.13
# nhel   = 0 # means sum over all helicity
# libmatrix = "matrix/"+nameX+"/"+str(mX)+"/"
# THDM.setParameters(nameX,mX,cuts,type,sba)

type   = THDM.model.type
sba    = THDM.model.sba
mX     = THDM.model.mX
nameX  = THDM.model.nameX
cuts   = THDM.model.cuts
mgpath = THDM.model.mgpath
alphaS = THDM.model.alphaS
nhel   = THDM.model.nhel
libmatrix = "matrix/"+nameX+"/"+str(mX)+"/"
THDM.setParameters(nameX,mX,cuts,type,sba)
###############
# [20] tanb=0.660000 sba=1.000000 cba=0.000000 wA=52.510675 wH=29.593605 YMT=261.363636 YMB=3.102000 YMC=2.151515 YMM=0.069736 YMTAU=1.172820
index=20 ######
THDM.setModules(os.getcwd(),libmatrix,nameX,len(THDM.parameters),"SM")
THDM.setModules(os.getcwd(),libmatrix,nameX,len(THDM.parameters),"X",index)
###############


### mass histos
hmSM     = TH1D("hmSM",     ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSMXrwt = TH1D("hmSMXrwt", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSMXgen = TH1D("hmSMXgen", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)

hmSM.Sumw2()
hmSM.SetLineWidth(2)
hmSM.SetLineColor(ROOT.kBlack)
hmSM.SetMarkerColor(ROOT.kBlack)
hmSM.SetMarkerStyle(24)

hmSMXrwt.Sumw2()
hmSMXrwt.SetLineWidth(2)
hmSMXrwt.SetLineColor(ROOT.kRed)
hmSMXrwt.SetMarkerColor(ROOT.kRed)
hmSMXrwt.SetMarkerStyle(20)

hmSMXgen.Sumw2()
hmSMXgen.SetLineWidth(2)
hmSMXgen.SetLineColor(ROOT.kAzure+9)
hmSMXgen.SetMarkerColor(ROOT.kAzure+9)
hmSMXgen.SetMarkerStyle(24)


fSM = TFile("tops.SM.500k.root","READ")
tSM = fSM.Get("SM")
pSM = ROOT.vector(TLorentzVector)()
iSM = std.vector(int)()
tSM.SetBranchAddress("p4",pSM);
tSM.SetBranchAddress("id",iSM);

fXX = TFile("tops.SMIH_valid.500k.root","READ")
# fXX = TFile("tops.SMIH_valid.50k.root","READ")
tXX = fXX.Get("SMIH_valid")
pXX = ROOT.vector(TLorentzVector)()
iXX = std.vector(int)()
tXX.SetBranchAddress("p4",pXX);
tXX.SetBranchAddress("id",iXX);


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
   ## the ME^2 and the weight
   me2SM = THDM.modules['matrix2SMpy'].get_me(P,alphaS,nhel)                   ### calculate the SM ME^2
   me2XX = THDM.modules['matrix2'+nameX+str(index)+'py'].get_me(P,alphaS,nhel) ### calculate the X ME^2
   weightX = me2XX/me2SM                                                       ### calculate the weight

   hmSM.Fill(mtt)
   hmSMXrwt.Fill(mtt,weightX)

   n+=1



n=1
for event in tXX:
   if(n%10000==0): print "processed |SM+X|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()
   hmSMXgen.Fill(mtt)
   n+=1


tanb = '%.2f' % THDM.parameters[index].get("tanb")
wXprcnt = '%.1f' % (THDM.parameters[index].get("w"+nameX)/mX*100.)
stanb = tanb.replace(".","")
plot(hmSM,hmSMXrwt,hmSMXgen,tanb,wXprcnt,"2HDM.reweighting."+nameX+"."+str(mX)+"GeV.tanb"+stanb+".pdf")
