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

hmX0rwt.Sumw2()
hmX0rwt.SetLineWidth(2)
hmX0rwt.SetLineColor(ROOT.kRed)
hmX0rwt.SetMarkerColor(ROOT.kRed)
hmX0rwt.SetMarkerStyle(20)

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

hmX0gen.Sumw2()
hmX0gen.SetLineWidth(2)
hmX0gen.SetLineColor(ROOT.kAzure+9)
hmX0gen.SetMarkerColor(ROOT.kAzure+9)
hmX0gen.SetMarkerStyle(24)

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
lumiData = 20.3 # fb-1
sqrts = 8 #or 13 TeV
docuts = False
#############################

fpath = "/Users/hod/MC/Pythia/pythia8215/examples/"
fSM = TFile(fpath+"tops.SM_nobmass_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.SM_nobmass.root","READ")
tSM = fSM.Get("SM_nobmass_8TeV") if(sqrts==8) else fSM.Get("SM_nobmass")
# pSM = ROOT.vector(TLorentzVector)()
# iSM = std.vector(int)()
# tSM.SetBranchAddress("p4",pSM);
# tSM.SetBranchAddress("id",iSM);

fXX = TFile(fpath+"tops.SMIA_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.SMIA.root","READ")
tXX = fXX.Get("SMIA_8TeV") if(sqrts==8) else fXX.Get("SMIA")
# pXX = ROOT.vector(TLorentzVector)()
# iXX = std.vector(int)()
# tXX.SetBranchAddress("p4",pXX);
# tXX.SetBranchAddress("id",iXX);

fX0 = TFile(fpath+"tops.A_8TeV.root","READ") if(sqrts==8) else TFile(fpath+"tops.A.root","READ")
tX0 = fX0.Get("A_8TeV") if(sqrts==8) else fXX.Get("A")
# pX0 = ROOT.vector(TLorentzVector)()
# iX0 = std.vector(int)()
# tX0.SetBranchAddress("p4",pX0);
# tX0.SetBranchAddress("id",iX0);

nkevents = tXX.GetEntries()/1000
legcme = "#sqrt{#it{s}} = 8 TeV" if("8TeV" in fSM.GetName()) else "#sqrt{#it{s}} = 13 TeV"
scme = "8TeV" if("8TeV" in fSM.GetName()) else "13TeV"
leglumi = "20.3 fb^{-1}"


n=1
sumofweights_X0 = 0
sumofweights_SM = 0
for event in tSM:
   if(n%10000==0): print "processed |SM|^2 and reweighting ", n
   g1=event.p4[0]
   g2=event.p4[1]
   t1=event.p4[2]
   t2=event.p4[3]
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

   hmSMgen.Fill(mtt)
   hmXXrwt.Fill(mtt,weightXX)
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
   t1=event.p4[2]
   t2=event.p4[3]
   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue
   mtt = (t1+t2).M()
   hmXXgen.Fill(mtt)
   n+=1

n=1
for event in tX0:
   if(n%10000==0): print "processed |X|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   if(docuts and (t1.Pt()<60 or t2.Pt()<60)):               continue
   if(docuts and (abs(t1.Eta())>2.5 or abs(t2.Eta())>2.5)): continue
   mtt = (t1+t2).M()
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
cnv.SaveAs(pdfname+")")
