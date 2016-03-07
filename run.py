#!/usr/bin/python

import ROOT
from ROOT import std, gROOT, gStyle, gPad, TCanvas, TH1, TH1D, TH2D, TLegend, TLine, TFile, TTree, TLorentzVector, TMath, TVirtualPad, TEventList
import matrix2SMpy
import matrix2SMIHpy
import matrix2Hpy
import matrix2Apy
import kinematics
import THDM


def invert_momenta(p):
   """ fortran/C-python do not order table in the same order"""
   new_p = []
   for i in range(len(p[0])):  new_p.append([0]*len(p))
   for i, onep in enumerate(p):
      for j, x in enumerate(onep):
         new_p[j][i] = x
   return new_p


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
   gStyle.SetOptTitle(0);
   gStyle.SetOptStat(0);
   gStyle.SetOptFit(0);
   gStyle.SetPadTickX(1);
   gStyle.SetPadTickY(1);

def plot(h1, h2, h3, fname, ymin=-1, ymax=-1):
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
	
   hs = th1n_tmp.Clone("subtr");
   hs.SetTitle(";"+sXtitle+";|SM+#it{A}|^{2}-|SM|^{2}");
   hs.Reset();
   for b in range (1, hs.GetNbinsX()+1):
      v1 = th1n_tmp.GetBinContent(b);
      v2 = th1d_tmp.GetBinContent(b);
      d1 = th1n_tmp.GetBinError(b);
      d2 = th1d_tmp.GetBinError(b);
      r  = v2-v1;
      dr = TMath.Sqrt(d1*d1 + d2*d2);
      hs.SetBinContent(b,r);
      hs.SetBinError(b,dr);

   rmin=-2350;
   rmax=+2350;
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
   hr.SetTitle(";"+sXtitle+";|SM+#it{A}|^{2}_{rwt}/|SM+#it{A}|^{2}_{gen}")
   hr.Divide(h3)

   # hr.Reset()
   # for b in range (1, hr.GetNbinsX()+1):
   #    v1 = th1n_tmp.GetBinContent(b);
   #    v2 = th1d_tmp.GetBinContent(b);
   #    d1 = th1n_tmp.GetBinError(b);
   #    d2 = th1d_tmp.GetBinError(b);
   #    if(v2==0. or v1==0.): continue
   #    r  = v2/v1;
   #    dr = (v2/v1)*TMath.Sqrt((d1/v1)*(d1/v1) + (d2/v2)*(d2/v2));
   #    hr.SetBinContent(b,r);
   #    hr.SetBinError(b,dr);

   rmin=+0.8;
   rmax=+1.2;
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
   
   leg = TLegend(0.5,0.55,0.87,0.9,"","brNDC");
   leg.SetFillStyle(4000); # will be transparent
   leg.SetFillColor(0);
   leg.SetTextFont(42);
   leg.SetBorderSize(0);
   leg.AddEntry(0, "MadGraph+Pythia8", "");
   leg.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} (500k events)", "");
   leg.AddEntry(0, "Truth level for now", "");
   leg.AddEntry(h1,"|SM|^{2}","ple");
   leg.AddEntry(h2,"|SM+#it{A}|^{2} reweighted","ple");
   leg.AddEntry(h3,"|SM+#it{A}|^{2} generated","ple");
   leg.AddEntry(0, "#it{m}_{#it{A}} = 500 GeV", "");
   
   tvp_hists.cd();
   if(ymin>-1): h3.SetMinimum(ymin);
   if(ymax>-1): h3.SetMaximum(ymax);
   h3.Draw();
   h2.Draw("same");
   h1.Draw("same");
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


def plotnorm(h1, h2, h3, fname, particle, swap=False):
   leg = TLegend(0.15,0.65,0.50,0.92,"","brNDC");
   leg.SetFillStyle(4000); # will be transparent
   leg.SetFillColor(0);
   leg.SetTextFont(42);
   leg.SetBorderSize(0);
   leg.AddEntry(0, "MadGraph+Pythia8", "");
   leg.AddEntry(0, "#it{gg}#rightarrow#it{t}#bar{#it{t}} (500k events)", "");
   leg.AddEntry(0, "Truth level for now", "");
   leg.AddEntry(h1,"|SM|^{2}","ple");
   leg.AddEntry(h2,"|#it{"+particle+"}|^{2} reweighted","ple");
   leg.AddEntry(h3,"|#it{"+particle+"}|^{2} generated","ple");
   leg.AddEntry(0, "#it{m}_{#it{"+particle+"}} = 500 GeV", "");
   # leg.AddEntry(0, "in |#it{m}_{#it{t#bar{t}}}-#it{m}_{#it{"+particle+"}}|<100 GeV", "");

   cnv = TCanvas("cnv","",600,600);
   cnv.Draw();
   cnv.SetTicks(1,1);
   h1.GetYaxis().SetTitle("Norm.")
   h2.GetYaxis().SetTitle("Norm.")
   h3.GetYaxis().SetTitle("Norm.")
   if(swap):
      h2.DrawNormalized()
      h1.DrawNormalized("same")
   else:
      h1.DrawNormalized()
      h2.DrawNormalized("same")
   h3.DrawNormalized("same")
   leg.Draw("same")

   cnv.Update()
   cnv.RedrawAxis()
   cnv.SaveAs(fname)



##################
### Begin code ###
##################


THDM.testTHDM()
quit()



setStyle()

matrix2SMpy.initialise('param_card_SM.dat')
matrix2SMIHpy.initialise('param_card_SMIH.dat')
matrix2Hpy.initialise('param_card_H.dat')
matrix2Apy.initialise('param_card_A.dat')
alphaS = 0.13
nhel = 0 # means sum over all helicity
#### dummy test
# p = [[   0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03],
#     [   0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03],
#     [   0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03],
#     [   0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03]]
# P=invert_momenta(p)
# me2SM = matrix2SMpy.get_me(P, alphaS, nhel)
# me2SMIH = matrix2SMIHpy.get_me(P, alphaS, nhel)
# me2H = matrix2Hpy.get_me(P, alphaS, nhel)
# me2A = matrix2Apy.get_me(P, alphaS, nhel)
# print me2SM
# print me2SMIH
# print me2H
# print me2A





### mass histos
hmSM    = TH1D("hmSM",    ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSMIH  = TH1D("hmSMIH",  ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSMIH0 = TH1D("hmSMIH0", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSMx   = TH1D("hmSMx",   ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",100,400,600)
hmH     = TH1D("hmH",     ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",100,400,600)
hmA     = TH1D("hmA",     ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",100,400,600)
hmH0    = TH1D("hmH0",    ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",100,400,600)
hmA0    = TH1D("hmA0",    ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",100,400,600)

hmSM.Sumw2()
hmSM.SetLineWidth(2)
hmSM.SetLineColor(ROOT.kBlack)
hmSM.SetMarkerColor(ROOT.kBlack)
hmSM.SetMarkerStyle(24)

hmSMx.Sumw2()
hmSMx.SetLineWidth(2)
hmSMx.SetLineColor(ROOT.kBlack)
hmSMx.SetMarkerColor(ROOT.kBlack)
hmSMx.SetMarkerStyle(24)

hmSMIH.Sumw2()
hmSMIH.SetLineWidth(2)
hmSMIH.SetLineColor(ROOT.kRed)
hmSMIH.SetMarkerColor(ROOT.kRed)
hmSMIH.SetMarkerStyle(20)

hmSMIH0.Sumw2()
hmSMIH0.SetLineWidth(2)
hmSMIH0.SetLineColor(ROOT.kAzure+9)
hmSMIH0.SetMarkerColor(ROOT.kAzure+9)
hmSMIH0.SetMarkerStyle(24)

hmH.Sumw2()
hmH.SetLineWidth(2)
hmH.SetLineColor(ROOT.kViolet+1)
hmH.SetMarkerColor(ROOT.kViolet+1)
hmH.SetMarkerStyle(22)

hmA.Sumw2()
hmA.SetLineWidth(2)
hmA.SetLineColor(ROOT.kViolet+1)
hmA.SetMarkerColor(ROOT.kViolet+1)
hmA.SetMarkerStyle(23)

hmA0.Sumw2()
hmA0.SetLineWidth(2)
hmA0.SetLineColor(ROOT.kGreen+1)
hmA0.SetMarkerColor(ROOT.kGreen+1)
hmA0.SetMarkerStyle(32)

hmH0.Sumw2()
hmH0.SetLineWidth(2)
hmH0.SetLineColor(ROOT.kGreen+1)
hmH0.SetMarkerColor(ROOT.kGreen+1)
hmH0.SetMarkerStyle(26)




### cos(theta*) Helicity histos
hctbSM    = TH1D("hctbSM",    ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)
hctbSMIH  = TH1D("hctbSMIH",  ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)
hctbSMIH0 = TH1D("hctbSMIH0", ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)
hctbH     = TH1D("hctbH",     ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)
hctbA     = TH1D("hctbA",     ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)
hctbA0    = TH1D("hctbA0",    ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)
hctbH0    = TH1D("hctbH0",    ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-1.,+1.)

hctbSM.Sumw2()
hctbSM.SetLineWidth(2)
hctbSM.SetLineColor(ROOT.kBlack)
hctbSM.SetMarkerColor(ROOT.kBlack)
hctbSM.SetMarkerStyle(24)

hctbSMIH.Sumw2()
hctbSMIH.SetLineWidth(2)
hctbSMIH.SetLineColor(ROOT.kRed)
hctbSMIH.SetMarkerColor(ROOT.kRed)
hctbSMIH.SetMarkerStyle(20)

hctbSMIH0.Sumw2()
hctbSMIH0.SetLineWidth(2)
hctbSMIH0.SetLineColor(ROOT.kAzure+9)
hctbSMIH0.SetMarkerColor(ROOT.kAzure+9)
hctbSMIH0.SetMarkerStyle(24)

hctbH.Sumw2()
hctbH.SetLineWidth(2)
hctbH.SetLineColor(ROOT.kViolet+1)
hctbH.SetMarkerColor(ROOT.kViolet+1)
hctbH.SetMarkerStyle(22)

hctbA.Sumw2()
hctbA.SetLineWidth(2)
hctbA.SetLineColor(ROOT.kViolet+1)
hctbA.SetMarkerColor(ROOT.kViolet+1)
hctbA.SetMarkerStyle(23)

hctbA0.Sumw2()
hctbA0.SetLineWidth(2)
hctbA0.SetLineColor(ROOT.kGreen+1)
hctbA0.SetMarkerColor(ROOT.kGreen+1)
hctbA0.SetMarkerStyle(32)

hctbH0.Sumw2()
hctbH0.SetLineWidth(2)
hctbH0.SetLineColor(ROOT.kGreen+1)
hctbH0.SetMarkerColor(ROOT.kGreen+1)
hctbH0.SetMarkerStyle(26)



### cos(theta*) Collins-Soper histos
hctcSM    = TH1D("hctcSM",    ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)
hctcSMIH  = TH1D("hctcSMIH",  ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)
hctcSMIH0 = TH1D("hctcSMIH0", ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)
hctcH     = TH1D("hctcH",     ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)
hctcA     = TH1D("hctcA",     ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)
hctcA0    = TH1D("hctcA0",    ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)
hctcH0    = TH1D("hctcH0",    ";Truth cos#it{#theta}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-1.,+1.)

hctcSM.Sumw2()
hctcSM.SetLineWidth(2)
hctcSM.SetLineColor(ROOT.kBlack)
hctcSM.SetMarkerColor(ROOT.kBlack)
hctcSM.SetMarkerStyle(24)

hctcSMIH.Sumw2()
hctcSMIH.SetLineWidth(2)
hctcSMIH.SetLineColor(ROOT.kRed)
hctcSMIH.SetMarkerColor(ROOT.kRed)
hctcSMIH.SetMarkerStyle(20)

hctcSMIH0.Sumw2()
hctcSMIH0.SetLineWidth(2)
hctcSMIH0.SetLineColor(ROOT.kAzure+9)
hctcSMIH0.SetMarkerColor(ROOT.kAzure+9)
hctcSMIH0.SetMarkerStyle(24)

hctcH.Sumw2()
hctcH.SetLineWidth(2)
hctcH.SetLineColor(ROOT.kViolet+1)
hctcH.SetMarkerColor(ROOT.kViolet+1)
hctcH.SetMarkerStyle(22)

hctcA.Sumw2()
hctcA.SetLineWidth(2)
hctcA.SetLineColor(ROOT.kViolet+1)
hctcA.SetMarkerColor(ROOT.kViolet+1)
hctcA.SetMarkerStyle(23)

hctcA0.Sumw2()
hctcA0.SetLineWidth(2)
hctcA0.SetLineColor(ROOT.kGreen+1)
hctcA0.SetMarkerColor(ROOT.kGreen+1)
hctcA0.SetMarkerStyle(32)

hctcH0.Sumw2()
hctcH0.SetLineWidth(2)
hctcH0.SetLineColor(ROOT.kGreen+1)
hctcH0.SetMarkerColor(ROOT.kGreen+1)
hctcH0.SetMarkerStyle(26)





### phi* Helicity histos
hpbSM    = TH1D("hpbSM",    ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpbSMIH  = TH1D("hpbSMIH",  ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpbSMIH0 = TH1D("hpbSMIH0", ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpbH     = TH1D("hpbH",     ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpbA     = TH1D("hpbA",     ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpbA0    = TH1D("hpbA0",    ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpbH0    = TH1D("hpbH0",    ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Helicity frame;Events",80,-TMath.Pi(),+TMath.Pi())

hpbSM.Sumw2()
hpbSM.SetLineWidth(2)
hpbSM.SetLineColor(ROOT.kBlack)
hpbSM.SetMarkerColor(ROOT.kBlack)
hpbSM.SetMarkerStyle(24)

hpbSMIH.Sumw2()
hpbSMIH.SetLineWidth(2)
hpbSMIH.SetLineColor(ROOT.kRed)
hpbSMIH.SetMarkerColor(ROOT.kRed)
hpbSMIH.SetMarkerStyle(20)

hpbSMIH0.Sumw2()
hpbSMIH0.SetLineWidth(2)
hpbSMIH0.SetLineColor(ROOT.kAzure+9)
hpbSMIH0.SetMarkerColor(ROOT.kAzure+9)
hpbSMIH0.SetMarkerStyle(24)

hpbH.Sumw2()
hpbH.SetLineWidth(2)
hpbH.SetLineColor(ROOT.kViolet+1)
hpbH.SetMarkerColor(ROOT.kViolet+1)
hpbH.SetMarkerStyle(22)

hpbA.Sumw2()
hpbA.SetLineWidth(2)
hpbA.SetLineColor(ROOT.kViolet+1)
hpbA.SetMarkerColor(ROOT.kViolet+1)
hpbA.SetMarkerStyle(23)

hpbA0.Sumw2()
hpbA0.SetLineWidth(2)
hpbA0.SetLineColor(ROOT.kGreen+1)
hpbA0.SetMarkerColor(ROOT.kGreen+1)
hpbA0.SetMarkerStyle(32)

hpbH0.Sumw2()
hpbH0.SetLineWidth(2)
hpbH0.SetLineColor(ROOT.kGreen+1)
hpbH0.SetMarkerColor(ROOT.kGreen+1)
hpbH0.SetMarkerStyle(26)




### phi* Collins-Soper histos
hpcSM    = TH1D("hpcSM",    ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpcSMIH  = TH1D("hpcSMIH",  ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpcSMIH0 = TH1D("hpcSMIH0", ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpcH     = TH1D("hpcH",     ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpcA     = TH1D("hpcA",     ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpcA0    = TH1D("hpcA0",    ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())
hpcH0    = TH1D("hpcH0",    ";Truth #it{#phi}_{#it{t}#bar{#it{t}}} Collins-Soper frame;Events",80,-TMath.Pi(),+TMath.Pi())

hpcSM.Sumw2()
hpcSM.SetLineWidth(2)
hpcSM.SetLineColor(ROOT.kBlack)
hpcSM.SetMarkerColor(ROOT.kBlack)
hpcSM.SetMarkerStyle(24)

hpcSMIH.Sumw2()
hpcSMIH.SetLineWidth(2)
hpcSMIH.SetLineColor(ROOT.kRed)
hpcSMIH.SetMarkerColor(ROOT.kRed)
hpcSMIH.SetMarkerStyle(20)

hpcSMIH0.Sumw2()
hpcSMIH0.SetLineWidth(2)
hpcSMIH0.SetLineColor(ROOT.kAzure+9)
hpcSMIH0.SetMarkerColor(ROOT.kAzure+9)
hpcSMIH0.SetMarkerStyle(24)

hpcH.Sumw2()
hpcH.SetLineWidth(2)
hpcH.SetLineColor(ROOT.kViolet+1)
hpcH.SetMarkerColor(ROOT.kViolet+1)
hpcH.SetMarkerStyle(22)

hpcA.Sumw2()
hpcA.SetLineWidth(2)
hpcA.SetLineColor(ROOT.kViolet+1)
hpcA.SetMarkerColor(ROOT.kViolet+1)
hpcA.SetMarkerStyle(23)

hpcA0.Sumw2()
hpcA0.SetLineWidth(2)
hpcA0.SetLineColor(ROOT.kGreen+1)
hpcA0.SetMarkerColor(ROOT.kGreen+1)
hpcA0.SetMarkerStyle(32)

hpcH0.Sumw2()
hpcH0.SetLineWidth(2)
hpcH0.SetLineColor(ROOT.kGreen+1)
hpcH0.SetMarkerColor(ROOT.kGreen+1)
hpcH0.SetMarkerStyle(26)




### pT histos
hpt1SM    = TH1D("hpt1SM",    ";Truth #it{p}_{T}^{#it{t}} [GeV];Events",25,0,500)
hpt1SMIH  = TH1D("hpt1SMIH",  ";Truth #it{p}_{T}^{#it{t}} [GeV];Events",25,0,500)
hpt1SMIH0 = TH1D("hpt1SMIH0", ";Truth #it{p}_{T}^{#it{t}} [GeV];Events",25,0,500)

hpt2SM    = TH1D("hpt2SM",    ";Truth #it{p}_{T}^{#bar{#it{t}}} [GeV];Events",25,0,500)
hpt2SMIH  = TH1D("hpt2SMIH",  ";Truth #it{p}_{T}^{#bar{#it{t}}} [GeV];Events",25,0,500)
hpt2SMIH0 = TH1D("hpt2SMIH0", ";Truth #it{p}_{T}^{#bar{#it{t}}} [GeV];Events",25,0,500)

hpt1SM.Sumw2()
hpt1SM.SetLineWidth(2)
hpt1SM.SetLineColor(ROOT.kBlack)
hpt1SM.SetMarkerColor(ROOT.kBlack)
hpt1SM.SetMarkerStyle(24)
hpt2SM.Sumw2()
hpt2SM.SetLineWidth(2)
hpt2SM.SetLineColor(ROOT.kBlack)
hpt2SM.SetMarkerColor(ROOT.kBlack)
hpt2SM.SetMarkerStyle(24)

hpt1SMIH.Sumw2()
hpt1SMIH.SetLineWidth(2)
hpt1SMIH.SetLineColor(ROOT.kRed)
hpt1SMIH.SetMarkerColor(ROOT.kRed)
hpt1SMIH.SetMarkerStyle(20)
hpt2SMIH.Sumw2()
hpt2SMIH.SetLineWidth(2)
hpt2SMIH.SetLineColor(ROOT.kRed)
hpt2SMIH.SetMarkerColor(ROOT.kRed)
hpt2SMIH.SetMarkerStyle(20)

hpt1SMIH0.Sumw2()
hpt1SMIH0.SetLineWidth(2)
hpt1SMIH0.SetLineColor(ROOT.kAzure+9)
hpt1SMIH0.SetMarkerColor(ROOT.kAzure+9)
hpt1SMIH0.SetMarkerStyle(24)
hpt2SMIH0.Sumw2()
hpt2SMIH0.SetLineWidth(2)
hpt2SMIH0.SetLineColor(ROOT.kAzure+9)
hpt2SMIH0.SetMarkerColor(ROOT.kAzure+9)
hpt2SMIH0.SetMarkerStyle(24)



### eta histos
heta1SM    = TH1D("heta1SM",    ";Truth #it{#eta}_{#it{t}};Events",25,-5,+5)
heta1SMIH  = TH1D("heta1SMIH",  ";Truth #it{#eta}_{#it{t}};Events",25,-5,+5)
heta1SMIH0 = TH1D("heta1SMIH0", ";Truth #it{#eta}_{#it{t}};Events",25,-5,+5)
                                                  
heta2SM    = TH1D("heta2SM",    ";Truth #it{#eta}_{#bar{#it{t}}};Events",25,-5,+5)
heta2SMIH  = TH1D("heta2SMIH",  ";Truth #it{#eta}_{#bar{#it{t}}};Events",25,-5,+5)
heta2SMIH0 = TH1D("heta2SMIH0", ";Truth #it{#eta}_{#bar{#it{t}}};Events",25,-5,+5)

heta1SM.Sumw2()
heta1SM.SetLineWidth(2)
heta1SM.SetLineColor(ROOT.kBlack)
heta1SM.SetMarkerColor(ROOT.kBlack)
heta1SM.SetMarkerStyle(24)
heta2SM.Sumw2()
heta2SM.SetLineWidth(2)
heta2SM.SetLineColor(ROOT.kBlack)
heta2SM.SetMarkerColor(ROOT.kBlack)
heta2SM.SetMarkerStyle(24)

heta1SMIH.Sumw2()
heta1SMIH.SetLineWidth(2)
heta1SMIH.SetLineColor(ROOT.kRed)
heta1SMIH.SetMarkerColor(ROOT.kRed)
heta1SMIH.SetMarkerStyle(20)
heta2SMIH.Sumw2()
heta2SMIH.SetLineWidth(2)
heta2SMIH.SetLineColor(ROOT.kRed)
heta2SMIH.SetMarkerColor(ROOT.kRed)
heta2SMIH.SetMarkerStyle(20)

heta1SMIH0.Sumw2()
heta1SMIH0.SetLineWidth(2)
heta1SMIH0.SetLineColor(ROOT.kAzure+9)
heta1SMIH0.SetMarkerColor(ROOT.kAzure+9)
heta1SMIH0.SetMarkerStyle(24)
heta2SMIH0.Sumw2()
heta2SMIH0.SetLineWidth(2)
heta2SMIH0.SetLineColor(ROOT.kAzure+9)
heta2SMIH0.SetMarkerColor(ROOT.kAzure+9)
heta2SMIH0.SetMarkerStyle(24)



### phi histos
hphi1SM    = TH1D("hphi1SM",    ";Truth #it{#phi}_{#it{t}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi1SMIH  = TH1D("hphi1SMIH",  ";Truth #it{#phi}_{#it{t}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi1SMIH0 = TH1D("hphi1SMIH0", ";Truth #it{#phi}_{#it{t}};Events",25,-TMath.Pi(),+TMath.Pi())
                                                  
hphi2SM    = TH1D("hphi2SM",    ";Truth #it{#phi}_{#bar{#it{t}}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi2SMIH  = TH1D("hphi2SMIH",  ";Truth #it{#phi}_{#bar{#it{t}}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi2SMIH0 = TH1D("hphi2SMIH0", ";Truth #it{#phi}_{#bar{#it{t}}};Events",25,-TMath.Pi(),+TMath.Pi())

hphi1SM.Sumw2()
hphi1SM.SetLineWidth(2)
hphi1SM.SetLineColor(ROOT.kBlack)
hphi1SM.SetMarkerColor(ROOT.kBlack)
hphi1SM.SetMarkerStyle(24)
hphi2SM.Sumw2()
hphi2SM.SetLineWidth(2)
hphi2SM.SetLineColor(ROOT.kBlack)
hphi2SM.SetMarkerColor(ROOT.kBlack)
hphi2SM.SetMarkerStyle(24)

hphi1SMIH.Sumw2()
hphi1SMIH.SetLineWidth(2)
hphi1SMIH.SetLineColor(ROOT.kRed)
hphi1SMIH.SetMarkerColor(ROOT.kRed)
hphi1SMIH.SetMarkerStyle(20)
hphi2SMIH.Sumw2()
hphi2SMIH.SetLineWidth(2)
hphi2SMIH.SetLineColor(ROOT.kRed)
hphi2SMIH.SetMarkerColor(ROOT.kRed)
hphi2SMIH.SetMarkerStyle(20)

hphi1SMIH0.Sumw2()
hphi1SMIH0.SetLineWidth(2)
hphi1SMIH0.SetLineColor(ROOT.kAzure+9)
hphi1SMIH0.SetMarkerColor(ROOT.kAzure+9)
hphi1SMIH0.SetMarkerStyle(24)
hphi2SMIH0.Sumw2()
hphi2SMIH0.SetLineWidth(2)
hphi2SMIH0.SetLineColor(ROOT.kAzure+9)
hphi2SMIH0.SetMarkerColor(ROOT.kAzure+9)
hphi2SMIH0.SetMarkerStyle(24)




fSM   = TFile("tops.SM.500k.root","READ")
tSM   = fSM.Get("SM")
fSMIH = TFile("tops.SMIH.root","READ")
tSMIH = fSMIH.Get("SMIH")
fA    = TFile("tops.A.root","READ")
tA    = fA.Get("A")
fH    = TFile("tops.H.root","READ")
tH    = fH.Get("H")

pSM   = ROOT.vector(TLorentzVector)()
iSM   = std.vector(int)()
pSMIH = ROOT.vector(TLorentzVector)()
iSMIH = std.vector(int)()
pA    = ROOT.vector(TLorentzVector)()
iA    = std.vector(int)()
pH    = ROOT.vector(TLorentzVector)()
iH    = std.vector(int)()

tSM.SetBranchAddress("p4",pSM);
tSM.SetBranchAddress("id",iSM);
tSMIH.SetBranchAddress("p4",pSMIH);
tSMIH.SetBranchAddress("id",iSMIH);
tA.SetBranchAddress("p4",pA);
tA.SetBranchAddress("id",iA);
tH.SetBranchAddress("p4",pH);
tH.SetBranchAddress("id",iH);


n=1
for event in tSM:
   if(n%50000==0): print "processed |SM|^2 and reweighting ", n
   g1=event.p4[0]
   g2=event.p4[1]
   t1=event.p4[2]
   t2=event.p4[3]
   id1=event.id[2]
   id2=event.id[3]
   mtt = (t1+t2).M()
   pT1 = 0.
   pT2 = 0.
   eta1 = -100.
   eta2 = -100.
   phi1 = -100.
   phi2 = -100.
   if(id1>0):
      pT1 = t1.Pt()
      pT2 = t2.Pt()
      eta1 = t1.Eta()
      eta2 = t2.Eta()
      phi1 = t1.Phi()
      phi2 = t2.Phi()
   else:
      pT1 = t2.Pt()
      pT2 = t1.Pt()
      eta1 = t2.Eta()
      eta2 = t1.Eta()
      phi1 = t2.Phi()
      phi2 = t1.Phi()

   p = [[ g1.E(), g1.Px(), g1.Py(), g1.Pz() ],
        [ g2.E(), g2.Px(), g2.Py(), g2.Pz() ],
        [ t1.E(), t1.Px(), t1.Py(), t1.Pz() ],
        [ t2.E(), t2.Px(), t2.Py(), t2.Pz() ]]
   P=invert_momenta(p)

   me2SM   = matrix2SMpy.get_me(P,alphaS,nhel)
   me2SMIH = matrix2SMIHpy.get_me(P,alphaS,nhel)
   me2H    = matrix2Hpy.get_me(P,alphaS,nhel)
   me2A    = matrix2Apy.get_me(P,alphaS,nhel)

   weightSMIH = me2SMIH/me2SM
   weightH    = me2H/me2SM
   weightA    = me2A/me2SM

   # print "mtt="+str(mtt)+": |SM|^2="+str(me2SM)+" |SM+A|^2="+str(me2SMIH)

   hmSM.Fill(mtt)
   hmSMx.Fill(mtt)
   hmSMIH.Fill(mtt,weightSMIH)
   hmH.Fill(mtt,weightH)
   hmA.Fill(mtt,weightA)
   # if(TMath.Abs(mtt-500.)<100.):
   #    charge1 = id1
   #    cosThHE = kinematics.CostHE(t1,charge1,t2)
   #    cosThCS = kinematics.CostCS(t1,charge1,t2)
   #    phiHE   = kinematics.PhiHE(t1,charge1,t2)
   #    phiCS   = kinematics.PhiCS(t1,charge1,t2)
   #    hctbSM.Fill(cosThHE)
   #    hctbSMIH.Fill(cosThHE,weightSMIH)
   #    hctbH.Fill(cosThHE,weightH)
   #    # hctbA.Fill(cosThHE,weightA)
   #    hctcSM.Fill(cosThCS)
   #    hctcSMIH.Fill(cosThCS,weightSMIH)
   #    hctcH.Fill(cosThCS,weightH)
   #    # hctcA.Fill(cosThCS,weightA)
   #    hpbSM.Fill(phiHE)
   #    hpbSMIH.Fill(phiHE,weightSMIH)
   #    hpbH.Fill(phiHE,weightH)
   #    # hpbA.Fill(phiHE,weightA)
   #    hpcSM.Fill(phiCS)
   #    hpcSMIH.Fill(phiCS,weightSMIH)
   #    hpcH.Fill(phiCS,weightH)
   #    # hpcA.Fill(phiCS,weightA)
   charge1 = id1
   cosThHE = kinematics.CostHE(t1,charge1,t2)
   cosThCS = kinematics.CostCS(t1,charge1,t2)
   phiHE   = kinematics.PhiHE(t1,charge1,t2)
   phiCS   = kinematics.PhiCS(t1,charge1,t2)
   hctbSM.Fill(cosThHE)
   hctbSMIH.Fill(cosThHE,weightSMIH)
   hctbH.Fill(cosThHE,weightH)
   hctbA.Fill(cosThHE,weightA)
   hctcSM.Fill(cosThCS)
   hctcSMIH.Fill(cosThCS,weightSMIH)
   hctcH.Fill(cosThCS,weightH)
   hctcA.Fill(cosThCS,weightA)
   hpbSM.Fill(phiHE)
   hpbSMIH.Fill(phiHE,weightSMIH)
   hpbH.Fill(phiHE,weightH)
   hpbA.Fill(phiHE,weightA)
   hpcSM.Fill(phiCS)
   hpcSMIH.Fill(phiCS,weightSMIH)
   hpcH.Fill(phiCS,weightH)
   hpcA.Fill(phiCS,weightA)

   hpt1SM.Fill(pT1)
   hpt1SMIH.Fill(pT1,weightSMIH)
   hpt2SM.Fill(pT2)
   hpt2SMIH.Fill(pT2,weightSMIH)
   heta1SM.Fill(eta1)
   heta1SMIH.Fill(eta1,weightSMIH)
   heta2SM.Fill(eta2)
   heta2SMIH.Fill(eta2,weightSMIH)
   hphi1SM.Fill(phi1)
   hphi1SMIH.Fill(phi1,weightSMIH)
   hphi2SM.Fill(phi2)
   hphi2SMIH.Fill(phi2,weightSMIH)
   n+=1

n=1
for event in tSMIH:
   if(n%50000==0): print "processed |SM+A|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   id1=event.id[2]
   id2=event.id[3]
   mtt = (t1+t2).M()
   pT1 = 0.
   pT2 = 0.
   eta1 = 0.
   eta2 = 0.
   phi1 = 0.
   phi2 = 0.
   if(id1>0):
      pT1 = t1.Pt()
      pT2 = t2.Pt()
      eta1 = t1.Eta()
      eta2 = t2.Eta()
      phi1 = t1.Phi()
      phi2 = t2.Phi()
   else:
      pT1 = t2.Pt()
      pT2 = t1.Pt()
      eta1 = t2.Eta()
      eta2 = t1.Eta()
      phi1 = t2.Phi()
      phi2 = t1.Phi()

   hmSMIH0.Fill(mtt)
   # if(TMath.Abs(mtt-500.)<100.):
   #    charge1 = id1
   #    cosThHE = kinematics.CostHE(t1,charge1,t2)
   #    cosThCS = kinematics.CostCS(t1,charge1,t2)
   #    phiHE = kinematics.PhiHE(t1,charge1,t2)
   #    phiCS = kinematics.PhiCS(t1,charge1,t2)
   #    hctbSMIH0.Fill(cosThHE)
   #    hctcSMIH0.Fill(cosThCS)
   #    hpbSMIH0.Fill(phiHE)
   #    hpcSMIH0.Fill(phiCS)
   charge1 = id1
   cosThHE = kinematics.CostHE(t1,charge1,t2)
   cosThCS = kinematics.CostCS(t1,charge1,t2)
   phiHE = kinematics.PhiHE(t1,charge1,t2)
   phiCS = kinematics.PhiCS(t1,charge1,t2)
   hctbSMIH0.Fill(cosThHE)
   hctcSMIH0.Fill(cosThCS)
   hpbSMIH0.Fill(phiHE)
   hpcSMIH0.Fill(phiCS)

   hpt1SMIH0.Fill(pT1)
   hpt2SMIH0.Fill(pT2)
   heta1SMIH0.Fill(eta1)
   heta2SMIH0.Fill(eta2)  
   hphi1SMIH0.Fill(phi1)
   hphi2SMIH0.Fill(phi2)
   n+=1

n=1
for event in tA:
   if(n%50000==0): print "processed |A|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   id1=event.id[2]
   id2=event.id[3]
   mtt = (t1+t2).M()
   charge1 = id1
   cosThHE = kinematics.CostHE(t1,charge1,t2)
   cosThCS = kinematics.CostCS(t1,charge1,t2)
   phiHE = kinematics.PhiHE(t1,charge1,t2)
   phiCS = kinematics.PhiCS(t1,charge1,t2)
   hmA0.Fill(mtt)
   # if(TMath.Abs(mtt-500)<100):
   #    hctbA0.Fill(cosThHE)
   #    hctcA0.Fill(cosThCS)
   #    hpbA0.Fill(phiHE)
   #    hpcA0.Fill(phiCS)
   hctbA0.Fill(cosThHE)
   hctcA0.Fill(cosThCS)
   hpbA0.Fill(phiHE)
   hpcA0.Fill(phiCS)
   n+=1

n=1
for event in tH:
   if(n%50000==0): print "processed |H|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   id1=event.id[2]
   id2=event.id[3]
   mtt = (t1+t2).M()
   hmH0.Fill(mtt)
   # if(TMath.Abs(mtt-500.)<100.):
   #    charge1 = id1
   #    cosThHE = kinematics.CostHE(t1,charge1,t2)
   #    cosThCS = kinematics.CostCS(t1,charge1,t2)
   #    phiHE = kinematics.PhiHE(t1,charge1,t2)
   #    phiCS = kinematics.PhiCS(t1,charge1,t2)
   #    hctbH0.Fill(cosThHE)
   #    hctcH0.Fill(cosThCS)
   #    hpbH0.Fill(phiHE)
   #    hpcH0.Fill(phiCS)
   charge1 = id1
   cosThHE = kinematics.CostHE(t1,charge1,t2)
   cosThCS = kinematics.CostCS(t1,charge1,t2)
   phiHE = kinematics.PhiHE(t1,charge1,t2)
   phiCS = kinematics.PhiCS(t1,charge1,t2)
   hctbH0.Fill(cosThHE)
   hctcH0.Fill(cosThCS)
   hpbH0.Fill(phiHE)
   hpcH0.Fill(phiCS)

   n+=1




plot(hmSM,hmSMIH,hmSMIH0,"2HDM.reweighting.pdf(")
# plot(hctbSM,hctbSMIH,hctbSMIH0,"2HDM.reweighting.pdf")
# plot(hctcSM,hctcSMIH,hctcSMIH0,"2HDM.reweighting.pdf",-1,45000)
# plot(hpbSM,hpbSMIH,hpbSMIH0,"2HDM.reweighting.pdf")
# plot(hpcSM,hpcSMIH,hpcSMIH0,"2HDM.reweighting.pdf")
plot(hpt1SM,hpt1SMIH,hpt1SMIH0,"2HDM.reweighting.pdf")
plot(hpt2SM,hpt2SMIH,hpt2SMIH0,"2HDM.reweighting.pdf")
plot(heta1SM,heta1SMIH,heta1SMIH0,"2HDM.reweighting.pdf",-1,60000)
plot(heta2SM,heta2SMIH,heta2SMIH0,"2HDM.reweighting.pdf",-1,60000)
plot(hphi1SM,hphi1SMIH,hphi1SMIH0,"2HDM.reweighting.pdf",10000,35000)
plot(hphi2SM,hphi2SMIH,hphi2SMIH0,"2HDM.reweighting.pdf",10000,35000)


plotnorm(hmSMx,hmH,hmH0,"2HDM.reweighting.pdf","H",True)
plotnorm(hctbSM,hctbH,hctbH0,"2HDM.reweighting.pdf","H")
plotnorm(hctcSM,hctcH,hctcH0,"2HDM.reweighting.pdf","H")
# plotnorm(hpbSM,hpbH,hpbH0,"2HDM.reweighting.pdf","H")
# plotnorm(hpcSM,hpcH,hpcH0,"2HDM.reweighting.pdf","H")
plotnorm(hmSMx,hmA,hmA0,"2HDM.reweighting.pdf","A",True)
plotnorm(hctbSM,hctbA,hctbA0,"2HDM.reweighting.pdf","A")
plotnorm(hctcSM,hctcA,hctcA0,"2HDM.reweighting.pdf)","A")
# plotnorm(hpbSM,hpbH,hpbH0,"2HDM.reweighting.pdf","A")
# plotnorm(hpcSM,hpcH,hpcH0,"2HDM.reweighting.pdf)","A")
