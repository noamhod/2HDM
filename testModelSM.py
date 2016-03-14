#!/usr/bin/python
import ROOT
from ROOT import std, gROOT, gStyle, gPad, TCanvas, TH1, TH1D, TH2D, TLegend, TLine, TFile, TTree, TLorentzVector, TMath, TVirtualPad, TEventList

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


setStyle()

fMG = TFile("tops.MG_SM.500k.root","READ")
tMG = fMG.Get("MG_SM")
pMG = ROOT.vector(TLorentzVector)()
iMG = std.vector(int)()
tMG.SetBranchAddress("p4",pMG);
tMG.SetBranchAddress("id",iMG);

fSM = TFile("tops.SM.500k.root","READ")
tSM = fSM.Get("SM")
pSM = ROOT.vector(TLorentzVector)()
iSM = std.vector(int)()
tSM.SetBranchAddress("p4",pSM);
tSM.SetBranchAddress("id",iSM);

hmMG = TH1D("hmMG", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmMG.Sumw2()
hmMG.SetLineWidth(2)
hmMG.SetLineColor(ROOT.kBlack)
hmMG.SetMarkerColor(ROOT.kBlack)
hmMG.SetMarkerStyle(20)

hmSM = TH1D("hmSMX", ";Truth #it{m}_{#it{t}#bar{#it{t}}} [GeV];Events",25,350,850)
hmSM.Sumw2()
hmSM.SetLineWidth(2)
hmSM.SetLineColor(ROOT.kAzure+9)
hmSM.SetMarkerColor(ROOT.kAzure+9)
hmSM.SetMarkerStyle(24)


### pT histos
hpt1MG  = TH1D("hpt1SM",    ";Truth #it{p}_{T}^{#it{t}} [GeV];Events",25,0,500)
hpt1SM  = TH1D("hpt1SMIH",  ";Truth #it{p}_{T}^{#it{t}} [GeV];Events",25,0,500)
hpt2MG  = TH1D("hpt2SM",    ";Truth #it{p}_{T}^{#bar{#it{t}}} [GeV];Events",25,0,500)
hpt2SM  = TH1D("hpt2SMIH",  ";Truth #it{p}_{T}^{#bar{#it{t}}} [GeV];Events",25,0,500)

hpt1MG.Sumw2()
hpt1MG.SetLineWidth(2)
hpt1MG.SetLineColor(ROOT.kBlack)
hpt1MG.SetMarkerColor(ROOT.kBlack)
hpt1MG.SetMarkerStyle(24)
hpt2MG.Sumw2()
hpt2MG.SetLineWidth(2)
hpt2MG.SetLineColor(ROOT.kBlack)
hpt2MG.SetMarkerColor(ROOT.kBlack)
hpt2MG.SetMarkerStyle(24)

hpt1SM.Sumw2()
hpt1SM.SetLineWidth(2)
hpt1SM.SetLineColor(ROOT.kAzure+9)
hpt1SM.SetMarkerColor(ROOT.kAzure+9)
hpt1SM.SetMarkerStyle(20)
hpt2SM.Sumw2()
hpt2SM.SetLineWidth(2)
hpt2SM.SetLineColor(ROOT.kAzure+9)
hpt2SM.SetMarkerColor(ROOT.kAzure+9)
hpt2SM.SetMarkerStyle(20)


### eta histos
heta1MG  = TH1D("heta1SM",    ";Truth #it{#eta}_{#it{t}};Events",25,-5,+5)
heta1SM  = TH1D("heta1SMIH",  ";Truth #it{#eta}_{#it{t}};Events",25,-5,+5)
heta2MG  = TH1D("heta2SM",    ";Truth #it{#eta}_{#bar{#it{t}}};Events",25,-5,+5)
heta2SM  = TH1D("heta2SMIH",  ";Truth #it{#eta}_{#bar{#it{t}}};Events",25,-5,+5)

heta1MG.Sumw2()
heta1MG.SetLineWidth(2)
heta1MG.SetLineColor(ROOT.kBlack)
heta1MG.SetMarkerColor(ROOT.kBlack)
heta1MG.SetMarkerStyle(24)
heta2MG.Sumw2()
heta2MG.SetLineWidth(2)
heta2MG.SetLineColor(ROOT.kBlack)
heta2MG.SetMarkerColor(ROOT.kBlack)
heta2MG.SetMarkerStyle(24)

heta1SM.Sumw2()
heta1SM.SetLineWidth(2)
heta1SM.SetLineColor(ROOT.kAzure+9)
heta1SM.SetMarkerColor(ROOT.kAzure+9)
heta1SM.SetMarkerStyle(20)
heta2SM.Sumw2()
heta2SM.SetLineWidth(2)
heta2SM.SetLineColor(ROOT.kAzure+9)
heta2SM.SetMarkerColor(ROOT.kAzure+9)
heta2SM.SetMarkerStyle(20)


### phi histos
hphi1MG  = TH1D("hphi1SM",    ";Truth #it{#phi}_{#it{t}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi1SM  = TH1D("hphi1SMIH",  ";Truth #it{#phi}_{#it{t}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi2MG  = TH1D("hphi2SM",    ";Truth #it{#phi}_{#bar{#it{t}}};Events",25,-TMath.Pi(),+TMath.Pi())
hphi2SM  = TH1D("hphi2SMIH",  ";Truth #it{#phi}_{#bar{#it{t}}};Events",25,-TMath.Pi(),+TMath.Pi())

hphi1MG.Sumw2()
hphi1MG.SetLineWidth(2)
hphi1MG.SetLineColor(ROOT.kBlack)
hphi1MG.SetMarkerColor(ROOT.kBlack)
hphi1MG.SetMarkerStyle(24)
hphi1MG.SetMinimum(0)
hphi2MG.Sumw2()
hphi2MG.SetLineWidth(2)
hphi2MG.SetLineColor(ROOT.kBlack)
hphi2MG.SetMarkerColor(ROOT.kBlack)
hphi2MG.SetMarkerStyle(24)
hphi2MG.SetMinimum(0)

hphi1SM.Sumw2()
hphi1SM.SetLineWidth(2)
hphi1SM.SetLineColor(ROOT.kAzure+9)
hphi1SM.SetMarkerColor(ROOT.kAzure+9)
hphi1SM.SetMarkerStyle(20)
hphi2SM.Sumw2()
hphi1SM.SetMinimum(0)
hphi2SM.SetLineWidth(2)
hphi2SM.SetLineColor(ROOT.kAzure+9)
hphi2SM.SetMarkerColor(ROOT.kAzure+9)
hphi2SM.SetMarkerStyle(20)
hphi2SM.SetMinimum(0)



n=1
for event in tMG:
   if(n%50000==0): print "processed MG |SM|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()
   hmMG.Fill(mtt)

   id1=event.id[2]
   id2=event.id[3]
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
   hpt1MG.Fill(pT1)
   hpt2MG.Fill(pT2)
   heta1MG.Fill(eta1)
   heta2MG.Fill(eta2)
   hphi1MG.Fill(phi1)
   hphi2MG.Fill(phi2)

   n+=1

n=1
for event in tSM:
   if(n%50000==0): print "processed 2HDM |SM|^2 generated ", n
   t1=event.p4[2]
   t2=event.p4[3]
   mtt = (t1+t2).M()
   hmSM.Fill(mtt)
   id1=event.id[2]
   id2=event.id[3]
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
   hpt1SM.Fill(pT1)
   hpt2SM.Fill(pT2)
   heta1SM.Fill(eta1)
   heta2SM.Fill(eta2)
   hphi1SM.Fill(phi1)
   hphi2SM.Fill(phi2)

   n+=1

cnv = TCanvas("cnv","",600,600)
cnv.Draw()
cnv.cd()
hmMG.Draw()
hmSM.Draw("same")
cnv.Update()
cnv.RedrawAxis()
cnv.SaveAs("testModelSM.pdf(")

cnv = TCanvas("cnv","",1200,600)
cnv.Divide(2,1)
cnv.Draw()
cnv.cd(1)
hpt1MG.Draw()
hpt1SM.Draw("same")
cnv.cd(2)
hpt2MG.Draw()
hpt2SM.Draw("same")
cnv.Update()
cnv.RedrawAxis()
cnv.SaveAs("testModelSM.pdf")

cnv = TCanvas("cnv","",1200,600)
cnv.Divide(2,1)
cnv.Draw()
cnv.cd(1)
heta1MG.Draw()
heta1SM.Draw("same")
cnv.cd(2)
heta2MG.Draw()
heta2SM.Draw("same")
cnv.Update()
cnv.RedrawAxis()
cnv.SaveAs("testModelSM.pdf")

cnv = TCanvas("cnv","",1200,600)
cnv.Divide(2,1)
cnv.Draw()
cnv.cd(1)
hphi1MG.Draw()
hphi1SM.Draw("same")
cnv.cd(2)
hphi2MG.Draw()
hphi2SM.Draw("same")
cnv.Update()
cnv.RedrawAxis()
cnv.SaveAs("testModelSM.pdf)")
