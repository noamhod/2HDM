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



fSM   = TFile("tops.SM.500k.root","READ")
tSM   = fSM.Get("SM")
pSM   = ROOT.vector(TLorentzVector)()
iSM   = std.vector(int)()
tSM.SetBranchAddress("p4",pSM);
tSM.SetBranchAddress("id",iSM);


#### get the model
type      = THDM.model.type
sba       = THDM.model.sba
mX        = THDM.model.mX
nameX     = THDM.model.nameX
cuts      = THDM.model.cuts
mgpath    = THDM.model.mgpath
alphaS    = THDM.model.alphaS
nhel      = THDM.model.nhel
libmatrix = "matrix/"+nameX+"/"+str(mX)+"/"+str(sba)+"/"
THDM.setParameters(nameX,mX,cuts,type,sba)




THDM.setModules(os.getcwd(),libmatrix,nameX,len(THDM.parameters),"All")
print THDM.modules

### output tree
Nk = tSM.GetEntries()
sNk = str(Nk/1000)+"k"
outFile = "2HDM."+nameX+"."+str(mX)+"GeV.tree."+sNk+".root" 
newFile = TFile(outFile,"RECREATE") 
tnew = tSM.CloneTree(0)
wgt  = std.vector(float)()
tnb  = std.vector(float)()
wdtA = std.vector(float)()
wdtH = std.vector(float)()
ymt   = std.vector(float)()
ymb   = std.vector(float)()
ymc   = std.vector(float)()
ymtau = std.vector(float)()
ymmu  = std.vector(float)()
tnew.Branch("wgt",wgt)
tnew.Branch("tnb",tnb)
tnew.Branch("wdtA",wdtA)
tnew.Branch("wdtH",wdtH)
tnew.Branch("ymt",ymt)
tnew.Branch("ymb",ymb)
tnew.Branch("ymc",ymc)
tnew.Branch("ymtau",ymtau)
tnew.Branch("ymmu",ymmu)
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
   me2SM = THDM.modules['matrix2SMpy'].get_me(P,alphaS,nhel)  ### calculate the SM ME^2

   ### clear decoration branches
   wgt.clear()
   tnb.clear()
   wdtA.clear()
   wdtH.clear()

   ymt.clear()
   ymb.clear()
   ymc.clear()
   ymtau.clear()
   ymmu.clear()

   ### loop over all model points
   for i in range(0,len(THDM.parameters)):

      ## the ME^2 and the weight
      me2XX = THDM.modules['matrix2'+nameX+str(i)+'py'].get_me(P,alphaS,nhel) ### calculate the X ME^2
      weightX = me2XX/me2SM                                                   ### calculate the weight

      # print "[event "+str(n)+"][point "+str(i)+"] m="+str(mtt)+" tanb="+str(THDM.parameters[i].get("tanb"))+" --> weight="+str(weightX)
      wgt.push_back(weightX)
      tnb.push_back(THDM.parameters[i].get("tanb"))
      wdtA.push_back(THDM.parameters[i].get("wA"))
      wdtH.push_back(THDM.parameters[i].get("wH"))
      ymt.push_back(THDM.parameters[i].get("YMT"))
      ymb.push_back(THDM.parameters[i].get("YMB"))
      ymc.push_back(THDM.parameters[i].get("YMC"))
      ymtau.push_back(THDM.parameters[i].get("YMTAU"))
      ymmu.push_back(THDM.parameters[i].get("YMM"))
      
   tnew.Fill()
   n+=1

# use GetCurrentFile just in case we went over the
# (customizable) maximum file size
tnew.GetCurrentFile().Write() 
tnew.GetCurrentFile().Close()