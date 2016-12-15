#!/usr/bin/python

import ROOT
from ROOT import std, gROOT, TLorentzVector, TMath, TVector3

def cosThetaDitop(pa, pb):
   paMag = TMath.Sqrt( pa.Px()*pa.Px() + pa.Py()*pa.Py() + pa.Pz()*pa.Pz() )
   pbMag = TMath.Sqrt( pb.Px()*pb.Px() + pb.Py()*pb.Py() + pb.Pz()*pb.Pz() )
   costh = (pa.Px()*pb.Px() + pa.Py()*pb.Py() + pa.Pz()*pb.Pz()) / (paMag*pbMag)
   return costh

def QT(pa, pb):
   pTmp = pa+pb
   return pTmp.Perp()

def ySystem(pa, pb):
   pTmp = pa+pb
   return pTmp.Rapidity()

def magBetaSystem(pa, pb):
   pTmp = pa+pb
   return pTmp.Beta()

def systemBoostVector(pa, pb):
   pTmp = pa+pb
   return pTmp.BoostVector()

def Boost(pboost, p):
   pBoosted = p.Clone("")
   pBoosted.Boost(-1.*pBoost.BoostVector())
   return pBoosted

def Boost(pa, pb, p):
   pTmp = pa+pb
   pBoosted = p.Clone("")
   pBoosted.Boost(-1.*pTmp.BoostVector())
   return pBoosted

def betaiSystem(pa, pb, i):
   pTmp = pa+pb
   E = pTmp.E()
   if(E<=0.):  return -99999999.
   if  (i==1): return pTmp.Px()/E;
   elif(i==2): return pTmp.Py()/E;
   elif(i==3): return pTmp.Pz()/E;
   return -99999999.

def betaTSystem(pa, pb):
   pTmp = pa+pb
   beta = 0
   if(pTmp.E()!=0.): return pTmp.Pt()/pTmp.E()
   return -999.

def magBetaSystem(pa, pb, ea, eb):
   if((ea+eb)>0.): return (pa+pb)/(ea+eb)
   return -999999.

def cosphi(reference1, target1, reference2, target2):
   targettmp1 = target1.Clone("tmp1")
   targettmp2 = target2.Clone("tmp2")
   vboost1 = reference1.BoostVector()
   targettmp1.Boost(-vboost1)
   vboost2 = reference2.BoostVector()
   targettmp2.Boost(-vboost2)
   cosp = targettmp1.Vect()*targettmp2.Vect() # scalar product
   cosp /= (targettmp1.Vect().Mag()*targettmp2.Vect().Mag())
   return cosp

def costheta(reference, target):
   targettmp = target.Clone("tmp")
   vboost = reference.BoostVector()
   targettmp.Boost(-vboost)
   cost = targettmp.Vect()*vboost # scalar product
   cost /= (targettmp.Vect().Mag()*vboost.Mag())
   return cost



def cosThetaTrue(pa, ida, pb, idb, doflip=False):
   #  http://xrootd.slac.stanford.edu/BFROOT/www/doc/workbook_backup_010108/analysis/analysis.html
   #  A useful quantity in many analyses is the helicity angle.
   #  In the reaction Y . X . a + b, the helicity angle of
   #  particle a is the angle measured in the rest frame of the
   #  deidaying parent particle, X, between the direction of the
   #  deiday daughter a and the direction of the grandparent particle Y.
   pTmp = pa+pb;                     # this is the mumu system (Z) 4vector
   ZboostVector = pTmp.BoostVector() # this is the 3vector of the Z
   p = TLorentzVector(0,0,-1,1)      # this is the muon 4vector
   if  (ida>0): p.SetPxPyPzE(pa.Px(),pa.Py(),pa.Pz(),pa.E())
   else:        p.SetPxPyPzE(pb.Px(),pb.Py(),pb.Pz(),pb.E())
   p.Boost(-ZboostVector) # boost p to the Ditop CM (rest) frame
   cosThetaB = p.Vect().Z()/p.Vect().Mag()
   if(ySystem(pa,pb)<0. and doflip): cosThetaB *= -1. # reclassify ???
   return cosThetaB

def cosThetaBoost(pa, ida, pb, idb, doflip=False):
   #  http://xrootd.slac.stanford.edu/BFROOT/www/doc/workbook_backup_010108/analysis/analysis.html
   #  A useful quantity in many analyses is the helicity angle.
   #  In the reaction Y . X . a + b, the helicity angle of
   #  particle a is the angle measured in the rest frame of the
   #  deidaying parent particle, X, between the direction of the
   #  deiday daughter a and the direction of the grandparent particle Y.
   pTmp = pa+pb;                     # this is the mumu system (Z) 4vector
   ZboostVector = pTmp.BoostVector() # this is the 3vector of the Z
   p = TLorentzVector(0,0,-1,1)      # this is the muon 4vector
   if  (ida>0): p.SetPxPyPzE(pa.Px(),pa.Py(),pa.Pz(),pa.E())
   else:        p.SetPxPyPzE(pb.Px(),pb.Py(),pb.Pz(),pb.E())
   p.Boost(-ZboostVector) # boost p to the Ditop CM (rest) frame
   cosThetaB = p.Vect()*pTmp.Vect()/(p.P()*pTmp.P())
   if(ySystem(pa,pb)<0. and doflip): cosThetaB *= -1. # reclassify ???
   return cosThetaB

def cosThetaCollinsSoper(pa, ida, pb, idb, doflip=True):
   # this will work only for leptons e, mu and tau
   # by default it is assumed that pa is the lepton
   # if instead pb is the lepton, then the result is
   # reclassified by a (-) sign - see line 4.
   mass = (pa+pb).M()
   mass2 = (pa+pb).M()*(pa+pb).M()
   QT2 = (pa+pb).Pt()*(pa+pb).Pt() # QT(pa,pb)*QT(pa,pb)
   # cosThetaCS = 2.*( pa.Plus()*pb.Minus() - pa.Minus()*pb.Plus() ) / TMath.Sqrt( mass2 * (mass2 + QT2) )
   cosThetaCS = 2.*( pa.Pz()*pb.E() - pa.E()*pb.Pz() ) / TMath.Sqrt( mass2 * (mass2 + QT2) )
   if(ida<0): cosThetaCS *= -1. # if pb is the lepton
   if(ySystem(pa,pb)<0 and doflip): cosThetaCS *= -1. # reclassify
   return cosThetaCS



def CostCS(p1, charge1, p2):
   #Cosine of the theta decay angle (top (Q=+2/3)) in the Collins-Soper frame
   pTop1CM = TLorentzVector(0,0,-1,1)  # In the CM. frame
   pTop2CM = TLorentzVector(0,0,-1,1)  # In the CM. frame
   pProjCM = TLorentzVector(0,0,-1,1) # In the CM. frame
   pTargCM = TLorentzVector(0,0,-1,1) # In the CM. frame
   pDitopCM = TLorentzVector(0,0,-1,1) # In the CM. frame
   pTop1Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTop2Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pProjDitop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTargDitop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   beta = TVector3(0,0,0)
   zaxisCS = TVector3(0,0,0)
   mp = 0.93827231
   ep = 6500.

   # Fill the Lorentz vector for projectile and target in the CM frame
   pProjCM.SetPxPyPzE(0.,0.,-ep,TMath.Sqrt(ep*ep+mp*mp))
   pTargCM.SetPxPyPzE(0.,0.,+ep,TMath.Sqrt(ep*ep+mp*mp))

   # Get the Topons parameters in the CM frame 
   pTop1CM.SetPxPyPzE(p1.Px(),p1.Py(),p1.Pz(),p1.E())
   pTop2CM.SetPxPyPzE(p2.Px(),p2.Py(),p2.Pz(),p2.E())

   # Obtain the Ditop parameters in the CM frame
   pDitopCM=pTop1CM+pTop2CM

   # Translate the Ditop parameters in the Ditop rest frame
   beta=(-1./pDitopCM.E())*pDitopCM.Vect()
   if(beta.Mag()>=1): return 666.
   pTop1Ditop=pTop1CM
   pTop2Ditop=pTop2CM
   pProjDitop=pProjCM
   pTargDitop=pTargCM
   pTop1Ditop.Boost(beta)
   pTop2Ditop.Boost(beta)
   pProjDitop.Boost(beta)
   pTargDitop.Boost(beta)

   # Determine the z axis for the CS angle 
   zaxisCS=(((pProjDitop.Vect()).Unit())-((pTargDitop.Vect()).Unit())).Unit();

   # Determine the CS angle (angle between Top+ and the z axis defined above)
   cost = -999

   if(charge1>0): cost = zaxisCS.Dot((pTop1Ditop.Vect()).Unit())
   else:          cost = zaxisCS.Dot((pTop2Ditop.Vect()).Unit())

   return cost


def CostHE(p1, charge1, p2):
   #Cosine of the theta decay angle (top (Q=+2/3)) in the Helicity frame
   pTop1CM = TLorentzVector(0,0,-1,1)   # In the CM frame 
   pTop2CM = TLorentzVector(0,0,-1,1)   # In the CM frame 
   pDitopCM = TLorentzVector(0,0,-1,1) # In the CM frame 
   pTop1Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTop2Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   beta = TVector3(0,0,0)
   zaxisCS = TVector3(0,0,0)

   # Get the muons parameters in the CM frame
   pTop1CM.SetPxPyPzE(p1.Px(),p1.Py(),p1.Pz(),p1.E())
   pTop2CM.SetPxPyPzE(p2.Px(),p2.Py(),p2.Pz(),p2.E())

   # Obtain the Ditop parameters in the CM frame
   pDitopCM=pTop1CM+pTop2CM

   # Translate the muon parameters in the Ditop rest frame
   beta=(-1./pDitopCM.E())*pDitopCM.Vect()
   if(beta.Mag()>=1): return 666.
   pTop1Ditop=pTop1CM
   pTop2Ditop=pTop2CM
   pTop1Ditop.Boost(beta)
   pTop2Ditop.Boost(beta)

   # Determine the z axis for the calculation of the polarization angle (i.e. the direction of the Ditop in the CM system)
   zaxis=(pDitopCM.Vect()).Unit()

   # Calculation of the polarization angle (angle between mu+ and the z axis defined above)
   cost = -999.
   if(charge1>0): cost = zaxis.Dot((pTop1Ditop.Vect()).Unit())
   else:          cost = zaxis.Dot((pTop2Ditop.Vect()).Unit())

   return cost


def PhiCS(p1, charge1, p2):
   # Phi decay angle (top (Q=+2/3)) in the Collins-Soper frame
   pTop1CM = TLorentzVector(0,0,-1,1)  # In the CM frame
   pTop2CM = TLorentzVector(0,0,-1,1)  # In the CM frame
   pProjCM = TLorentzVector(0,0,-1,1)  # In the CM frame
   pTargCM = TLorentzVector(0,0,-1,1)  # In the CM frame
   pDitopCM = TLorentzVector(0,0,-1,1) # In the CM frame
   pTop1Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTop2Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pProjDitop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTargDitop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   beta = TVector3(0,0,0)
   yaxisCS = TVector3(0,0,0)
   xaxisCS = TVector3(0,0,0)
   zaxisCS = TVector3(0,0,0)
   mp = 0.93827231
   ep = 6500.

   # Fill the Lorentz vector for projectile and target in the CM frame
   pProjCM.SetPxPyPzE(0.,0.,-ep,TMath.Sqrt(ep*ep+mp*mp))
   pTargCM.SetPxPyPzE(0.,0.,+ep,TMath.Sqrt(ep*ep+mp*mp))

   # Get the muons parameters in the CM frame 
   pTop1CM.SetPxPyPzE(p1.Px(),p1.Py(),p1.Pz(),p1.E())
   pTop2CM.SetPxPyPzE(p2.Px(),p2.Py(),p2.Pz(),p2.E())

   # Obtain the Ditop parameters in the CM frame
   pDitopCM=pTop1CM+pTop2CM

   # Translate the Ditop parameters in the Ditop rest frame
   beta=(-1./pDitopCM.E())*pDitopCM.Vect()
   if(beta.Mag()>=1): return 666.
   pTop1Ditop=pTop1CM
   pTop2Ditop=pTop2CM
   pProjDitop=pProjCM
   pTargDitop=pTargCM
   pTop1Ditop.Boost(beta)
   pTop2Ditop.Boost(beta)
   pProjDitop.Boost(beta)
   pTargDitop.Boost(beta)

   # Determine the z axis for the CS angle 
   zaxisCS=(((pProjDitop.Vect()).Unit())-((pTargDitop.Vect()).Unit())).Unit()
   yaxisCS=(((pProjDitop.Vect()).Unit()).Cross((pTargDitop.Vect()).Unit())).Unit()
   xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit()

   phi = -999.
   if(charge1>0.): phi = TMath.ATan2((pTop1Ditop.Vect()).Dot(yaxisCS),((pTop1Ditop.Vect()).Dot(xaxisCS)))
   else:           phi = TMath.ATan2((pTop2Ditop.Vect()).Dot(yaxisCS),((pTop2Ditop.Vect()).Dot(xaxisCS)))
   if(phi>TMath.Pi()): phi = phi-TMath.Pi()

   return phi

def PhiHE(p1, charge1, p2):
   # Phi decay angle (top (Q=+2/3)) in the Helicity frame
   pTop1Lab = TLorentzVector(0,0,-1,1)  # In the lab. frame 
   pTop2Lab = TLorentzVector(0,0,-1,1)  # In the lab. frame 
   pProjLab = TLorentzVector(0,0,-1,1)  # In the lab. frame 
   pTargLab = TLorentzVector(0,0,-1,1)  # In the lab. frame 
   pDitopLab = TLorentzVector(0,0,-1,1) # In the lab. frame 
   pTop1Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTop2Ditop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pProjDitop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   pTargDitop = TLorentzVector(0,0,-1,1) # In the Ditop rest frame
   beta = TVector3(0,0,0)
   xaxis = TVector3(0,0,0)
   yaxis = TVector3(0,0,0)
   zaxis = TVector3(0,0,0)
   mp = 0.93827231
   ep = 6500.

   # Get the muons parameters in the LAB frame
   pTop1Lab.SetPxPyPzE(p1.Px(),p1.Py(),p1.Pz(),p1.E())
   pTop2Lab.SetPxPyPzE(p2.Px(),p2.Py(),p2.Pz(),p2.E())
   
   # Obtain the Ditop parameters in the LAB frame
   pDitopLab=pTop1Lab+pTop2Lab
   zaxis=(pDitopLab.Vect()).Unit()
   
   # Translate the muon parameters in the Ditop rest frame
   beta=(-1./pDitopLab.E())*pDitopLab.Vect()
   if(beta.Mag()>=1.): return 666.
   
   pProjLab.SetPxPyPzE(0.,0.,-ep,TMath.Sqrt(ep*ep+mp*mp))
   pTargLab.SetPxPyPzE(0.,0.,+ep,TMath.Sqrt(ep*ep+mp*mp))
   
   pProjDitop=pProjLab
   pTargDitop=pTargLab
   
   pProjDitop.Boost(beta)
   pTargDitop.Boost(beta)
   
   yaxis=((pProjDitop.Vect()).Cross(pTargDitop.Vect())).Unit()
   xaxis=(yaxis.Cross(zaxis)).Unit()
   
   pTop1Ditop=pTop1Lab
   pTop2Ditop=pTop2Lab
   pTop1Ditop.Boost(beta)
   pTop2Ditop.Boost(beta)

   phi = -999.
   if(charge1>0.): phi = TMath.ATan2((pTop1Ditop.Vect()).Dot(yaxis),(pTop1Ditop.Vect()).Dot(xaxis))
   else:           phi = TMath.ATan2((pTop2Ditop.Vect()).Dot(yaxis),(pTop2Ditop.Vect()).Dot(xaxis))
   
   return phi