#!/usr/bin/python

import ROOT
from ROOT import gROOT, gDirectory, TFile, TTree, TMath, TEventList, TEntryList, TDirectory, TCut
import numpy as NUMPY
from enum import Enum
import subprocess # just to call an arbitrary command e.g. 'ls'
from contextlib import contextmanager
import os
import sys
from pipes import quote
from pprint import pprint
import imp


class t2HDM:
   """The 2HDM definitions"""
   def __init__(self, nameX="A", mX=500, type=2, sba=1, mintanb=0.3):
      self.nameX   = nameX
      self.mX      = mX
      self.type    = type
      self.sba     = sba
      self.mintanb = mintanb
      self.cuts  = "m"+nameX+"=="+str(mX)+" && sba=="+str(sba)+" && tanb>="+str(mintanb)
      self.cuts += " && TMath::ATan(tanb)>0. && TMath::ATan(tanb)<TMath::Pi()/2."
      self.cuts += " && TMath::Abs(cba)<=1."
      self.cuts += " && type=="+str(type)
      self.cuts += " && (status&3)==0" 
   def constrainWidth(wmin,wmax):
      self.cuts += " && (width_"+nameX+"/m"+nameX+">"+str(wmin)+" && width_"+nameX+"/m"+nameX+"<"+str(wmax)+")"

   mgpath = "/Users/hod/MC/MadGraph/MG5_aMC_v2_3_3_tests/"
   alphaS = 0.13 # For the ME^2 calculations
   nhel   = 0    # means sum over all helicity
##########################
### load the default model
model = t2HDM() ##########
##########################



class Fermions(Enum):
   u = 1
   d = 2
   s = 3
   c = 4
   b = 5
   t = 6
   e = 11
   mu = 13
   tau = 15

MC = 1.42
MB = 4.7
MT = 172.5
MM = 0.10566
MTAU = 1.777

### the data structure
parameters = []

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def invert_momenta(p):
   """ fortran/C-python do not order table in the same order"""
   new_p = []
   for i in range(len(p[0])):  new_p.append([0]*len(p))
   for i, onep in enumerate(p):
      for j, x in enumerate(onep):
         new_p[j][i] = x
   return new_p


def couplings(type,nameX,fermion,tanb,sba,cba):
   g = 0.
   if(type==1):
      if(nameX=="h"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = sba+cba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = sba+cba/tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = sba+cba/tanb
      if(nameX=="H"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = cba-sba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = cba-sba/tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = cba-sba/tanb
      if(nameX=="A"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = +1./tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = -1./tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = -1./tanb
   if(type==2):
      if(nameX=="h"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = sba+cba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = sba-cba*tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = sba-cba*tanb
      if(nameX=="H"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = cba-sba/tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = cba+sba*tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = cba+sba*tanb
      if(nameX=="A"):
         if(fermion==Fermions.u or fermion==Fermions.c  or fermion==Fermions.t):   g = +1./tanb
         if(fermion==Fermions.d or fermion==Fermions.s  or fermion==Fermions.b):   g = tanb
         if(fermion==Fermions.e or fermion==Fermions.mu or fermion==Fermions.tau): g = tanb
   return g

def setParameters(nameX,mX,cuts="",type=2,sba=1):
   f = TFile("thdm_grid_v163_13TeV.root","READ")
   t = f.Get("thdm")
   b_tb  = NUMPY.zeros(1, dtype=float)
   b_sba = NUMPY.zeros(1, dtype=float)
   b_cba = NUMPY.zeros(1, dtype=float)
   b_wA  = NUMPY.zeros(1, dtype=float)
   b_wH  = NUMPY.zeros(1, dtype=float)
   b_mA  = NUMPY.zeros(1, dtype=float)
   b_mH  = NUMPY.zeros(1, dtype=float)
   t.SetBranchAddress("tanb",b_tb);
   t.SetBranchAddress("sba",b_sba);
   t.SetBranchAddress("cba",b_cba);
   t.SetBranchAddress("width_A",b_wA);
   t.SetBranchAddress("width_H",b_wH);
   t.SetBranchAddress("mA",b_mA);
   t.SetBranchAddress("mH",b_mH);

   print "Before cuts = "+str(t.GetEntries())
   print "cuts: "+cuts
   t.Draw(">>elist",cuts,"entrylist");
   elist = gDirectory.Get("elist");
   t.SetEntryList(elist);
   Nlist = elist.GetN()
   print "After cuts = "+str(Nlist)

   global parameters
   for i in range(0,Nlist):
      n = elist.Next()
      t.GetEntry(n)
      tanb = b_tb[0]
      sba  = b_sba[0]
      cba  = b_cba[0]
      wA   = b_wA[0]
      wH   = b_wH[0]
      mA   = b_mA[0]
      mH   = b_mH[0]
      YMT   = couplings(type,nameX,Fermions.t,tanb,sba,cba)*MT
      YMB   = couplings(type,nameX,Fermions.b,tanb,sba,cba)*MB
      YMC   = couplings(type,nameX,Fermions.c,tanb,sba,cba)*MC
      YMM   = couplings(type,nameX,Fermions.mu,tanb,sba,cba)*MM
      YMTAU = couplings(type,nameX,Fermions.tau,tanb,sba,cba)*MTAU
      print "["+str(i)+"] tanb="+'%.6f' % tanb+" sba="+'%.6f' % sba+" cba="+'%.6f' % cba+" wA="+'%.6f' % wA+" wH="+'%.6f' % wH+" YMT="+'%.6f' % YMT+" YMB="+'%.6f' % YMB+" YMC="+'%.6f' % YMC+" YMM="+'%.6f' % YMM+" YMTAU="+'%.6f' % YMTAU
      adict = {}
      adict = {'tanb':tanb, 'sba':sba, 'cba':cba, 'mA':mA, 'wA':wA, 'mH':mH, 'wH':wH, 'YMT':YMT, 'YMB':YMB, 'YMC':YMC, 'YMM':YMM, 'YMTAU':YMTAU}
      parameters.append(adict.copy())
   print "N="+str(len(parameters))+" parameters set !"


def compileSM(mgpath,nameX,mX):	
   X       = "matrix/"+nameX+"/"+str(mX)+"/"
   command = "./bin/mg5_aMC noam/proc_card_mg5_SM_partonlevel.minimal.dat"
   procdir = "ggtt-SM-partonlevel/SubProcesses/P1_gg_ttx_no_hh1/"
   matxdir = "/Users/hod/MC/Pythia/pythia8215/examples/"+X
   libdir  = matxdir+"SM/"

   global parameters

   p = subprocess.Popen("rm -rf "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("mkdir -p "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()

   # enter the directory like this:
   with cd(mgpath):
      # make sure the old directory is removed
      p = subprocess.Popen("rm -rf "+procdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      
      # make sure the original SM card file is used
      forig = "models/Higgs_Effective_Couplings_FormFact/parameters.py_ORIG"
      fused = "models/Higgs_Effective_Couplings_FormFact/parameters.py"
      p = subprocess.Popen("/bin/cp -f "+forig+" "+fused, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      # execute the generation of the process
      p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()

      # go to make the library
      with cd(procdir):
         p = subprocess.Popen("make", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()

         ### cahnge the makefile, make and copy
         p = subprocess.Popen('sed -i -e "s/MENUM)py/MENUM)SMpy/g" '+mgpath+procdir+'makefile', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         p = subprocess.Popen("make matrix2SMpy.so", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         p = subprocess.Popen("cp matrix2SMpy.so "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         p = subprocess.Popen("cp ../../Cards/*.dat "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()


def compileX(index,mgpath,nameX,mX):
   X       = "matrix/"+nameX+"/"+str(mX)+"/"
   command = ""
   procdir = ""
   if(nameX=="A"):
      command = "./bin/mg5_aMC noam/proc_card_mg5_A_partonlevel.minimal.dat"
      procdir = "ggtt-A-partonlevel/SubProcesses/P1_gg_ttx_no_h/"
   if(nameX=="H"):
      command = "./bin/mg5_aMC noam/proc_card_mg5_H_partonlevel.minimal.dat"
      procdir = "ggtt-H-partonlevel/SubProcesses/P1_gg_ttx_no_h1/"
   if(nameX=="AH"):
      command = "./bin/mg5_aMC noam/proc_card_mg5_AH_partonlevel.minimal.dat"
      procdir = "ggtt-AH-partonlevel/SubProcesses/P1_gg_ttx/"
   
   thisdir = "/Users/hod/MC/Pythia/pythia8215/examples/"
   matxdir = thisdir+X
   libdir  = matxdir+str(index)+"/"

   global parameters

   if(index==0):
      p = subprocess.Popen("rm -rf "+matxdir+"*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      # print err

   p = subprocess.Popen("mkdir -p "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   # print err

   # enter the directory like this:
   with cd(mgpath):
      # make sure the old directory is removed
      p = subprocess.Popen("rm -rf "+procdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      # print err

      # modify the main parameters card
      ifname = "models/Higgs_Effective_Couplings_FormFact/parameters.py_ORIG"
      ofname = ifname.replace("_ORIG","")
      p = subprocess.Popen("/bin/cp -f "+ifname+" "+ofname, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
      out, err = p.communicate()
      # print err

      replacements = {}
      replacements.update({' = 111.,':          ' = '+str(parameters[index].get("mA"))+','})
      replacements.update({' = 0.11111,':       ' = '+str(parameters[index].get("wA"))+','})
      replacements.update({' = 120.,':          ' = '+str(parameters[index].get("mH"))+','})
      replacements.update({' = 0.00575308848,': ' = '+str(parameters[index].get("wH"))+','})
      replacements.update({' = 1.27,':          ' = '+str(parameters[index].get("YMC"))+','})
      replacements.update({' = 4.2,':           ' = '+str(parameters[index].get("YMB"))+','})
      replacements.update({' = 164.5,':         ' = '+str(parameters[index].get("YMT"))+','})
      replacements.update({' = 0.10567,':       ' = '+str(parameters[index].get("YMM"))+','})
      replacements.update({' = 1.778,':         ' = '+str(parameters[index].get("YMTAU"))+','})
      replacements.update({' = 172.,':          ' = 172.5,'})
      for sold in replacements.keys():
         snew = replacements[sold]
         p = subprocess.Popen('sed -i -- "s/'+sold+'/'+snew+'/g" '+ofname, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)     

      # execute the generation of the process
      p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      # print err

      # go to make the library
      with cd(procdir):
         ### cahnge the makefile, make and copy
         p = subprocess.Popen('sed -i -e "s/MENUM)py/MENUM)'+nameX+str(index)+'py/g" '+mgpath+procdir+'makefile', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         p = subprocess.Popen("make matrix2"+nameX+str(index)+"py.so", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         p = subprocess.Popen("cp matrix2"+nameX+str(index)+"py.so "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()
         p = subprocess.Popen("cp ../../Cards/*.dat "+libdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         out, err = p.communicate()


modules = {}
def setModules(basepath,libmatrix,nameX,nX,libs="All",index=-1):
   if(libs=="All" or libs=="AllX"):
      for i in range (0,nX):
         name = 'matrix2'+nameX+str(i)+'py'
         sindex = str(i)
         print "changing dir to: "+libmatrix+sindex+"/"
         with cd(libmatrix+sindex+"/"):
            print "in "+os.getcwd()+", trying to import ",name
            module_info = imp.find_module(name,[basepath+"/"+libmatrix+str(i)+"/"])
            modules.update({name:imp.load_module(name, *module_info)})
            modules[name].initialise("param_card.dat")
            print "Successfully initialised ",name
   if(libs=="X" and index!=-1):
      name = 'matrix2'+nameX+str(index)+'py'
      sindex = str(index)
      with cd(libmatrix+sindex+"/"):
         print "in "+os.getcwd()+", trying to import ",name
         module_info = imp.find_module(name,[basepath+"/"+libmatrix+str(index)+"/"])
         modules.update({name:imp.load_module(name, *module_info)})
         modules[name].initialise("param_card.dat")
         print "Successfully initialised ",name
   if(libs=="All" or libs=="SM"):
      name = 'matrix2SMpy'
      with cd(libmatrix+"SM/"):
         print "in "+os.getcwd()+", trying to import "+name
         module_info = imp.find_module(name,[basepath+"/"+libmatrix+"SM/"])
         modules.update({name:imp.load_module(name, *module_info)})
         modules[name].initialise("param_card.dat")
         print "Successfully initialised ",name


def testImport(nameX,mX,index=-1):
   libmatrix = "matrix/"+nameX+"/"+str(mX)+"/"
   if(index<0): setModules(os.getcwd(),libmatrix,nameX,len(parameters),"SM",index)
   else:        setModules(os.getcwd(),libmatrix,nameX,len(parameters),"X",index)
   alphaS = 0.13
   nhel = 0 # means sum over all helicity
   p = [[   0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03],
        [   0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03],
        [   0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03],
        [   0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03]]
   P=invert_momenta(p)
   ## the ME^2 and the weight
   me2 = -1
   if(index<0): me2 = modules['matrix2SMpy'].get_me(P,alphaS,nhel)                   ### calculate the SM ME^2
   else:        me2 = modules['matrix2'+nameX+str(index)+'py'].get_me(P,alphaS,nhel) ### calculate the X ME^2
   return me2


def makeSM(mgpath,nameX,mX,test=False):
   compileSM(mgpath,nameX,mX)
   if(test):
      me2 = testImport(nameX,mX)
      print "Done making library SM -> test ME^2="+str(me2)
   else:
      print "Done making library SM"


def make2HDM(mgpath,nameX,mX,test=False):
   for i in range(0,len(parameters)):
      compileX(i,mgpath,nameX,mX)
      if(test):
        me2 = testImport(nameX,mX,i)
        print "Done making library "+str(i)+" -> test ME^2="+str(me2)
      else:
        print "Done making library "+str(i)


def makeAll(mgpath,nameX,mX,test=False):	
   make2HDM(mgpath,nameX,mX,test)
   makeSM(mgpath,nameX,mX,test)


def testTHDM(mgpath,nameX,mX):
   type=model.type
   sba=model.sba
   cuts = model.cuts
   setParameters(nameX,mX,cuts,type,sba)
   print parameters
   makeAll(mgpath,nameX,mX,True)