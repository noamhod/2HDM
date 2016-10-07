#!/usr/bin/python

import kinematics
import THDM

##################
### Begin code ###
##################

### make the definitions
mX     = THDM.model.mX
sba    = THDM.model.sba
nameX  = THDM.model.nameX
mgpath = THDM.model.mgpath
### make the libraries
THDM.testTHDM(mgpath)