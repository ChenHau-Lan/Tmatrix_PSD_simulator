# The subrountine to calculate the normalized DSD based on Scaling Normalization (Lee et al. 2004)
import numpy as np
from scipy.special import gamma, gammainc
#--------------------------------------------------------------------
# the list of subroutines:
# 1. 
#--------------------------------------------------------------------
def calc_moments_fc(ND,aD,dD,pn):
    # ND is array
    if len(ND.shape)==1:
        Mn=np.sum(ND*(aD**pn)*dD)
    else:
        Mn=np.sum(ND*(aD**pn)*dD,axis=-1)
    return Mn

def trunc_gamma(x,x_max=10,x_min=0):
    return gamma(x)*gammainc(x,x_max)
def gn(n,mu,c,Dmax):
    return trunc_gamma(mu+n/c,Dmax)