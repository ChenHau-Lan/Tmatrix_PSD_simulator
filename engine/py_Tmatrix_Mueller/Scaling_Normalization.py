# The subrountine to calculate the normalized DSD based on Scaling Normalization (Lee et al. 2004)
import numpy as np
from scipy.special import gamma
#--------------------------------------------------------------------
# the list of subroutines:
# 1. Dc_NC_2mom_fc
#--------------------------------------------------------------------
def calc_moments_fc(ND,aD,dD,pn):
    # ND is array
    if len(ND.shape)==1:
        Mn=np.sum(ND*(aD**pn)*dD)
    else:
        Mn=np.sum(ND*(aD**pn)*dD,axis=-1)
    return Mn
######################################################################
# 1. hxgg_fc
# the generalized gamma distribution of Normalized DSD h(x) 
def hxgg_fc(x,mu,c,mi,mj):
    # The input
    #  x  : the normalized DSD diameter: D/Dm
    #  mu : the shape parameter of generalized gamma distribution
    #  c  : the scale parameter of generalized gamma distribution
    #  mi : the smaller order of moment
    #  mj : the larger order of moment
    x[x==0]=np.nan
    gi=gamma(mu+mi/c)
    gj=gamma(mu+mj/c)
    hx=c*gi**((mj+c*mu)/(mi-mj))*gj**((-mi-c*mu)/(mi-mj))*x**(c*mu-1)*np.exp(-(gi/gj)**(c/(mi-mj))*x**c)
    return hx
######################################################################
# 2. normalized_ND_2mom_fc
# the normalized DSD based on 2 moments Scaling Normalization (Lee et al. 2004)
def normalized_ND_2mom_fc(ND,mi,mj,aD,dD):
    # The input
    #  ND : DSD bins
    #  mi : the smaller order of moment
    #  mj : the larger order of moment
    #  aD : DSD bin size
    #  dD : the interval of DSD bins
    # The output
    #  hx : the normalized DSD
    #  DDm: the normalized diameter
    Mi=DSD_micro.calc_moments_fc(ND,aD,dD,mi)
    Mj=DSD_micro.calc_moments_fc(ND,aD,dD,mj)  
    N0=Mi**((mj+1)/(mj-mi))*Mj**((mi+1)/(mi-mj))
    Dm=(Mj/Mi)**(1/(mj-mi))
    # calculate normalized ND
    # DDm is the normalized diameter DDm[n*m]=aD[n]/Dm[m]
    DDm=(np.asmatrix(1/Dm).T.dot(np.asmatrix(aD)))
    hx=(np.asmatrix(ND)/np.asmatrix(N0).T)
    return hx,DDm
######################################################################
# 3 ND_from_hxgg_fc
# calculate the normal DSD [N(D)] from moments Mi and Mj , and the normalized DSD h(x) from GG model
def ND_from_hxgg_fc(Mi,Mj,aD,mi,mj,mu_assumed,c_assumed):
    N0=Mi**((mj+1)/(mj-mi))*Mj**((mi+1)/(mi-mj))
    Dm=(Mj/Mi)**(1/(mj-mi))
    dx=aD/Dm
    hx_est=hxgg_fc(dx,mu_assumed,c_assumed,mi,mj)
    ND=hx_est*N0
    return ND
######################################################################
######################################################################
######################################################################
######################################################################
# 
def calc_DcNc_EVEN_mon_fc(ND,aD,dD,npar):
    # The input
    #  ND   : DSD bins
    #  aD   : DSD bin size
    #  dD   : the interval of DSD bins
    #  npar : the list of moments, first element is P1.
    # The output
    #  Dc : the mass-weighted mean diameter
    #  Nc : the normalized intercept parameter
    #------------------------------------------------------------------
    # The moments
    if (len(ND.shape)==1):
        Mn=np.zeros((1,len(npar)))
    else:
        Mn=np.zeros((ND.shape[0],len(npar)))
    for i in range(len(npar)):
        Mn[:,i]=calc_moments_fc(ND,aD,dD,npar[i])
    # calculate the Dc and Nc
    pn=npar[::2].sum()-npar[1::2].sum()
    p1=npar[0]
    Nc=Mn[:,0]**(1-(p1+1)/pn)
    Dc=Mn[:,0]**(1/pn)
    for i in range(1,len(npar)):
        Nc=Nc*Mn[:,i]**((-1)**(i-1)*(p1+1)/pn)
        Dc=Dc*Mn[:,i]**((-1)**(i)*1/pn)
    return Dc,Nc

def calc_DcNc_ODD_mon_fc(ND,aD,dD,npar):
    # The input
    #  ND   : DSD bins
    #  aD   : DSD bin size
    #  dD   : the interval of DSD bins
    #  npar : the list of moments, first element is P1.
    # The output
    #  Dc : the mass-weighted mean diameter
    #  Nc : the normalized intercept parameter
    #------------------------------------------------------------------
    # The moments
    if (len(ND.shape)==1):
        Mn=np.zeros((1,len(npar)))
    else:
        Mn=np.zeros((ND.shape[0],len(npar)))
    for i in range(len(npar)):
        Mn[:,i]=calc_moments_fc(ND,aD,dD,npar[i])
    # calculate the Dc and Nc
    pn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
    p1=npar[0]
    Nc=Mn[:,0]**(1-(p1+1)/pn)
    Dc=Mn[:,0]**(1/pn)
    for i in range(1,len(npar)):
        Nc=Nc*Mn[:,i]**((-1)**(i-1)*(p1+1)/pn)
        Dc=Dc*Mn[:,i]**((-1)**(i)*1/pn)
    Nc=Nc*Mn[:,-2]**((p1+1)/pn)
    Dc=Dc*Mn[:,-2]**(-1/pn)
    return Dc,Nc

def calc_DcNc_EVEN_Mn_fc(Mn,npar):
    # The input
    #------------------------------------------------------------------
    # The moments already calculated
    # calculate the Dc and Nc
    pn=npar[::2].sum()-npar[1::2].sum()
    p1=npar[0]
    Nc=Mn[:,0]**(1-(p1+1)/pn)
    Dc=Mn[:,0]**(1/pn)
    for i in range(1,len(npar)):
        Nc=Nc*Mn[:,i]**((-1)**(i-1)*(p1+1)/pn)
        Dc=Dc*Mn[:,i]**((-1)**(i)*1/pn)
    return Dc,Nc
def calc_DcNc_ODD_Mn_fc(Mn,npar):
    #------------------------------------------------------------------
    # The moments already calculated
    # calculate the Dc and Nc
    pn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
    p1=npar[0]
    Nc=Mn[:,0]**(1-(p1+1)/pn)
    Dc=Mn[:,0]**(1/pn)
    for i in range(1,len(npar)):
        Nc=Nc*Mn[:,i]**((-1)**(i-1)*(p1+1)/pn)
        Dc=Dc*Mn[:,i]**((-1)**(i)*1/pn)
    Nc=Nc*Mn[:,-2]**((p1+1)/pn)
    Dc=Dc*Mn[:,-2]**(-1/pn)
    return Dc,Nc
def calc_DcNc_from_N_Mn_fc(Mn,npar):
    #------------------------------------------------------------------
    # The moments
    num_npar=npar.size
    # get num_npar is even or odd
    odd_even=num_npar%2
    if num_npar<2:
        print('The number of moments should be larger than 1')
        Dc,Nc=np.nan,np.nan
    elif odd_even==0:
        # print('The number of moments should be even')
        pn=npar[::2].sum()-npar[1::2].sum()
        if pn==0:
            print('This combination of moments is not valid (=0)')
            Dc,Nc=np.nan,np.nan
        else:
            Dc,Nc=calc_DcNc_EVEN_Mn_fc(Mn,npar)
    else: # odd
        # print('The number of moments should be odd')
        pn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
        if pn==0:
            print('This combination of moments is not valid (=0)')
            Dc,Nc=np.nan,np.nan
        else:
            Dc,Nc=calc_DcNc_ODD_Mn_fc(Mn,npar)
    return Dc,Nc

######################################################################
# 10. calc_DcNc_Nmom_fc
# The subroutine to calculate the Dc and Nc from DSD using N moment scaling normalization
def calc_DcNc_Nmon_fc(ND,aD,dD,npar):
    # The input
    #  ND   : DSD bins
    #  aD   : DSD bin size
    #  dD   : the interval of DSD bins
    #  npar : the list of moments, first element is P1.
    # The output
    #  Dc : the mass-weighted mean diameter
    #  Nc : the normalized intercept parameter
    #------------------------------------------------------------------
    # The moments
    num_npar=npar.size
    # get num_npar is even or odd
    odd_even=num_npar%2
    if num_npar<2:
        print('The number of moments should be larger than 1')
        Dc,Nc=np.nan,np.nan
    elif odd_even==0:
        # print('The number of moments should be even')
        pn=npar[::2].sum()-npar[1::2].sum()
        if pn==0:
            print('This combination of moments is not valid (=0)')
            Dc,Nc=np.nan,np.nan
        else:
            Dc,Nc=calc_DcNc_EVEN_mon_fc(ND,aD,dD,npar)
    else: # odd
        # print('The number of moments should be odd')
        pn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
        if pn==0:
            print('This combination of moments is not valid (=0)')
            Dc,Nc=np.nan,np.nan
        else:
            Dc,Nc=calc_DcNc_ODD_mon_fc(ND,aD,dD,npar)
    return Dc,Nc
#----------------------------------------------------------
def calc_normalized_ND_Nmom_fc(ND,aD,dD,npar):
    Dc,Nc=calc_DcNc_Nmon_fc(ND,aD,dD,npar)
    aDc=(np.asmatrix(1/Dc).T.dot(np.asmatrix(aD)))
    dDc=(np.asmatrix(1/Dc).T.dot(np.asmatrix(dD)))
    hx=(np.asmatrix(ND)/np.asmatrix(Nc).T)
    hx=np.array(hx)
    aDc=np.array(aDc)
    dDc=np.array(dDc)
    return hx,aDc,dDc
def calc_normalized_ND_from_DcNc_fc(ND,aD,dD,Dc,Nc):
    aDc=(np.asmatrix(1/Dc).T.dot(np.asmatrix(aD)))
    dDc=(np.asmatrix(1/Dc).T.dot(np.asmatrix(dD)))
    hx=(np.asmatrix(ND)/np.asmatrix(Nc).T)
    hx=np.array(hx)
    aDc=np.array(aDc)
    dDc=np.array(dDc)
    return hx,aDc,dDc
#====================================================================
def hxgg_Nmon_EVEN_fc(x,mu,c,npar):
    Cn=npar[::2].sum()-npar[1::2].sum()
    p1=npar[0]
    gp1=gamma(mu+p1/c)
    gcn=1
    for i in range(len(npar)):
        gcn=gcn*gamma(mu+npar[i]/c)**((-1)**i)
    hx=c*(gp1**-1)*(gcn**((p1+c*mu)/Cn))*x**(c*mu-1)*np.exp(-(gcn)**(c/Cn)*x**c)
    return hx
def hxgg_Nmon_ODD_fc(x,mu,c,npar):
    Cn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
    p1=npar[0]
    gp1=gamma(mu+p1/c)
    gcn=1
    for i in range(len(npar)):
        gcn=gcn*gamma(mu+npar[i]/c)**((-1)**i)
    gcn=gcn/gamma(mu+npar[-2]/c)
    hx=c*(gp1**-1)*(gcn**((p1+c*mu)/Cn))*x**(c*mu-1)*np.exp(-(gcn)**(c/Cn)*x**c)
    return hx
#----------------------------------------------------------
def hxgg_Nmon_fc(x,mu,c,npar):
    x[x==0]=np.nan
    # The moments
    num_npar=npar.size
    # get num_npar is even or odd
    odd_even=num_npar%2
    hx=np.ones_like(x)*np.nan
    if num_npar<2:
        print('The number of moments should be larger than 1')
    elif odd_even==0: 
        # print('The number of moments should be even')
        Cn=npar[::2].sum()-npar[1::2].sum()
        if Cn==0:
            print('This combination of moments is not valid (=0)')
        else:
            hx=hxgg_Nmon_EVEN_fc(x,mu,c,npar)
    else: # odd
        # print('The number of moments should be odd')
        Cn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
        if Cn==0:
            print('This combination of moments is not valid (=0)')
        else:
            hx=hxgg_Nmon_ODD_fc(x,mu,c,npar)
    return hx

def calc_Mn_from_NcDc_fc(Dc,Nc,pn):
    # The input
    #  Dc : the mass-weighted mean diameter
    #  Nc : the normalized intercept parameter
    #------------------------------------------------------------------
    # calculate the Mi and Mj
    Mpn=Nc*Dc**(pn+1)
    return Mpn
#====================================================================
# Mp_from_NcDc_fc
def calc_Mz_from_Cn_NcDc_fc(Dc,Nc,Cn,z):
    # The input
    #  Dc : the mass-weighted mean diameter
    #  Nc : the normalized intercept parameter
    #------------------------------------------------------------------
    # calculate the Mz
    Mz=Cn*Nc*Dc**(z+1)
    return Mz

def calc_Mz_from_NcDc_hxgg_fc(Dc,Nc,mu,c,npar,z,aD,dD):
    ND=calc_ND_from_DcNc_hxgg_fc(Dc,Nc,mu,c,npar,aD)
    Mz=calc_moments_fc(ND,aD,dD,z)
    return Mz

#====================================================================
# calc_ND_from_DcNc_hxgg_fc
def calc_ND_from_DcNc_hxgg_fc(Dc,Nc,mu,c,npar,aD):
    # The input
    #  Dc : the mass-weighted mean diameter
    #  Nc : the normalized intercept parameter
    #  mu : the shape parameter of generalized gamma distribution
    #  c  : the scale parameter of generalized gamma distribution
    #  npar : the list of moments, first element is P1.
    # The output
    #  ND : the normalized DSD
    #------------------------------------------------------------------
    # calculate the Mi and Mj
    x=aD/Dc
    hx=hxgg_Nmon_fc(x,mu,c,npar)
    ND=hx*Nc
    return ND

#====================================================================
def fitting_GG_hx(hxar,DDar,mi,mj):
    import scipy.optimize as optimize
    from scipy.special import gamma
    def model(p,x,mi,mj):
        mu,c =p     
        gi=gamma(mu+mi/c)
        gj=gamma(mu+mj/c)  
        if gi==np.inf:
            gi=gamma(0.001)
        if gj==np.inf:
            gj=gamma(0.001)        
        return c*gi**((mj+c*mu)/(mi-mj))*gj**((-mi-c*mu)/(mi-mj))*x**(c*mu-1)*np.exp(-(gi/gj)**(c/(mi-mj))*x**c)
    def fun(p,x,mi,mj,y):
        return np.log10(model(p,x,mi,mj))-np.log10(y)

    ijk=np.where(hxar>0)[0]
    hxar=hxar[ijk]
    DDar=DDar[ijk]  
    # fitting
    x0=np.array([2,1])
    # print(DDar)
    # print(hxar)
    fit_res=optimize.least_squares(fun,x0,args=(DDar,mi,mj,hxar),bounds=((-5,0),(np.inf,np.inf)),jac='cs',verbose=0)
    mu=fit_res.x[0]
    c=fit_res.x[1]   
    return mu,c

#====================================================================
def Cn_fc(npar):
    npar=np.array(npar)
    num_npar=npar.size
    odd_even=num_npar%2
    if odd_even==0:
        Cn=npar[::2].sum()-npar[1::2].sum()
    else:
        Cn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
    return Cn
def gcn_fc(npar,mu,c):
    from scipy.special import gamma
    npar=np.array(npar)
    num_npar=npar.size
    odd_even=num_npar%2
    Cn=Cn_fc(npar)
    if odd_even==0:
        gcn=1
        for i in range(len(npar)):
            gcn=gcn*gamma(mu+npar[i]/c)**((-1)**i)
    else:
        gcn=1
        for i in range(len(npar)):
            gcn=gcn*gamma(mu+npar[i]/c)**((-1)**i)
        gcn=gcn/gamma(mu+npar[-2]/c)
    return gcn

def alpha_pz_fc(pnar,pz,mu,c):
    from scipy.special import gamma
    p1=pnar[0]
    Cn=Cn_fc(pnar)
    gcn=gcn_fc(pnar,mu,c)
    gp1=gamma(mu+p1/c)
    gpz=gamma(mu+pz/c)
    alpha_pn=gpz/gp1*(gcn)**((p1-pz)/Cn)
    return alpha_pn