"""
The Subroutine of class is used to define the particle size distribution (PSD)
the type of PSD is included
- ExponentialPSD
- UnnormalizedGammaPSD
- GammaPSD
- BinnedPSD
"""
from .Scaling_Normalization import hxgg_Nmon_fc
from datetime import datetime
# try:
#     import cPickle as pickle
# except ImportError:
#     import pickle
import warnings
import numpy as np
# from scipy.integrate import trapz
from scipy.special import gamma, gammainc



class PSD(object):
    def __call__(self, D):
        if np.shape(D) == ():
            return 0.0
        else:
            return np.zeros_like(D)

    def __eq__(self, other):
        return False


class ExponentialPSD(PSD):
    """Exponential particle size distribution (PSD).
    
    Callable class to provide an exponential PSD with the given 
    parameters. The attributes can also be given as arguments to the 
    constructor.

    The PSD form is:
    N(D) = N0 * exp(-Lambda*D)

    Attributes:
        N0: the intercept parameter.
        Lambda: the inverse scale parameter        
        D_max: the maximum diameter to consider (defaults to 11/Lambda,
            i.e. approx. 3*D0, if None)

    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters larger than D_max.
    """

    def __init__(self, N0=np.array([1.0]), Lambda=np.array([1.0]), D_max=None):
        # check the arguments is a ndarray or not
        if isinstance(N0, np.ndarray) and isinstance(Lambda, np.ndarray):
            if N0.size != Lambda.size:
                raise ValueError("The number of N0 and Lambda must be the same")
            self.N0 = N0.astype(float)
            self.Lambda = Lambda.astype(float)
            self.num = len(N0)
            D0 = (3.67+0)/Lambda
            self.D0 = D0
            if D_max is None:
                self.D_max = 3.0*D0
            else:
                # if D_max is a scalar, make it an array
                if isinstance(D_max, (int, float)):
                    self.D_max = np.ones(self.num)*D_max
                elif isinstance(D_max, np.ndarray):
                    self.D_max = D_max
            if self.D_max.size != self.N0.size:
                raise ValueError("The number of D_max and N0 must be the same")
        else: 
            raise ValueError("The N0 and Lambda must be ndarray")
        


    def __call__(self, D):
        psd = np.zeros([self.num,D.size], dtype=float)
        for i in range(self.num):
            psd[i,:] = self.N0[i] * np.exp(-self.Lambda[i]*D)
            if np.shape(D) == ():
                if D > self.D_max[i]:
                    psd[i,:] = np.nan
            else:
                psd[i, D > self.D_max[i]] = np.nan
        return psd



class GammaPSD(PSD):
    """Gamma particle size distribution (PSD).
    
    Callable class to provide an gamma PSD with the given 
    parameters. The attributes can also be given as arguments to the 
    constructor.

    The PSD form is:
    N(D) = N0 * D**mu * exp(-Lambda*D)

    Attributes:
        N0: the intercept parameter.
        Lambda: the inverse scale parameter
        mu: the shape parameter
        D_max: the maximum diameter to consider (defaults to 11/Lambda,
            i.e. approx. 3*D0, if None)

    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters larger than D_max.
    """

    def __init__(self, N0=np.array([1.0]), Lambda=np.array([1.0]),mu=np.array([0.0]), D_max=None):
        # check the arguments is a ndarray or not
        if isinstance(N0, np.ndarray) and isinstance(Lambda, np.ndarray):
            if N0.size != Lambda.size:
                raise ValueError("The number of N0 and Lambda must be the same")
            self.N0 = N0.astype(float)
            self.Lambda = Lambda.astype(float)
            self.mu = mu.astype(float)
            self.num = len(N0)
            D0 = (3.67+mu)/Lambda
            self.D0 = D0
            if D_max is None:
                self.D_max = 3.0*D0
            else:
                # if D_max is a scalar, make it an array
                if isinstance(D_max, (int, float)):
                    self.D_max = np.ones(self.num)*D_max
                elif isinstance(D_max, np.ndarray):
                    self.D_max = D_max
            if self.D_max.size != self.N0.size:
                raise ValueError("The number of D_max and N0 must be the same")
        else: 
            raise ValueError("The N0 and Lambda must be ndarray")

    def __call__(self, D):
        psd = np.zeros([self.num,D.size], dtype=float)
        for i in range(self.num):
            psd[i,:] = self.N0[i] * D**self.mu[i] * np.exp(-self.Lambda[i]*D)
            if np.shape(D) == ():
                if D > self.D_max[i]:
                    psd[i,:] = np.nan
            else:
                psd[i, D > self.D_max[i]] = np.nan
        return psd
                
class NormalizedGammaPSD(PSD):
    """Normalized gamma particle size distribution (PSD).
    
    Callable class to provide a normalized gamma PSD with the given 
    parameters. The attributes can also be given as arguments to the 
    constructor.

    The PSD form is:
    N(D) = Nw * f_D0(mu) * (D/D0)**mu * exp(-(3.67+mu)*D/D0)
    or N(D) = Nw * f_Dm(mu) * (D/Dm)**mu * exp(-(4+mu)*D/Dm)
    f_D0(mu) = 6/(3.67**4) * (3.67+mu)**(mu+4)/Gamma(mu+4)
    f_Dm(mu) = 6/(4**4) * (4+mu)**(mu+4)/Gamma(mu+4)
    here we only consider the Dm 
    change the D0 to Dm from the relationship Dm = (4+mu)/(3.67+mu)*D0
    Attributes:
        D0: the median volume diameter.
        Nw: the intercept parameter.
        mu: the shape parameter.
        D_max: the maximum diameter to consider (defaults to 3*D0 when
            if None)

    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters larger than D_max.
    """       
    def __init__(self, Dm=None, D0=None, Nw=np.array([1.0]), mu=np.array([0.0]), D_max=None):
        # check the variable for input (not None)
        if Dm is None and D0 is None:
            raise ValueError("Either Dm or D0 must be given.")
        elif Dm is not None and D0 is not None:
            raise ValueError("Only one of Dm or D0 should be given.")
        elif Dm is not None and D0 is None:
            D0 = Dm*(3.67+mu)/(4+mu)
        elif Dm is None and D0 is not None:
            Dm = D0*(4+mu)/(3.67+mu)
        # check the arguments is a ndarray or not
        if isinstance(Dm, np.ndarray) and isinstance(mu, np.ndarray) and isinstance(Nw, np.ndarray):
            if Dm.size != mu.size or Dm.size != Nw.size or mu.size != Nw.size:
                raise ValueError("The number of Dm, mu and Nw must be the same")
            self.Dm = Dm.astype(float)
            self.mu = mu.astype(float)
            self.num = len(Dm)
            D0 = Dm*(3.67+mu)/(4+mu)
            self.D0 = D0
            self.Nw = Nw
            self.nf = 6/(4**4) * (4+mu)**(mu+4)/gamma(mu+4)
            if D_max is None:
                self.D_max = 3.0*D0
            if self.D_max.size != self.Dm.size:
                raise ValueError("The number of D_max and Dm must be the same")
        else: 
            raise ValueError("The input variable must be ndarray")
        
    def __call__(self, D):
        psd = np.zeros([self.num,D.size], dtype=float)
        for i in range(self.num):
            d = (D/self.Dm[i])
            psd[i,:] = self.Nw[i] * self.nf[i] *  d ** self.mu[i] * np.exp(-(4+self.mu[i])*d)
            if np.shape(D) == ():
                if D > self.D_max[i]:
                    psd[i,:] = np.nan
            else:
                psd[i, D > self.D_max[i]] = np.nan
        return psd

class BinnedPSD(PSD):
    """Binned gamma particle size distribution (PSD).
    
    Callable class to provide a binned PSD with the given bin edges and PSD
    values.

    Args (constructor):
        The first argument to the constructor should specify n+1 bin edges, 
        and the second should specify n bin_psd values.        
        
    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters outside the bins.
    """
    
    def __init__(self, bin_D, bin_dD, bin_psd):
        if len(bin_D) != len(bin_dD) or len(bin_D) != len(bin_psd) or len(bin_dD) != len(bin_psd):
            raise ValueError("There input variable must be the same length")
        self.D=bin_D
        self.dD=bin_dD
        self.psd=bin_psd
        num= bin_psd.shape[0]
        self.num = len(bin_D)
        #---------------------------------------------------------
        # there is some variable that need to be calculated
        Dmax = np.zeros(num)
        for i in range(num):
            ind=np.where(bin_psd[i,:]>0)[-1]
            Dmax[i]=bin_D[ind[-1]]
        self.D_max = Dmax

    # no __call_ now, because it is not needed
    # def __call__(self, D):
    #     psd = np.zeros([self.num,D.size], dtype=float)
    #     for i in range(self.num):
    #         d = (D/self.Dm[i])
    #         psd[i,:] = self.Nw[i] * self.nf[i] *  d ** self.mu[i] * np.exp(-(4+self.mu[i])*d)
    #         if np.shape(D) == ():
    #             if D > self.D_max[i]:
    #                 psd[i,:] = np.nan
    #         else:
    #             psd[i, D > self.D_max[i]] = np.nan
    #     return psd

class GeneralizedGammaPSD(PSD):
    """Generalized gamma particle size distribution (PSD).
    
    Callable class to provide a generalized gamma PSD with the given 
    parameters. The attributes can also be given as arguments to the 
    constructor.

    The PSD form is:
    N(D) = M0 * (c * Lambda )/ gamma(mu) * D ** (c*mu-1) * exp(-(Lambda*D)**c))

    Attributes:
        M0: the intercept parameter. also the 0th moment of the PSD
        Lambda: the slope parameter
        mu: the shape parameter
        c: the shape parameter    
        D_max: the maximum diameter to consider (defaults to 3*D0 when
            if None)

    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters larger than D_max.
    """ 
    def __init__(self, M0=np.array([1.0]), Lambda=np.array([1.0]), mu=np.array([0.0]), c=np.array([1.0]), D_max=None):
        # check the variable for input (not None)
        if isinstance(M0, np.ndarray) and isinstance(Lambda, np.ndarray) and isinstance(mu, np.ndarray) and isinstance(c, np.ndarray):
            if M0.size != Lambda.size or M0.size != mu.size or M0.size != c.size or Lambda.size != mu.size or Lambda.size != c.size or mu.size != c.size:
                raise ValueError("The number of M0, Lambda, mu and c must be the same")
            self.M0 = M0.astype(float)
            self.Lambda = Lambda.astype(float)
            self.mu = mu.astype(float)
            self.c = c.astype(float)
            self.num = len(M0)
            D0 = (3.67+mu)/Lambda
            if D_max is None:
                self.D_max = 3.0*D0
            if self.D_max.size != self.M0.size:
                raise ValueError("The number of D_max and M0 must be the same")
        else: 
            raise ValueError("The input variable must be ndarray")
            
    def __call__(self, D):
        psd = np.zeros([self.num,D.size], dtype=float)
        for i in range(self.num):
            psd[i,:] = self.M0[i] * (self.c[i] * self.Lambda[i])/ gamma(self.mu[i]) * D ** (self.c[i]*self.mu[i]-1) * np.exp(-(self.Lambda[i]*D)**self.c[i])    
            if np.shape(D) == ():
                if D > self.D_max[i]:
                    psd[i,:] = np.nan
            else:
                psd[i, D > self.D_max[i]] = np.nan
        return psd




class NormalizedGeneralizedGammaPSD(PSD):
    """Normalized Generalized gamma particle size distribution (PSD).
    
    Callable class to provide a Normalized Generalized Gamma DSD
    parameters. The attributes can also be given as arguments to the 
    constructor.

    The PSD form is:
    N(D) = Nc * h(D/Dc)
    Attributes:
        Nc: the characteristic number concentration
        Dc: the characteristic diameter
        h(D/Dc): the normalized generalized gamma function
    Nc and Dc are calculated from the N-moment of the PSD

    the normalized generalized gamma function is:
    h(x) = c*(gp1**-1)*(gcn**((p1+c*mu)/Cn))*x**(c*mu-1)*np.exp(-(gcn)**(c/Cn)*x**c)
    x is the normalized diameter D/Dc
    p1 is the P1 moment of the PSD
    Cn is the constant calculated from the array of moments (pnar)

    gp1 = gamma(1+p1/Cn)

    hx=c*(gp1**-1)*(gcn**((p1+c*mu)/Cn))*x**(c*mu-1)*np.exp(-(gcn)**(c/Cn)*x**c)

    Attributes:
        M0: the intercept parameter. also the 0th moment of the PSD
        Lambda: the slope parameter
        mu: the shape parameter
        c: the shape parameter    
        D_max: the maximum diameter to consider (defaults to 3*D0 when
            if None)

    Args (call):
        D: the particle diameter.

    Returns (call):
        The PSD value for the given diameter.    
        Returns 0 for all diameters larger than D_max.
    """
    def __init__(self, Dc=np.array([1.0]), Nc=np.array([1.0]), mu=np.array([0.0]), c=np.array([1.0]), pnar=np.array([3,6]), D_max=None):
        # check the variable for input (not None)
        if isinstance(Dc, np.ndarray) and isinstance(Nc, np.ndarray) and isinstance(mu, np.ndarray) and isinstance(c, np.ndarray):
            if Dc.size != Nc.size or Dc.size != mu.size or Dc.size != c.size or Nc.size != mu.size or Nc.size != c.size or mu.size != c.size:
                raise ValueError("The number of variables must be the same")
            self.Dc = Dc.astype(float)
            self.Nc = Nc.astype(float)
            self.mu = mu.astype(float)
            self.c = c.astype(float)
            self.num = len(Dc)
            self.pnar = pnar
            if D_max is None:
                self.D_max = 3.0*Dc # there is no D0, so we use Dc
            if self.D_max.size != self.Dc.size:
                raise ValueError("The number of D_max and M0 must be the same")
        else: 
            raise ValueError("The input variable must be ndarray")

    
    def __call__(self, D):
        psd = np.zeros([self.num,D.size], dtype=float)
        for i in range(self.num):
            x = D/self.Dc[i]
            hx = hxgg_Nmon_fc(x, self.mu[i], self.c[i], self.pnar)
            psd[i,:] = self.Nc[i] * hx
            if np.shape(D) == ():
                if D > self.D_max[i]:
                    psd[i,:] = np.nan
            else:
                psd[i, D > self.D_max[i]] = np.nan
        return psd
##=========================================================
#---------------------------------------------------------
def trunc_gamma(x,x_max):
    return gamma(x)*gammainc(x,x_max)
def gn(n,mu,c,Dmax):
    return trunc_gamma(mu+n/c,Dmax)
#---------------------------------------------------------
def Mn_from_DcNc(Dc,Nc,pn,npar,mu,c,Dmax):
    p1=npar[0]
    gp1=gn(p1,mu,c,Dmax)
    gpn=gn(pn,mu,c,Dmax)
    num_npar=npar.size
    odd_even=num_npar%2
    if odd_even==0: #even
        Cn=npar[::2].sum()-npar[1::2].sum()
        gcn=1
        for i in range(len(npar)):
            gcn=gcn*gn(npar[i],mu,c,Dmax)**((-1)**i)
    else: # odd
        Cn=npar[::2].sum()-npar[1::2].sum()-npar[-2]
        gcn=1
        for i in range(len(npar)):
            gcn=gcn*gn(npar[i],mu,c,Dmax)**((-1)**i)
        gcn=gcn/gn(npar[-2],mu,c,Dmax)
    Mn=gpn/gp1*(gcn)**((p1-pn)/Cn)*Nc*Dc**(pn+1)
    return Mn
##=========================================================
    #---------------------------------------------------------
    # function to get mu and c from the Nc and Dc
def CF_shape_from_Nc_Dc( DSD_shape, Dc, Nc, pnar,Dmax):
    mu,c=DSD_shape
    p1 = pnar[0]
    # ---------------------------------------------------------
    # 1. p1 is folloewed the self - consistency 
    Mp1= Nc*Dc**(p1+1)
    # 2. use guess mu and c to calculate other moments
    Mn = np.zeros_like(pnar, dtype=float)
    Mn[0] = Mp1
    for i in range(1,len(pnar)):
        Mn[i] =  Mn_from_DcNc(Dc,Nc,pnar[i],pnar,mu,c,Dmax)
    # 3. calculate Mp1 from lambda from other moments
    ld = np.zeros_like(pnar, dtype=float)
    for i in range(len(pnar)-1):
        p2=pnar[i]
        p3=pnar[i+1]
        # ld[i] = 
        ld[i]=((Mn[i+1]/Mn[i])*(gn(p2,mu,c,Dmax)/gn(p3,mu,c,Dmax)))**(1/(p2-p3))
    # 4. calculate Mp1 from lambda from other moments
    M0 = Mn_from_DcNc(Dc,Nc,0,pnar,mu,c,Dmax)
    Mn_est = np.zeros_like(pnar, dtype=float)
    for i in range(len(pnar)):
        Mn_est[i] = M0*ld[i]**(-p1)*trunc_gamma(mu+p1/c,Dmax)/gamma(mu)
    # ignore the first and last lambda calculation
    Mn_est[0] = Mp1
    Mn_est[-1] = Mp1
    # print(Mn_est)
    CF = np.sum((Mn_est-Mp1)**2)
    return CF

def get_shape_from_Nc_Dc(Dc,Nc,pnar,Dmax, mu_guess=1.0, c_guess=1.):
    """
    Get the mu and c from the Nc and Dc
    here DSD only can do in one case.
    """
    from scipy.optimize import minimize
    #---------------------------------------------------------
    # 1. guess the mu and c
    DSD_shape = (mu_guess, c_guess)
    #---------------------------------------------------------
    # 2. minimize the function to get mu and c
    res = minimize(CF_shape_from_Nc_Dc, DSD_shape, args=(Dc, Nc, pnar,Dmax), method='Nelder-Mead')
    mu, c = res.x

    return mu, c


        
