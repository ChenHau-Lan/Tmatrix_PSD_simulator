# this subroutines are used to calculate the radar parameters
# for the given T-matrix

import numpy as np

# load the scatter.py module --> use scatter object
from .scatter import Scatterer as scat
from scipy import integrate
# # define the function to calculate the radar parameters
# def calc_dual(scat,ND,aD,dD):
# 	"""
# 	Here is to calculate the dual-polarization radar parameters
# 	ZHH: horizontal reflectivity factor (in dB)
# 	ZVV: vertical reflectivity factor (in dB)
# 	ZDR: differential reflectivity (in dB)
# 	LDR: linear depolarization ratio (in dB)
# 	KDP: specific differential phase (in deg/km)
# 	RHOHV: copolar correlation coefficient (unitless)
# 	ATTH: horizontal attenuation (in dB/km)
# 	ATTV: vertical attenuation (in dB/km)
# 	"""
# 	#ZHH
# 	Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
# 	Zh =integrate.simpson(NFAC*Fhh*ND,x=aD,dx=dD)
# 	ZHH=10*np.log10(Zh)
# 	#ZVV
# 	Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
# 	Zv =NFAC*integrate.simpson(Fvv*ND,x=aD,dx=dD)
# 	#ZDR
# 	ZDR=10*np.log10(Zh/Zv)
#     #KDP
# 	RADDEG = 180./PI
# 	KDP=integrate.simpson(RADDEG*K[:,2,3]*0.1*ND,x=aD,dx=dD)
# 	#ZDP
# 	return ZHH,ZDR,KDP  


def calc_Zh(scat,ND,aD,dD):
	#ZHH: horizontal reflectivity factor (in dB)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
	# Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	Zh =integrate.simpson(NFAC*Fhh*ND,x=aD,dx=dD)
	return Zh

def calc_ZH(scat,ND,aD,dD):
	#ZHH: horizontal reflectivity factor (in dB)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
	# Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	Zh =integrate.simpson(NFAC*Fhh*ND,x=aD,dx=dD)
	ZHH=10*np.log10(Zh)
	return ZHH

def calc_Zv(scat,ND,aD,dD):
	#ZVV: vertical reflectivity factor (in dB)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	# Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
	Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	Zv =integrate.simpson(NFAC*Fvv*ND,x=aD,dx=dD)
	return Zv

def calc_ZV(scat,ND,aD,dD):
	#ZVV: vertical reflectivity factor (in dB)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	# Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
	Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	Zv =integrate.simpson(NFAC*Fvv*ND,x=aD,dx=dD)
	ZVV=10*np.log10(Zv)
	return ZVV

def calc_Zdr(scat,ND,aD,dD):
	#ZDR: differential reflectivity (in dB)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
	Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	Zh =integrate.simpson(Fhh*ND,x=aD,dx=dD)
	Zv =integrate.simpson(Fvv*ND,x=aD,dx=dD)
	return Zh/Zv

# def calc_ZDP(scat,ND,aD,dD):


def calc_ZDR(scat,ND,aD,dD):
	#ZDR: differential reflectivity (in dB)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])	
	Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	Zh =integrate.simpson(Fhh*ND,x=aD,dx=dD)
	Zv =integrate.simpson(Fvv*ND,x=aD,dx=dD)
	ZDR=10*np.log10(Zh/Zv)
	return ZDR

def calc_KDP(scat,ND,aD,dD):
	#KDP: specific differential phase (in deg/km)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	RADDEG = 180./PI
	KDP=integrate.simpson(RADDEG*K[:,2,3]*0.1*ND,x=aD,dx=dD)
	return KDP

def calc_RHOHV(scat,ND,aD,dD):
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	Fhh=(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])
	MAG1=integrate.simpson(Fhh*ND,x=aD,dx=dD)
	Fvv=(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])
	MAG2=integrate.simpson(Fvv*ND,x=aD,dx=dD)
	DENOM = (MAG1)**(1/2)*(MAG2)**(1/2) 
	RHORE = integrate.simpson((S[:,2,2] + S[:,3,3])*ND,x=aD,dx=dD)
	RHORE = -RHORE/DENOM
	RHOIM = integrate.simpson((S[:,2,3] - S[:,3,2])*ND,x=aD,dx=dD)
	RHOIM = -RHOIM/DENOM
	RHOHV = (RHORE**2 + RHOIM**2)**(1/2)
	return RHOHV

def calc_ATTH(scat,ND,aD,dD):
	#ATTH: horizontal attenuation (in dB/km)
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	ATTH = 0.4343*integrate.simpson((K[:,0,0]-K[:,0,1]) *ND,x=aD,dx=dD)
	return ATTH

def calc_DATT(scat,ND,aD,dD):
	#DATT: DIFFERENTIAL ATTENUATION DB/KM
	S=scat.S
	K=scat.K
	wl=scat.wavelength
	MAGKSQ=0.92
	PI=np.pi
	NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
	NFAC = NFAC*4.*PI
	DATT = -0.4343*2.*integrate.simpson((K[:,0,1]) *ND,x=aD,dx=dD)
	return DATT
