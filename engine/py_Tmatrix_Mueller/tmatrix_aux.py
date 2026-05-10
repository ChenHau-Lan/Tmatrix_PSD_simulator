#current version
VERSION = "0.1.1"

#typical wavelengths [mm] at different bands
wl_S = 10.7
wl_C = 5.35
wl_X = 3.33
wl_Ku = 2.20
wl_Ka = 0.843
wl_W = 0.319

#typical values of K_w_sqr at different bands
K_w_sqr = {wl_S: 0.93, wl_C: 0.93, wl_X: 0.93, wl_Ku: 0.93, wl_Ka: 0.92, 
  wl_W: 0.75}

#Drop Shape Relationship Functions
def dsr_thurai_2007(D_eq):
	"""
	Drop shape relationship function from Thurai2007
    (http://dx.doi.org/10.1175/JTECH2051.1) paper.
    Arguments:
        D_eq: Drop volume-equivalent diameter (mm)
    Returns:
        r: The vertical-to-horizontal drop axis ratio. Note: the Scatterer class
        expects horizontal to vertical, so you should pass 1/dsr_thurai_2007
    """
	AVORB = 1.065 - 6.25e-2*D_eq - 3.99e-3*D_eq**2 + 7.66e-4*D_eq**3 - \
			4.095e-5*D_eq**4 # D_eq >= 1.5
	AVORB[D_eq < 1.5] = 1.173 - 0.5165*D_eq[D_eq < 1.5] + 0.4698*D_eq[D_eq < 1.5]**2 - 0.1317*D_eq[D_eq < 1.5]**3 - \
			8.5e-3*D_eq[D_eq < 1.5]**4
	AVORB[D_eq < 0.7] = 1.0
	return AVORB
def dsr_green_1975(D_eq):
	"""
	Drop shape relationship function from Green1975 [eq can see on Zhang et al. 2001]
    (DOI: 10.1109/36.917906) paper.
    Arguments:
        D_eq: Drop volume-equivalent diameter (mm)
    Returns:
        r: The vertical-to-horizontal drop axis ratio. Note: the Scatterer class
        expects horizontal to vertical, so you should pass 1/dsr_thurai_2007
    """
	AVORB= 1.0148 -0.020465*D_eq -0.020048*D_eq**2 +0.003095*D_eq**3 -0.0001453*D_eq**4
	return AVORB
def dsr_bradnes_2003(D_eq):
	"""
	Drop shape relationship function from Bradnes2003
    (https://doi.org/10.1175/1520-0450(2003)042<0652:AEOADD>2.0.CO;2) paper.
    Arguments:
        D_eq: Drop volume-equivalent diameter (mm)
    Returns:
        r: The vertical-to-horizontal drop axis ratio. Note: the Scatterer class
        expects horizontal to vertical, so you should pass 1/dsr_thurai_2007
    """
	AVORB= 0.9951 +0.0251*D_eq -0.03644*D_eq**2+0.0050303*D_eq**3 -0.0002492*D_eq**4
	return AVORB
def dsr_chang_2009(D_eq):
	"""
	Drop shape relationship function from Chang2009
    (https://doi.org/10.1175/1520-0450(2003)042<0652:AEOADD>2.0.CO;2) paper.
    Arguments:
        D_eq: Drop volume-equivalent diameter (mm)
    Returns:
        r: The vertical-to-horizontal drop axis ratio. Note: the Scatterer class
        expects horizontal to vertical, so you should pass 1/dsr_thurai_2007
    """
	AVORB= 0.98287 +0.042514*D_eq -0.033439*D_eq**2+0.0043402*D_eq**3 -0.00019233*D_eq**4
	return AVORB