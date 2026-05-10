# %%
# Import modules    
import numpy as np
import os
import subprocess
import shutil
import warnings

folder = os.path.dirname(os.path.abspath(__file__))
inputfile = folder+'/fortran_tm/FOLDERNAME'
with open(inputfile, 'w') as f:
	f.write("'"+folder+'/fortran_tm/'+"'"+'\n')

class Scatterer(object):
	"""T-Matrix scattering from nonspherical particles.
	Class for simulating scattering from nonspherical particles with the 
	T-Matrix method. Uses a wrapper to the Fortran code by M. Mishchenko.
	Attributes: 
	radius :Equivalent radius (mm). [array]
	wavelength :Wavelength (cm).
	tempt : droplet temperature (C)
	NM : the number of integration points
	NRANK : the number of RANK.
	axis_ratio :Axis ratio of the particle.
	----------------------------------------------------------------
	!!!! Here is the parameters not used in this version. !!!!!
	shape: Particle shape. 
	alpha, beta: The Euler angles of the particle orientation (degrees).
	thet0, thet: The zenith angles of incident and scattered radiation (degrees).
	phi0, phi: The azimuth angles of incident and scattered radiation (degrees).
	m :Complex refractive index.
	orient: Orientation of the particle.
	or_pdf: Particle orientation PDF for orientational averaging.
	psd: Particle size distribution.
	psd_integrator: Set this to a PSDIntegrator instance to enable size
			distribution integration. If this is None (default), 
	"""
	#--------------------------------------------------------------	
	_attr_list = set(["radius", "axis_ratio", "wavelength", "IB",
				   "IOPt","tempt","NM","NRANK","Kw_sqr","anginc","ddelt","ndgs",
				   "alpha","beta","thet0","thet","phi0","phi","eastrt","eastop",
				   "eainc","ntheta","nphi","dclass","distyp","dmin","dmax",
				   "dpar1","dpar2","dpar3","m","orient","scatter","or_pdf",
				   "psd_integrator","psd"])
	
	_deprecated_aliases = {"axi": "radius",
		"lam": "wavelength",
		"eps": "axis_ratio",
		"rat": "radius_type",
		"np": "shape",
		"scatter": "orient"}
	
	def __init__(self, **kwargs):
		"""Initialize the Scatterer instance.           """
		self.radius = np.array([0.1])
		self.axis_ratio = np.array([1])
		self.wavelength = 1.0
		self.IB =8
		self.IOPT=1
		self.tempt = 25.0
		self.NM = 4
		self.NRANK = 8
		self.Kw_sqr = 0.93
		self.anginc = 90.0
		self.ddelt = 1e-3
		self.ndgs = 2
		self.alpha = 0.0
		self.beta = 0.0
		self.thet0 = 0.0
		self.thet = 0.0
		self.phi0 = 0.0
		self.phi = 180.0
		self.eastrt = 0.0
		self.eainc = 0.1
		self.ntheta = 9
		self.nphi = 7
		self.dclass = 1
		self.distyp=2
		self.dmin=0.0
		self.dmax=10.0
		self.dpar1=0.0
		self.dpar2=0.0
		self.dpar3=0.0
		self.m = complex(2,0)
		self.orient = None
		self.or_pdf = None
		self.psd = None
		self.psd_integrator = None
		# filedir = '/NAS1/mchlan/code/Tmatrix_Mueller/04_Tmatrix/'
		#--------------------------------------------------------------
		self.suppress_warning = kwargs["suppress_warning"] if \
		"suppress_warning" in kwargs else False

		for attr in self.__class__._deprecated_aliases:            
			if attr in kwargs:
				self._warn_deprecation(attr)
				self.__dict__[self._deprecated_aliases[attr]] = kwargs[attr]
				
		for attr in self._attr_list:
			if attr in kwargs:
				self.__dict__[attr] = kwargs[attr]

	def set_axisratio(self, fun):
		"""Set the axis ratio of the particle.

		Args:
			fun: A function of one variable that returns the axis ratio
				for a given equivalent volume diameter.
		"""
		self.axis_ratio = fun(self.radius)

	def set_geometry(self, geom):
		"""A convenience function to set the geometry variables.

		Args:
			geom: A tuple containing (thet0, thet, phi0, phi, alpha, beta).
			See the Scatterer class documentation for a description of these
			angles.
		"""
		(self.thet0, self.thet, self.phi0, self.phi, self.alpha, 
			self.beta) = geom

	def get_geometry(self):
		"""A convenience function to get the geometry variables.

		Returns:
			A tuple containing (thet0, thet, phi0, phi, alpha, beta).
			See the Scatterer class documentation for a description of these
			angles.
		"""
		return (self.thet0, self.thet, self.phi0, self.phi, self.alpha, 
			self.beta)	
	def _warn_deprecation(self, attr):
		if not self.suppress_warning:
			replacement = self._deprecated_aliases[attr]           
			warnings.simplefilter("always")
			warnings.warn(("The attribute '{attr}' is deprecated and may " + \
				"be removed in a future version. It has been renamed to " + \
				"'{replacement}'.").format(attr=attr, 
				replacement=replacement), DeprecationWarning)
			warnings.filters.pop(0)
	
	def __getattr__(self, name):
		if name == "_aliases":
			raise AttributeError
		if name in self._deprecated_aliases:
			self._warn_deprecation(name)
		name = self._deprecated_aliases.get(name, name)  
		return object.__getattribute__(self, name)


	def __setattr__(self, name, value):
		if name in self._deprecated_aliases:
			self._warn_deprecation(name)
		name = self._deprecated_aliases.get(name, name)
		object.__setattr__(self, name, value)  

	def write_input(self):
		# filedir = '/NAS1/mchlan/code/Tmatrix_Mueller/04_Tmatrix/'
		# make the inputfile : input, raininfo2.dat
		inputfile = folder+'/fortran_tm/input'
		with open(inputfile, 'w') as f:
			f.write(str(self.IB)+'\n')
			f.write(str(self.tempt)+'\n')
			f.write(str(self.wavelength)+'\n')

	def write_raininfo(self):
		# filedir = '/NAS1/mchlan/code/Tmatrix_Mueller/04_Tmatrix/'
		ML=len(self.radius)
		D_eq = self.radius/10 # input radius in cm
		AVORB = self.axis_ratio
		# print(AVORB)
		inputfile = folder+'/fortran_tm/raininfo.dat'
		with open(inputfile, 'w') as f:
			f.write(str(ML)+' '+str(self.IOPT)+'  # OF DATA SETS   IOPT=2 FOR SPONGY GRAUPEL'+'\n')
			for n in range(ML):
				f.write(str(self.NM).zfill(2)+' '
						+str(self.NRANK).zfill(2)+' '
						+"%6.4f" %D_eq[n]+' '  
						+"%7.5f" %AVORB[n]+'  '
						+"%5.2f" %self.anginc+' '
						+'\n')
			f.write('\n')
			f.write(' m'+' NR'+'    Deq'+'   AVORB'+'  anginc'+'\n')

	def write_tmat_filelist(self):
		# write the filelist for tmat_py.exe : input, raininfo.dat
		inputfile = folder+'/fortran_tm/filelist_tmat'
		with open(inputfile, 'w') as f:
			f.write("'"+folder+'/fortran_tm/input'+"'"+'\n')
			f.write("'"+folder+'/fortran_tm/raininfo.dat'+"'"+'\n')

	def write_mueller(self):
		# filedor
		ML=len(self.radius)
		inputfile = folder+'/fortran_tm/mueller.inp'
		icheck=0
		npart=1
		Dmin=self.dmin/10 # in cm
		Dmax=self.dmax/10 # in cm
		Dstep=0.01   # not used
		Distype=self.distyp
		Dpar1=self.dpar1
		Dpar2=self.dpar2
		if Dpar2<0.005: 
			Dpar2=0.005
		Dpar3=self.dpar3
		classname=['RAINDROPS',
			 		'RAINDROPS(FROZEN)',
					 'HAIL(DRY)',
					 'HAIL(SPONGY)',
					 'HAIL(WET)',
					 'GRAUPEL(DRY)',
					 'GRAUPEL(SPONGY)',
					 'GRAUPEL(WET)',
					 'CRYSTALS(NEEDLES)',
					 'CRYSTALS(PLATES)',
					 'AGGREGATES(SPONGY)',
					 'AGGREGATES(0.2)',
					 'AGGREGATES(0.5)',
					 'AGGREGATES(0.8)',
					]
		classar=np.zeros(len(classname))
		classar[self.dclass-1]=1
		classar = classar==1
		with open(inputfile, 'w') as f:
			f.write(str(icheck)+'                             CHECK (=0 , =1 CREATES DIAGNOSTIC FILE)\n')	
			f.write("%3.1f"%self.eastrt+'                          ELEVATION ANGLE START (DEGREES)\n')
			f.write("%3.1f"%self.eastrt+'                          ELEVATION ANGLE STOP (DEGREES)\n')
			f.write("%3.1f"%self.eainc+'                          ELEVATION ANGLE INCREMENT (DEGREES)\n')
			f.write("%2d"%self.ntheta+'                             ENTER ORDER OF INTEGRATION OVER THETA\n')
			f.write("%2d"%self.nphi+'                               ENTER ORDER OF INTEGRATION OVER PHI\n')
			f.write(str(ML)+'                          # OF DSD (DROP SIZE DISTRIBUTION PARAMETERS)\n')
			f.write(str(npart)+'                             # OF DIFFERENT PARTICLE TYPES (SPECIES)\n')
			for n in range(len(classname)):
				f.write(str(classar[n])[0]+'  '
						+"%5.3f"%Dmin+'  '
						+"%5.3f"%Dmax+'  '
						+"%5.3f"%Dstep+'  '
						+"%2d"%Distype+'  '
						+"%5.3f"%Dpar1+'  '
						+"%5.3f"%Dpar2+'  '
						+"%5.3f"%Dpar3+'     '
						+classname[n]+'\n'
				)
			f.write('\n')
			f.write('T/F DMIN DMAX DSTEP  DTYP DPAR1  DPAR2  DPAR3\n')
			f.write("   (ALL DIA'S IN CMS)\n")
			f.write('                      DTYP = DISTRIBUTION TYPE FOR CANTING ANGLE\n')
			f.write('                         0 = RANDOM\n')
			f.write('                         1 = SIMPLE HARMONIC (DPAR1=MEAN, DPAR2=THETAM)\n')
			f.write('                         2 = GAUSSIAN (DPAR1=MEAN, DPAR2=SIGMA)\n')
			f.write('                         3 = LANGEVIN (DPAR1=KAPPA, DPAR2=THETA0, DPAR3=PHI0)\n')
			f.write('                         4 = FISHER(MEAN THETA=0) (DPAR1=KAPPA)\n')

	def write_mueller_filelist(self):
		# write the filelist for tmat_py.exe : mueller.inp, out1_tmat
		inputfile = folder+'/fortran_tm/filelist_mueller'
		with open(inputfile, 'w') as f:
			f.write("'"+folder+'/fortran_tm/mueller.inp'+"'"+'\n')
			f.write("'"+folder+'/fortran_tm/out1_tmat'+"'"+'\n')


	def read_scat(self):
		""" 
		Here is read the scattering outputfile : scat_parall, scat_perpn
		and mueller outputfile : out1_Smatrix, out1_Mmatrix
		to get the fa. fb and matrix S and M.
		"""
		ML=len(self.radius)
		# filedir = '/NAS1/mchlan/code/Tmatrix_Mueller/'
		# read the outputfile : scat
		inputfile = folder+'/fortran_tm/scat_parall'
		self.SCAT_PAR,self.EXTN_PAR,self.RADA_PAR, \
		self.REFO_PAR,self.IMFO_PAR,self.REBA_PAR,self.IMBA_PAR=np.loadtxt(inputfile,unpack=True)
		inputfile = folder+'/fortran_tm/scat_perpn'
		self.SCAT_PER,self.EXTN_PER,self.RADA_PER, \
		self.REFO_PER,self.IMFO_PER,self.REBA_PER,self.IMBA_PER=np.loadtxt(inputfile,unpack=True)	
		# read the outputfile : mueller
		inputfile = folder+'/fortran_tm/out1_Smatrix'
		S=np.zeros([ML,4,4])
		data=np.loadtxt(inputfile)
		for n in range(ML):
			S[n,0,:]=data[n*4+0,:]
			S[n,1,:]=data[n*4+1,:]
			S[n,2,:]=data[n*4+2,:]
			S[n,3,:]=data[n*4+3,:]
		inputfile = folder+'/fortran_tm/out1_Kmatrix'
		K=np.zeros([ML,4,4])
		data=np.loadtxt(inputfile)
		for n in range(ML):
			K[n,0,:]=data[n*4+0,:]
			K[n,1,:]=data[n*4+1,:]
			K[n,2,:]=data[n*4+2,:]
			K[n,3,:]=data[n*4+3,:]
		self.S=S
		self.K=K
	
	def calc_scat(self):
		"""Calculate the scattering properties of the particle.
		"""
		wl=self.wavelength
		S=self.S
		K=self.K
		MAGKSQ=0.92
		PI=np.pi
		NFAC =1.E6*(wl**4)/((PI**5)*(MAGKSQ))
		NFAC = NFAC*4.*PI
		Fhh=0.5*(S[:,0,0] - S[:,0,1] - S[:,1,0]+ S[:,1,1])*1e6	
		Fvv=0.5*(S[:,0,0] + S[:,0,1] + S[:,1,0]+ S[:,1,1])*1e6
		Fvh=0.5*(S[:,0,0] - S[:,0,1] + S[:,1,0]- S[:,1,1])*1e6	
		# for Rhohv
		RHORE =-1.*((S[2,2] + S[3,3]))
		RHOIM =-1.*((S[3,2] - S[2,3]))

		self.Fhh=Fhh
		self.Fvv=Fvv
		self.Fvh=Fvh
		self.RHORE=RHORE
		self.RHOIM=RHOIM
		#--------------------------------------------------------------
		# forward scatter
		#Re[fa(0)-fb(0)] --> for kdp
		DPHA = K[2,3]*0.1 
		ATTH = 0.4343*(K[0,0] - K[0,1] ) 
		ATTV = 0.4343*(K[0,0] + K[0,1] ) 
		DATT = -0.4343*2.* K[0,1]

		self.DPHA=DPHA
		self.ATTH=ATTH
		self.ATTV=ATTV
		self.DATT=DATT

	def get_table(self):
		"""Get the the table from fortran : 04-Tmatrix.
		"""
		# filedir = '/NAS1/mchlan/code/Tmatrix_Mueller/04_Tmatrix/'
		# make the inputfile : input, raininfo2.dat
		self.write_input()
		self.write_raininfo()
		# self.write_tmat_filelist()
		self.write_mueller()
		# self.write_mueller_filelist()
		with open(f'{folder}/fortran_tm/FOLDERNAME', 'r') as stdin, open(f'{folder}/fortran_tm/log_tmat', 'w') as stdout:
			subprocess.run([f'{folder}/fortran_tm/tmat_py.exe'], stdin=stdin, stdout=stdout, check=True, cwd=f'{folder}/fortran_tm')
		# The legacy Mueller executable opens hard-coded relative names for
		# non-rain species. Keep those aliases synced to the table just built.
		for alias in ('inp', 'inp1', 'inp2'):
			shutil.copyfile(f'{folder}/fortran_tm/out1_tmat', f'{folder}/fortran_tm/{alias}')
		with open(f'{folder}/fortran_tm/FOLDERNAME', 'r') as stdin, open(f'{folder}/fortran_tm/log_mueller', 'w') as stdout:
			subprocess.run([f'{folder}/fortran_tm/mueller_py.exe'], stdin=stdin, stdout=stdout, check=True, cwd=f'{folder}/fortran_tm')
		self.read_scat()
		# self.calc_scat()

	# def calc_(self)
