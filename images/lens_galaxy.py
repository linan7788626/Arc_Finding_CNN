#!/home/wtluo/anaconda/bin/python2.7

import pyfits as pf
import numpy as np
from optparse import OptionParser
from astropy.cosmology import Planck13
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#import galsim as gs

nMgyCount_r=0.004760406   # nanomaggies per count for SDSS detector.
pixsize    =0.396         # pixel size of SDSS detector.
sky_r      =5.98          # SDSS typical r band sky
softbias   =1000.0        # SDSS softbias
Mgy2nanoMgy=10e+9         # nanoMaggy to Maggy

#----Fundamental Plain------------------------------
def Brightness(Re,Vd):
	a       =1.49
	b       =0.2
	c       =-8.778
	mag_e   =((np.log10(Re)-a*np.log10(Vd)-c)/b)+20.09 # Bernardi et al 2003
	nanoMgy =Mgy2nanoMgy*10.0**(-(mag_e-22.5)/2.5)

	counts  =nanoMgy/nMgyCount_r

	return counts
#----de Vaucouleurs profile-------------------------
def deVaucouleurs(R,Re,Vd,e,phi,Npix):
	count   =Brightness(R,Vd)
	x,y     =np.mgrid[:Npix,:Npix]
	Rpix    =R/pixsize
	theta   =phi*np.pi/180.
	xc      =int(Npix/2)
	yc      =int(Npix/2)

	xx      =x-xc
	yy      =y-yc
	rx      =xx*np.cos(theta)+yy*np.sin(theta)
	ry      =-xx*np.sin(theta)+yy*np.cos(theta)
	rr      =np.sqrt(rx*rx/(1.0-e)+ry*ry*(1.0-e))
	image   =count*np.exp(-7.669*((rr/Rpix)**0.25-1.0))
	soften  =count*np.exp(-7.669*((0.02)**0.25-1.0))
	ix      =np.where(image>=soften)
	image[ix]=soften


	return image

#----End of all funcs-----------------------

parser=OptionParser()
parser.add_option("--R_eff",dest="R_eff",default=2.9918,
		    help="Effective radius of the galaxy in arcseconds",metavar="value",
		    type="float")
parser.add_option("--zlens",dest="zlens",default=0.298,
		    help="the redshift of galaxy",metavar="value",
		    type="float")
parser.add_option("--VelDis",dest="VelDis",default=300.0,
		    help="Velocity Dispersion ",metavar="value",
		    type="float")
parser.add_option("--Npix",dest="Npix",default=128,
		    help="Number of Pixels for one dimension ",metavar="value",
		    type="float")
parser.add_option("--E",dest="E",default=0.0,
		    help="Ellipticity of galaxy ",metavar="value",
		    type="float")
parser.add_option("--Phi",dest="Phi",default=0,
		    help="Orientation of galaxy wrt North",metavar="value",
		    type="float")
parser.add_option("--NoiseVar",dest="NoiseVar",default=0.0,
                    help="Noise variance",metavar="value",
                    type="float")
parser.add_option("--NoiseType",dest="NoiseType",default="Gaussian",
                    help="Gaussian or Poisson",metavar="value",
                    type="string")

(o,args)=parser.parse_args()
R       =o.R_eff     #vd is velocity dispersion.
zl      =o.zlens     #zl is the redshift of the lens galaxy.
Vd      =o.VelDis    #Velocity Dispersion.
Npix    =o.Npix      #Image dimension
e       =o.E         #Ellipticity of simulated galaxy
phi     =o.Phi       #Orientation angle wrt North
nstd    =np.sqrt(o.NoiseVar)

dA      =Planck13.comoving_distance(zl).value*1000./(1+zl)
Re      =dA*np.sin(R*np.pi/180./3600.)
print Re
imgal   =deVaucouleurs(R,Re,Vd,e,phi,Npix)

skycount=sky_r/(nMgyCount_r)
if o.NoiseType=='Poisson':
        print '## Hey! You choose Poisson noise. Which is not good for now.'
        noise=np.random.poisson(nstd,(Npix,Npix))-nstd
if o.NoiseType=='Gaussian':
        print '## Hey! You choose Gaussian noise.'
        noise=nstd*np.random.normal(0.0,1.0,(Npix,Npix))


import matplotlib.pyplot as plt
plt.imshow((imgal+noise+skycount+softbias),interpolation='Nearest',cmap=cm.gray)
plt.colorbar()
plt.title('R_eff=10')
plt.show()

