

import numpy as np
import pyfits as pf

dlambda      =702.4   # ACS F814w band width
lam_f814w    =8057.0  # Wavelenght of the band center

zp_cosABmag  =25.947  # Zero point of ACS F814w in AB mag
expt_cos     =1000.0  # Exposure time of ACS
gain_acs     =1.0     # Gain value of F814w in ACS
pixel_acs    =15      # in microns with 1 micron=1e-6 meters

zp_lsstABmag =24.000  # Lsst zero point in AB mag
nvisit_lsst  =180.0   # Number of visit for the same area
expt_lsst    =30.0    # Exposure time of each visit
gaini_lsst   =2.0     # Temporarily to set 2.0???????? need to be further confirmed
Alsst        =319.0/9.6

zp_desABmag =24.000  # DES zero point in AB mag
nvisit_des  =10.0    # Number of visit for the same area
expt_des    =90.0    # Exposure time of each visit
gaini_dex   =2.0     # Temporarily to set 2.0???????? need to be further confirmed
Ades        =319.0/9.6

#def photon2ABmag(dlambda,lam_f814w,image):
def photon2ABmag(image):
	idx       =np.where(image<=0.0)
	image[idx]=1e-7
	image_Jy  =image/1.51e+7/(dlambda/lam_f814w)
	image_AB  =-2.5*np.log10(image_Jy/3631.0)

	return image_AB

	
def ABmag2lsstccd(image_AB):	
	exptime   =expt_lsst*nvisit_lsst
	s0        =Alsst*1.681
	image_lsst=s0*10.0**(-0.4*(image_AB-zp_lsstABmag))*exptime
	#assuming gain lsst =1
	return image_lsst

def ABmag2desccd(image_AB):	
	exptime   =expt_des*nvisit_des
	s0        =Ades*1.681
	image_des =s0*10.0**(-0.4*(image_AB-zp_desABmag))*exptime
	return image_des

def ABmag2sdssccd(image_AB):	
	exptime   =54
	s0        =Ades*1.681
	image_sdss=s0*10.0**(-0.4*(image_AB-zp_desABmag))*exptime
	return image_sdss


	
