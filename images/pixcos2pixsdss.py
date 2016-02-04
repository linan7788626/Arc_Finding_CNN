import pyfits as pf
import numpy as np

#The conversion between the value of ACS pix to SDSS pixel is simple.
#The ACS pix value after drizzle is in unit of electrons/seconds which
#is physical and SDSS should have the same electrons/seconds from the same
# source. But the exposure time is only 53 seconds for SDSS, so the conversion
# is just ACS pixel value \time 53/gain. Here gain is the parameter for SDSS
# CCD relating electrons to DN(Data Number from CCD).
gain    =4.7
expsdss =53.0
aa_sdss =-24.149
aa_cos  =25.523
kk      =0.156347
airmass =1.201824

#Convert CCD value of cosmos in cps to counts in SDSS CCD.
def pixcos2pixsdss(image):
	image=image*expsdss/gain	
	
	return image*10**(aa_sdss+aa_cos)


	return im_mag
#Convert ST magnitude to CCD value(SDSS)	
def mag2sdssccd(image):
	factor= 10.0**(0.4*(aa_sdss+kk*airmass))
	frat  = 10.0**(image/(-2.5))
	im_ccd= (frat/factor)*expsdss/gain
	
	return im_ccd	
#Convert ACS ccd value to ST magnitude
def cosccd2mag(image):
	im_mag=-2.5*np.log10(image)+aa_cos
	return im_mag

#conver ST magnitude to cosmos CCD value
def mag2cosccd(image):
	im_ccd =10.0**((image-aa_cos)/(-2.5))

	return im_ccd
	
