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
aa_cos  =-25.523
kk      =0.156347
airmass =1.201824

#Convert CCD value of cosmos in cps to counts in SDSS CCD.
def pixcos2pixsdss(image):
	image=image*expsdss/gain	
	
	return image*10**(aa_sdss-aa_cos)

#Convert CCD value in counts(SDSS) to Pogson magnitude.
def sdssccd2mag(image):
	factor= 10.0**(0.4*(aa_sdss+kk*airmass))
	ff0   = (image*gain/expsdss)*factor
	im_mag= -2.5*np.log10(factor*(image*gain/expsdss))

	return im_mag
#Convert Pogson magnitude to CCD value(SDSS)	
def mag2sdssccd(image):
	factor= 10.0**(0.4*(aa_sdss+kk*airmass))
	ff0   = 10.0**(image/(-2.5))
	im_ccd= (ff0/factor)*expsdss/gain
	
	return im_ccd	

#convert CCD value in cps(count per second with gain=1) to Pogson magnitude
#, but given the SDSS zeropoint so that we calculate the magnitude should be
# SDSS-like.
def cosccd2mag(image):
	factor=10.0**(0.4*(aa_sdss+kk*airmass))
	ff0   =image*factor      # gain=1 and its alreay count per second.
	im_mag=-2.5*np.log10(ff0)

	return im_mag

#conver Pogson magnitude to cosmos CCD value
def mag2cosccd(image):
	factor=10.0**(0.4*(aa_sdss+kk*airmass))
	ff0   =10.0**(image/(-2.5))
	im_ccd=ff0/factor

	return im_ccd
	
