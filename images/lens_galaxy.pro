;+
; NAME:
;    lens_galaxy
;
; PURPOSE:
;    Simulate SDSS-like lens image simply using
;    de Vaucouleure profile given effective radius
;    in arcseconds, velocity dispersion, redshift of
;    the lens, image size in unit of pixels by pixels,
;    ellipticity of the galaxy and orientation of
;    the galaxy. You could also add noise to the image
;    , either Gaussian or Poisson wtih the variance of
;    the noise.
;
; CALLING SEQUENCE:
;    lens_galaxy,[R_eff=,zlens=,VelDis=,Npix=,E=,Phi=,
;    NoiseVar=, NoiseType= ]
;
; OPTIONAL KEYWORD PARAMETERS:
;    R_eff         - Effective radius of the lens galaxy
;    zlens         - The redshift of the lens galaxy
;    VelDis        - Velocity dispersion of the halo hosting the lens galaxy
;    Npix          - The size of the postage stamp image
;    E             - Ellipticity of the lens galaxy
;    Phi           - Orientation of the lens galaxy wrt the North.
;    NoiseVar      - The variance of the noise
;    NoiseType     - Gaussian or Poisson
;
; PROCEDURES CALLED(all included in idlutils,except for imdisp.pro)
;    readfits
;    writefits
;    angdidis
;    deVaucouleurs
;    brightness
;    *imdisp         -If you do not want to use this procedure, please comment this line.
;
; HISTORY
;    Written by Wentao Luo, Feb. 1st 2016 from
;    python version 'lens_galaxy.py'
;---------------------------------------------------
;
;-- Functions called in main program---------------
;-- brightness function---------------------------
FUNCTION brightness,Re,Vd
	COMMON share,Mgy2nanoMgy,nMgyCount_r,pi
	a       =1.49
	b       =0.2
	c       =-8.778
        ; These parameters are from Bernardi et al 2003
	mag_e   =(ALOG10(Re)-a*ALOG10(Vd)-c)/b+20.09
	nanoMgy =Mgy2nanoMgy*10.0^(-(mag_e-22.5)/2.5)
	counts  =nanoMgy/nMgyCount_r

	RETURN,counts
END
;-------------------
FUNCTION deVaucouleurs,Re,R_eff,Vd,e,phi,Npix
	COMMON share,Mgy2nanoMgy,nMgyCount_r,pi
	count   =BRIGHTNESS(Re,Vd)
	x       =FINDGEN(Npix)
	y       =FINDGEN(Npix)
	x       =x#REPLICATE(1,Npix)
	y       =REPLICATE(1,Npix)#y
	xc      =FIX(Npix/2.0)
	yc      =FIX(Npix/2.0)
	pixsize =0.396

	xx      =x-xc
    yy      =y-yc
	theta   =phi*pi/180.
	Rpix    =R_eff/pixsize
	rx      =xx*COS(theta)+yy*SIN(theta)
	ry      =-xx*SIN(theta)+yy*COS(theta)
	rr      =SQRT(rx*rx/(1.0-e)+ry*ry*(1.0-e))
	image   =count*EXP(-7.669*((rr/Rpix)^0.25-1.0))

	soften  =count*EXP(-7.669*((0.2)^0.25-1.0))
	ind     =WHERE(image ge soften)
	image[ind]=soften

	RETURN,image
END
;-- End of call functions-------------------------


PRO lens_galaxy, R_eff    =R_eff,      $
				 zlens    =zlens,      $
				 VelDis   =VelDis,     $
				 Npix     =Npix,       $
				 E        =E,          $
				 Phi      =phi,        $
				 NoiseVar =NoiseVar,   $
				 NoiseType=NoiseType

IF NOT KEYWORD_SET(R_eff)     THEN R_eff    =3.0
IF NOT KEYWORD_SET(zlens)     THEN zlens    =0.3
IF NOT KEYWORD_SET(VelDis)    THEN VelDis   =200.0
IF NOT KEYWORD_SET(Npix)      THEN Npix     =128
IF NOT KEYWORD_SET(E)         THEN E        =0.0
IF NOT KEYWORD_SET(Phi)       THEN Phi      =0.0
IF NOT KEYWORD_SET(NoiseVar)  THEN NoiseVar =58.0
IF NOT KEYWORD_SET(NoiseType) THEN NoiseType='Gaussian'

COMMON share,Mgy2nanoMgy,nMgyCount_r,pi
;-- Cosmological and Observational(SDSS) Parameters
pi         =3.14159265
H0         =69.6    ; The Hubble constant at z=0 from Plank cosmology
Omeg_m     =0.286
Omeg_l     =0.714

sky_r      =5.98
Mgy2nanoMgy=10e+9
nMgyCount_r=0.004760406
softbias   =1000.0
;-- End of Parameters needed

nstd       =SQRT(NoiseVar)
dA         =ANGDIDIS(zlens,Omeg_m,Omeg_l)*3000.0*1000.0
Re         =dA*SIN(R_eff*pi/180./3600.)

imgal      =deVaucouleurs(Re,R_eff,VelDis,E,Phi,Npix)
skycount   =sky_r/nMgyCount_r

IF NoiseType EQ 'Poisson'  THEN $
     noise = RANDOMN(seed,Npix,Npix,POISSON=SQRT(NoiseVar))-SQRT(NoiseVar)
IF NoiseType EQ 'Gaussian' THEN $
     noise = SQRT(NoiseVar)*RANDOMN(seed,Npix,Npix)

IMDISP,imgal+noise+softbias+skycount
;writefits,'test_im.fits',imgal+noise+softbias+skycount


END
