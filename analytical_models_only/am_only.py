import numpy as np
import pylab as pl
import pyfits
#import scipy.ndimage.filters as snf
import scipy.signal as ss
import mycosmology as mm
#import congrid
#import rebin_array_2d
from astropy.cosmology import Planck13


nMgyCount_r=0.004760406   # nanomaggies per count for SDSS detector.
sky_r      =5.98          # SDSS typical r band sky
softbias   =1000.0        # SDSS softbias
Mgy2nanoMgy=10e+9         # nanoMaggy to Maggy
skycount=sky_r/(nMgyCount_r)


def noise_map(nx1,nx2,nstd,NoiseType):
    if NoiseType=='Poisson':
    	noise=np.random.poisson(nstd,(nx1,nx2))-nstd
    if NoiseType=='Gaussian':
    	noise=nstd*np.random.normal(0.0,1.0,(nx1,nx2))
    return noise
#--------------------------------------------------------------------
def make_r_coor(nc,dsx):

    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0

    x2,x1 = np.meshgrid(x1,x2)
    return x1,x2
def make_c_coor(nc,dsx):

    bsz = nc*dsx
    x1,x2 = np.mgrid[0:(bsz-dsx):nc*1j,0:(bsz-dsx):nc*1j]-bsz/2.0+dsx/2.0
    return x1,x2

#--------------------------------------------------------------------
def lens_equation_sie(x1,x2,lpar):
    xc1 = lpar[0]   #x coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]   #y coordinate of the center of lens (in units of Einstein radius).
    q   = lpar[2]   #Ellipticity of lens.
    rc  = lpar[3]   #Core size of lens (in units of Einstein radius).
    re  = lpar[4]   #Einstein radius of lens.
    pha = lpar[5]   #Orintation of lens.

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1-xc1)*cosa+(x2-xc2)*sina
    xt2 = (x2-xc2)*cosa-(x1-xc1)*sina

    phi = np.sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc)
    sq = np.sqrt(1.0-q*q)
    pd1 = phi+rc/q
    pd2 = phi+rc*q
    fx1 = sq*xt1/pd1
    fx2 = sq*xt2/pd2
    qs = np.sqrt(q)

    a1 = qs/sq*np.arctan(fx1)
    a2 = qs/sq*np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 =-sina

    fx11 = xt11/pd1-xt1*(xt1*q*q*xt11+xt2*xt21)/(phi*pd1*pd1)
    fx22 = xt22/pd2-xt2*(xt1*q*q*xt12+xt2*xt22)/(phi*pd2*pd2)
    fx12 = xt12/pd1-xt1*(xt1*q*q*xt12+xt2*xt22)/(phi*pd1*pd1)
    fx21 = xt21/pd2-xt2*(xt1*q*q*xt11+xt2*xt21)/(phi*pd2*pd2)

    a11 = qs/(1.0+fx1*fx1)*fx11
    a22 = qs/(1.0-fx2*fx2)*fx22
    a12 = qs/(1.0+fx1*fx1)*fx12
    a21 = qs/(1.0-fx2*fx2)*fx21

    rea11 = (a11*cosa-a21*sina)*re
    rea22 = (a22*cosa+a12*sina)*re
    rea12 = (a12*cosa-a22*sina)*re
    rea21 = (a21*cosa+a11*sina)*re

    y11 = 1.0-rea11
    y22 = 1.0-rea22
    y12 = 0.0-rea12
    y21 = 0.0-rea21

    jacobian = y11*y22-y12*y21
    mu = 1.0/jacobian

    res1 = (a1*cosa-a2*sina)*re
    res2 = (a2*cosa+a1*sina)*re
    return res1,res2,mu
#--------------------------------------------------------------------
def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x-xcen)*np.cos(phirad)+(y-ycen)*np.sin(phirad)
    ynew = (y-ycen)*np.cos(phirad)-(x-xcen)*np.sin(phirad)
    return (xnew,ynew)

def gauss_2d(x, y, par):
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt(((xnew**2)*par[4]+(ynew**2)/par[4]))/np.abs(par[1])
    res = par[0]*np.exp(-res0**2.0)
    return res

def re_sv(sv,z1,z2):
    res = 4.0*np.pi*(sv**2.0/mm.vc**2.0)*mm.Da2(z1,z2)/mm.Da(z2)*mm.apr
    return res

#----Fundamental Plain------------------------------
def Brightness(Re,Vd):
    a       =1.49
    b       =0.2
    c       =-8.778
    mag_e   =((np.log10(Re)-a*np.log10(Vd)-c)/b)+20.09 # Bernardi et al 2003
    nanoMgy =Mgy2nanoMgy*10.0**(-(mag_e-22.5)/2.5)
    counts  =nanoMgy/nMgyCount_r
    return counts

def de_vaucouleurs_2d(x,y,par):
    #[I0, Re, xc1,xc2,q,pha]
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt((xnew**2)*par[4]+(ynew**2)/par[4])/par[1]
    #res = par[0]*np.exp(-par[1]*res0**0.25)
    res = par[0]*np.exp(-7.669*(res0**0.25-1.0))
    return res

##----de Vaucouleurs profile-------------------------
#def deVaucouleurs(x,y,xc,yc,counts,R,e,phi):
    #theta   =phi*np.pi/180.

    #xx      =x-xc
    #yy      =y-yc
    #rx      =xx*np.cos(theta)+yy*np.sin(theta)
    #ry      =-xx*np.sin(theta)+yy*np.cos(theta)
    #rr      =np.sqrt(rx*rx/(1.0-e)+ry*ry*(1.0-e))
    #image   =counts*np.exp(-7.669*((rr/R)**0.25-1.0))
    #soften  =counts*np.exp(-7.669*((0.02)**0.25-1.0))
    #ix      =np.where(image>=soften)
    #image[ix]=soften


    #return image
#--------------------------------------------------------------------
def main():

    #dsx     = 0.396         # pixel size of SDSS detector.
    dsx = 0.05         # pixel size of SDSS detector.
    R  = 2.9918     #vd is velocity dispersion.
    zl = 0.2     #zl is the redshift of the lens galaxy.
    zs = 1.0
    vd = 220    #Velocity Dispersion.
    nnn = 128      #Image dimension
    bsz = dsx*nnn
    nstd = 59

    xx01 = np.linspace(-bsz/2.0,bsz/2.0,nnn)+0.5*dsx
    xx02 = np.linspace(-bsz/2.0,bsz/2.0,nnn)+0.5*dsx
    xi2,xi1 = np.meshgrid(xx01,xx02)
    #----------------------------------------------------------------------
    g_amp = 10.0       # peak brightness value
    g_sig = 0.02    # Gaussian "sigma" (i.e., size)
    g_xcen = 0.06   # x position of center (also try (0.0,0.14)
    g_ycen = 0.11    # y position of center
    g_axrat = 1.0   # minor-to-major axis ratio
    g_pa = 0.0      # major-axis position angle (degrees) c.c.w. from x axis
    gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
    #----------------------------------------------------------------------
    #g_source = 0.0*xi1
    #g_source = gauss_2d(xi1,xi2,gpar) # modeling source as 2d Gaussian with input parameters.
    #----------------------------------------------------------------------
    xc1 = 0.0       #x coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0       #y coordinate of the center of lens (in units of Einstein radius).
    q   = 0.7       #Ellipticity of lens.
    rc  = 0.0       #Core size of lens (in units of Einstein radius).
    re  = re_sv(vd,zl,zs)       #Einstein radius of lens.
    pha = 45.0      #Orintation of lens.
    lpar = np.asarray([xc1,xc2,q,rc,re,pha])
    #----------------------------------------------------------------------
    ai1,ai2,mua = lens_equation_sie(xi1,xi2,lpar)
    yi1 = xi1-ai1
    yi2 = xi2-ai2
    #----------------------------------------------------------------------
    gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
    g_limage = gauss_2d(yi1,yi2,gpar)

    #print np.sum(g_limage)
    #pl.figure()
    #pl.contourf(g_limage)
    #pl.colorbar()

    #g_limage = rebin_array_2d.rebin(g_limage,[128,128])
    #-------------------------------------------------------------
    dA = Planck13.comoving_distance(zl).value*1000./(1+zl)
    Re = dA*np.sin(R*np.pi/180./3600.)
    counts  =Brightness(R,vd)
    vpar = np.asarray([counts,Re,xc1,xc2,q,pha])
    #g_lens = deVaucouleurs(xi1,xi2,xc1,xc2,counts,R,1.0-q,pha)
    g_lens = de_vaucouleurs_2d(xi1,xi2,vpar)

    g_lens = g_lens/1e4

    file_psf = "../PSF_and_noise/sdsspsf.fits"
    g_psf = pyfits.getdata(file_psf)-1000.0
    g_psf = g_psf/np.sum(g_psf)
    g_limage = ss.fftconvolve(g_limage+g_lens,g_psf,mode="same")

    g_noise = noise_map(nnn,nnn,nstd,"Gaussian")


    g_limage = g_limage+(g_noise+skycount)/200.0

    pl.figure()
    pl.contourf(g_limage)
    pl.colorbar()

    output_filename = "../output_fits/test.fits"
    pyfits.writeto(output_filename,g_limage,clobber=True)

    pl.show()

    return 0
#------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
