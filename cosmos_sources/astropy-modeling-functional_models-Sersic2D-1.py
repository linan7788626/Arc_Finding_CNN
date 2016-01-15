import numpy as np
from astropy.modeling.models import Sersic2D

def build_source(x1,x2,xc1,xc2,ell,pha,amp,reff,nsersic,nnn,boxsize):
    pha = np.rad2deg(pha)
    mod = Sersic2D(amp, reff, nsersic, xc1, xc2,ell,pha)
    res = mod(x1,x2)
    return res

def readin_morphology_catalog(filename="/Users/uranus/Desktop/working/cosmos_catalogs/astrometry.cosmos_morph_cassata_1_116795.tbl"):
    cgal_ids,cgal_mags,cgal_reff,cgal_aratio,cgal_class,cgal_classweight = np.loadtxt(filename,dtype=str, comments='/', skiprows=43, usecols=(0,3,5,10,11,12), unpack=True)

    idx1 = cgal_mags != 'null'
    idx2 = cgal_reff != 'null'
    idx3 = cgal_aratio != 'null'
    idx4 = cgal_class != 'null'
    idx5 = cgal_classweight != 'null'

    idx = idx1&idx2&idx3&idx4&idx5

    cgal_ids_tmp = cgal_ids[idx].astype('int')
    cgal_mags_tmp = cgal_mags[idx].astype('float32')
    cgal_reff_tmp = cgal_reff[idx].astype('float32')
    cgal_aratio_tmp = cgal_aratio[idx].astype('float32')
    cgal_class_tmp = cgal_class[idx].astype('int')
    cgal_classweight_tmp = cgal_classweight[idx].astype('float32')

    return cgal_ids_tmp,cgal_mags_tmp,cgal_reff_tmp,cgal_aratio_tmp,cgal_class_tmp,cgal_classweight_tmp

def readin_redshifts_catalog(filename="/Users/uranus/Desktop/working/cosmos_catalogs/astrometry.cosmos_zphot_mag253922.tbl"):
    cgal_ids,cgal_zs,cgal_mags = np.loadtxt(filename,dtype=str,
                                            comments='/', skiprows=117,
                                            usecols=(0,6,46), unpack=True)
    idx1 = cgal_zs != 'null'
    idx2 = cgal_mags != 'null'
    idx = idx1&idx2

    cgal_ids_tmp = cgal_ids[idx].astype('int')
    cgal_zs_tmp = cgal_zs[idx].astype('float32')
    cgal_mags_tmp = cgal_mags[idx].astype('float32')

    idz = cgal_zs_tmp > 1.0

    return cgal_ids_tmp[idz],cgal_zs_tmp[idz],cgal_mags_tmp[idz]

if __name__ == '__main__':
    bsz = 3.0
    nnn = 256
    dsx = bsz/nnn
    xc1 = 0.0
    xc2 = 0.0

    x1 = np.linspace(0.0,bsz-dsx,nnn-1)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0.0,bsz-dsx,nnn-1)-bsz/2.0+dsx/2.0
    x1,x2 = np.meshgrid(x1,x2)


    #print id_morph_list
    #print mags_morph_list
    #print reff_list
    #print aratio_list
    #print class_list
    #print classweight_list

    id_list,zs_list,mags_list = readin_redshifts_catalog()
    id_morph_list,mags_morph_list,reff_list,aratio_list,class_list,classweight_list = readin_morphology_catalog()


    id_overlap = np.intersect1d(id_morph_list,id_list)
    #id_overlap = np.in1d(id_morph_list,id_list).nonzero()
    idx = np.random.randint(0,len(id_overlap))
    #print len(id_overlap)
    #print type(id_overlap[0])
    #print type(id_list[0])
    #print type(id_morph_list[0])
    #print id_list[id_list == id_overlap[idx]]
    #print id_morph_list[id_morph_list == id_overlap[idx]]

    idg = id_morph_list == id_overlap[idx]

    source_img = build_source(x1,x2,xc1,xc2,aratio_list[idg],45.0,mags_morph_list[idg],reff_list[idg],class_list[idg],nnn,bsz)

    import pylab as pl
    pl.contourf(np.log10(source_img))
    pl.colorbar()
    pl.show()

    #print len(id_list)
    #print id_list
    #print zs_list
    #print mags_list


#plt.figure()
#plt.imshow(log_img, origin='lower', interpolation='nearest',
           #vmin=-1, vmax=2)
#plt.xlabel('x')
#plt.ylabel('y')
#cbar = plt.colorbar()
#cbar.set_label('Log Brightness', rotation=270, labelpad=25)
#cbar.set_ticks([-1, 0, 1, 2], update_ticks=True)
#plt.show()
