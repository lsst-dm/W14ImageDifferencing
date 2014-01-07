import os
import os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass

def integrate(sed, bp):
    wavelen = sed.wavelen
    fnu = sed.fnu
    if sed.needResample(wavelen=wavelen, wavelen_match=bp.wavelen):
        wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
    return np.sum(fnu * bp.phi)

filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
uBand  = Bandpass()
gBand  = Bandpass()
rBand  = Bandpass()
iBand  = Bandpass()
zBand  = Bandpass()
yBand  = Bandpass()
uBand.readThroughput(os.path.join(filtDir, "total_u.dat"))
gBand.readThroughput(os.path.join(filtDir, "total_g.dat"))
rBand.readThroughput(os.path.join(filtDir, "total_r.dat"))
iBand.readThroughput(os.path.join(filtDir, "total_i.dat"))
zBand.readThroughput(os.path.join(filtDir, "total_z.dat"))
yBand.readThroughput(os.path.join(filtDir, "total_y.dat"))

catDir = os.environ["CAT_SHARE_DATA"]
aStar  = Sed() 
gStar  = Sed() 
mStar  = Sed() 
aStar.readSED_flambda(os.path.join(catDir, "data/starSED/kurucz", "kp01_9750.fits_g45_9830.gz"))
gStar.readSED_flambda(os.path.join(catDir, "data/starSED/kurucz", "km20_6000.fits_g30_6020.gz"))
mStar.readSED_flambda(os.path.join(catDir, "data/starSED/mlt", "m2.0Full.dat"))

for sname, sed in zip(("A", "G", "M"), (aStar, gStar, mStar)):
    sed.flambdaTofnu() 
    rBand.sbTophi()
    rflux = integrate(sed, rBand)
    for bpname, bp in zip(("u", "g", "i", "z", "y"),
                  (uBand, gBand, iBand, zBand, yBand)):
        bp.sbTophi()
        flux = integrate(sed, bp)
        dmag = -2.5 * np.log10(flux/rflux)
        print sname, "%s-r"%bpname, dmag
