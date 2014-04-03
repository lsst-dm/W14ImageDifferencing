import os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
fluxFormatter = FormatStrFormatter('%d')
dcrFormatter = FormatStrFormatter('%+0.2f')

# THIS ONE LOOKS AT THE AIRMASS-DEPENDENT FILTER PROFILE
# And only looks at 30,50 deg to compare to phosim outputs.

#def refract(wavelength, zd, docorr=True, P=600.0, T=7.0, f=8.0):
def refract(wavelength, zd, docorr=True, P=520.0, T=20.0, f=8.0):
    # ZD = zenith distance in radians
    # Wavelength = microns
    xnm1  = 10**-6 * ( 64.328 + 29498.1 / (146.0 - 1/wavelength**2) + 255.4 / (41.0 - 1/wavelength**2) )

    if docorr:
        # Correct for temperature and pressure
        # T in C
        # P in mm Hg
        xnm1 *= P * (1 + (1.049 - 0.0157 * T) * 10**-6 * P) / (720.883 * (1 + 0.003661 * T))

        # Correct for water vapor in atmosphere
        # f = water vapor pressure in mm of Hg
        # T = air temperature in C
        xnm1 -= 10**-6 * f * (0.0624 - 0.000680/wavelength**2) / (1 + 0.003661 * T)

    xn = xnm1 + 1
    r0 = (xn**2 - 1) / (2 * xn**2)
    of  = r0 * np.tan(zd)
    return of

if __name__ == "__main__":
    tp_dir = os.environ["THROUGHPUTS_DIR"]
    common_components = map(lambda x: os.path.join(tp_dir, "baseline", x), 
                            ["detector.dat", "m1.dat", "m2.dat", "m3.dat", 
                             "lens1.dat", "lens2.dat", "lens3.dat"])

    seds = []
    catDir = os.environ["CAT_SHARE_DATA"]
    oStar  = Sed() 
    gStar  = Sed() 
    mStar  = Sed() 
    oStar.readSED_flambda(os.path.join(catDir, "data/starSED/kurucz", "kp01_9750.fits_g45_9830.gz"))
    gStar.readSED_flambda(os.path.join(catDir, "data/starSED/kurucz", "km20_6000.fits_g30_6020.gz"))
    mStar.readSED_flambda(os.path.join(catDir, "data/starSED/mlt", "m2.0Full.dat"))
    seds.append(oStar)
    seds.append(gStar)
    seds.append(mStar)

    for sed in seds:
        sed.flambdaTofnu() 

    filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
    for airmass, zd in zip((15, 11), (50, 30)):
        for sidx, sed in enumerate(seds):
            for bpname in ("g", "r", "i"):
                bp = Bandpass(wavelen_max=1200) 
                components = common_components + [os.path.join(tp_dir, "baseline", "filter_%s.dat" % (bpname)), 
                                                  os.path.join(tp_dir, "atmos", "atmos_%d.dat" % (airmass))]
                #print components
                bp.readThroughputList(components)
                bp.sbTophi()

                wavelen           = sed.wavelen
                fnu               = sed.fnu
                wavelen, fnu      = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
                wavelenf, flambda = sed.fnuToflambda(wavelen, fnu)

                waveleng            = gStar.wavelen
                fnug                = gStar.fnu
                waveleng, fnug      = gStar.resampleSED(waveleng, fnug, wavelen_match=bp.wavelen)
                wavelenfg, flambdag = gStar.fnuToflambda(waveleng, fnug)

                flux    = fnu * bp.phi
                fluxg   = fnug * bp.phi
            
                leff1   = np.exp(np.sum(fnu * bp.phi * np.log(wavelen)) / np.sum(fnu * bp.phi))
                leff2   = np.exp(np.sum(flambda * bp.phi * np.log(wavelenf)) / np.sum(flambda * bp.phi))

                leff1g  = np.exp(np.sum(fnug * bp.phi * np.log(waveleng)) / np.sum(fnug * bp.phi))
                leff2g  = np.exp(np.sum(flambdag * bp.phi * np.log(wavelenfg)) / np.sum(flambdag * bp.phi))

                off1    = refract(leff1*10**-3,    zd * np.pi / 180.) * 180. / np.pi * 3600.
                off2    = refract(leff2*10**-3,    zd * np.pi / 180.) * 180. / np.pi * 3600.
                off     = refract(wavelen*10**-3,  zd * np.pi / 180.) * 180. / np.pi * 3600.
                off     = np.sum(off * flux) / np.sum(flux)

                off1g   = refract(leff1g*10**-3,   zd * np.pi / 180.) * 180. / np.pi * 3600.
                off2g   = refract(leff2g*10**-3,   zd * np.pi / 180.) * 180. / np.pi * 3600.
                offg    = refract(waveleng*10**-3, zd * np.pi / 180.) * 180. / np.pi * 3600.
                offg    = np.sum(offg * fluxg) / np.sum(fluxg)

                print "%s %d %d : %+.4f %+.4f %+.4f" % (
                    bpname, sidx, zd, off1-off1g, off2-off2g, off-offg)
        print
