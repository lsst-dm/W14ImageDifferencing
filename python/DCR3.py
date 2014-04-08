import os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
fluxFormatter = FormatStrFormatter('%d')
dcrFormatter = FormatStrFormatter('%+0.2f')
majorLocator   = MultipleLocator(3)

# THIS ONE LOOKS AT THE AIRMASS-DEPENDENT FILTER PROFILE
# And makes a heat map of delta r in airmass vs. delta_theta

def getOffset(wavelen, flux, zd):
    off = refract(wavelen*10**-3,  zd * np.pi / 180.) * 180. / np.pi * 3600.
    return np.sum(off * flux) / np.sum(flux)

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
    aStar  = Sed() 
    gStar  = Sed() 
    mStar  = Sed() 
    aStar.readSED_flambda(os.path.join(catDir, "data/starSED/kurucz", "kp01_9750.fits_g45_9830.gz"))
    gStar.readSED_flambda(os.path.join(catDir, "data/starSED/kurucz", "km20_6000.fits_g30_6020.gz"))
    mStar.readSED_flambda(os.path.join(catDir, "data/starSED/mlt", "m2.0Full.dat"))
    seds.append(aStar)
    seds.append(gStar)
    seds.append(mStar)

    for sed in seds:
        sed.flambdaTofnu() 

    filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
    for bpname in ("g", "r", "i"):
        bp = Bandpass(wavelen_max=1200) 
        components = common_components + [os.path.join(tp_dir, "baseline", "filter_%s.dat" % (bpname)), 
                                          os.path.join(tp_dir, "atmos", "atmos_%d.dat" % (12))]
        bp.readThroughputList(components)
        bp.sbTophi()

        wavelenA            = aStar.wavelen
        fnuA                = aStar.fnu
        wavelenA, fnuA      = aStar.resampleSED(wavelenA, fnuA, wavelen_match=bp.wavelen)

        wavelenG            = gStar.wavelen
        fnuG                = gStar.fnu
        wavelenG, fnuG      = gStar.resampleSED(wavelenG, fnuG, wavelen_match=bp.wavelen)

        wavelenM            = mStar.wavelen
        fnuM                = mStar.fnu
        wavelenM, fnuM      = mStar.resampleSED(wavelenM, fnuM, wavelen_match=bp.wavelen)

        fluxA   = fnuA * bp.phi
        fluxG   = fnuG * bp.phi
        fluxM   = fnuM * bp.phi

        airmass1   = 1.25
        dairmasses = np.arange(-0.25, 0.26, 0.05)
        #airmass1   = 1.00
        #dairmasses = np.arange(0., 0.51, 0.05)
        #airmass1   = 1.50
        #dairmasses = np.arange(-0.5, 0.1, 0.05)

        dthetas    = np.arange(0, 181, 5)
        drA        = np.empty((len(dairmasses)-1, len(dthetas)-1))
        drM       =  np.empty((len(dairmasses)-1, len(dthetas)-1))

        zd1       = np.arctan(airmass1) * 180 / np.pi
        dA1       = getOffset(wavelenA, fluxA, zd1) - getOffset(wavelenG, fluxG, zd1) 
        dM1       = getOffset(wavelenM, fluxM, zd1) - getOffset(wavelenG, fluxG, zd1) 
        for ia, dairmass in enumerate(dairmasses):
            zd2       = np.arctan(airmass1 + dairmass) * 180 / np.pi
            dA2       = getOffset(wavelenA, fluxA, zd2) - getOffset(wavelenG, fluxG, zd2) 
            dM2       = getOffset(wavelenM, fluxM, zd2) - getOffset(wavelenG, fluxG, zd2) 

            for it, dtheta in enumerate(dthetas):
                try:
                    drA[ia, it] = dA1**2 + dA2**2 - 2 * dA1 * dA2 * np.cos(dtheta * np.pi / 180)
                    drM[ia, it] = dM1**2 + dM2**2 - 2 * dM1 * dM2 * np.cos(dtheta * np.pi / 180)
                except:
                    pass

        drA = np.log10(drA)
        drM = np.log10(drM)

        fig = plt.figure()
        sp1 = fig.add_subplot(211)
        sp2 = fig.add_subplot(212)
        hmap1 = sp1.pcolor(drA, cmap=plt.cm.Greys, shading="faceted", vmin=-4, vmax=-1)
        sp1.set_xticklabels(dthetas[::2], minor=False, weight="bold", fontsize=12)
        sp1.set_yticklabels(airmass1+dairmasses[::2], minor=False, weight="bold", fontsize=12)
        sp1.set_title("A star, %s-band" % (bpname), weight="bold", fontsize=14)
        sp1.contour(drA, levels=(np.log10(0.00114), np.log10(0.00204)), colors=("b", "r"), linewidths=(3,), linestyles=("solid",))
        hmap2 = sp2.pcolor(drM, cmap=plt.cm.Greys, shading="faceted", vmin=-4, vmax=-1)
        sp2.xaxis.set_major_locator(majorLocator)
        sp2.set_xticklabels(dthetas[::3], minor=False, weight="bold", fontsize=12)
        sp2.set_yticklabels(airmass1+dairmasses[::2], minor=False, weight="bold", fontsize=12)
        sp2.set_title("M star, %s-band" % (bpname), weight="bold", fontsize=14)
        sp2.contour(drM, levels=(np.log10(0.00114), np.log10(0.00204)), colors=("b", "r"), linewidths=(3,), linestyles=("solid",))
        cb1 = plt.colorbar(hmap1, ax=sp1, fraction=0.1)
        cb1.set_label("Log10 Offset('')", weight="bold", fontsize=10)
        cb2 = plt.colorbar(hmap2, ax=sp2, fraction=0.1)
        cb2.set_label("Log10 Offset('')", weight="bold", fontsize=10)
        sp1.set_ylabel("Airmass", weight="bold", fontsize=13)
        sp2.set_ylabel("Airmass", weight="bold", fontsize=13)
        plt.setp(sp1.get_xticklabels(), visible=False)
        sp2.set_xlabel("Delta Theta (degrees)", weight="bold", fontsize=13)

    plt.show()
    import pdb; pdb.set_trace()
