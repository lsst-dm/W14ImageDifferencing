import os
import numpy as np
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
fluxFormatter = FormatStrFormatter('%d')
dcrFormatter = FormatStrFormatter('%+0.2f')

def dcr(wavelength, zd, docorr=True, P=600.0, T=7.0, f=8.0):
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

bps = []
filtDir = os.environ["LSST_THROUGHPUTS_BASELINE"]
gBand  = Bandpass()
rBand  = Bandpass()
iBand  = Bandpass()
gBand.readThroughput(os.path.join(filtDir, "total_g.dat"))
rBand.readThroughput(os.path.join(filtDir, "total_r.dat"))
iBand.readThroughput(os.path.join(filtDir, "total_i.dat"))
bps.append(gBand)
bps.append(rBand)
bps.append(iBand)
bNames = ["g-band", "r-band", "i-band"]
bColors = ["g", "r", "m"]


# Note on file names.  
# "k" = kurucz
# "m/p" = minus or plus metallicity
# _temp_ = effective temperature bin
# _g50_ = log g
# _temp.gz = actual temperature
#
# G-star:
# http://www.stsci.edu/hst/observatory/cdbs/k93models.html
# G0V       6030     +4.39       kp00_6000[g45] 
#
# Lets find the most common G0V-star used in the sims:
# grep 6000 /tmp/Star_SEDS.dat | grep -v berg | awk '{split($0, f, ","); print f[1],f[2] }' | sort -n -k 2 | tail
# > km20_6000.fits_g30_6020 1246
#
# Most common M-dwarf
# grep Full /tmp/Star_SEDS.dat | grep -v berg | awk '{split($0, f, ","); print f[1],f[2] }' | sort -n -k 2 | tail
# > m2.0Full.dat 37498
# 
# Hottest temp in sims:
# grep -v berg /tmp/Star_SEDS.dat | grep -v dat | awk '{split($0, f, "_"); print f[2] }' | awk '{split($0, f, "."); print f[1] }' | sort -n -k 1 | tail
# > 9750
# This is like an A0V star
# A0V       9520     +4.14       kp00_9500[g40] 
#
# kp01_9750.fits_g45_9830

# Most common hot star (not used)
# O5V      44500     +4.04      kp00_45000[g50] 

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
for redshift in (0, 0.5):
    agn    = Sed()
    agn.readSED_flambda(os.path.join(catDir, "data/agnSED/agn.spec.gz"))
    agn.redshiftSED(redshift)
    seds.append(agn)
sNames = ["A0V Star", "G0V Star", "M2.0V Star", "AGN(z=0)", "AGN(z=0.5)"]

fig1 = plt.figure(figsize=(12,6))
fig2 = plt.figure(figsize=(12,6))
fig3 = plt.figure(figsize=(12,6))
fig4 = plt.figure(figsize=(12,6))
gs   = gridspec.GridSpec(4, 6)
sp01 = fig1.add_subplot(gs[3, 0])
sp02 = fig2.add_subplot(gs[3, 0])
spz2 = fig2.add_subplot(gs[3, 5])
sp03 = fig3.add_subplot(gs[3, 0])
spz3 = fig3.add_subplot(gs[3, 5])
sp04 = fig4.add_subplot(gs[3, 0])
spz4 = fig4.add_subplot(gs[3, 5])
spz2.yaxis.tick_right()
spz2.yaxis.set_label_position("right")
spz3.yaxis.tick_right()
spz3.yaxis.set_label_position("right")
spz4.yaxis.tick_right()
spz4.yaxis.set_label_position("right")

for sidx, sed in enumerate(seds):
    sed.flambdaTofnu() 

    ax = fig1.add_subplot(gs[0, sidx+1], sharex=sp01)
    ax.plot(sed.wavelen, sed.fnu, lw=2)
    plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.set_title(sNames[sidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )
    if sidx == 0:
        ax.set_ylabel(r"$f_\nu(\lambda)$", weight="bold", fontsize=16)

    ax = fig2.add_subplot(gs[0, sidx+1], sharex=sp02)
    ax.plot(sed.wavelen, sed.fnu, lw=2)
    plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.set_title(sNames[sidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )

    ax = fig3.add_subplot(gs[0, sidx+1], sharex=sp02)
    ax.plot(sed.wavelen, sed.fnu, lw=2)
    plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.set_title(sNames[sidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )

    ax = fig4.add_subplot(gs[0, sidx+1], sharex=sp02)
    ax.plot(sed.wavelen, sed.fnu, lw=2)
    plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.set_title(sNames[sidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )

for bidx, bp in enumerate(bps):
    bp.sbTophi()

    if bidx == 3:
        ax = sp01
    else:
        ax = fig1.add_subplot(gs[bidx+1, 0], sharex=sp01)
        plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.plot(bp.wavelen, bp.phi, "%s-" % (bColors[bidx]), lw=2)
    ax.set_ylabel(bNames[bidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )

    if bidx == 3:
        ax = sp02
    else:
        ax = fig2.add_subplot(gs[bidx+1, 0], sharex=sp02)
        plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.plot(bp.wavelen, bp.phi, "%s-" % (bColors[bidx]), lw=2)
    ax.set_ylabel(bNames[bidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )

    if bidx == 3:
        ax = sp02
    else:
        ax = fig3.add_subplot(gs[bidx+1, 0], sharex=sp02)
        plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.plot(bp.wavelen, bp.phi, "%s-" % (bColors[bidx]), lw=2)
    ax.set_ylabel(bNames[bidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )

    if bidx == 3:
        ax = sp02
    else:
        ax = fig4.add_subplot(gs[bidx+1, 0], sharex=sp02)
        plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
    ax.plot(bp.wavelen, bp.phi, "%s-" % (bColors[bidx]), lw=2)
    ax.set_ylabel(bNames[bidx], weight="bold", fontsize=16)
    ax.yaxis.set_major_formatter( fluxFormatter )


for bidx, bp in enumerate(bps):
    for sidx, sed in enumerate(seds):
        wavelen = sed.wavelen
        fnu = sed.fnu
        if sed.needResample(wavelen=wavelen, wavelen_match=bp.wavelen):
            wavelen, fnu = sed.resampleSED(wavelen, fnu, wavelen_match=bp.wavelen)
        flux = fnu * bp.phi
        #print "FLUX", bidx, sidx, np.sum(flux)

        ax = fig1.add_subplot(gs[bidx+1, sidx+1], sharex=sp01)
        ax.plot(wavelen, flux, "k-", lw=2)
        plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
        ax.yaxis.set_major_formatter( fluxFormatter )

        # Print out effective wavelengths
        leff = np.exp(np.sum(fnu * bp.phi * np.log(wavelen)) / np.sum(fnu * bp.phi))
        print "Effective wavelength bp %d, sed %d : %.2f -> %.2f" % (bidx, sidx, bp.calcEffWavelen()[0]*10, leff*10)

        # Absolute offsets
        offsets = []
        zds = (0, 10, 20, 30, 40, 50)
        for zd in zds:
            offset = dcr(wavelen*10**-3, zd * np.pi / 180.) * 180. / np.pi * 3600.
            offsets.append(np.sum(offset * flux) / np.sum(flux))
        if bidx == 2 and sidx == 4:
            ax = spz2
        else:
            ax = fig2.add_subplot(gs[bidx+1, sidx+1], sharey=spz2)
            plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
        ax.plot(zds, offsets, "k-", lw=2)
        # Without STP correction
        offsets = []
        for zd in zds:
            offset = dcr(wavelen*10**-3, zd * np.pi / 180., docorr=False) * 180. / np.pi * 3600.
            offsets.append(np.sum(offset * flux) / np.sum(flux))
        ax.plot(zds, offsets, "r--", lw=2)

        

        # Relative to g-band
        sedg = seds[1]
        waveleng = sedg.wavelen
        fnug = sedg.fnu
        if sedg.needResample(wavelen=waveleng, wavelen_match=bp.wavelen):
            waveleng, fnug = sed.resampleSED(waveleng, fnug, wavelen_match=bp.wavelen)
        fluxg = fnug * bp.phi
        offsets = []
        for zd in zds:
            offset  = dcr(wavelen*10**-3, zd * np.pi / 180.) * 180. / np.pi * 3600.
            offsetg = dcr(waveleng*10**-3, zd * np.pi / 180.) * 180. / np.pi * 3600.
            offsets.append(np.sum(offset * flux) / np.sum(flux) - np.sum(offsetg * fluxg) / np.sum(fluxg))
        if bidx == 2 and sidx == 4:
            ax = spz3
        else:
            ax = fig3.add_subplot(gs[bidx+1, sidx+1], sharey=spz3)
            plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
        ax.plot(zds, offsets, "k-", lw=2)
        ax.axhline(y=0, c='k', linestyle=':')
        # Without STP correction
        offsets = []
        for zd in zds:
            offset  = dcr(wavelen*10**-3, zd * np.pi / 180., docorr=False) * 180. / np.pi * 3600.
            offsetg = dcr(waveleng*10**-3, zd * np.pi / 180., docorr=False) * 180. / np.pi * 3600.
            offsets.append(np.sum(offset * flux) / np.sum(flux) - np.sum(offsetg * fluxg) / np.sum(fluxg))
        ax.plot(zds, offsets, "r--", lw=2)

        # LLC relative to URC
        offsets = []
        zds1 = np.array(zds)
        zds2 = zds1 + np.sqrt(2) * 1024 * 4 * 0.2 / 3600.
        for i in range(len(zds)):
            offset1 = dcr(wavelen*10**-3, zds1[i] * np.pi / 180.) * 180. / np.pi * 3600.
            offset2 = dcr(wavelen*10**-3, zds2[i] * np.pi / 180.) * 180. / np.pi * 3600.
            o1      = np.sum(offset1 * flux) / np.sum(flux)
            o2      = np.sum(offset2 * flux) / np.sum(flux)
            offsets.append(o1 - o2)
        if bidx == 2 and sidx == 4:
            ax = spz4
        else:
            ax = fig4.add_subplot(gs[bidx+1, sidx+1], sharey=spz4)
            plt.setp(ax.get_xticklabels()+ax.get_yticklabels(), visible=False)
        ax.plot(zds, offsets, "k-", lw=2)
        ax.axhline(y=0, c='k', linestyle='--')
        # Without STP correction
        offsets = []
        for i in range(len(zds)):
            offset1 = dcr(wavelen*10**-3, zds1[i] * np.pi / 180., docorr=False) * 180. / np.pi * 3600.
            offset2 = dcr(wavelen*10**-3, zds2[i] * np.pi / 180., docorr=False) * 180. / np.pi * 3600.
            o1      = np.sum(offset1 * flux) / np.sum(flux)
            o2      = np.sum(offset2 * flux) / np.sum(flux)
            offsets.append(o1 - o2)
        ax.plot(zds, offsets, "r--", lw=2)

sp01.set_xlim(350, 900)
plt.setp(sp01.get_yticklabels(), visible=False)
plt.setp(sp01.get_xticklabels(), weight="bold", fontsize=11, rotation=45)
sp01.set_xlabel("Wavelength (nm)", weight="bold", fontsize=13)

sp02.set_xlim(350, 900)
plt.setp(sp02.get_yticklabels(), visible=False)
plt.setp(sp02.get_xticklabels(), weight="bold", fontsize=11, rotation=45)
sp02.set_xlabel("Wavelength (nm)", weight="bold", fontsize=13)
spz2.set_xlim(0, 50)
plt.setp(spz2.get_xticklabels(), weight="bold", fontsize=11) #, rotation=45)
plt.setp(spz2.get_yticklabels(), weight="bold", fontsize=11)
spz2.set_xlabel("ZD (deg)", weight="bold", fontsize=13)
spz2.set_ylabel(r'$\Delta x$"', weight="bold", fontsize=15, rotation=360)
fig2.text(0.075, 0.825, r"$\frac{\int\ R(\lambda) f(\lambda) d\lambda}{\int\ f(\lambda) d\lambda}$",
          weight="bold", fontsize=25)

sp03.set_xlim(350, 900)
plt.setp(sp03.get_yticklabels(), visible=False)
plt.setp(sp03.get_xticklabels(), weight="bold", fontsize=11, rotation=45)
sp03.set_xlabel("Wavelength (nm)", weight="bold", fontsize=13)
spz3.set_xlim(0, 50)
plt.setp(spz3.get_xticklabels(), weight="bold", fontsize=11) #, rotation=45)
plt.setp(spz3.get_yticklabels(), weight="bold", fontsize=11)
spz3.set_xlabel("ZD (deg)", weight="bold", fontsize=13)
spz3.set_ylabel(r'$\Delta x $"', weight="bold", fontsize=15, rotation=360)
fig3.text(0.01, 0.825, r"$\frac{\int\ R(\lambda) f(\lambda) d\lambda}{\int\ f(\lambda) d\lambda}\ -\ \frac{\int\ R(\lambda) f(\lambda)_{G0V}\ d\lambda}{\int\ f(\lambda)_{G0V}\ d\lambda}$",
          weight="bold", fontsize=19)
spz3.set_ylim(-0.15, 0.15)

sp04.set_xlim(350, 900)
plt.setp(sp04.get_yticklabels(), visible=False)
plt.setp(sp04.get_xticklabels(), weight="bold", fontsize=11, rotation=45)
sp04.set_xlabel("Wavelength (nm)", weight="bold", fontsize=13)
spz4.set_xlim(0, 50)
plt.setp(spz4.get_xticklabels(), weight="bold", fontsize=11) #, rotation=45)
plt.setp(spz4.get_yticklabels(), weight="bold", fontsize=11)
spz4.set_xlabel("ZD (deg)", weight="bold", fontsize=13)
spz4.set_ylabel(r'$\Delta x $"', weight="bold", fontsize=15, rotation=360)
fig4.text(0.01, 0.825, r"$\frac{\int\ R(\lambda)_{LLC} f(\lambda) d\lambda}{\int\ f(\lambda) d\lambda}\ -\ \frac{\int\ R(\lambda)_{RHC} f(\lambda)\ d\lambda}{\int\ f(\lambda)\ d\lambda}$",
          weight="bold", fontsize=17)
spz4.set_ylim(-1.01, .11)

fig1.subplots_adjust(hspace=0.025, wspace=0.025)
fig2.subplots_adjust(hspace=0.025, wspace=0.025)
fig3.subplots_adjust(hspace=0.025, wspace=0.025)
fig4.subplots_adjust(hspace=0.025, wspace=0.025)
plt.show()


######## NOTE: make these bigger.  also calculate the variation of this effect across the focal plane and a chip
