import scipy
import lsst.sims.catalogs.measures.example_utils as utils
import numpy as np
import sys, re


# Come from colors.py
colors_vs_r = {}
for sed in ("kp01_9750.fits_g45_9830.gz",
            "km20_6000.fits_g30_6020.gz",
            "m2.0Full.dat.gz"):
    colors_vs_r[sed] = {}
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["u"] =  0.60680054126
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["g"] = -0.23920305940
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["r"] =  0.0
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["i"] =  0.23364024549
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["z"] =  0.41330200452
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["y"] =  0.45125140946

colors_vs_r["km20_6000.fits_g30_6020.gz"]["u"] =  1.21602658889
colors_vs_r["km20_6000.fits_g30_6020.gz"]["g"] =  0.30637498242
colors_vs_r["km20_6000.fits_g30_6020.gz"]["r"] =  0.0
colors_vs_r["km20_6000.fits_g30_6020.gz"]["i"] = -0.10157478397
colors_vs_r["km20_6000.fits_g30_6020.gz"]["z"] = -0.11823111346
colors_vs_r["km20_6000.fits_g30_6020.gz"]["y"] = -0.11561614594

colors_vs_r["m2.0Full.dat.gz"]["u"] =  3.6508634144
colors_vs_r["m2.0Full.dat.gz"]["g"] =  1.3542089873
colors_vs_r["m2.0Full.dat.gz"]["r"] =  0.0
colors_vs_r["m2.0Full.dat.gz"]["i"] = -1.0173818086
colors_vs_r["m2.0Full.dat.gz"]["z"] = -1.4889633446 
colors_vs_r["m2.0Full.dat.gz"]["y"] = -1.6751067151


fMapper = {"g": 1, "r": 2, "i": 3}

raRad     = 79.6892650466 * np.pi / 180
#decRad    = -9.70229566883 * np.pi / 180 # Original data 
decRad    = -0.517781017                  # Putting the star field through zenith
rotSkyRad = 251.068956606 * np.pi / 180   # Put stars in same orientation on the chip
mjd0      = 51130.2617188                 # Original MJD from trim files
#mjds      = [51130.0783855, 51130.1443577, 51130.2728299, 51130.4006077, 51130.4665799] # these put things at airmass 1.3, 2.0
mjds      = [51130.1117188, 51130.1763021, 51130.2728299, 51130.3686632, 51130.4332466]

if 0:
    for minute in np.arange(-360, 360, 1.0):
        mjd    = mjd0 + minute / 24. / 60.
        params = utils.makeObsParamsRaDecSky(raRad, decRad, mjd, "i", rotSkyRad=rotSkyRad)
        print mjd, params["Unrefracted_Altitude"][0] * 180 / np.pi
    sys.exit(1)



utilsKeys = ['Opsim_azimuth', 'Opsim_moondec', 'Opsim_rottelpos', 'Opsim_moonalt', 'Opsim_altitude', 'Opsim_rotskypos', 'Opsim_moonra', 'Opsim_sunalt', 'Unrefracted_RA', 'Opsim_dist2moon', 'Opsim_filter', 'Unrefracted_Dec', "Opsim_obshistid", "Opsim_expmjd", "Unrefracted_Altitude", "Unrefracted_Azimuth"]

ddecl     = -19.96437131727
trimfile  = sys.argv[1]
newlines  = []
header    = []
for line in open(trimfile).readlines():
    if not line.startswith("object"):
        fields = line.split()
        if not fields[0] in utilsKeys:
            header.append(line)
        if fields[0] == "Opsim_obshistid":
            oid = fields[1][-1]
            if oid == "6":
                oid = "0"
        continue

    fields = line.split()
    decl   = fields[3]
    newlines.append(re.sub(decl, str(float(decl)+ddecl), line))
    #print line,



for imjd, mjd in enumerate(mjds):
    for ifilt, filt in enumerate(("g", "r", "i")):
        obshistid = "%d00%d00%s" % (ifilt+1, imjd, oid)
        outfile = "%s.trim" % (obshistid)
        buff = open(outfile, "w")
        buff.write("%s %s\n" % ("Opsim_obshistid", obshistid))
        buff.write("%s %f\n" % ("Opsim_expmjd", mjd))
        for h in header:
            buff.write(h)
        params = utils.makeObsParamsRaDecSky(raRad, decRad, mjd, filt, rotSkyRad=rotSkyRad)
        for key, value in params.iteritems():
            if key == "Opsim_filter":
                value = fMapper[value[0]]
            else:
                value = value[0] * 180. / np.pi  # turn from rad into degrees
            buff.write("%s %s\n" % (key, value))
        for line in newlines:
            fields = line.split()
            sedKey = fields[5].split("/")[-1]
            rmag   = fields[4]
            line   = re.sub(rmag, str(float(rmag)-colors_vs_r[sedKey][filt]), line)
            buff.write(line)
        buff.close()
        #print obshistid, params


