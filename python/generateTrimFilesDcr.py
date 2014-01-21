import scipy
import lsst.sims.catalogs.measures.example_utils as utils
import numpy as np
import sys, re

fMapper = {"g": 1, "r": 2, "i": 3}

raRad     = 79.6892650466 * np.pi / 180
#decRad    = -9.70229566883 * np.pi / 180 # Original data 
decRad    = -0.517781017                  # Putting the star field through zenith
rotSkyRad = 251.068956606 * np.pi / 180   # Put stars in same orientation on the chip
mjd0      = 51130.2617188                 # Original MJD from trim files
mjds      = [51130.0783855, 51130.1443577, 51130.2728299, 51130.4006077, 51130.4665799]

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
            buff.write(line)
        buff.close()
        #print obshistid, params


if 0:
    for minute in np.arange(-360, 360, 1.0):
        mjd    = mjd0 + minute / 24. / 60.
        params = utils.makeObsParamsRaDecSky(raRad, decRad, mjd, "i", rotSkyRad=rotSkyRad)
        print mjd, params["Opsim_altitude"][0] * 180 / np.pi
