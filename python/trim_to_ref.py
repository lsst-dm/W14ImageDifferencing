import sys
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom

# Come from colors.py
colors_vs_r = {}
for sed in ("kp01_9750.fits_g45_9830.gz",
            "km20_6000.fits_g30_6020.gz",
            "m2.0Full.dat.gz"):
    colors_vs_r[sed] = {}
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["u"] =  0.60680054126
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["g"] = -0.23920305940
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["i"] =  0.23364024549
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["z"] =  0.41330200452
colors_vs_r["kp01_9750.fits_g45_9830.gz"]["y"] =  0.45125140946

colors_vs_r["km20_6000.fits_g30_6020.gz"]["u"] =  1.21602658889
colors_vs_r["km20_6000.fits_g30_6020.gz"]["g"] =  0.30637498242
colors_vs_r["km20_6000.fits_g30_6020.gz"]["i"] = -0.10157478397
colors_vs_r["km20_6000.fits_g30_6020.gz"]["z"] = -0.11823111346
colors_vs_r["km20_6000.fits_g30_6020.gz"]["y"] = -0.11561614594

colors_vs_r["m2.0Full.dat.gz"]["u"] =  3.6508634144
colors_vs_r["m2.0Full.dat.gz"]["g"] =  1.3542089873
colors_vs_r["m2.0Full.dat.gz"]["i"] = -1.0173818086
colors_vs_r["m2.0Full.dat.gz"]["z"] = -1.4889633446 
colors_vs_r["m2.0Full.dat.gz"]["y"] = -1.6751067151

print "# id,dbid,raJ2000,decJ2000,glat,glon,uDerived,gDerived,rDerived,iDerived,zDerived,yDerived,isStar,varClass,objClass,redshift,semiMajorBulge,semiMinorBulge,semiMajorDisk,semiMinorDisk,sedFilename,properMotionRa,properMotionDec,parallax,radialVelocity"

for line in open(sys.argv[1]).readlines():
    if not line.startswith("object"):
        continue

    obj, oid, ra, decl, rmag, sed, a, b, c, d, e, f, star, g, h = line.split()
    c    = afwCoord.makeCoord(afwCoord.FK5, float(ra) * afwGeom.degrees, float(decl) * afwGeom.degrees)
    glat = c.toGalactic()[0].asDegrees()
    glon = c.toGalactic()[1].asDegrees()

    sedKey = sed.split("/")[-1]
    umag = float(rmag) + colors_vs_r[sedKey]["u"]
    gmag = float(rmag) + colors_vs_r[sedKey]["g"]
    imag = float(rmag) + colors_vs_r[sedKey]["i"]
    zmag = float(rmag) + colors_vs_r[sedKey]["z"]
    ymag = float(rmag) + colors_vs_r[sedKey]["y"]

    print ",".join(map(str, (oid, oid, ra, decl, glat, glon,
                             umag, gmag, rmag, imag, zmag, ymag,
                             1, 0, 5, None, None, None, None, None,
                             sed, 0., 0., 0., 0.)))

