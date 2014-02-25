import os
import numpy as np
from lsst.obs.lsstSim import LsstSimMapper as Mapper 
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.meas.astrom as measAstrom
import lsst.afw.display.ds9 as ds9
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom

def getDcrAmp(zenith, filterName, gr):
    # Derived from DCR.py
    dcr50 = {}
    dcr50["g"] = {}
    dcr50["r"] = {}
    dcr50["i"] = {}
    dcr50["g"]["A"] = 55.74
    dcr50["r"]["A"] = 55.69
    dcr50["i"]["A"] = 55.58
    dcr50["g"]["G"] = 55.12
    dcr50["r"]["G"] = 55.10
    dcr50["i"]["G"] = 55.08
    dcr50["g"]["M"] = 54.83
    dcr50["r"]["M"] = 54.83
    dcr50["i"]["M"] = 54.81
    if zenith < 1:
        return 0.0

    if gr < 0: star="A"
    elif gr < 1: star="G"
    else: star="M"

    return dcr50[filterName][star] / 0.2 # arcseconds to pixels

# Airmass of the template image
suffixT = "C"
airmassIdT = 2
astrometer = measAstrom.Astrometry(measAstrom.MeasAstromConfig())
indir = "/nfs/lsst/home/becker/Winter2014"

# long, lat (long is +ve E of Greenwich), alt (I'm assuming its meters?
lsst = afwCoord.Observatory(-70.7494167 * afwGeom.degrees, -30.2446389 * afwGeom.degrees, 2662)

for filterName, filterId in zip(("g", "r", "i"), (1, 2, 3)):
    template = "%d00%d000" % (filterId, airmassIdT)
    
    # Airmass of the science image
    for suffixI, airmassIdI, frameId, mjdI in zip(("A", "C", "E"), (0, 2, 4), (1, 2, 3), (51130.11171, 51130.27283, 51130.43324)):
        visit   = int("%d00%d002" % (filterId, airmassIdI))

        mapper  = Mapper(root=os.path.join(indir, "outputs8%s%s_doPreConvolveFalse" % (suffixT, suffixI)), calibRoot = None, outputRoot = None)
        butler  = dafPersist.ButlerFactory(mapper = mapper).create()
        dataId  = {"filter": filterName, "visit": visit, "raft": "2,2", "sensor": "1,1"}
        calexp  = butler.get(datasetType="calexp", dataId=dataId)
        rotTel  = calexp.getMetadata().get("SPIDANG")
        zenith  = calexp.getMetadata().get("ZENITH")
        diffim  = butler.get(datasetType="deepDiff_differenceExp", dataId=dataId)
        src     = butler.get(datasetType="src", dataId=dataId)
        matches = astrometer.useKnownWcs(src, exposure=calexp).matches

        if filterId == 1 and airmassIdI == 0:
            # Make sure ds9 has all the frames initialized
            for i in range(1, 10): ds9.mtv(calexp, frame=i)

        ds9.mtv(calexp, frame=frameId)
        ds9.ds9Cmd("mask transparency 80")
        if suffixI != "C":
            for match in matches:
                refObj = match.first
                gr     = -2.5 * np.log10(refObj["g"]) - -2.5 * np.log10(refObj["r"])
                amp    = getDcrAmp(zenith, filterName, gr)
                if amp == 0:
                    continue
                if gr < 0: color="blue"
                elif gr < 1: color="green"
                else: color="red"
                src    = match.second
                x1, y1 = src.getCentroid() # End of the Dcr vector
                dx     = amp * -np.sin(np.pi * -rotTel / 180.) # Minus since we are rotating clockwise from y-axis
                dy     = amp *  np.cos(np.pi * -rotTel / 180.) 
                x0     = x1 - dx
                y0     = y1 - dy
                cmd    = "regions command {vector %f %f %f %f # vector=1 color=%s width=3}" % (x0, y0, amp, -rotTel+90, color)
                ds9.ds9Cmd(cmd, silent=False)
    
        ds9.mtv(calexp, frame=frameId+3)
        ds9.ds9Cmd("mask transparency 80")
        if suffixI != "C":
            for match in matches:
                refObj = match.first
                gr     = -2.5 * np.log10(refObj["g"]) - -2.5 * np.log10(refObj["r"])
                amp    = getDcrAmp(zenith, filterName, gr)
                gamp   = getDcrAmp(zenith, filterName, 0.306) # Refernce star km20_6000.fits_g30_6020.gz
                amp   -= gamp
                if amp == 0:
                    continue
                amp   *= 10.0
                if gr < 0: color="blue"
                elif gr < 1: color="green"
                else: color="red"
                src    = match.second
                x1, y1 = src.getCentroid() # End of the Dcr vector
                dx     = amp * -np.sin(np.pi * -rotTel / 180.)
                dy     = amp *  np.cos(np.pi * -rotTel / 180.)
                x0     = x1 - dx
                y0     = y1 - dy
                if amp < 0: ang = -rotTel - 90
                else:       ang = -rotTel + 90
                cmd    = "regions command {vector %f %f %f %f # vector=1 color=%s width=3}" % (x1, y1, amp, ang, color)
                ds9.ds9Cmd(cmd, silent=False)

        ds9.mtv(diffim, frame=frameId+6)
        ds9.ds9Cmd("mask transparency 80")

        # Skip the at zenith one
        if suffixI != "C":
            # Epoch of image
            obsDate  = dafBase.DateTime(mjdI, dafBase.DateTime.MJD, dafBase.DateTime.TAI)
            epoch    = obsDate.get(dafBase.DateTime.EPOCH)
            for match in matches:
                refObj = match.first

                #ra       = 79.6892650466
                #decl     = -29.6666669861
                ra       = refObj.getCoord()[0].asDegrees() + 0.003 # Offset a bit
                decl     = refObj.getCoord()[1].asDegrees() + 0.003
                wcs      = diffim.getWcs()
                center   = afwCoord.makeCoord(afwCoord.FK5, ra*afwGeom.degrees, decl*afwGeom.degrees, epoch)
                xyc      = wcs.skyToPixel(center)
                raoff    = afwCoord.makeCoord(afwCoord.FK5, (ra+0.003)*afwGeom.degrees, decl*afwGeom.degrees, epoch)
                xydra    = wcs.skyToPixel(raoff)
                decoff   = afwCoord.makeCoord(afwCoord.FK5, ra*afwGeom.degrees, (decl+0.003)*afwGeom.degrees, epoch)
                xyddec   = wcs.skyToPixel(decoff)
                cmd      = "regions command {line %f %f %f %f # color=magenta width=3}" % (xyc[0], xyc[1], xydra[0], xydra[1])  # East
                ds9.ds9Cmd(cmd, silent=False)
                cmd      = "regions command {line %f %f %f %f # color=yellow width=3}" % (xyc[0], xyc[1], xyddec[0], xyddec[1]) # North
                ds9.ds9Cmd(cmd, silent=False)
                
                centeraa = center.toTopocentric(lsst, obsDate) # az, alt
                azoff    = afwCoord.TopocentricCoord((centeraa.getAzimuth().asDegrees() + 0.003)*afwGeom.degrees,
                                                     centeraa.getAltitude().asDegrees()*afwGeom.degrees, epoch, lsst)
                altoff   = afwCoord.TopocentricCoord(centeraa.getAzimuth().asDegrees()*afwGeom.degrees,
                                                     (centeraa.getAltitude().asDegrees() + 0.003)*afwGeom.degrees, epoch, lsst)
        
                off1   = azoff.toFk5(epoch)
                xy1    = wcs.skyToPixel(off1)
                off2   = altoff.toFk5(epoch)
                xy2    = wcs.skyToPixel(off2)
                cmd      = "regions command {line %f %f %f %f # color=green width=3}" % (xyc[0], xyc[1], xy1[0], xy1[1]) # positive Az
                ds9.ds9Cmd(cmd, silent=False)
                cmd      = "regions command {line %f %f %f %f # color=cyan width=3}" % (xyc[0], xyc[1], xy2[0], xy2[1])  # positive Alt
                ds9.ds9Cmd(cmd, silent=False)
    
            # Angle of DCR!
            dx     = xyc[0] - xy2[0]
            dy     = xyc[1] - xy2[1]
            print "ANGLE TO ZENITH", dx, dy
        
    import pdb; pdb.set_trace()
    
    
    
            
