import sys, os
import matplotlib.pyplot as plt
import numpy as np
import lsst.pipe.base as pipeBase
from lsst.obs.lsstSim import LsstSimMapper as Mapper 
import lsst.daf.persistence as dafPersist
import lsst.afw.display.ds9 as ds9

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please
pexLog.Trace_setVerbosity("", -100) # No logging info please

def checkDipoles(diaSrc):
    dxdyp = []
    dxdyn = []
    for src in diaSrc:
        centPosNaive     = src.get("centroid.dipole.naive.pos")
        centPosNaiveErr  = src.get("centroid.dipole.naive.pos.err")
        centPosNaiveFlag = src.get("centroid.dipole.naive.pos.flags")
        fluxPosNaive     = src.get("flux.dipole.naive.pos")
        fluxPosNaiveErr  = src.get("flux.dipole.naive.pos.err")
        fluxPosNaiveFlag = src.get("flux.dipole.naive.pos.flags")

        centNegNaive     = src.get("centroid.dipole.naive.neg")
        centNegNaiveErr  = src.get("centroid.dipole.naive.neg.err")
        centNegNaiveFlag = src.get("centroid.dipole.naive.neg.flags")
        fluxNegNaive     = src.get("flux.dipole.naive.neg")
        fluxNegNaiveErr  = src.get("flux.dipole.naive.neg.err")
        fluxNegNaiveFlag = src.get("flux.dipole.naive.neg.flags")

        ######

        centPosPsf     = src.get("flux.dipole.psf.pos.centroid")
        centPosPsfErr  = src.get("flux.dipole.psf.pos.centroid.err")
        centPosPsfFlag = src.get("flux.dipole.psf.pos.centroid.flags")
        fluxPosPsf     = src.get("flux.dipole.psf.pos")
        fluxPosPsfErr  = src.get("flux.dipole.psf.pos.err")
        fluxPosPsfFlag = src.get("flux.dipole.psf.pos.flags")

        centNegPsf     = src.get("flux.dipole.psf.neg.centroid")
        centNegPsfErr  = src.get("flux.dipole.psf.neg.centroid.err")
        centNegPsfFlag = src.get("flux.dipole.psf.neg.centroid.flags")
        fluxNegPsf     = src.get("flux.dipole.psf.neg")
        fluxNegPsfErr  = src.get("flux.dipole.psf.neg.err")
        fluxNegPsfFlag = src.get("flux.dipole.psf.neg.flags")

        chi2dof        = src.get("flux.dipole.psf.chi2dof")

        #print centPosPsf, centNegPsf, fluxPosPsf, fluxNegPsf, centPosPsfFlag, fluxPosPsfFlag, centNegPsfFlag, fluxNegPsfFlag
        if not (centPosPsfFlag or fluxPosPsfFlag or centNegPsfFlag or fluxNegPsfFlag):
            dx = centPosPsf[0] - centNegPsf[0]
            dy = centPosPsf[1] - centNegPsf[1]
            totalFlux = fluxPosPsf + fluxNegPsf
            sigmaFlux = totalFlux / np.sqrt(fluxPosPsfErr**2 + fluxNegPsfErr**2)
            if False:
                print centNegPsf, centPosPsf, dx, dy, totalFlux, sigmaFlux, chi2dof
                cmd = "regions command {line %f %f %f %f # color=magenta width=3}" % (centNegPsf[0], centNegPsf[1], centPosPsf[0], centPosPsf[1])
                ds9.ds9Cmd(cmd, silent=False)

            # Lets look at dy, since dx~0 for visits D and E
            if dy < 0.0:
                dxdyn.append((dx, dy))
            else:
                dxdyp.append((dx, dy))

    dxn = np.array([x[0] for x in dxdyn])
    dyn = np.array([x[1] for x in dxdyn])
    dxp = np.array([x[0] for x in dxdyp])
    dyp = np.array([x[1] for x in dxdyp])

    return dxn,dyn,dxp,dyp

if __name__ == "__main__":
    indir = "/nfs/lsst/home/becker/Winter2014"
    for doPreConvolve in (False,):
        fig = plt.figure()
        for visitId in (1, 2, 3):
            figx  = visitId

            for filterName, filterId in zip(("g", "r", "i"), (1, 2, 3)):
                figy  = filterId
                sp    = fig.add_subplot(3, 3, 3*(figy-1)+figx)
                fps   = np.zeros((5, 5))

                # Just look at dipoles from zenith template
                for suffixT, airmassIdT in zip(("C",), (2,)):

                    for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"),
                                                   (0, 1, 2, 3, 4)):

                        visit = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                        mapper  = Mapper(root=os.path.join(indir, "outputs8b%s%s_doPreConvolve%s" % (suffixT, suffixI, doPreConvolve)), calibRoot = None, outputRoot = None)
                        butler  = dafPersist.ButlerFactory(mapper = mapper).create()

                        nTotal  = 0
                        frameId = 1
                        dxns = np.array(())
                        dyns = np.array(())
                        dxps = np.array(())
                        dyps = np.array(())
                        for sx in range(3):
                            for sy in range(3):
                                sensor    = "%d,%d" % (sx, sy)
                                dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                                if not butler.datasetExists(datasetType="deepDiff_differenceExp", dataId=dataId):
                                    continue
                                src       = butler.get(datasetType="deepDiff_diaSrc", dataId = dataId)
                                if False:
                                    diffim    = butler.get(datasetType="deepDiff_differenceExp", dataId=dataId)
                                    ds9.mtv(diffim, frame=frameId); frameId += 1
                                    ds9.ds9Cmd("mask transparency 80")

                                #print "#", suffixT, suffixI, dataId
                                dxn, dyn, dxp, dyp = checkDipoles(src)
                                dxns = np.append(dxns, dxn)
                                dyns = np.append(dyns, dyn)
                                dxps = np.append(dxps, dxp)
                                dyps = np.append(dyps, dyp)

                        distn = np.sqrt(dxns**2 + dyns**2)
                        angn  = np.arctan(dyns/dxns) * 180 / np.pi
                        distp = np.sqrt(dxps**2 + dyps**2)
                        angp  = np.arctan(dyps/dxps) * 180 / np.pi
                        if len(dxns)>0 and len(dxps)>0:
                            print filterName, "%s%s" % (suffixT, suffixI), visit, "dy<0 : dx=%+.3f dy=%+.3f d=%.3f ang=%.3f n=%d" % (np.median(dxns), np.median(dyns), np.median(distn), np.median(angn), len(dxns))
                            print filterName, "%s%s" % (suffixT, suffixI), visit, "dy>0 : dx=%+.3f dy=%+.3f d=%.3f ang=%.3f n=%d" % (np.median(dxps), np.median(dyps), np.median(distp), np.median(angp), len(dxps))

