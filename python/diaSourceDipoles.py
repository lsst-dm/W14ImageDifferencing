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
            print centNegPsf, centPosPsf, dx, dy, totalFlux, sigmaFlux, chi2dof
            cmd = "regions command {line %f %f %f %f # color=magenta width=3}" % (centNegPsf[0], centNegPsf[1], centPosPsf[0], centPosPsf[1])
            ds9.ds9Cmd(cmd, silent=False)
              
            if dx < 0:
                dxdyn.append((dx, dy))
            else:
                dxdyp.append((dx, dy))
    dxn = np.array([x[0] for x in dxdyn])
    dyn = np.array([x[1] for x in dxdyn])
    dxp = np.array([x[0] for x in dxdyp])
    dyp = np.array([x[1] for x in dxdyp])

    distn = np.sqrt(dxn**2 + dyn**2)
    angn  = np.arctan(dyn/dxn) * 180 / np.pi
    distp = np.sqrt(dxp**2 + dyp**2)
    angp  = np.arctan(dyp/dxp) * 180 / np.pi

    print "DX NEGAGIVE:", np.median(dxn), np.median(dyn), np.median(distn), np.median(angn)
    print "DX POSITIVE:", np.median(dxp), np.median(dyp), np.median(distp), np.median(angp)

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

                #for suffixT, airmassIdT in zip(("A", "B", "C", "D", "E"),
                #                               (0, 1, 2, 3, 4)):

                # Just look at dipoles from zenith template
                for suffixT, airmassIdT in zip(("C",), (2,)):

                    #for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"),
                    #                               (0, 1, 2, 3, 4)):

                    for suffixI, airmassIdI in zip(("A",), (0,)):

                        visit = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                        mapper  = Mapper(root=os.path.join(indir, "outputs8%s%s_doPreConvolve%s_test2" % (suffixT, suffixI, doPreConvolve)), calibRoot = None, outputRoot = None)
                        butler  = dafPersist.ButlerFactory(mapper = mapper).create()

                        nTotal  = 0
                        frameId = 1
                        for sx in range(3):
                            for sy in range(3):
                                sensor    = "%d,%d" % (sx, sy)
                                dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                                if not butler.datasetExists(datasetType="deepDiff_differenceExp", dataId=dataId):
                                    continue
                                src       = butler.get(datasetType="deepDiff_diaSrc", dataId = dataId)
                                diffim    = butler.get(datasetType="deepDiff_differenceExp", dataId=dataId)
                                ds9.mtv(diffim, frame=frameId); frameId += 1
                                ds9.ds9Cmd("mask transparency 80")

                                print "#", suffixT, suffixI, dataId
                                checkDipoles(src)
                        import pdb; pdb.set_trace()
