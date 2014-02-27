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

        if (not centPosNaiveFlag) and (not centNegNaiveFlag):
            dx = centPosNaive[0] - centNegNaive[0]
            dy = centPosNaive[1] - centNegNaive[1]
            print centNegNaive, centPosNaive, dx, dy
            cmd = "regions command {line %f %f %f %f # color=magenta width=3}" % (centNegNaive[0], centNegNaive[1], centPosNaive[0], centPosNaive[1])
            ds9.ds9Cmd(cmd, silent=False)
                

if __name__ == "__main__":
    indir = "/nfs/lsst/home/becker/Winter2014"
    for doPreConvolve in (True, False):
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

                    for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"),
                                                   (0, 1, 2, 3, 4)):

                        visit = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                        mapper  = Mapper(root=os.path.join(indir, "outputs8%s%s_doPreConvolve%s" % (suffixT, suffixI, doPreConvolve)), calibRoot = None, outputRoot = None)
                        butler  = dafPersist.ButlerFactory(mapper = mapper).create()

                        nTotal  = 0
                        frameId = 1
                        for sx in range(3):
                            for sy in range(3):
                                sensor    = "%d,%d" % (sx, sy)
                                dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                                src       = butler.get(datasetType="deepDiff_diaSrc", dataId = dataId)
                                diffim    = butler.get(datasetType="deepDiff_differenceExp", dataId=dataId)
                                ds9.mtv(diffim, frame=frameId); frameId += 1
                                ds9.ds9Cmd("mask transparency 80")

                                print "#", suffixT, suffixI, dataId
                                checkDipoles(src)
                        import pdb; pdb.set_trace()
