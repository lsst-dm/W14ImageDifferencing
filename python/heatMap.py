import sys, os
import matplotlib.pyplot as plt
import numpy as np
import lsst.pipe.base as pipeBase
from lsst.obs.lsstSim import LsstSimMapper as Mapper 
import lsst.daf.persistence as dafPersist

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please
pexLog.Trace_setVerbosity("", -100) # No logging info please

def countSources(diaSrc):
    refMatchId = diaSrc["refMatchId"]
    flags      = diaSrc["centroid.sdss.flags"]
    x          = diaSrc["centroid.sdss.x"]
    y          = diaSrc["centroid.sdss.y"]
    edge       = diaSrc["flags.pixel.edge"]
    matchIdx   = (flags == False) & (np.isfinite(x)) & (np.isfinite(y)) & (edge == False) & (refMatchId > 0)
    orphanIdx  = (flags == False) & (np.isfinite(x)) & (np.isfinite(y)) & (edge == False) & (refMatchId == 0)
    matchSrc   = diaSrc[matchIdx]
    orphanSrc  = diaSrc[orphanIdx]
    return len(matchSrc), len(orphanSrc)

fontref = 10
if __name__ == "__main__":
    indir = "/nfs/lsst/home/becker/Winter2014"
    for doPreConvolve in (True, False):

        fig = plt.figure(figsize=(9.95, 9.56), dpi=100)
        for visitId in (1, 2, 3):
            figx  = visitId

            for filterName, filterId in zip(("g", "r", "i"), (1, 2, 3)):
                figy  = filterId
                sp    = fig.add_subplot(3, 3, 3*(figy-1)+figx)
                fps   = np.zeros((5, 5))

                for suffixT, airmassIdT in zip(("A", "B", "C", "D", "E"),
                                               (0, 1, 2, 3, 4)):

                    for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"),
                                                   (0, 1, 2, 3, 4)):

                        visit = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                        mapper  = Mapper(root=os.path.join(indir, "outputs8%s%s_doPreConvolve%s" % (suffixT, suffixI, doPreConvolve)), calibRoot = None, outputRoot = None)
                        butler  = dafPersist.ButlerFactory(mapper = mapper).create()

                        nTotal = []
                        for sx in range(3):
                            for sy in range(3):
                                sensor    = "%d,%d" % (sx, sy)
                                dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                                src       = butler.get(datasetType="deepDiff_diaSrc", dataId = dataId)
                                sourceM, sourceO = countSources(src)
                                nTotal.append(sourceM + sourceO)
                        fps[airmassIdT,airmassIdI] = np.median(nTotal)
                        print filterName, suffixT, suffixI, nTotal

                print filterName, visitId, fps
                sp.pcolor(fps, cmap=plt.cm.Greys, shading="faceted", vmin=0, vmax=300)
                sp.set_xticks(np.arange(fps.shape[0])+0.5, minor=False)
                sp.set_yticks(np.arange(fps.shape[1])+0.5, minor=False)
                sp.invert_yaxis()
                sp.set_xticklabels(("A", "B", "C", "D", "E"), minor=False, weight="bold", fontsize=fontref+1)
                sp.set_yticklabels(("A", "B", "C", "D", "E"), minor=False, weight="bold", fontsize=fontref+1)
                for i in range(5):
                    for j in range(5):
                        sp.text(i+0.5, j+0.5, "%d" % (int(fps[j,i])), weight="bold", color="r", fontsize=fontref, horizontalalignment="center", verticalalignment="center")
                if figy == 1:
                    sp.set_title("Seeing %d" % (visitId), weight="bold", fontsize=fontref+4)
                if figx == 1:
                    sp.set_ylabel("Image Id; %s-band" % (filterName), weight="bold", fontsize=fontref+2) 
                if figy == 3:
                    sp.set_xlabel("Template Id", weight="bold", fontsize=fontref+2)
                    
        fig.suptitle("Prefilter=%s" % (doPreConvolve), weight="bold", fontsize=fontref+6)
        import pdb; pdb.set_trace()
