import os, sys
from lsst.obs.lsstSim import LsstSimMapper as Mapper
from lsst.pipe.tasks.repositoryIterator import RepositoryIterator
import lsst.meas.algorithms as measAlg
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
import numpy as np
import matplotlib.pyplot as plt
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please

if __name__ == "__main__":
    repoIter = RepositoryIterator(
        formatStr = "outputs5_doPreConvolve%(doPreConvolve)s",
        doPreConvolve = [True, False],
    )
    prefix = "/nfs/lsst/home/becker/Winter2014/"
    
    for ind, repoInfo in enumerate(repoIter):
        outDir = os.path.join(prefix, repoInfo.name)
        print "#", outDir
        mapper = Mapper(root = outDir, calibRoot = None, outputRoot = None)
        butler = dafPersist.ButlerFactory(mapper = mapper).create()
    
        fig = plt.figure()
        spx = 0
        for filt, visit in zip(("g", "r", "i"), (4866601, 5866601, 6866601)):
            for sx in range(3):
                for sy in range(3):
                    sensor = "%d,%d" % (sx, sy)
                    dataId   = {"visit": visit, "raft": "2,2", "sensor": sensor}
                    metadata = butler.get(datasetType="deepDiff_metadata", dataId = dataId)
                    print filt, visit, sx, sy,
                    for key in ("ALBasisDegGauss", "ALBasisNGauss", "ALBasisSigGauss", "ALKernelSize"):
                        print metadata.get("winter2013ImageDifference:subtract").get(key),
                    print
