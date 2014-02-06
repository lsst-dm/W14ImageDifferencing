import os, sys
from lsst.obs.lsstSim import LsstSimMapper as Mapper
from lsst.pipe.tasks.repositoryIterator import RepositoryIterator
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please

repoIter = RepositoryIterator(
    formatStr = "outputs5_doPreConvolve%(doPreConvolve)s",
    doPreConvolve = [True, False],
)
prefix = "/astro/store/shared-scratch1/acbecker/LSST/Winter2014/"

for ind, repoInfo in enumerate(repoIter):
    outDir = os.path.join(prefix, repoInfo.name)
    print "#", outDir
    mapper = Mapper(root = outDir, calibRoot = None, outputRoot = None)
    butler = dafPersist.ButlerFactory(mapper = mapper).create()

    for filt, visits in zip(("g", "r", "i"), ((4866601, 4866602, 4866603), 
                                              (5866601, 5866602, 5866603),
                                              (6866601, 6866602, 6866603))):
        for visit in visits:
            for sx in range(3):
                for sy in range(3):
                    sensor = "%d,%d" % (sx, sy)
                    dataId   = {"visit": visit, "raft": "2,2", "sensor": sensor}
                    metadata = butler.get(datasetType="deepDiff_metadata", dataId = dataId)
                    diffexp  = butler.get(datasetType="deepDiff_differenceExp", dataId = dataId)
                    import pdb; pdb.set_trace()
