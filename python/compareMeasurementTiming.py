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

if __name__ == "__main__":
    prefix = "/nfs/lsst/home/becker/Winter2014/"
    suffixT = "C" 
    for filterName, filterId in zip(("g", "r", "i"), (1, 2, 3)):
        for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"), (0, 1, 2, 3, 4)):
            for visitId in (1, 2, 3):
                visit   = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                mapperA = Mapper(root=os.path.join(prefix, "outputs8%s%s_doPreConvolveFalse" % (suffixT, suffixI)), calibRoot = None, outputRoot = None)
                butlerA = dafPersist.ButlerFactory(mapper = mapperA).create()
                mapperB = Mapper(root=os.path.join(prefix, "outputs8b%s%s_doPreConvolveFalse" % (suffixT, suffixI)), calibRoot = None, outputRoot = None)
                butlerB = dafPersist.ButlerFactory(mapper = mapperB).create()
                
                for sx in range(3):
                    for sy in range(3):
                        sensor    = "%d,%d" % (sx, sy)
                        dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                        
                        metadataA = butlerA.get(datasetType="deepDiff_metadata", dataId = dataId)
                        metadataB = butlerB.get(datasetType="deepDiff_metadata", dataId = dataId)

                        # NOTE here; we had thigns set to up only
                        # measure up to 200 false positives; so we
                        # need to only look at dia sources with less
                        # than this.  Convenient, because the task is
                        # not called dipoleMeasurement.

                        if not "winter2013ImageDifference:dipoleMeasurement" in metadataA.names():
                            continue
                        if not "NFalsePositivesTotal" in metadataA.get("winter2013ImageDifference").names():
                            continue
                        timeA     = metadataA.get("winter2013ImageDifference:dipoleMeasurement").get("measureEndCpuTime")-metadataA.get("winter2013ImageDifference:dipoleMeasurement").get("measureStartCpuTime")
                        timeB     = metadataB.get("winter2013ImageDifference:dipoleMeasurement").get("measureEndCpuTime")-metadataB.get("winter2013ImageDifference:dipoleMeasurement").get("measureStartCpuTime")
                        nFpA      = metadataA.get("winter2013ImageDifference").get("NFalsePositivesTotal")
                        nFpB      = metadataB.get("winter2013ImageDifference").get("NFalsePositivesTotal")

                        if timeA >= 1.0 and timeB >= 1.0:
                            print filterName, "%s%s" % (suffixT, suffixI), dataId, "%d in %d s (%.1f/s) vs %d in %d s (%.1f/s)  %.1f" % (nFpA, timeA, nFpA/timeA, nFpB, timeB, nFpB/timeB, timeA/timeB)

