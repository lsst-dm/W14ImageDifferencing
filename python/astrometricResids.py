import os
import numpy as np
from lsst.obs.lsstSim import LsstSimMapper as Mapper
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) 

if __name__ == "__main__":
    prefix = "/nfs/lsst/home/becker/Winter2014/"
    suffixT = "C" 
    for filterName, filterId in zip(("g", "r", "i"), (1, 2, 3)):
        for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"), (0, 1, 2, 3, 4)):
            bDlatMed = []; bDlatStd = []; bDlonMed = []; bDlonStd = []
            gDlatMed = []; gDlatStd = []; gDlonMed = []; gDlonStd = []
            rDlatMed = []; rDlatStd = []; rDlonMed = []; rDlonStd = []
            for visitId in (1, 2, 3):
                visit  = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                mapper = Mapper(root=os.path.join(prefix, "outputs8b%s%s_doPreConvolveFalse" % (suffixT, suffixI)), calibRoot = None, outputRoot = None)
                butler = dafPersist.ButlerFactory(mapper = mapper).create()
                for sx in range(3):
                    for sy in range(3):
                        sensor = "%d,%d" % (sx, sy)
                        dataId   = {"visit": visit, "raft": "2,2", "sensor": sensor}
                        metadata = butler.get(datasetType="deepDiff_metadata", dataId = dataId).get("winter2013ImageDifference")
                        bDlatMed.append(metadata.get("RegisterBlueLatOffsetMedian"))
                        bDlatStd.append(metadata.get("RegisterBlueLatOffsetStd"))
                        bDlonMed.append(metadata.get("RegisterBlueLongOffsetMedian"))
                        bDlonStd.append(metadata.get("RegisterBlueLongOffsetStd"))
                        gDlatMed.append(metadata.get("RegisterGreenLatOffsetMedian"))
                        gDlatStd.append(metadata.get("RegisterGreenLatOffsetStd"))
                        gDlonMed.append(metadata.get("RegisterGreenLongOffsetMedian"))
                        gDlonStd.append(metadata.get("RegisterGreenLongOffsetStd"))
                        rDlatMed.append(metadata.get("RegisterRedLatOffsetMedian"))
                        rDlatStd.append(metadata.get("RegisterRedLatOffsetStd"))
                        rDlonMed.append(metadata.get("RegisterRedLongOffsetMedian"))
                        rDlonStd.append(metadata.get("RegisterRedLongOffsetStd"))
                

            db = np.sqrt(np.array(bDlatMed)**2+np.array(bDlonMed)**2)
            dg = np.sqrt(np.array(gDlatMed)**2+np.array(gDlonMed)**2)
            dr = np.sqrt(np.array(rDlatMed)**2+np.array(rDlonMed)**2)
            # Use median
            bDist = np.median(db)
            gDist = np.median(dg)
            rDist = np.median(dr)

            # Std estimate
            bStd = 0.741 * (np.percentile(db, 75) - np.percentile(db, 25))
            gStd = 0.741 * (np.percentile(dg, 75) - np.percentile(dg, 25))
            rStd = 0.741 * (np.percentile(dr, 75) - np.percentile(dr, 25))
            print filterName, "%s%s" % (suffixT, suffixI), "%.4f %.4f, %.4f %.4f, %.4f %.4f" % (bDist, bStd, gDist, gStd, rDist, rStd)
    
