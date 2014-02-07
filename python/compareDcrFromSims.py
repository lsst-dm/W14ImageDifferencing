import numpy as np
import lsst.afw.image as afwImage
import lsst.daf.persistence as dafPersist
import lsst.meas.astrom as measAstrom
from lsst.obs.lsstSim import LsstSimMapper as Mapper 
import matplotlib.pyplot as plt

# Globals
mapperA = Mapper(root="/nfs/lsst/home/becker/Winter2014/inputs3", calibRoot = None, outputRoot = None)
butlerA = dafPersist.ButlerFactory(mapper = mapperA).create()
mapperB = Mapper(root="/nfs/lsst/home/becker/Winter2014/inputs5", calibRoot = None, outputRoot = None)
butlerB = dafPersist.ButlerFactory(mapper = mapperB).create()
astrometer = measAstrom.Astrometry(measAstrom.MeasAstromConfig())

# Necessary to get the colors out from astrometry.net
from lsst.meas.photocal.colorterms import Colorterm
class MyColorTerm(object):
    def get(self):
        return True
Colorterm._activeColorterms = {}
Colorterm._activeColorterms["i"] = MyColorTerm()

class RefMatchClass(object):
    def __init__(self, ref):
        self.ref   = ref
        self.refId = ref.getId()
        self.srcA  = None
        self.srcB  = None

if __name__ == "__main__":
    raft    = "2,2"
    filtA   = "i"
    for filtB, visitsB in zip(("g", "r", "i"), ((4866601, 4866602, 4866603, 48868666), 
                                                (5866601, 5866602, 5866603, 58868666),
                                                (6866601, 6866602, 6866603, 68868666))):
        for visitB in visitsB:
            if visitB/1e7 < 1:
                visitA = 6866600 + visitB%10
            else:
                visitB = 68868666
            for sx in range(3):
                for sy in range(3):
                    sensor = "%d,%d" % (sx, sy)
                    
                    dataIdA = {"filter": filtA, "visit": visitA, "raft": raft, "sensor": sensor}
                    dataIdB = {"filter": filtB, "visit": visitB, "raft": raft, "sensor": sensor}
                    
                    expA = butlerA.get(datasetType="calexp", dataId=dataIdA)
                    expB = butlerB.get(datasetType="calexp", dataId=dataIdB)
        
                    srcA = butlerA.get(datasetType="src", dataId=dataIdA)
                    srcB = butlerB.get(datasetType="src", dataId=dataIdB)

                    # I need to be able to get all the magnitudes here, not just i-band
                    matchA = astrometer.useKnownWcs(srcA, exposure=expA).matches
                    matchB = astrometer.useKnownWcs(srcB, exposure=expB).matches
               
                    offsetA = np.array([x[2] for x in matchA]) * 180 / np.pi * 3600
                    offsetB = np.array([x[2] for x in matchB]) * 180 / np.pi * 3600
                   
                    print "WCS astrom quality: %.2f (%.2f) and %.2f (%.2f)" % (np.mean(offsetA), np.std(offsetA),
                                                                               np.mean(offsetB), np.std(offsetB))
        
                    refMatches = {}
                    for match in matchA:
                        refId = match[0].getId()
                        if not refMatches.has_key(refId):
                            refMatches[refId] = RefMatchClass(match[0])
                        refMatches[refId].srcA = match[1]
    
                    # Just look for things in both images
                    for match in matchB:
                        refId = match[0].getId()
                        if refMatches.has_key(refId):
                            refMatches[refId].srcB = match[1]
                    badKeys = [key for key, value in refMatches.iteritems() if value.srcB is None]
                    for key in badKeys:
                        del refMatches[key]
    
                    fig = plt.figure()
                    sp  = fig.add_subplot(111)
                    for key, value in refMatches.iteritems():
                        gRef   = -2.5 * np.log10(value.ref.get("g"))
                        rRef   = -2.5 * np.log10(value.ref.get("r"))
                        if abs((gRef-rRef)-0.306) < 0.001:
                            sType = "G"
                        elif abs((gRef-rRef)-1.354) < 0.001:
                            sType = "M"
                        elif abs((gRef-rRef)+0.239) < 0.001:
                            sType = "A"
                        else:
                            sType = None
                        xA, yA = value.srcA.getCentroid()
                        xB, yB = value.srcB.getCentroid()
                        rA     = value.srcA.getCoord()[0]
                        dA     = value.srcA.getCoord()[1]
                        rB     = value.srcB.getCoord()[0]
                        dB     = value.srcB.getCoord()[1]

                        import pdb; pdb.set_trace()
                        dx, dy = (xA - xB), (yA - yB)
                        dist   = np.sqrt(dx**2 + dy**2)
                        theta  = np.arctan2(dy, dx)
                        if sType == "A": color = "b"
                        elif sType == "G": color="g"
                        elif sType == "M": color="r"
                        else: color="k"
    
                        sp.quiver(xA, yA, dx, dy, scale=0.1, angles='xy', color=color, pivot='middle', width=0.5, headwidth=1, headlength=2, units="xy")
                    plt.show()
    