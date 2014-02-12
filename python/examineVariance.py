import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.obs.lsstSim import LsstSimMapper as Mapper
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please

def convolveMi(mi, psf):
    # From detection.py
    shape = psf.computeShape()
    sigma = shape.getDeterminantRadius()
    kWidth = (int(sigma * 7 + 0.5) / 2) * 2 + 1 # make sure it is odd
    gaussFunc = afwMath.GaussianFunction1D(sigma)
    gaussKernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)

    cmi = mi.Factory(mi.getBBox(afwImage.PARENT))
    afwMath.convolve(cmi, mi, gaussKernel, afwMath.ConvolutionControl())
    return cmi

def miStats(mi):
    im  = mi.getImage()
    ma  = mi.getMask()
    va  = mi.getVariance()
    
    idx   = np.where(ma.getArray() == 0)
    medim = np.median(im.getArray()[idx])
    medva = np.median(va.getArray()[idx])
    iqrim = np.percentile(im.getArray()[idx], 75) - np.percentile(im.getArray()[idx], 25)
    
    return medim, medva, (0.741*iqrim)**2

if __name__ == "__main__":
    prefix = "/nfs/lsst/home/becker/Winter2014/"
    mapperT = Mapper(root = os.path.join(prefix, "outputs5b_doPreConvolveTrue"), calibRoot = None, outputRoot = None)
    mapperF = Mapper(root = os.path.join(prefix, "outputs5b_doPreConvolveFalse"), calibRoot = None, outputRoot = None)
    butlerT = dafPersist.ButlerFactory(mapper = mapperT).create()
    butlerF = dafPersist.ButlerFactory(mapper = mapperF).create()
    for filt, visits in zip(("g", "r", "i"), ((4866601, 4866602, 4866603), 
                                              (5866601, 5866602, 5866603),
                                              (6866601, 6866602, 6866603))):
        for visit in visits:
            for sx in range(3):
                for sy in range(3):
                    sensor    = "%d,%d" % (sx, sy)
                    dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                    calExp    = butlerT.get(datasetType="calexp", dataId = dataId) # Same for T and F
                    diffExpT  = butlerT.get(datasetType="deepDiff_differenceExp", dataId = dataId)
                    diffExpF  = butlerF.get(datasetType="deepDiff_differenceExp", dataId = dataId)
                    psfF      = diffExpF.getPsf()
                    cdiffMiF  = convolveMi(diffExpF.getMaskedImage(), psfF)

                    # First look at the empirical variance in the image, vs. median variance of Variance plane
                    medimT, medvaT, medvaIqrT = miStats(diffExpT.getMaskedImage())
                    medimF, medvaF, medvaIqrF = miStats(diffExpF.getMaskedImage())
                    
                    # Then look at difference of: standard difference image filtered with Psf vs. pre-filtered difference image
                    cdiffMiF -= diffExpT.getMaskedImage()
                    medimD, medvaD, medvaIqrD = miStats(cdiffMiF)

                    # Orginal noise properties of calexp:
                    medimO, medvaO, medvaIqrO = miStats(calExp.getMaskedImage())

                    print visit, sx, sy, "%.3f %.3f %.3f:" % (medvaO, medvaIqrO, medvaIqrO/medvaO), 
                    print "T %.3f %.3f %.3f %.3f :" % (medimT, medvaT, medvaIqrT, medvaIqrT/medvaT),
                    print "F %.3f %.3f %.3f %.3f :" % (medimF, medvaF, medvaIqrF, medvaIqrF/medvaF), 
                    print "D %.3f %.3f %.3f %.3f"   % (medimD, medvaD, medvaIqrD, medvaIqrD/medvaD)
                    #import pdb; pdb.set_trace()

                        
