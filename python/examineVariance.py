import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.obs.lsstSim import LsstSimMapper as Mapper
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please

def convolveMi(mi, psf):
    if 0:
        # From detection.py
        shape = psf.computeShape()
        sigma = shape.getDeterminantRadius()
        kWidth = (int(sigma * 7 + 0.5) / 2) * 2 + 1 # make sure it is odd
        gaussFunc = afwMath.GaussianFunction1D(sigma)
        kernel = afwMath.SeparableKernel(kWidth, kWidth, gaussFunc, gaussFunc)
    else:
        # Do this the exact same way as imageDifference.py prefilters
        kernel = psf.getLocalKernel() # This does this for the average position

    cmi = mi.Factory(mi.getBBox(afwImage.PARENT))
    afwMath.convolve(cmi, mi, kernel, afwMath.ConvolutionControl())
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
            medimTs = []
            medvaTs = []
            iqrvaTs = []
            medimFs = []
            medvaFs = []
            iqrvaFs = []
            medimCs = []
            medvaCs = []
            iqrvaCs = []
            medimOs = []
            medvaOs = []
            iqrvaOs = []

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
                    medimT, medvaT, iqrvaT = miStats(diffExpT.getMaskedImage())
                    medimF, medvaF, iqrvaF = miStats(diffExpF.getMaskedImage())
                    
                    # Then look at this in the standard diffim,
                    # filtered by its Psf.  I need to make sure the
                    # filtering is done in the exact same way as in
                    # prefiltering                    
                    medimC, medvaC, iqrvaC = miStats(cdiffMiF)

                    # Orginal noise properties of calexp.  These get
                    # renormalized when filtering by Psf, which ends
                    # up being a multiplication by somethign other
                    # than 1
                    medimO, medvaO, iqrvaO = miStats(calExp.getMaskedImage())

                    medimTs.append(medimT)
                    medvaTs.append(medvaT)
                    iqrvaTs.append(iqrvaT)
                    medimFs.append(medimF)
                    medvaFs.append(medvaF)
                    iqrvaFs.append(iqrvaF)
                    medimCs.append(medimC)
                    medvaCs.append(medvaC)
                    iqrvaCs.append(iqrvaC)
                    medimOs.append(medimO)
                    medvaOs.append(medvaO)
                    iqrvaOs.append(iqrvaO)
                    
                    
            print visit, "%.3f %.3f %.3f:" % (np.median(medvaO), np.median(iqrvaO), np.median(np.array(iqrvaO)/np.array(medvaO))), 
            print "T %.3f %.3f %.3f %.3f :" % (np.median(medimT), np.median(medvaT), np.median(iqrvaT), np.median(np.array(iqrvaT)/np.array(medvaT))),
            print "F %.3f %.3f %.3f %.3f :" % (np.median(medimF), np.median(medvaF), np.median(iqrvaF), np.median(np.array(iqrvaF)/np.array(medvaF))), 
            print "C %.3f %.3f %.3f %.3f : %.3f %.3f"   % (np.median(medimC), np.median(medvaC), np.median(iqrvaC), 
                                                           np.median(np.array(iqrvaC)/np.array(medvaC)), 
                                                           np.median(np.array(iqrvaC)/np.array(iqrvaT)), 
                                                           np.median(np.array(medvaC)/np.array(medvaT)))

                        
