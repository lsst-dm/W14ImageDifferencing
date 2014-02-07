import os, sys
from lsst.obs.lsstSim import LsstSimMapper as Mapper
import lsst.meas.algorithms as measAlg
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
import numpy as np
import matplotlib.pyplot as plt
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please

if __name__ == "__main__":
    prefix = "/nfs/lsst/home/becker/Winter2014/inputs5"
    mapper = Mapper(root = prefix, calibRoot = None, outputRoot = None)
    butler = dafPersist.ButlerFactory(mapper = mapper).create()
    
    for filt, visits in zip(("g", "r", "i"), ((4866601, 4866602, 4866603, 48868666), 
                                              (5866601, 5866602, 5866603, 58868666),
                                              (6866601, 6866602, 6866603, 68868666))):
            for visit in visits:
                sigmas = []
                for sx in range(3):
                    for sy in range(3):
                        sensor   = "%d,%d" % (sx, sy)
                        dataId   = {"visit": visit, "raft": "2,2", "sensor": sensor}
                        calexp   = butler.get(datasetType="calexp", dataId = dataId)
                        psf = calexp.getPsf()
                        width, height = psf.computeKernelImage().getDimensions()
                        psfAttr = measAlg.PsfAttributes(psf, width//2, height//2)
                        psfSigma = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT) # gaussian sigma in pixels
                        sigmas.append(psfSigma)
                        
                print filt, visit, "%.2f" % (np.mean(sigmas))
