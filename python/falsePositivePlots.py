import os, sys
from lsst.obs.lsstSim import LsstSimMapper as Mapper
from lsst.pipe.tasks.repositoryIterator import RepositoryIterator
import lsst.meas.algorithms as measAlg
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
import numpy as np
import matplotlib.pyplot as plt
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please

def nVsThreshold(threshold):
    return threshold * np.exp(-0.5 * threshold**2) / (2**2.5 * np.pi**1.5)
 
def nVsSeeing(seeing, nx=4000, ny=4072, threshold=5.0):
    # Seeing is Gaussian sigma
    density = nVsThreshold(threshold)
    # One-sided fluctuations
    ntot    = density * nx * ny / seeing**2
    # Since we are looking at both positive and negative going fluctuations
    return 2 * ntot

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
    return matchSrc, orphanSrc

if __name__ == "__main__":
    repoIter = RepositoryIterator(
        formatStr = "outputs5b_doPreConvolve%(doPreConvolve)s",
        doPreConvolve = [True, False],
    )
    prefix = "/nfs/lsst/home/becker/Winter2014/"
    
    for ind, repoInfo in enumerate(repoIter):
        outDir = os.path.join(prefix, repoInfo.name)
        if not os.path.isdir(outDir):
            print "# Missing", outDir
            continue
        print "#", outDir
        mapper = Mapper(root = outDir, calibRoot = None, outputRoot = None)
        butler = dafPersist.ButlerFactory(mapper = mapper).create()
    
        fig = plt.figure()
        spx = 0
        for filt, visits in zip(("g", "r", "i"), ((4866601, 4866602, 4866603), 
                                                  (5866601, 5866602, 5866603),
                                                  (6866601, 6866602, 6866603))):
            spy = 0
            for visit in visits:
                sp      = fig.add_subplot(3, 3, 3 * spx + spy + 1)
                sigmas  = []
                orphans = []
                totals  = []
                for sx in range(3):
                    for sy in range(3):
                        sensor = "%d,%d" % (sx, sy)
                        dataId   = {"visit": visit, "raft": "2,2", "sensor": sensor}
                        metadata = butler.get(datasetType="deepDiff_metadata", dataId = dataId)
                        diffExp  = butler.get(datasetType="deepDiff_differenceExp", dataId = dataId)
                        diaSrc   = butler.get(datasetType="deepDiff_diaSrc", dataId = dataId)
                        try:
                            foo = metadata.toString()
                        except:
                            print "# Missing", visit, sx, sy
                            continue
                        matchSrc, orphanSrc = countSources(diaSrc)
                        psf = diffExp.getPsf()
                        width, height = psf.computeKernelImage().getDimensions()
                        psfAttr = measAlg.PsfAttributes(psf, width//2, height//2)
                        psfSigma = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT) # gaussian sigma in pixels
                        print filt, visit, sx, sy, "%.2f %.2f" % (psfSigma, nVsSeeing(psfSigma)), len(diaSrc), len(matchSrc), len(orphanSrc)
                        #import pdb; pdb.set_trace()
                        
                        sigmas.append(psfSigma)
                        orphans.append(len(orphanSrc))
                        totals.append(len(diaSrc))
                meanSeeing = np.mean(sigmas)
                xt         = np.arange(4.8, 5.2, 0.1)
                yt         = np.zeros(len(xt))
                for i, t in enumerate(xt): yt[i] = nVsSeeing(meanSeeing, threshold=t)
                sp.plot(xt, yt, "k-")

                # One way to do it
                sp.errorbar(5.0, np.mean(orphans), np.std(orphans), fmt="go")
                sp.errorbar(5.0, np.mean(totals), np.std(totals), fmt="rs", alpha=0.5)

                # Try boxplots
                #boxp = sp.boxplot(orphans, positions=(5,), whis=100.0, widths=0.15)
                #void = [x.set_linestyle("-") for x in boxp["whiskers"]]

                # Try just plotting the points
                #sp.plot(4.98*np.ones(len(orphans)), orphans, "go")
                #sp.plot(5.02*np.ones(len(totals)), totals, "rs", alpha=0.5)

                plt.setp(sp.get_xticklabels()+sp.get_yticklabels(), fontsize=10, weight="bold")
                #sp.semilogy()
                if spx == 0:
                    sp.set_title("Visit %s" % (str(visit)[-1]), fontsize=14, weight="bold")
                if spx == 2:
                    sp.set_xlabel("Sigma", fontsize=12, weight="bold")
                if spy == 0:
                    sp.set_ylabel("Nfp", fontsize=12, weight="bold")
                if spy == 2:
                    sp2 = sp.twinx()
                    sp2.set_ylabel("%s-band" % (filt), fontsize=14, weight="bold", rotation="horizontal")
                    plt.setp(sp2.get_xticklabels()+sp2.get_yticklabels(), visible=False)
                    
                sp.set_ylim(0, 25)
                sp.set_xlim(4.8, 5.2)

                spy += 1
            spx += 1

        fig.suptitle("Prefilter %s" % (repoInfo.valTuple[0]), fontsize=16, weight="bold")
    plt.show()
