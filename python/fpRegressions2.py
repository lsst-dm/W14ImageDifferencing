import sys, os
import matplotlib.pyplot as plt
import numpy as np
import lsst.pipe.base as pipeBase
from lsst.obs.lsstSim import LsstSimMapper as Mapper 
import lsst.daf.persistence as dafPersist
import lsst.meas.algorithms as measAlg
import lsst.meas.astrom as measAstrom
import scipy.stats as stats
import cPickle
import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please
pexLog.Trace_setVerbosity("", -100) # No logging info please

# THis one normalizes by the number of stars in the field.  So its a likelihood.

# From DCR.py
dcr = {}
for pb in ["g", "r", "i"]:
    dcr[pb] = {}
    for sed in ["A", "G", "M"]:
        dcr[pb][sed] = {}
        dcr[pb][sed][0] = 0.0

# Note: these are the realized offsets!
dcr["g"]["A"][30] = 0.0205
dcr["g"]["A"][50] = 0.0421
dcr["r"]["A"][30] = 0.0062
dcr["r"]["A"][50] = 0.0127
dcr["i"]["A"][30] = 0.0025
dcr["i"]["A"][50] = 0.0049

dcr["g"]["G"][30] = 0.0
dcr["g"]["G"][50] = 0.0
dcr["r"]["G"][30] = 0.0 
dcr["r"]["G"][50] = 0.0
dcr["i"]["G"][30] = 0.0
dcr["i"]["G"][50] = 0.0

dcr["g"]["M"][30] = 0.0385
dcr["g"]["M"][50] = 0.0771
dcr["r"]["M"][30] = 0.0097
dcr["r"]["M"][50] = 0.0197
dcr["i"]["M"][30] = 0.0057
dcr["i"]["M"][50] = 0.0115


def nVsThreshold(threshold):
    return threshold * np.exp(-0.5 * threshold**2) / (2**2.5 * np.pi**1.5)
 
def nVsSeeing(seeing, nx=4000, ny=4072, threshold=5.0):
    density = nVsThreshold(threshold)
    ntot    = density * nx * ny / seeing**2
    return 2 * ntot

def countSources(diaSrc):
    global matches
    refMatchId = diaSrc["refMatchId"]
    flags      = diaSrc["centroid.sdss.flags"]
    x          = diaSrc["centroid.sdss.x"]
    y          = diaSrc["centroid.sdss.y"]
    edge       = diaSrc["flags.pixel.edge"]
    grs        = np.array([matches[x] if x in matches.keys() else 0.0 for x in refMatchId])
    AmatchIdx  = (flags == False) & (np.isfinite(x)) & (np.isfinite(y)) & (edge == False) & (refMatchId > 0) & (grs < 0)
    GmatchIdx  = (flags == False) & (np.isfinite(x)) & (np.isfinite(y)) & (edge == False) & (refMatchId > 0) & (grs > 0) & (grs < 1)
    MmatchIdx  = (flags == False) & (np.isfinite(x)) & (np.isfinite(y)) & (edge == False) & (refMatchId > 0) & (grs > 1)
    orphanIdx  = (flags == False) & (np.isfinite(x)) & (np.isfinite(y)) & (edge == False) & (refMatchId == 0)
    return np.sum(AmatchIdx), np.sum(GmatchIdx), np.sum(MmatchIdx), np.sum(orphanIdx)

if __name__ == "__main__":
    fontref = 10
    sigma2fwhm = 2. * np.sqrt(2. * np.log(2.))
    astrometer = measAstrom.Astrometry(measAstrom.MeasAstromConfig())
    indir = "/nfs/lsst/home/becker/Winter2014"

    if os.path.isfile("SAVEME2.pickle"):
        buff = open("SAVEME2.pickle", "rb")
        allResults = cPickle.load(buff)
        buff.close()
    else:
        allResults = {}
        for doPreConvolve in (True, False):
            results = {}
            for visitId in (1, 2, 3):
    
                #fig = plt.figure()
                AstarAll = []
                GstarAll = []
                MstarAll = []
                OrphAll  = []
                for filterName, filterId in zip(("g", "r", "i"), (1, 2, 3)):
                    for suffixT, airmassIdT in zip(("C",), (2,)):
                        for suffixI, airmassIdI, airmassI in zip(("A", "B", "C", "D", "E"),
                                                                 (0, 1, 2, 3, 4),
                                                                 (50, 30, 0, 30, 50)):
    
                            visit = int("%d00%d00%d" % (filterId, airmassIdI, visitId))
                            mapper  = Mapper(root=os.path.join(indir, "outputs8%s%s_doPreConvolve%s" % (suffixT, suffixI, doPreConvolve)), calibRoot = None, outputRoot = None)
                            butler  = dafPersist.ButlerFactory(mapper = mapper).create()
    
                            nA = []
                            nG = []
                            nM = []
                            nO = []
                            for sx in range(3):
                                for sy in range(3):
                                    sensor    = "%d,%d" % (sx, sy)
                                    dataId    = {"visit": visit, "raft": "2,2", "sensor": sensor}
                                    diaSrc    = butler.get(datasetType="deepDiff_diaSrc", dataId = dataId)
                                    diffExp   = butler.get(datasetType="deepDiff_differenceExp", dataId = dataId)
    
                                    # FWHM and number vs. seeing
                                    psf       = diffExp.getPsf()
                                    width, height = psf.computeKernelImage().getDimensions()
                                    psfAttr   = measAlg.PsfAttributes(psf, width//2, height//2)
                                    psfSigma  = psfAttr.computeGaussianWidth(psfAttr.ADAPTIVE_MOMENT) # gaussian sigma in pixels
                                    psfFwhm   = sigma2fwhm * psfSigma
                                    nExpect   = nVsSeeing(psfSigma)
    
                                    # Matches for source colors
                                    calexp    = butler.get(datasetType="calexp", dataId=dataId)
                                    src       = butler.get(datasetType="src", dataId=dataId)
                                    mat       = astrometer.useKnownWcs(src, exposure=calexp).matches
                                    matches   = {}
                                    nAcalexp  = 0.
                                    nGcalexp  = 0.
                                    nMcalexp  = 0.
                                    for m in mat:
                                        refId = m.first.getId()
                                        gr    = -2.5 * np.log10(m.first["g"]) - -2.5 * np.log10(m.first["r"])
                                        matches[refId] = gr
                                        if gr < 0: nAcalexp += 1
                                        elif (gr > 0) and (gr < 1): nGcalexp += 1
                                        else: nMcalexp += 1
    
                                    # Count the realized false positives, normalize by the number of stars in there
                                    Amatch, Gmatch, Mmatch, orphan = countSources(diaSrc)
                                    Aamp = dcr[filterName]["A"][airmassI]
                                    Gamp = dcr[filterName]["G"][airmassI]
                                    Mamp = dcr[filterName]["M"][airmassI]
                                    nA.append(Amatch/nAcalexp)
                                    nG.append(Gmatch/nGcalexp)
                                    nM.append(Mmatch/nMcalexp)
                                    nO.append(orphan)
    
                            AstarAll.append((Aamp, Aamp/psfFwhm, np.median(nA), nExpect, airmassI, filterName))
                            GstarAll.append((Gamp, Gamp/psfFwhm, np.median(nG), nExpect, airmassI, filterName))
                            MstarAll.append((Mamp, Mamp/psfFwhm, np.median(nM), nExpect, airmassI, filterName))
                            OrphAll.append((np.median(nO), nExpect, airmassI, filterName))
    
                            #plt.plot(Aamp, np.median(nA), "bo")
                            #plt.plot(Gamp, np.median(nG), "go")
                            #plt.plot(Mamp, np.median(nM), "ro")
    
                            print filterName, suffixT, suffixI, np.median(nA), np.median(nG), np.median(nM), np.median(nO)
                results[visitId] = ((AstarAll, GstarAll, MstarAll, OrphAll))
    
            fig = plt.figure()
            for visitId, shape in zip((1,2,3), ("o", "s", "^")): plt.plot([x[1] if x[1]>0 else 1e-4 for x in results[visitId][0]], [x[2]/x[3] for x in results[visitId][0]], "b%s" % (shape))
            for visitId, shape in zip((1,2,3), ("o", "s", "^")): plt.plot([x[1] if x[1]>0 else 1e-4 for x in results[visitId][1]], [x[2]/x[3] for x in results[visitId][1]], "g%s" % (shape))
            for visitId, shape in zip((1,2,3), ("o", "s", "^")): plt.plot([x[1] if x[1]>0 else 1e-4 for x in results[visitId][2]], [x[2]/x[3] for x in results[visitId][2]], "r%s" % (shape))
            plt.semilogx()
            #plt.semilogy()
            plt.xlabel("DCR Amplitude / FWHM")
            plt.ylabel("Nfp / Nexpect")
            plt.title("Prefilter %s" % (doPreConvolve))
            plt.axhline(y=1, c='k', linestyle='--')
    
            allResults[doPreConvolve] = results

        buff = open("SAVEME2.pickle", "wb")
        import cPickle
        cPickle.dump(allResults, buff)
        buff.close()

    import pdb; pdb.set_trace()
    for doPrefilter in (True, False):
        results = allResults[doPrefilter]
    
        fig = plt.figure()
    
        # A
        regX = np.array(())
        regY = np.array(())
        for visitId in (1,2,3): regX = np.concatenate((regX, np.array([x[1] if x[1]>0 else 1e-4 for x in results[visitId][0]])))
        for visitId in (1,2,3): regY = np.concatenate((regY, np.array([x[2] for x in results[visitId][0]])))
        idx = np.isfinite(np.log10(regY))
        regA = stats.linregress(np.log10(regX)[idx], np.log10(regY)[idx])
        print "A-star %s: logy = %.4f + %.4f logx.  logy=-2 at x=%.4f" % (doPrefilter, regA[1], regA[0], 10**((-regA[1]-2)/regA[0]))
        xmodA = np.arange(regX.min(), regX.max(), 0.0005)
        ymodA = 10**(regA[1] + regA[0] * np.log10(xmodA))
        plt.plot(xmodA, ymodA, "b-", alpha=0.75)
        plt.axvline(x=10**((-regA[1]-2)/regA[0]), c='b', linestyle='--', alpha=0.5)
    
        # M
        regX = np.array(())
        regY = np.array(())
        for visitId in (1,2,3): regX = np.concatenate((regX, np.array([x[1] if x[1]>0 else 1e-4 for x in results[visitId][2]])))
        for visitId in (1,2,3): regY = np.concatenate((regY, np.array([x[2] for x in results[visitId][2]])))
        idx = np.isfinite(np.log10(regY))
        regM = stats.linregress(np.log10(regX)[idx], np.log10(regY)[idx])
        print "M-star %s: logy = %.4f + %.4f logx.  logy=-2 at x=%.4f" % (doPrefilter, regM[1], regM[0], 10**((-regM[1]-2)/regM[0]))
        xmodM = np.arange(regX.min(), regX.max(), 0.0005)
        ymodM = 10**(regM[1] + regM[0] * np.log10(xmodM))
        plt.plot(xmodM, ymodM, "r-", alpha=0.75)
        plt.axvline(x=10**((-regM[1]-2)/regM[0]), c='r', linestyle='--', alpha=0.5)
        
        yfoo = np.linspace(-3, 1, 10)
        xrat = 10**((yfoo-regM[1])/regM[0]) / 10**((yfoo-regA[1])/regA[0])
        print "RATIOS", yfoo, xrat
    
        for visitId, shape in zip((1,2,3), ("o", "s", "^")): plt.plot([x[1] if x[1]>0 else 1e-4 for x in results[visitId][0]], [x[2] for x in results[visitId][0]], "b%s" % (shape))
        for visitId, shape in zip((1,2,3), ("o", "s", "^")): plt.plot([x[1] if x[1]>0 else 1e-4 for x in results[visitId][1]], [x[2] for x in results[visitId][1]], "g%s" % (shape))
        for visitId, shape in zip((1,2,3), ("o", "s", "^")): plt.plot([x[1] if x[1]>0 else 1e-4 for x in results[visitId][2]], [x[2] for x in results[visitId][2]], "r%s" % (shape))
        plt.axhline(y=0.01, c='k', linestyle='--', alpha=0.5)
        plt.semilogx()
        #plt.semilogy()
    
        plt.xlabel("DCR Amplitude / FWHM", weight="bold", fontsize=15)
        plt.ylabel("Likelihood of false positive", weight="bold", fontsize=15)
        plt.title("Prefilter %s" % (doPrefilter), weight="bold", fontsize=16)
        plt.setp(fig.gca().get_xticklabels()+fig.gca().get_yticklabels(), weight="bold", fontsize=14)
        plt.xlim(9e-5, 2e-2)
        plt.ylim(-0.1, 1.0)
        plt.show()
    
