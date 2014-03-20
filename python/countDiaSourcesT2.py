from operator import itemgetter
import sys
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.getRepositoryData import GetRepositoryDataTask
from lsst.pipe.tasks.repositoryIterator import RepositoryIterator, SourceData

import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("CameraMapper", -100) # No logging info please
pexLog.Trace_setVerbosity("", -100) # No logging info please


MaxRepos = 1e5

class GetSourcesTask(GetRepositoryDataTask):
    """Obtain sources from a repository
    """
    datasetType = None

    @classmethod
    def _makeArgumentParser(cls):
        # set dataSetType
        parser = pipeBase.ArgumentParser(name=cls._DefaultName) #, datasetType="deepDiff_diaSrc")
        parser.add_id_argument(name="--id", datasetType="deepDiff_diaSrc", 
                               help="data ID, e.g. --id visit=6866601 raft=2,2 sensor=1,1", level="sensor")
        #parser.add_id_argument(name="deepDiff_kernelSrc", datasetType="deepDiff_kernelSrc", 
        #                       help="deepDiff_kernelSrc", level="sensor")
        return parser

    @pipeBase.timeMethod
    def run(self, sensorRefList):
        """Get sources from a repository

        @return a Struct containing:
        - sourceDict: a dict of sensorRef: source table
        """
        idListStruct = self.getIdList(sensorRefList)
        sourceDict = dict()
        for datasetType in ("deepDiff_diaSrc", "deepDiff_kernelSrc"):
            sourceDict[datasetType] = [sensorRef.get(datasetType) for sensorRef in sensorRefList]
        return pipeBase.Struct(
            idValList = idListStruct.idValList,
            idKeyTuple = idListStruct.idKeyTuple,
            sourceDict = sourceDict,
        )


    # Override these methods to be no-ops.
    def writeConfig(self, butler, clobber=False):
        pass

    def writeMetadata(self, dataRef):
        pass

def countSources(sourceArr, dims=4000, buff=320):
    # Try this for a centroid first
    flags1         = sourceArr["centroid.sdss.flags"]
    x1             = sourceArr["centroid.sdss.x"]
    y1             = sourceArr["centroid.sdss.y"]

    # And then this
    flags2a        = sourceArr["flux.dipole.psf.neg.centroid.flags"]
    flags2b        = sourceArr["flux.dipole.psf.pos.centroid.flags"]
    x2             = sourceArr["flux.dipole.psf.centroid.x"]
    y2             = sourceArr["flux.dipole.psf.centroid.y"]

    # Other issues
    badcentroid    = sourceArr["flags.badcentroid"]
    edge           = sourceArr["flags.pixel.edge"]

    # This tells us if its even in this dataId
    visit          = sourceArr["visit"]

    # Match or orphan
    #srcMatchId    = sourceArr["srcMatchId"]
    # DUH!  False positives (vs. false negatives) match to sources
    srcMatchId    = sourceArr["refMatchId"]

    # Slicing:
    # First look at the things that could be naturally centroided
    if False:
        # This is a tighter set of cuts
        idx1M = (edge == False) & (badcentroid == False) & (flags1 == False) & \
            (visit > 0) & np.isfinite(x1) & np.isfinite(y1) & (x1>0.) & (y1>0.) & (srcMatchId > 0)
        idx1O = (edge == False) & (badcentroid == False) & (flags1 == False) & \
            (visit > 0) & np.isfinite(x1) & np.isfinite(y1) & (x1>0.) & (y1>0.) & (srcMatchId == 0)
    else:
        # Looser set of cuts (no badcentroid or sdss.flags)
        idx1M = (edge == False) & \
            (visit > 0) & np.isfinite(x1) & np.isfinite(y1) & (x1>0.) & (y1>0.) & (srcMatchId > 0)
        idx1O = (edge == False) & \
            (visit > 0) & np.isfinite(x1) & np.isfinite(y1) & (x1>0.) & (y1>0.) & (srcMatchId == 0)

    # Raw bulk false positives (outside of the edges)
    idx2 = (edge == False) & (visit > 0)

    # Slice and dice...
    sourceM = sourceArr[idx1M]
    sourceO = sourceArr[idx1O]
    nTot    = len(sourceArr[idx2])
    return nTot, sourceM, sourceO


repoIter = RepositoryIterator(
    formatStr = "outputs5b_doPreConvolve%(doPreConvolve)s",
    doPreConvolve = [False, True],
)


if __name__ == "__main__":
    for C in (1, 2, 3):
        for f in (4, 5, 6):
            visit = "%d86660%d" % (f, C)
            raft  = "2,2"
            for ind, repoInfo in enumerate(repoIter):
                outDir = repoInfo.name
                
                diaSrcData = SourceData(
                    datasetType = "deepDiff_diaSrc",
                    sourceKeyTuple = (
                        "flags.badcentroid",
                        "flags.pixel.edge",
                        "centroid.sdss.x",
                        "centroid.sdss.y",
                        "centroid.sdss.flags",
                        "centroid.gaussian.x",
                        "centroid.gaussian.y",
                        "centroid.gaussian.flags",
                        "centroid.naive.x",
                        "centroid.naive.y",
                        "centroid.naive.flags",
                        "centroid.dipole.naive.neg.x",
                        "centroid.dipole.naive.neg.y",
                        "centroid.dipole.naive.neg.flags",
                        "centroid.dipole.naive.pos.x",
                        "centroid.dipole.naive.pos.y",
                        "centroid.dipole.naive.pos.flags",
                        "flux.dipole.psf.centroid.x",
                        "flux.dipole.psf.centroid.y",
                        "flux.dipole.psf.centroid.flags",
                        "flux.dipole.psf.neg.centroid.x",
                        "flux.dipole.psf.neg.centroid.y",
                        "flux.dipole.psf.neg.centroid.flags",
                        "flux.dipole.psf.pos.centroid.x",
                        "flux.dipole.psf.pos.centroid.y",
                        "flux.dipole.psf.pos.centroid.flags",
                        "classification.dipole",
                        "srcMatchId",
                        "refMatchId",
                        ),
                    )
            
                for sx in (0, 1, 2):
                    for sy in (0, 1, 2):
                        sensor = "%d,%d" % (sx, sy)
                        task_args = [repoInfo.name, "--id", "visit=%s"%(visit), "raft=%s"%(raft), "sensor=%s"%(sensor)]
                        #print outDir, visit, raft, sensor, ":",
                        try:
                            fullResults = GetSourcesTask.parseAndRun(task_args, doReturnResults=True)
                        except Exception, e:
                            #print e
                            continue
                        else:
                            #print "OK", 
                            pass
                        taskResult = fullResults.resultList[0].result
                        numDiaSources = diaSrcData.addSourceMetrics(repoInfo = repoInfo, 
                                                                    idKeyTuple = taskResult.idKeyTuple,
                                                                    idValList = taskResult.idValList,
                                                                    sourceTableList = taskResult.sourceDict["deepDiff_diaSrc"])
                        #print numDiaSources

                try:
                    diaSrcData.finalize()
                except:
                    print "FAIL"
                    continue
                sourcesT = []
                sourcesM = []
                sourcesO = []
                for i in range(len(diaSrcData.repoArr)):
                    key = diaSrcData.repoArr[i]
                    nTot, sourceM, sourceO = countSources(diaSrcData.sourceArr[i])
                    sourcesT.append(nTot)
                    sourcesM.append(len(sourceM))
                    sourcesO.append(len(sourceO))
                print repoInfo.name, visit, diaSrcData.repoArr[0], len(sourcesT), int(np.mean(sourcesT)+0.5), int(np.mean(sourcesM)+0.5), int(np.mean(sourcesO)+0.5)
    
