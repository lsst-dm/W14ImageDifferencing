
# Airmass of the template image
for suffixT, airmassIdT in zip(("A", "B", "C", "D", "E"),
                               (0, 1, 2, 3, 4)):

    for filterName, filterId in zip(("g", "r", "i"),
                                    (1, 2, 3)):
        
        template = "%d00%d000" % (filterId, airmassIdT)

        # Airmass of the science image
        for suffixI, airmassIdI in zip(("A", "B", "C", "D", "E"),
                                       (0, 1, 2, 3, 4)):

            visits = "%d00%d001^%d00%d002^%d00%d003" % (filterId, airmassIdI,filterId, airmassIdI,filterId, airmassIdI)
            
            print '$PIPE_TASKS_DIR/bin/imageDifferenceWinter2013.py inputs8 --id visit=%s --output=/nfs/lsst/home/becker/Winter2014/outputs8%s%s_doPreConvolveTrue --configfile imageDifferenceConfigGstar.py --config diaSourceMatchRadius=2.0 doAddCalexpBackground=False winter2013TemplateId=%s doPreConvolve=True subtract.kernel.active.singleKernelClipping=False subtract.kernel.active.kernelSumClipping=False subtract.kernel.active.spatialKernelClipping=False subtract.kernel.active.scaleByFwhm=True subtract.kernel.active.spatialKernelOrder=5 subtract.kernel.active.alardNGauss=3 subtract.kernel.active.alardDegGauss="[5,2,2]" subtract.kernel.active.alardSigGauss="[1.0,1.0,1.0]" doSelectDcrCatalog=True  --clobber-config -j 27 >& outputs8%s%s_doPreConvolveTrue.out' % (visits, suffixT, suffixI, template, suffixT, suffixI)
            print '$PIPE_TASKS_DIR/bin/imageDifferenceWinter2013.py inputs8 --id visit=%s --output=/nfs/lsst/home/becker/Winter2014/outputs8%s%s_doPreConvolveFalse --configfile imageDifferenceConfigGstar.py --config diaSourceMatchRadius=2.0 doAddCalexpBackground=False winter2013TemplateId=%s doPreConvolve=False subtract.kernel.active.singleKernelClipping=False subtract.kernel.active.kernelSumClipping=False subtract.kernel.active.spatialKernelClipping=False subtract.kernel.active.scaleByFwhm=True subtract.kernel.active.spatialKernelOrder=5 subtract.kernel.active.alardNGauss=3 subtract.kernel.active.alardDegGauss="[5,2,2]" subtract.kernel.active.alardSigGauss="[1.0,1.0,1.0]" doSelectDcrCatalog=True  --clobber-config -j 27 >& outputs8%s%s_doPreConvolveFalse.out' % (visits, suffixT, suffixI, template, suffixT, suffixI)

    
