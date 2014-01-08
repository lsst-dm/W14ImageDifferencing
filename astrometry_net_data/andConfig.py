root.starGalaxyColumn = "starnotgal"
filters = ('u', 'g', 'r', 'i', 'z', 'y')
root.magColumnMap = dict([(f,f) for f in filters])
root.magErrorColumnMap = dict([(f, f + '_err') for f in filters])
root.indexFiles = [
    'index-140108000.fits',  
    'index-140108001.fits',  
    'index-140108002.fits',  
    'index-140108003.fits',
    'index-140108004.fits',
    ]
