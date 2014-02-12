root.starGalaxyColumn = "starnotgal"
filters = ('u', 'g', 'r', 'i', 'z', 'y')
root.magColumnMap = dict([(f,f) for f in filters])
root.magErrorColumnMap = dict([(f, f + '_err') for f in filters])
root.indexFiles = [
    'index-140212000.fits',  
    'index-140212001.fits',  
    'index-140212002.fits',  
    'index-140212003.fits',
    'index-140212004.fits',
    ]
