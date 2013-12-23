import sys, re

A0V = "starSED/kurucz/kp01_9750.fits_g45_9830.gz"
M2V = "starSED/mlt/m2.0Full.dat"
G0V = "starSED/kurucz/km20_6000.fits_g30_6020.gz"

trimIn = sys.argv[1]
oid = 0
for line in open(trimIn).readlines():
    if not line.startswith("object"):
        print line,
        continue
    
    if (oid % 10) == 0:
        star = A0V
    elif (oid % 10) == 5:
        star = M2V
    else:
        star = G0V

    print re.sub("starSED/kurucz/km50_5000.fits_g20_5140.gz", star, line), 
    oid += 1
