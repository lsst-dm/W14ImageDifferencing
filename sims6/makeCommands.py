import re, os

for trim in ["randomStars.2.seeing.0.6.trim", "randomStars.2.seeing.0.88.trim", 
             "randomStars.2.seeing.1.2.trim", "randomStars.2.template.trim"]:
   workdir = re.sub("trim", "work", trim)
   outdir = re.sub("trim", "out", trim)
   for sensor in ["R22_S00", "R22_S01", "R22_S02",
                  "R22_S10", "R22_S11", "R22_S12",
                  "R22_S20", "R22_S21", "R22_S22"]:
      if not os.path.isdir(os.path.join(workdir, sensor)):
         os.makedirs(os.path.join(workdir,sensor))
      if not os.path.isdir(os.path.join(outdir, sensor)):
         os.makedirs(os.path.join(outdir, sensor))
      print "python $PHOSIM_DIR/phosim.py $PWD/%s -w $PWD/%s/%s -o $PWD/%s/%s -c $PWD/clean.params.beta -s %s&" % (trim, 
              workdir, sensor, outdir, sensor, sensor)
