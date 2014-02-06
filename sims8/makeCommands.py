import re, os

for trim in os.listdir("."):
   if not trim.endswith(".trim"):
      continue

   workdir = re.sub("trim", "work", trim)
   outdir = re.sub("trim", "out", trim)
   sensors ="R22_S00|R22_S01|R22_S02|R22_S10|R22_S11|R22_S12|R22_S20|R22_S21|R22_S22"

   if not os.path.isdir(workdir):
      os.makedirs(workdir)
   if not os.path.isdir(outdir):
      os.makedirs(outdir)
   print 'python $PHOSIM_DIR/phosim.py $PWD/%s -w $PWD/%s -o $PWD/%s -c $PWD/clean.params.beta -s "%s" -p 16' % (trim, workdir, outdir, sensors)
