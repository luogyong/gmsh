#-*-coding:Utf-8-*-

# read the arguments
# 1st absolute path of the cm3py script
# 2nd name of the python script
# 3rd list of files to clean after the test separated by a semicolumn

import sys
import os

nargs = len(sys.argv)

if len(sys.argv)<3:
    print("Some arguments are missing cannot launch the test")
    sys.exit(1)
cm3py = sys.argv[1]
pyscript = sys.argv[2]
# launch the test
#import subprocess
#rcode = subprocess.call(['%s'%cm3py,' %s'%pyscript])
rcode = os.system( '%s %s'%(cm3py,pyscript))
# clean
delcode = 0
for f in sys.argv[3:]:
    dcode = os.system("rm -r %s"%f)
    if dcode != 0:
        delcode = dcode

# exit with the error code of the test
if rcode != 0 or delcode != 0:
    sys.exit(1)
