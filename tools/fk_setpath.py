# boilerplate to find fieldkit
import sys
import os
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2]) # remove script name and current directory
sys.path.append(libpath)
# end boilerplate


