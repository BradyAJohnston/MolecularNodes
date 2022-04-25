#For Blender 2.92:
#You have to have admininistor priviledges so for Windows users
#JUST THIS ONCE open Blender by R-clicking on it in the start menu and select"Run as Administrator"

#Open up the System Console (under the Window menu)
#Then paste and run the following into the scripting workspace and run it
# After it has run, then close and reopen Blender normally
# You can then then use the module as normal with 
# import MODULE

#--- FROM HERE ---#

import subprocess
import sys
import os
 
# path to python.exe
python_exe = os.path.realpath(sys.executable)
 
# upgrade pip
subprocess.call([python_exe, "-m", "ensurepip"])
subprocess.call([python_exe, "-m", "pip", "install", "--upgrade", "pip"])
 
# install required packages
subprocess.call([python_exe, "-m", "pip", "install", "atomium"])

print("Atomium Successfully Installed")



#def install(package):
#    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

#install("atomium")