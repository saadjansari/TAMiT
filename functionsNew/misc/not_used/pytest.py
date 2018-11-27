# Expected usage:
# 
# python * script_name * arg1 arg2
# script_name: pytest.py
# arg1: a or a-b where a, b are integers corresponding to cell_a - cell_b 
# arg2: a or a-b where a, b are integers corresponding to frame_a - frame_b
# 
# Note: Run in specific path that contains subfolders -> cell_x -> frame_x -> *pictures to open are over here*
# 

import os
from os.path import expanduser
import subprocess
import sys
pwd = os.getcwd()

# specify which file to open
pictureToOpen = "mt_fit_comparison.png"

# load in arguments corresponding to interesting cells and frames
cellRange = sys.argv[1]
frameRange = sys.argv[2]

cellInfo = cellRange.split('-')
frameInfo = frameRange.split('-')

cellStart = int(cellInfo[0])
if len(cellInfo) == 1:
	cellEnd = cellStart
else:
	cellEnd = int(cellInfo[1])

frameStart = int(frameInfo[0])
if len(frameInfo) == 1:
        frameEnd = frameStart
else:
        frameEnd = int(frameInfo[1])

for currCell in range(cellStart, 1+cellEnd, 1):

	for currFrame in range(frameStart, 1+frameEnd, 1):
		
		filepath = (pwd + "/cell_%d/frame_%d/" + pictureToOpen) %(currCell, currFrame)
		subprocess.call(['open', filepath])




