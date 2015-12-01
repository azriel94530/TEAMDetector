#!/usr/bin/python

###################################################################################################
# Take a path to a group of files and run the ReadTEAMData.py script on said files.               #
###################################################################################################
# Header, import statements etc.
import time
import sys
import glob
import os

####################################
#  BEGIN MAIN BODY OF THE CODE!!!  #
####################################

# Get the start time of this calculation
StartTime = time.time()

# Check for the appropriate number of arguments, and proceed if everything looks OK.
if(len(sys.argv) != 3):
  print "\tUSAGE: python BatchReadTeamData.py \"/path/to/the/list/of/TEAM/Detector/imges\" ModeNumber"
  print "\t       (Don\'t use a \'~\' in the path because it doesn't work with the glob package.)"
  print "\t       \"ModeNumber\" is either a 2 or a 3, corresponding to the trigger mode."
  exit()

TriggerMode = int(sys.argv[2])
if(TriggerMode == 2):
  print "\tUsing Mode 2 trigger..."
elif(TriggerMode == 3):
  print "\tUsing Mode 3 trigger..."
else:
  print "\tTrigger mode must be either 2 or 3."
  exit()

# Pull the path to the TEAM images in from the command line argument and create a list of file
# names from it.
PathToImageFiles = sys.argv[1]
print "\tReading in", PathToImageFiles
FileNameList = glob.glob(PathToImageFiles)
print "\t...Found", len(FileNameList), "files."

# Now loop over all those file names and run the read code on each one.
for filename in FileNameList:
  if(TriggerMode == 2): thisCommand = "python ReadTEAMData_mode2.py " + filename
  if(TriggerMode == 3): thisCommand = "python ReadTEAMData_mode3.py " + filename
  os.system(thisCommand)

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for the batch read to finish."
exit()
