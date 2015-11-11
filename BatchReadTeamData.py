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
if(len(sys.argv) != 2):
  print "\tUSAGE: python BatchReadTeamData.py \"/path/to/the/list/of/TEAM/Detector/imges\""
  print "\t         (Don\'t use a \'~\' because it doesn't work with the glob package.)"
  exit()

# Pull the path to the TEAM images in from the command line argument and create a list of file
# names from it.
PathToImageFiles = sys.argv[1]
print "\tReading in", PathToImageFiles
FileNameList = glob.glob(PathToImageFiles)
print "\t...Found", len(FileNameList), "files."

# Now loop over all those file names and run the read code on each one.
for filename in FileNameList:
  thisCommand = "python ReadTEAMData_mode2.py " + filename
  os.system(thisCommand)

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for the batch read to finish."
exit()
