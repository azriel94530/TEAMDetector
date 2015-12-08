#!/usr/bin/python

###################################################################################################
# Support functions for TEAM detector file analysis.                                              #
###################################################################################################

import sys
import numpy as np
import matplotlib.pyplot as plt
import os

def DarkCorrectMode2(rawimages):
  nImages = len(rawimages)
  DarkCorrectedImages = []
  nPixelsX = rawimages[0].shape[1]
  nPixelsY = rawimages[0].shape[0]
  # Construct the median dark image.
  MedianDarkImage = []
  print "\tConstructing dark image from median pixel values from all", nImages, "frames."
  irow = 0
  nStatusBarSteps = 100
  for row in range(nPixelsY):
    irow += 1
    #print irow, int(float(nPixelsY) / float(nStatusBarSteps))
    if((nPixelsX >= nStatusBarSteps) and (irow % int(float(nPixelsY) / float(nStatusBarSteps)) == 0)):
      StatusBar(irow, nPixelsY, nStatusBarSteps)
    thisDarkRow = [0] * nPixelsX
    thisDarkRow = np.array(thisDarkRow)
    for col in range(nPixelsX):
      thisPixelValueList = []
      for imageNumber in range(nImages):
        thisPixelValueList.append(rawimages[imageNumber][row][col])
      thisDarkRow[col] = np.median(thisPixelValueList)
    MedianDarkImage.append(thisDarkRow)
  print 
  MedianDarkImage = np.array(MedianDarkImage)
  # Actually do the dark correction.
  imageNumber = 0
  for rawimage in rawimages:
    imageNumber += 1
    print "\tDark-correcting exposure", imageNumber, "out of", str(nImages) + "..."
    thisDCImage = np.subtract(rawimage, MedianDarkImage)
    # Add the dark corrected image to the list of them. 
    DarkCorrectedImages.append(thisDCImage)
  return DarkCorrectedImages, MedianDarkImage

def DarkCorrectMode3(rawimages):
  nImages = len(rawimages)
  DarkCorrectedImages = []
  for imageNumber in range(0, nImages, 2):
    print "\tDark-correcting exposure", len(DarkCorrectedImages) + 1, "out of", str(nImages / 2) + "..."
    # Actually do the dark correction.
    thisDCImage = np.subtract(rawimages[imageNumber + 1], rawimages[imageNumber])
    # Add the dark corrected image to the list of them, so we can do things with them later in the code. 
    DarkCorrectedImages.append(thisDCImage)
  return DarkCorrectedImages

def StatusBar(current, total, steps):
  # Compute how far along we are...
  FractionComplete = float(current) / float(total)
  # Construct the string for the status bar.
  thisStatusBarString = '['
  for i in range(steps):
    #print float(i) / float(steps), FractionComplete * float(steps)
    if(float(i) / float(steps) <= FractionComplete):
      thisStatusBarString += '*'
    else:
      thisStatusBarString += ' '
  thisStatusBarString += ']'
  # If it's not the first one, erase the old status bar...
  if((float(current) / float(total)) > (1. / float(steps))): 
    backspaceString = '\b' * len(thisStatusBarString)
    sys.stdout.write(backspaceString)
  # Now write the new status bar.
  sys.stdout.write(thisStatusBarString)

def PlotHistogram(data, binedges, color, plottitle, xtitle, ytitle, plotfilepath):
  HistBinValues = plt.hist(data, bins=binedges, facecolor=color, alpha=0.75)[0]
  plt.xlabel(xtitle)
  plt.ylabel(ytitle)
  plt.title(plottitle)
  plt.grid(True)
  plt.yscale('log')
  if os.path.exists(plotfilepath):
    print "Deleting old version of", plotfilepath
    os.system("rm " + plotfilepath)
  plt.savefig(plotfilepath)
  plt.clf()
  return HistBinValues
