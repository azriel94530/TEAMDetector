#!/usr/bin/python

###################################################################################################
# Support functions for TEAM detector file analysis.                                              #
###################################################################################################

import sys
import numpy as np
import matplotlib.pyplot as plt
import os

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
