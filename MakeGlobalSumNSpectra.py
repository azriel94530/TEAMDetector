#!/usr/bin/python

####################################################################################################
# Read in the "summed Sum(N)" data from teh TEAM1K detector and add them up so that we can see     #
# more subtle features in these spectra.
####################################################################################################

# Header, import statements etc.
import time
import sys
import glob
import h5py
import matplotlib.pyplot as plt
import numpy as np

####################################
#  BEGIN MAIN BODY OF THE CODE!!!  #
####################################

# Get the start time of this calculation
StartTime = time.time()

# Check for the appropriate number of arguments, and proceed if everything looks OK.
if(len(sys.argv) != 2):
  print "\tUSAGE: python MakeGlobalSumNSpectra.py \"/path/to/the/list/of/TEAM/Detector/imges\""
  print "\t         (Don\'t use a \'~\' because it doesn't work with the glob package.)"
  exit()

# Get a list of the .dat files that went into this analysis.
PathToImageFiles = sys.argv[1]
print "\tReading in", PathToImageFiles + "..."
FileNameList = glob.glob(PathToImageFiles)
print "\t...Found", len(FileNameList), " image files.  Extracting summed Sum(N) spectra from hdf5 files."

# Open up the hdf5 files from which we're going to get the histogram values.
InputFiles = []
for name in FileNameList:
  HDF5Name = name.replace(".dat", "") + "/" + name.split("/")[-1].replace(".dat", ".hdf5")
  InputFiles.append(h5py.File(HDF5Name, "r"))

# Extract the summed Sum(N) spectra bin values.
Sum1SpectraVals, Sum9SpectraVals, Sum25SpectraVals = [], [], []
for thisFile in InputFiles:
  Sum1SpectraVals.append(np.array(thisFile['SummedSum01Vals']))
  Sum9SpectraVals.append(np.array(thisFile['SummedSum09Vals']))
  Sum25SpectraVals.append(np.array(thisFile['SummedSum25Vals']))

# Now sum each of the Sum(N) spectra to get the global ones...
GlobalSum1Spectrum  = np.zeros(len(Sum1SpectraVals[0]) )
GlobalSum9Spectrum  = np.zeros(len(Sum9SpectraVals[0]) )
GlobalSum25Spectrum = np.zeros(len(Sum25SpectraVals[0]))
for i in range(len(InputFiles)):
  GlobalSum1Spectrum  += Sum1SpectraVals[i]
  GlobalSum9Spectrum  += Sum9SpectraVals[i]
  GlobalSum25Spectrum += Sum25SpectraVals[i]

# Set up the x axis bins for these histograms.  These should really get read out of the hdf5 file...
xLo, xHi, xStep = -500., 9500., 10.
xBins = np.arange(xLo, xHi + xStep, xStep)
xBinCenters = [x + (0.5 * xStep) for x in xBins[:-1]]

# Now, plot each of the global Sum(N) histograms as a function of the value of the bin center.
plt.figure(num=None, figsize=(16, 9), dpi=80, facecolor='w', edgecolor='k')
plt.plot(xBinCenters, GlobalSum1Spectrum,  '-', color='k', linewidth=2.0, label='Sum(1)')
plt.plot(xBinCenters, GlobalSum9Spectrum,  '-', color='b', linewidth=2.0, label='Sum(9)')
plt.plot(xBinCenters, GlobalSum25Spectrum, '-', color='r', linewidth=2.0, label='Sum(25)')
plt.axis([xLo, xHi, 0.5, 1.05 * max(list(GlobalSum1Spectrum) + list(GlobalSum9Spectrum) + list(GlobalSum25Spectrum))])
plt.xlabel(InputFiles[0]['xHistoAxisTitle'][0])
plt.ylabel(InputFiles[0]['yHistoAxisTitle'][0])
plt.title('Global Sum(N) Histograms')
plt.grid(True)
plt.yscale('log')
plt.legend()
plt.savefig(FileNameList[0].replace(FileNameList[0].split("/")[-1], "") + "GobalSumHistograms.pdf")


# We need to do some fits to Landau functions, let's pull in ROOT and do that since it's easy to do that way...
import ROOT
import PythonTools_ROOT
import RootPlotLibs
ROOT.gStyle.SetOptFit(1)
Sum1RootHisto  = PythonTools_ROOT.MakePixValHisto("Sum1RootHisto",  "Sum(1) Histogram",  len(xBinCenters), xLo, xHi, ROOT.kBlack)
Sum9RootHisto  = PythonTools_ROOT.MakePixValHisto("Sum9RootHisto",  "Sum(9) Histogram",  len(xBinCenters), xLo, xHi, ROOT.kBlue)
Sum25RootHisto = PythonTools_ROOT.MakePixValHisto("Sum25RootHisto", "Sum(25) Histogram", len(xBinCenters), xLo, xHi, ROOT.kRed)
for i in range(len(xBinCenters)):
  Sum1RootHisto.SetBinContent(  Sum1RootHisto.FindBin(xBinCenters[i]), GlobalSum1Spectrum[i] )
  Sum9RootHisto.SetBinContent(  Sum9RootHisto.FindBin(xBinCenters[i]), GlobalSum9Spectrum[i] )
  Sum25RootHisto.SetBinContent(Sum25RootHisto.FindBin(xBinCenters[i]), GlobalSum25Spectrum[i])
#LandauFunc1 = PythonTools_ROOT.GetLandauFunction("LandauFunc1", 100., 9000., Sum1RootHisto.GetLineColor(), 2, 1000., 525., 50.)#Mode3
LandauFunc1 = PythonTools_ROOT.GetLandauFunction("LandauFunc1", 100., 9000., Sum1RootHisto.GetLineColor(), 2, 1000., 475., 50.)#Mode2
Sum1RootHisto.Fit("LandauFunc1",   "LEM", "", 400., 800.)
Sum1FitAnnot = PythonTools_ROOT.MakeFitAnnotationLandau(LandauFunc1)
Sum1RootHisto.GetXaxis().SetRangeUser(0., 1500.)
Sum1RootHisto.GetYaxis().SetRangeUser(0., 300.)
#LandauFunc9 = PythonTools_ROOT.GetLandauFunction("LandauFunc9", 100., 9000., Sum9RootHisto.GetLineColor(), 2, 1000., 550., 140.)# Mode 3
LandauFunc9 = PythonTools_ROOT.GetLandauFunction("LandauFunc9", 100., 9000., Sum9RootHisto.GetLineColor(), 2, 1400., 425., 120.)# Mode 2
Sum9RootHisto.Fit("LandauFunc9",   "LEM", "", 300., 1000.)
Sum9FitAnnot = PythonTools_ROOT.MakeFitAnnotationLandau(LandauFunc9)
Sum9RootHisto.GetXaxis().SetRangeUser(0., 1500.)
Sum9RootHisto.GetYaxis().SetRangeUser(0., 300.)
#LandauFunc25 = PythonTools_ROOT.GetLandauFunction("LandauFunc25", 100., 9000., Sum25RootHisto.GetLineColor(), 2, 1000., 900., 140.)# Mode 3
LandauFunc25 = PythonTools_ROOT.GetLandauFunction("LandauFunc25", 100., 9000., Sum25RootHisto.GetLineColor(), 2, 1000., 450., 150.)# Mode 2
Sum25RootHisto.Fit("LandauFunc25", "LEM", "", 400., 1400.)
Sum25FitAnnot = PythonTools_ROOT.MakeFitAnnotationLandau(LandauFunc25)
Sum25RootHisto.GetXaxis().SetRangeUser(0., 3000.)
Sum25RootHisto.GetYaxis().SetRangeUser(0., 300.)
aCanvas, aPad = RootPlotLibs.GetReadyToPlot()
aCanvas.Draw()
aCanvas.cd()
aPad.Draw()
aPad.cd()
Histos = [Sum1RootHisto, Sum9RootHisto, Sum25RootHisto]
Annots = [Sum1FitAnnot, Sum9FitAnnot, Sum25FitAnnot]
for i in range(len(Histos)):
  Histos[i].Draw()
  Annots[i].Draw()
  aCanvas.Update()
  aCanvas.SaveAs(FileNameList[0].replace(FileNameList[0].split("/")[-1], "") + Histos[i].GetName() + ".Fits.pdf")

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for this code to run."
exit()
