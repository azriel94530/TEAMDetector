#!/usr/bin/python

####################################################################################################
# Take the Sum(N) spectra that are created by ReadTEAMData.py and sum up a bunch of them for better#
# statistics.  Then do a fit to the 55Fe Peak
####################################################################################################

# Header, import statements etc.
import time
import sys
import glob
import h5py
import ROOT
import numpy as np
import array
import RootPlotLibs
import PythonTools
import PythonTools_ROOT

####################################
#  BEGIN MAIN BODY OF THE CODE!!!  #
####################################

# Get the start time of this calculation
StartTime = time.time()

# ROOT housekeeping...
ROOT.gROOT.Reset()
#ROOT.gROOT.ProcessLine(".L ../../CompiledTools.C+")

# Set some flags for how verbose our input and output are going to be.
Debugging = False
VerboseProcessing = True

# Check for the appropriate number of arguments, and proceed if everything looks OK.
if(len(sys.argv) != 2):
  print "\tUSAGE: python MakeGlobalSumNSpectra.py \"/path/to/the/list/of/TEAM/Detector/imges\""
  print "\t         (Don\'t use a \'~\' because it doesn't work with the glob package.)"
  exit()

# Get a list of the .dat files that went into this analysis.
PathToImageFiles = sys.argv[1]
print "\tReading in", PathToImageFiles + "..."
FileNameList = glob.glob(PathToImageFiles)
print "\t...Found", len(FileNameList), " image files.  Extracting summed Sum(N) spectra from files."

# Pick whether we're looking at root files or hdf5 files.
#FileType = "root"
FileType = "hdf5"

# Open up the files from which we're going to get the histogram values.
InputFiles = []
for name in FileNameList:
  if(FileType == "hdf5"):
    HDF5Name = name.replace(".dat", "") + "/" + name.split("/")[-1].replace(".dat", ".hdf5")
    InputFiles.append(h5py.File(HDF5Name, "r"))
  elif(FileType == "root"):
    RootName = name.replace(".dat", "") + "/" + name.split("/")[-1].replace(".dat", ".root")
    InputFiles.append(ROOT.TFile(RootName, "r"))

# Crack open the files, get the Sum(N) histograms...
Sum01Histos, Sum09Histos, Sum25Histos = [], [], []
for inputfile in InputFiles:
  #print inputfile.keys()
  if(FileType == "hdf5"):
    if(VerboseProcessing): print "\tProcessing file:", inputfile.filename + "..."
    Sum01Histos.append(np.array(inputfile['SummedSum01Vals']))
    Sum09Histos.append(np.array(inputfile['SummedSum09Vals']))
    Sum25Histos.append(np.array(inputfile['SummedSum25Vals']))
  elif(FileType == "root"):
    Sum01Histos.append(inputfile.Get("SummedSum01Spectrum"))
    Sum09Histos.append(inputfile.Get("SummedSum09Spectrum"))
    Sum25Histos.append(inputfile.Get("SummedSum25Spectrum"))

# Now create the global Sum(N) histograms.
xLo, xHi, xStep = -500., 2000., 10.
xBins = np.arange(xLo, xHi + xStep, xStep)
xBinCenters = [x + (0.5 * xStep) for x in xBins[:-1]]
GlobalSum1Spectrum  = PythonTools_ROOT.MakePixValHisto("GlobalSum01Spectrum", "Global Sum(1) Spectrum",  len(xBins), xLo, xHi, ROOT.kBlack)
GlobalSum9Spectrum  = PythonTools_ROOT.MakePixValHisto("GlobalSum09Spectrum", "Global Sum(9) Spectrum",  len(xBins), xLo, xHi, ROOT.kBlue)
GlobalSum25Spectrum = PythonTools_ROOT.MakePixValHisto("GlobalSum25Spectrum", "Global Sum(25) Spectrum", len(xBins), xLo, xHi, ROOT.kRed)
if(FileType == "hdf5"):
  GlobalSum1Vals  = np.zeros(len(Sum01Histos[0]) )
  GlobalSum9Vals  = np.zeros(len(Sum09Histos[0]) )
  GlobalSum25Vals = np.zeros(len(Sum25Histos[0]))
  for i in range(len(InputFiles)):
    GlobalSum1Vals  += Sum01Histos[i]
    GlobalSum9Vals  += Sum09Histos[i]
    GlobalSum25Vals += Sum25Histos[i]
  for i in range(len(xBinCenters)):
    GlobalSum1Spectrum.SetBinContent( GlobalSum1Spectrum.FindBin( xBinCenters[i]), GlobalSum1Vals[i] )
    GlobalSum9Spectrum.SetBinContent( GlobalSum9Spectrum.FindBin( xBinCenters[i]), GlobalSum9Vals[i] )
    GlobalSum25Spectrum.SetBinContent(GlobalSum25Spectrum.FindBin(xBinCenters[i]), GlobalSum25Vals[i])
elif(FileType == "root"):
  for i in range(len(InputFiles)):
    GlobalSum01Spectrum.Add(Sum01Histos[i], 1.)
    GlobalSum09Spectrum.Add(Sum09Histos[i], 1.)
    GlobalSum25Spectrum.Add(Sum25Histos[i], 1.)
GlobalSumSpectra = [GlobalSum1Spectrum, GlobalSum9Spectrum, GlobalSum25Spectrum]

# Get set up to plot some things.
aCanvas, aPad = RootPlotLibs.GetReadyToPlot()
aCanvas.Draw()
aCanvas.cd()
aPad.SetLeftMargin(0.08)
aPad.SetRightMargin(0.015)
aPad.SetBottomMargin(0.09)
aPad.SetLogx(0)
aPad.SetLogy(0)
aPad.Draw()
aPad.cd()

# Plot the global sum spectra and then start doing some fits.
FitModel_Sum01 = PythonTools_ROOT.GetRWFitModel("FitModel_Sum01", GlobalSum1Spectrum,  400., 50.)
FitModel_Sum09 = PythonTools_ROOT.GetRWFitModel("FitModel_Sum09", GlobalSum9Spectrum,  400., 100.)
FitModel_Sum25 = PythonTools_ROOT.GetRWFitModel("FitModel_Sum25", GlobalSum25Spectrum, 400., 100.)
FitModels = [FitModel_Sum01, FitModel_Sum09, FitModel_Sum25]
FitLos = [325., 250., 300.]
FitHis = [800., 700., 850.]
for i in range(len(GlobalSumSpectra)):
  GlobalSumSpectra[i].Draw()
  GlobalSumSpectra[i].Fit(FitModels[i], "LEM", "", FitLos[i], FitHis[i])
  FitComponents = PythonTools_ROOT.GetRWFitModelComponents(FitModels[i], FitLos[i], FitHis[i])
  for fitcomp in FitComponents:
    fitcomp.Draw("same")
  FitAnnotation = PythonTools_ROOT.MakeFitAnnotationRW(FitModels[i])
  FitAnnotation.Draw()
  aCanvas.Update()
  aCanvas.SaveAs(PathToImageFiles.replace(PathToImageFiles.split("/")[-1], GlobalSumSpectra[i].GetName() + ".png"))

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for this code to run."
exit()
