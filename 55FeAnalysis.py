#!/usr/bin/python

####################################################################################################
# Take the Sum(N) spectra that are created by ReadTEAMData.py and sum up a bunch of them for better#
# statistics.  Then do a fit to the 55Fe Peak
####################################################################################################

# Header, import statements etc.
import time
import sys
import numpy
import array
import ROOT
import RootPlotLibs
import PythonTools

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

# Pull in the path to the root files we're going to look at:
TEAMDetectorDataPath = "/Users/vmgehman/Documents/Detectorstan/TEAMDetector/Data/scan_0000000186/"
InputFilePaths = []
for i in range(10):
  InputFilePaths.append(TEAMDetectorDataPath + "team1k_0000000186_000000000" + str(i) + "/team1k_0000000186_000000000" + str(i) + ".root")
for i in range(10, 100):
  InputFilePaths.append(TEAMDetectorDataPath + "team1k_0000000186_00000000" + str(i) + "/team1k_0000000186_00000000" + str(i) + ".root")
InputFilePaths.append(TEAMDetectorDataPath + "team1k_0000000186_0000000100/team1k_0000000186_0000000100.root")

# Crack open the files, get the Sum(N) histograms...
InputFiles = []
for filepath in InputFilePaths:
  InputFiles.append(ROOT.TFile(filepath))
Sum01Histos = []
Sum09Histos = []
Sum25Histos = []
for inputfile in InputFiles:
  Sum01Histos.append(inputfile.Get("SummedSum01Spectrum"))
  Sum09Histos.append(inputfile.Get("SummedSum09Spectrum"))
  Sum25Histos.append(inputfile.Get("SummedSum25Spectrum"))

# Now create the global Sum(N) histograms.
xLo = -500.
xHi = 1500.
nBins = 200
GlobalSum01Spectrum = PythonTools.MakePixValHisto("GlobalSum01Spectrum", "Global Sum(1) Spectrum for Run: "  + TEAMDetectorDataPath.split("/")[-2], nBins, xLo, xHi, ROOT.kRed)
GlobalSum09Spectrum = PythonTools.MakePixValHisto("GlobalSum09Spectrum", "Global Sum(9) Spectrum for Run: "  + TEAMDetectorDataPath.split("/")[-2], nBins, xLo, xHi, ROOT.kBlue)
GlobalSum25Spectrum = PythonTools.MakePixValHisto("GlobalSum25Spectrum", "Global Sum(25) Spectrum for Run: " + TEAMDetectorDataPath.split("/")[-2], nBins, xLo, xHi, ROOT.kGreen)

# Actually add the histograms we just pulled in to the global sum histograms
for sum1histo in Sum01Histos:
  GlobalSum01Spectrum.Add(sum1histo, 1.)
for sum9histo in Sum09Histos:
  GlobalSum09Spectrum.Add(sum9histo, 1.)
for sum25histo in Sum25Histos:
  GlobalSum25Spectrum.Add(sum25histo, 1.)
GlobalSumSpectra = [GlobalSum01Spectrum, GlobalSum09Spectrum, GlobalSum25Spectrum]

# Get set up to plot some things.
aCanvas, aPad = RootPlotLibs.GetReadyToPlot()
aCanvas.Draw()
aCanvas.cd()
aPad.SetLeftMargin(0.08)
aPad.SetRightMargin(0.015)
aPad.SetBottomMargin(0.09)
aPad.SetLogy(1)
aPad.Draw()
aPad.cd()

# Plot the global sum spectra and then start doing some fits.
FitModel_Sum01 = PythonTools.GetRWFitModel("FitModel_Sum01", GlobalSum01Spectrum, 400., 50.)
FitModel_Sum09 = PythonTools.GetRWFitModel("FitModel_Sum09", GlobalSum09Spectrum, 400., 100.)
FitModel_Sum25 = PythonTools.GetRWFitModel("FitModel_Sum25", GlobalSum25Spectrum, 400., 100.)
FitModels = [FitModel_Sum01, FitModel_Sum09, FitModel_Sum25]
FitLos = [325., 250., 300.]
FitHis = [800., 700., 850.]
for i in range(len(GlobalSumSpectra)):
  GlobalSumSpectra[i].Draw()
  GlobalSumSpectra[i].Fit(FitModels[i], "QLLEM", "", FitLos[i], FitHis[i])
  FitComponents = PythonTools.GetRWFitModelComponents(FitModels[i], FitLos[i], FitHis[i])
  for fitcomp in FitComponents:
    fitcomp.Draw("same")
  FitAnnotation = PythonTools.MakeFitAnnotationRW(FitModels[i])
  FitAnnotation.Draw()
  aCanvas.Update()
  aCanvas.SaveAs(TEAMDetectorDataPath + "/" + GlobalSumSpectra[i].GetName() + ".pdf")

# Save the global sum(N) spectra to a root file so that we can do stuff with them later.
OutputFile = ROOT.TFile(TEAMDetectorDataPath + "/GlobalSumSpectra.root", "recreate")
for spec in GlobalSumSpectra:
  spec.Write()
for fitmodel in FitModels:
  fitmodel.Write()
OutputFile.Close()

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for this code to run."
exit()
