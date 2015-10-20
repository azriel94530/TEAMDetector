#!/usr/bin/python

####################################################################################################
# Open a binary file file from the TEAM detector, read and decode it using the struct package, and #
# then push it into a series of root TH2D objects.  Once we have all these nice TH2D objects, we   #
# then look at each pixel in all 32 images in the file we're examining, pick the median value of   #
# those 32 values, and construct an average dark image that we can then subtract off of each image #
# in the file.                                                                                     #
####################################################################################################

# Header, import statements etc.
import time
import sys
import os
import struct
import ROOT
import numpy
import RootPlotLibs
import PythonTools

####################################
#  BEGIN MAIN BODY OF THE CODE!!!  #
####################################

# Get the start time of this calculation
StartTime = time.time()

# Set some flags for how verbose our input and output are going to be.
Debugging = False
VerboseProcessing = True

if(len(sys.argv) != 2):
  print "Usage: [python] ReadTEAMData.py path/to/binary/image/file"
  exit()

# Pull in the path to the binary file we're going to look at
InputFilePath = sys.argv[1]
if(VerboseProcessing):
  print "\tReading in: '" + InputFilePath + "' for analysis."

# Make a directory to hold the output files as well as the raw and dark-corrected images:
OutputDir = InputFilePath.replace(".dat", "")
if(not(os.path.isdir(OutputDir))):
  os.system("mkdir " + OutputDir)
RawOutputDir = InputFilePath.replace(".dat", "") + "/" + "Raw"
if(not(os.path.isdir(RawOutputDir))):
  os.system("mkdir " + RawOutputDir)
DCOutputDir = InputFilePath.replace(".dat", "") + "/" + "DarkCorr"
if(not(os.path.isdir(DCOutputDir))):
  os.system("mkdir " + DCOutputDir)
SSOutputDir = InputFilePath.replace(".dat", "") + "/" + "SumSpec"
if(not(os.path.isdir(SSOutputDir))):
  os.system("mkdir " + SSOutputDir)

# Make a place to plot all these things...
aCanvas, aPad = RootPlotLibs.GetReadyToPlot()
aCanvas.Draw()
aCanvas.cd()
aPad.SetLeftMargin(0.06)
aPad.SetBottomMargin(0.08)
aPad.Draw()
aPad.cd()

# Open up the file with the images we're after, and define the file structure.
thisFile = open(InputFilePath, "rb")
# Header information:
HeaderLength = 9#32 bit integer words
HeaderWordLength = 4#The header words are all 4 bytes/32 bits
FirstHeaderFormatCode     = "<L"#First word is a little-endian, unsigned long
RestOfTheHeaderFormatCode = "<l"#The rest are little-endian, signed long
# Image information
nImagesPerFile = 32
nPixelsX = 1024
nPixelsY = 1024
lPixelX = 10. / float(nPixelsX)#Sensor is 10 mm x 10 mm
lPixelY = 10. / float(nPixelsY)
PixelWordLength = 2#The pixels are 16 bit/2 byte words
PixelFormatCode = "<H"

# Define the name for the animated gif file we're about to make of all the frames in this file.
RawGifFileName = RawOutputDir + "/" + InputFilePath.split("/")[-1].replace("dat", "gif")
if os.path.isfile(RawGifFileName):
  os.system("rm " + RawGifFileName)
DCGifFileName = DCOutputDir + "/" + InputFilePath.split("/")[-1].replace("dat", "gif")
if os.path.isfile(DCGifFileName):
  os.system("rm " + DCGifFileName)

# Loop over the file...  Starting with the header.
for i in range(HeaderLength):
  if(i == 0):
    thisCode = FirstHeaderFormatCode
  else:
    thisCode = RestOfTheHeaderFormatCode
  thisValue = struct.unpack(thisCode, thisFile.read(HeaderWordLength))
  #print thisValue[0]
# Now, the raw images...
RawImageTH2Ds = []
for imageNumber in range(nImagesPerFile):
  print "\tProcessing image", imageNumber + 1, "out of", str(nImagesPerFile) + "..."
  RawImageTH2Ds.append(PythonTools.MakeTH2DObject(nPixelsX, nPixelsY, lPixelX, lPixelY))
  for yPixel in range(nPixelsY):
    for xPixel in range(nPixelsX):
      thisPixelValue = struct.unpack(PixelFormatCode, thisFile.read(PixelWordLength))
      RawImageTH2Ds[-1].SetBinContent(xPixel, yPixel, thisPixelValue[0])
  RawImageTH2Ds[-1].SetName("TEAMimage_" + str(imageNumber))
  RawImageTH2Ds[-1].SetTitle("Image from TEAM Detector at Time Stamp: " + str(int(RawImageTH2Ds[-1].GetBinContent(1, 0))))
  RawImageTH2Ds[-1].Draw("colz")
  aCanvas.Update()
  ImagePlotFilePath = RawOutputDir + "/" + RawImageTH2Ds[-1].GetName() + ".png"
  if os.path.exists(ImagePlotFilePath):
    print "\tDeleting old version of", ImagePlotFilePath
    os.system("rm " + ImagePlotFilePath)
  aCanvas.SaveAs(ImagePlotFilePath)
  aPad.Print(RawGifFileName + "+", "gif+01")

# Close up the binary file.
thisFile.close()

# Now construct the median dark image.
MedianDarkImage = PythonTools.MakeTH2DObject(nPixelsX, nPixelsY, lPixelX, lPixelY)
MedianDarkImage.SetName("MedianDarkImage")
MedianDarkImage.SetTitle("Median Dark Image for TEAM Detector File: " + InputFilePath)
for yPixel in range(nPixelsY):
  for xPixel in range(nPixelsX):
    thisPixelValueList = []
    for imageNumber in range(nImagesPerFile):
      thisPixelValueList.append(RawImageTH2Ds[imageNumber].GetBinContent(xPixel, yPixel))
    MedianDarkImage.SetBinContent(xPixel, yPixel, numpy.median(thisPixelValueList))
DarkImageFilePath = OutputDir + "/" + MedianDarkImage.GetName() + ".png"
MedianDarkImage.Draw("colz")
aCanvas.Update()
aCanvas.SaveAs(DarkImageFilePath)

# And now, let's dark correct the raw images.
DCImageTH2Ds = []
for imageNumber in range(nImagesPerFile):
  print "\tDark-correcting image", imageNumber + 1, "out of", str(nImagesPerFile) + "..."
  DCImageTH2Ds.append(PythonTools.MakeTH2DObject(nPixelsX, nPixelsY, lPixelX, lPixelY))
  DCImageTH2Ds[imageNumber].SetName("DC_TEAMimage_" + str(imageNumber))
  DCImageTH2Ds[imageNumber].SetTitle("Dark-Corrected Image from TEAM Detector at Time Stamp: " + str(int(RawImageTH2Ds[-1].GetBinContent(1, 0))))
  DCImageTH2Ds[imageNumber].Add(RawImageTH2Ds[imageNumber], 1.)
  DCImageTH2Ds[imageNumber].Add(MedianDarkImage, -1.)
  DCImageTH2Ds[imageNumber].Draw("colz")
  aCanvas.Update()
  aCanvas.SaveAs(DCOutputDir + "/" + DCImageTH2Ds[imageNumber].GetName() + ".png")
  aPad.Print(DCGifFileName + "+", "gif+01")

# Now that we have dark-corrected images, let's create the calorimetric spectra: Sum01, Sum09, and
# Sum25.
aPad.SetLeftMargin(0.08)
aPad.SetRightMargin(0.02)
aPad.SetLogy(1)
Sum01Spectra = []
Sum09Spectra = []
Sum25Spectra = []
LocalMaxThreshold = 100.#[ADC Counts]
SumSpectrumThreshold = 25.#[ADC Counts]
xLo = -500.
xHi = 1500.
nBins = 200
for imageNumber in range(nImagesPerFile):
  # First, add the SumXX spectra for this image.
  Sum01Spectra.append(PythonTools.MakePixValHisto("Sum01Spectra" + str(imageNumber), "Sum(1) Spectrum for TEAM Detector at Time Stamp: " + str(int(RawImageTH2Ds[-1].GetBinContent(1, 0))), nBins, xLo, xHi, ROOT.kRed))
  Sum09Spectra.append(PythonTools.MakePixValHisto("Sum09Spectra" + str(imageNumber), "Sum(9) Spectrum for TEAM Detector at Time Stamp: " + str(int(RawImageTH2Ds[-1].GetBinContent(1, 0))), nBins, xLo, xHi, ROOT.kBlue))
  Sum25Spectra.append(PythonTools.MakePixValHisto("Sum25Spectra" + str(imageNumber), "Sum(25) Spectrum for TEAM Detector at Time Stamp: " + str(int(RawImageTH2Ds[-1].GetBinContent(1, 0))), nBins, xLo, xHi, ROOT.kBlack))
  # Next, need to make a list of the local maxima in each image.
  for yPixel in range(2, nPixelsY - 2):  #Looping over the entire image less the outer two pixesls so
    for xPixel in range(2, nPixelsX - 2):#that we can unambiguously construct the Sum09 and Sum25 spectra.
      LocalMax = True
      thisPixelValue = DCImageTH2Ds[imageNumber].GetBinContent(xPixel, yPixel)
      NeighboringPixelValues = []
      for i in range(-1, 2):
        for j in range(-1, 2):
          if((i != 0) or (j != 0)): NeighboringPixelValues.append(DCImageTH2Ds[imageNumber].GetBinContent(xPixel + i, yPixel + j))
      #print len(NeighboringPixelValues)
      NextNeighboringPixelValues = []
      for i in range(-2, 3):
        for j in range(-2, 3):
          if((i != 0) or (j != 0)): NextNeighboringPixelValues.append(DCImageTH2Ds[imageNumber].GetBinContent(xPixel + i, yPixel + j))
      #print len(NextNeighboringPixelValues)
      for npv in NextNeighboringPixelValues:
        if(thisPixelValue < npv): LocalMax = False
      if(thisPixelValue < LocalMaxThreshold): LocalMax = False
      # If we're still a local maximum, fill the Sum spectra...
      if(LocalMax):
        Sum01Spectra[imageNumber].Fill(thisPixelValue)
        NeighboringPixelSum = 0.
        for npv in NeighboringPixelValues:
          if(npv > SumSpectrumThreshold): NeighboringPixelSum += npv
        Sum09Spectra[imageNumber].Fill(thisPixelValue + NeighboringPixelSum)
        NeighboringPixelSum = 0.
        for npv in NextNeighboringPixelValues:
          if(npv > SumSpectrumThreshold): NeighboringPixelSum += npv
        Sum25Spectra[imageNumber].Fill(thisPixelValue + NeighboringPixelSum)
  # Save the sum spectra to disk...
  Sum01Spectra[imageNumber].Draw()
  aCanvas.SaveAs(SSOutputDir + "/" + Sum01Spectra[imageNumber].GetName() + ".pdf")
  Sum09Spectra[imageNumber].Draw()
  aCanvas.SaveAs(SSOutputDir + "/" + Sum09Spectra[imageNumber].GetName() + ".pdf")
  Sum25Spectra[imageNumber].Draw()
  aCanvas.SaveAs(SSOutputDir + "/" + Sum25Spectra[imageNumber].GetName() + ".pdf")

# It's probably a good idea to sum up the different sum spectra over each frame.
SummedSum01Spectrum = PythonTools.MakePixValHisto("SummedSum01Spectrum", "Sum(1)  Spectrum for TEAM Detector Summed Over All Frames in " + InputFilePath, nBins, xLo, xHi, ROOT.kRed)
SummedSum09Spectrum = PythonTools.MakePixValHisto("SummedSum09Spectrum", "Sum(9)  Spectrum for TEAM Detector Summed Over All Frames in " + InputFilePath, nBins, xLo, xHi, ROOT.kBlue)
SummedSum25Spectrum = PythonTools.MakePixValHisto("SummedSum25Spectrum", "Sum(25) Spectrum for TEAM Detector Summed Over All Frames in " + InputFilePath, nBins, xLo, xHi, ROOT.kBlack)
for spectrum in Sum01Spectra:
  SummedSum01Spectrum.Add(spectrum, 1.)
for spectrum in Sum09Spectra:
  SummedSum09Spectrum.Add(spectrum, 1.)
for spectrum in Sum25Spectra:
  SummedSum25Spectrum.Add(spectrum, 1.)
for summedspec in [SummedSum01Spectrum, SummedSum09Spectrum, SummedSum25Spectrum]:
  summedspec.Draw()
  SummedSpecFilePath = OutputDir + "/" + summedspec.GetName() + ".pdf"
  aCanvas.Update()
  aCanvas.SaveAs(SummedSpecFilePath)

# Now save everything we just did to a root file...
OutputFilePath = OutputDir + "/" + InputFilePath.split("/")[-1].replace("dat", "root")
if os.path.exists(OutputFilePath):
  print "\tDeleting old version of", OutputFilePath
  os.system("rm " + OutputFilePath)
aRootFile = ROOT.TFile(OutputFilePath, "recreate")
ListOfLists = [RawImageTH2Ds, DCImageTH2Ds, Sum01Spectra, Sum09Spectra, Sum25Spectra]
for thislist in ListOfLists:
  for thisobject in thislist:
    thisobject.Write()
for thisobject in [SummedSum01Spectrum, SummedSum09Spectrum, SummedSum25Spectrum, MedianDarkImage]:
  thisobject.Write()
aRootFile.Close()

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for this code to run."
exit()
