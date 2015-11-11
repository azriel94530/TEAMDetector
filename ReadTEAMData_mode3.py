#!/usr/bin/python

####################################################################################################
# Open a binary file file from the TEAM detector running in "mode 3,"" read and decode it using the#
# struct package, and then push it into a series of 2D histograms using imshow.  Once we have all  #
#these images in a more python-friendly format, then subtract the pre-trigger dark image (numbers: #
# 0, 2, 4, 6...) from the post trigger exposure (numbers: 1, 3, 5, 7...).  Once that's done, we'll #
# go ahead and save a png image of the raw and dark corrected images, and then write the image     #
# array to an HDF5 file.                                                                           #
####################################################################################################

# Header, import statements etc.
import time
import sys
import os
import struct
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filters
import scipy.ndimage as ndimage
import h5py
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
  print "Usage: [python] ReadTEAMData_mode3.py path/to/binary/image/file"
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
lSensorX = 10.#[mm]
lSensorY = 10.#[mm]
lPixelX = 10. / float(nPixelsX)#Sensor is 10 mm x 10 mm
lPixelY = 10. / float(nPixelsY)
PixelWordLength = 2#The pixels are 16 bit/2 byte words
PixelFormatCode = "<H"

# Parameters to set up the plots of all these images we're reading in...
xFigSize, yFigSize = 12., 9.
xAxisTitle, yAxisTitle = 'x Position [mm]', 'y Position [mm]'

# Let's open up an hdf5 file to save our output to.
OutputFilePath = OutputDir + "/" + InputFilePath.split("/")[-1].replace("dat", "hdf5")
OutputFile = h5py.File(OutputFilePath,'w')
OutputFile.create_dataset('FigSize', data=[xFigSize, yFigSize])
OutputFile.create_dataset('xImageAxisTitle', data=[xAxisTitle])
OutputFile.create_dataset('yImageAxisTitle', data=[yAxisTitle])

# Loop over the file...  Starting with the header.
for i in range(HeaderLength):
  if(i == 0):
    thisCode = FirstHeaderFormatCode
  else:
    thisCode = RestOfTheHeaderFormatCode
  thisValue = struct.unpack(thisCode, thisFile.read(HeaderWordLength))
  #print thisValue[0]
# Now, read the raw images...
ImagesInThisFile = []
plt.figure(num=None, figsize=(xFigSize, yFigSize), dpi=80, facecolor='w', edgecolor='k')
for imageNumber in range(nImagesPerFile):
  print "\tProcessing image", imageNumber + 1, "out of", str(nImagesPerFile) + "..."
  thisImage = []
  for row in range(nPixelsY):
    thisRow = [0] * nPixelsX
    thisRow = np.array(thisRow)
    for col in range(nPixelsX):
      thisPixelValue = int(struct.unpack(PixelFormatCode, thisFile.read(PixelWordLength))[0])
      thisRow[col] = thisPixelValue
    thisImage.append(thisRow)
  thisImage = np.array(thisImage)
  ImagesInThisFile.append(thisImage)
  # Now that we've got the most recent raw image, let's save it to the output file and plot it.
  #OutputFile.create_dataset('RawImage_' + str(imageNumber), data=thisImage)
  plt.imshow(thisImage, alpha=0.75, aspect='auto', origin='lower', extent=[0.,lSensorX, 0.,lSensorY], interpolation='none')
  plt.colorbar()
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.title("Image from TEAM Detector at Time Stamp: " + str(thisImage[0][1]))
  plt.grid(True)
  # Now, save the newly reconstructed raw image to a png file.
  ImagePlotFilePath = RawOutputDir + "/TEAMimage_" + str(imageNumber) + ".png"
  if os.path.exists(ImagePlotFilePath):
    print "Deleting old version of", ImagePlotFilePath
    os.system("rm " + ImagePlotFilePath)
  plt.savefig(ImagePlotFilePath)
  plt.clf()
# Now that we're done with it, close up the binary file.
thisFile.close()
# And convert ImagesInThisFile from a list of np.arrays to a three-dimensional np.array
ImagesInThisFile = np.array(ImagesInThisFile)

# And now, let's dark correct the raw images.
DCImages = []
for imageNumber in range(0, nImagesPerFile, 2):
  print "\tDark-correcting exposure", len(DCImages) + 1, "out of", str(nImagesPerFile / 2) + "..."
  # Actually do the dark correction.
  thisDCImage = np.subtract(ImagesInThisFile[imageNumber + 1], ImagesInThisFile[imageNumber])
  # And save it to the hdf5 file.
  OutputFile.create_dataset('DCImage_' + str(imageNumber), data=thisDCImage)
  # Add the dark corrected image to the list of them, so we can do things with them later in the code. 
  DCImages.append(thisDCImage)

# Now that we have dark-corrected images, let's create the calorimetric spectra: Sum01, Sum09, and
# Sum25.
Sum01Spectra = []
Sum09Spectra = []
Sum25Spectra = []
LocalMaxThreshold = 200.#[ADC Counts]
LocalMaxNeighborhood = 10
SumNThreshold = -1000.
# Step over the dark-corrected images...
imageNumber = -1
for dcimage in DCImages:
  imageNumber += 1
  print "\tCreating Sum(N) spectra for dark corrected image", imageNumber + 1, "out of", str(len(DCImages)) + "..."
  # First, make a list of the local maxima in each image.
  frameMaxima = filters.maximum_filter(dcimage, LocalMaxNeighborhood)
  # This is a mask with the same shape as the image that has Boolean values for each pixel to
  # denote the presence of a local mximum or not.
  LocalMaxMask = (dcimage == frameMaxima)
  # This is a similar image mask to make the threshold cut.
  ThresholdMask = ((dcimage) > LocalMaxThreshold)
  # Now we apply the threshold cut to the local max mask image.
  LocalMaxMask[ThresholdMask == False] = False
  # Now apply a unique label to every feature in LocalMaxMask, so that we can grab their
  # coordinates in the next step.
  LabeledMask, nFeatures = ndimage.label(LocalMaxMask)
  # Get the coordinates of each local maximum pixel.
  LocalMaxCoords = np.array(ndimage.center_of_mass(dcimage, LabeledMask, range(1, nFeatures + 1)))
  # Step over the local maximum coordinates, and get rid of the ones that are too close to the edge of the image.
  EdgeBoundary = 1 # Just cut away the outer 1 pixel.  We can fiddle with this later if we have to...
  EdgeRows = []
  for i in range(len(LocalMaxCoords)):
    EdgeHit = False
    if((LocalMaxCoords[i, 0] < EdgeBoundary) or (LocalMaxCoords[i, 0] >= (nPixelsY - EdgeBoundary))): EdgeHit = True
    if((LocalMaxCoords[i, 1] < EdgeBoundary) or (LocalMaxCoords[i, 1] >= (nPixelsY - EdgeBoundary))): EdgeHit = True
    if(EdgeHit): EdgeRows.append(i)
  LocalMaxCoords = np.delete(LocalMaxCoords, EdgeRows, 0)
  # Plot the dark corrected image with the coordinates of the local maxima marked on them.
  thisDCImageName = "DC_TEAMimage_" + str(imageNumber + 1)
  thisDCImageTitle = "Dark-Corrected Image from TEAM Detector at Time Stamp: " + str(ImagesInThisFile[2*imageNumber][0][1])
  plt.imshow(dcimage, alpha=0.75, aspect='auto', origin='lower', extent=[0.,lSensorX, 0.,lSensorY], interpolation='none')
  plt.colorbar()
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.title(thisDCImageTitle)
  plt.grid(True)
  ImagePlotFilePath = DCOutputDir + "/" + thisDCImageName + ".png"
  if os.path.exists(ImagePlotFilePath):
    print "Deleting old version of", ImagePlotFilePath
    os.system("rm " + ImagePlotFilePath)
  plt.savefig(ImagePlotFilePath)
  # Scale the list of x,y coordinates so that they can be plotted over the existing image, and then do so.
  LocalMaxDisplayX = LocalMaxCoords[:, 1] * ((lSensorX) / dcimage.shape[1]) + (0.5 * lPixelX)
  LocalMaxDisplayY = LocalMaxCoords[:, 0] * ((lSensorY) / dcimage.shape[0]) + (0.5 * lPixelY)
  plt.plot(LocalMaxDisplayX, LocalMaxDisplayY, 'ro')
  # Now, save the newly reconstructed dark corredted image to a png file.
  ImagePlotFilePath = DCOutputDir + "/" + thisDCImageName + ".Annotated.png"
  if os.path.exists(ImagePlotFilePath):
    print "Deleting old version of", ImagePlotFilePath
    os.system("rm " + ImagePlotFilePath)
  plt.savefig(ImagePlotFilePath)
  plt.clf()  
  # Get the values of the local maxima and the surrounding pixels so that we can build up the Sum(N) spectra.
  thisSum1Spectrum = []
  thisSum9Spectrum = []
  thisSum25Spectrum = []
  for coord in LocalMaxCoords:
    if(dcimage[coord[0]][coord[1]] > LocalMaxThreshold):
      thisSum1Spectrum.append(dcimage[coord[0]][coord[1]])
      # Build up the Sum(9) spectrum by collecting the values of the neighboring pixels.
      thisSum9Val = 0.
      for i in range(-1, 2):
        for j in range(-1, 2):
          # Make sure we haven't stepped off the edge of the image...
          if(((coord[0] + i) >= 0) and ((coord[0] + i) < nPixelsY) and 
             ((coord[1] + j) >= 0) and ((coord[1] + j) < nPixelsX) and 
             (dcimage[coord[0] + i][coord[1] + j] > SumNThreshold)): thisSum9Val += dcimage[coord[0] + i][coord[1] + j]
      thisSum9Spectrum.append(thisSum9Val)
      # And now the Sum(25) spectrum and the next-to-neighboring pixels.
      thisSum25Val = 0.
      for i in range(-2, 3):
        for j in range(-2, 3):
          # Make sure we haven't stepped off the edge of the image...
          if(((coord[0] + i) >= 0) and ((coord[0] + i) < nPixelsY) and 
             ((coord[1] + j) >= 0) and ((coord[1] + j) < nPixelsX) and 
             (dcimage[coord[0] + i][coord[1] + j] > SumNThreshold)): thisSum25Val += dcimage[coord[0] + i][coord[1] + j]
      thisSum25Spectrum.append(thisSum25Val)
  Sum01Spectra.append(thisSum1Spectrum)
  Sum09Spectra.append(thisSum9Spectrum)
  Sum25Spectra.append(thisSum25Spectrum)

# Now, make some histograms
xLo, xHi, xStep = -500., 9500., 10.
xBins = np.arange(xLo, xHi + xStep, xStep)
xBinCenters = [x + (0.5 * xStep) for x in xBins[:-1]]
xAxisTitle, yAxisTitle = 'Background Corrected ADC Value', 'Counts per ' + str(xStep) + ' ADC Unit Bin'
OutputFile.create_dataset('xHistoAxisTitle', data=[xAxisTitle])
OutputFile.create_dataset('yHistoAxisTitle', data=[yAxisTitle])
for imageNumber in range(nImagesPerFile / 2):
  # Sum(1)
  thisPlotTitle = 'Sum(1) Spectrum from TEAM Detector at Time Stamp: ' + str(ImagesInThisFile[imageNumber][0][1])
  ImagePlotFilePath = SSOutputDir + "/Sum01Spectrum" + str(imageNumber) + ".pdf"
  Sum01Vals = PythonTools.PlotHistogram(Sum01Spectra[imageNumber], xBins, 'r', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
  OutputFile.create_dataset('Sum01HistoVals_' + str(imageNumber), data=Sum01Vals)
  # Sum(9)
  thisPlotTitle = 'Sum(9) Spectrum from TEAM Detector at Time Stamp: ' + str(ImagesInThisFile[imageNumber][0][1])
  ImagePlotFilePath = SSOutputDir + "/Sum09Spectrum" + str(imageNumber) + ".pdf"
  Sum09Vals = PythonTools.PlotHistogram(Sum09Spectra[imageNumber], xBins, 'g', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
  OutputFile.create_dataset('Sum09HistoVals_' + str(imageNumber), data=Sum09Vals)
  # Sum(25)
  thisPlotTitle = 'Sum(25) Spectrum from TEAM Detector at Time Stamp: ' + str(ImagesInThisFile[imageNumber][0][1])
  ImagePlotFilePath = SSOutputDir + "/Sum25Spectrum" + str(imageNumber) + ".pdf"
  Sum25Vals = PythonTools.PlotHistogram(Sum25Spectra[imageNumber], xBins, 'b', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
  OutputFile.create_dataset('Sum25HistoVals_' + str(imageNumber), data=Sum25Vals)

# It's probably a good idea to sum up the different sum spectra over each frame.
SummedSum01Spectrum, SummedSum09Spectrum, SummedSum25Spectrum = [], [], []
for sum1spec  in Sum01Spectra: SummedSum01Spectrum += sum1spec
for sum9spec  in Sum09Spectra: SummedSum09Spectrum += sum9spec
for sum25spec in Sum25Spectra: SummedSum25Spectrum += sum25spec
# Summed Sum(1)
thisPlotTitle = 'Sum(1) Spectrum from TEAM Detector Summed Over all frames'
ImagePlotFilePath = OutputDir + "/SummedSum01Spectrum.pdf"
SummedSum01Vals = PythonTools.PlotHistogram(SummedSum01Spectrum, xBins, 'r', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
OutputFile.create_dataset('SummedSum01Vals', data=SummedSum01Vals)
# Summed Sum(9)
thisPlotTitle = 'Sum(9) Spectrum from TEAM Detector Summed Over all frames'
ImagePlotFilePath = OutputDir + "/SummedSum09Spectrum.pdf"
SummedSum09Vals = PythonTools.PlotHistogram(SummedSum09Spectrum, xBins, 'g', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
OutputFile.create_dataset('SummedSum09Vals', data=SummedSum09Vals)
# Summed Sum(25)
thisPlotTitle = 'Sum(25) Spectrum from TEAM Detector Summed Over all frames'
ImagePlotFilePath = OutputDir + "/SummedSum25Spectrum.pdf"
SummedSum25Vals = PythonTools.PlotHistogram(SummedSum25Spectrum, xBins, 'b', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
OutputFile.create_dataset('SummedSum25Vals', data=SummedSum25Vals)

# Close the hdf5 file...
OutputFile.close()

# Get the end time and report how long this calculation took
StopTime = time.time()
print "It took", StopTime - StartTime, "seconds for this code to run."
exit()
