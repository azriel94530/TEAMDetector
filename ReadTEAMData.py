#!/usr/bin/python

####################################################################################################
# Open a binary file file from the TEAM detector running in either trigger mode (2 or 3, selected  #
# by command line argument).  Then read and decode it using the the data file using the struct     #
# package, and then push it into a series of 2D histograms using imshow.  Once we have all these   #
# images in a more python-friendly format, then perform a dark correction based on which trigger   #
# mode was selected.  Once that's done, we'll go ahead and save a png image of the raw and dark    #
# corrected images, and then write the image array to an HDF5 file.                                #
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
from sklearn.cluster import MeanShift, estimate_bandwidth
from itertools import cycle
import PythonTools

####################################
#  BEGIN MAIN BODY OF THE CODE!!!  #
####################################

# Get the start time of this calculation
StartTime = time.time()

# Set some flags for how verbose our input and output are going to be.
Debugging = False
VerboseProcessing = True

if(len(sys.argv) != 3):
  print "Usage: [python] ReadTEAMData_mode3.py path/to/binary/image/file TriggerMode"
  exit()

# Pull in the path to the binary file we're going to look at.
InputFilePath = sys.argv[1]
if(VerboseProcessing):
  print "\tReading in: '" + InputFilePath + "' for analysis."

# Check that the trigger mode we read in makes sense.
TriggerMode = int(sys.argv[2])
print "\tTrigger mode set to", TriggerMode
if(TriggerMode == 2):
  print "\tWill construct dark image from the median of each pixel value across the whole data set."
elif(TriggerMode == 3):
  print "\tWill use even image numbers (starting from zero) to dark correct odd numbered exposures."
else:
  print "\tTrigger mode needs to be either a \'2\' or a \'3\'.  It seems to be something else..."
  exit()

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
nADCchannels = 16
xPixelsPerReadout = nPixelsX / nADCchannels
lSensorX = 10.#[mm]
lSensorY = 10.#[mm]
lPixelX = lSensorX / float(nPixelsX)#Sensor is 10 mm x 10 mm
lPixelY = lSensorY / float(nPixelsY)
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
  print "\tReading in image", imageNumber + 1, "out of", str(nImagesPerFile) + "..."
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
if(TriggerMode == 2): 
  DCImages, MedianDarkImage = PythonTools.DarkCorrectMode2(ImagesInThisFile)
  MedianDarkImageName = "MedianDarkImage"
  MedianDarkImageTitle = "Median Dark Image for " + InputFilePath
  OutputFile.create_dataset('DarkImage', data=MedianDarkImage)
  plt.imshow(MedianDarkImage, alpha=0.75, aspect='auto', origin='lower', extent=[0.,lSensorX, 0.,lSensorY], interpolation='none')
  plt.colorbar()
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  plt.title(MedianDarkImageTitle)
  plt.grid(True)
  DarkImageFilePath = OutputDir + "/" + MedianDarkImageName + ".png"
  plt.savefig(DarkImageFilePath)
  plt.clf()

if(TriggerMode == 3): DCImages = PythonTools.DarkCorrectMode3(ImagesInThisFile)

# Now that we have dark-corrected images, let's create the calorimetric spectra: Sum01, Sum09, and
# Sum25.
Sum01Spectra = []
Sum09Spectra = []
Sum25Spectra = []
LocalMaxThreshold = 80.#[ADC Counts]
LocalMaxNeighborhood = 2
SumNThreshold = 40.
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
  #################################################################################################
  # DATA QUALITY CUT NUMBER 1:                                                                    #
  # Step over the local maximum coordinates, and get rid of the ones that are too close to the    #
  # edge of the image.                                                                            #
  #################################################################################################
  EdgeBoundary = 2 # Just cut away the outer 2 pixels.  We can fiddle with this later if we have to...
  EdgePixelIndecies = []
  for i in range(len(LocalMaxCoords)):
    EdgeHit = False
    if((LocalMaxCoords[i, 0] < EdgeBoundary) or (LocalMaxCoords[i, 0] >= (nPixelsY - EdgeBoundary))): EdgeHit = True
    if((LocalMaxCoords[i, 1] < EdgeBoundary) or (LocalMaxCoords[i, 1] >= (nPixelsY - EdgeBoundary))): EdgeHit = True
    if(EdgeHit): EdgePixelIndecies.append(i)
  print "\tDeleted", len(EdgePixelIndecies), "hit pixels with the edge cut."
  LocalMaxCoords = np.delete(LocalMaxCoords, EdgePixelIndecies, 0)
  #################################################################################################
  # DATA QUALITY CUT NUMBER 2:                                                                    #
  # Step over the remaining local maximum coordinates and see if there is a hit above threshold   #
  # in in a different ADC channel at the same time to remove noise events.  This is done by       #
  # in every pixel that is a multiple of xPixelsPerReadout (defined in line 70 at time of writing)#
  # and checking if that is above threshold.                                                      #
  #################################################################################################
  NoiseCutPixelIndecies = []
  for i in range(len(LocalMaxCoords)):
    NoiseCut = False
    for NoiseCheckPixel in range(xPixelsPerReadout, nPixelsX, xPixelsPerReadout):
      CheckColumn = LocalMaxCoords[i, 1] + NoiseCheckPixel
      if(CheckColumn >= nPixelsX): CheckColumn -= nPixelsX
      if(dcimage[LocalMaxCoords[i, 0], CheckColumn] > (0.5 * LocalMaxThreshold)): NoiseCut = True
    if(NoiseCut): NoiseCutPixelIndecies.append(i)
  print "\tDeleted", len(NoiseCutPixelIndecies), "hit pixels with the nosie cut."
  LocalMaxCoords = np.delete(LocalMaxCoords, NoiseCutPixelIndecies, 0)
  # Plot the dark corrected image with the coordinates of the local maxima marked on them.
  thisDCImageName = "DC_TEAMimage_" + str(imageNumber + 1)
  thisDCImageTitle = "Dark-Corrected Image from TEAM Detector at Time Stamp: " + str(dcimage[0][1])
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
  #print LocalMaxCoords
  LocalMaxDisplayX = LocalMaxCoords[:, 1] * ((lSensorX) / dcimage.shape[1]) + (0.5 * lPixelX)
  LocalMaxDisplayY = LocalMaxCoords[:, 0] * ((lSensorY) / dcimage.shape[0]) + (0.5 * lPixelY)
  plt.plot(LocalMaxDisplayX, LocalMaxDisplayY, 'ro')
  plt.title("Found " + str(len(LocalMaxDisplayX)) + " Local Maxima")
  # Now, save the newly reconstructed dark corredted image to a png file.
  ImagePlotFilePath = DCOutputDir + "/" + thisDCImageName + ".Annotated.png"
  if os.path.exists(ImagePlotFilePath):
    print "Deleting old version of", ImagePlotFilePath
    os.system("rm " + ImagePlotFilePath)
  plt.savefig(ImagePlotFilePath)
  plt.clf()
  #################################################################################################
  # See if we can get the sklearn clustering algorithm to do something good here...               #
  #################################################################################################
  # First, let's redraw the image we're working with
  plt.imshow(dcimage, alpha=0.75, aspect='auto', origin='lower', interpolation='none', extent=[0.,lSensorX, 0.,lSensorY])
  plt.colorbar()
  plt.xlabel(xAxisTitle)
  plt.ylabel(yAxisTitle)
  #plt.title(thisDCImageTitle)
  plt.grid(True)
  # Set some thresholds to decide which clusters we want to keep...
  ClusterIntegralThreshold = 1e4
  ClusterSizeThreshold = 100
  # Estimate the bandwidth of this image.  The "quantile" argument changes how sensitive the
  # clustering algorithm is to fainter clusters.  Leave this around 0.3-0.5 to find the beam and
  # diffraction spots.  Dial it down to 0.1 to pick up a bunch of single electron scatters.  0.7
  # and higher seems to drop the diffraction spots, but they are still bright enough to drag the 
  # remaining cluster off of the main beam spot.  We should stick with something around 0.4 for
  # now.  "n_samples" sets how long the bandwidth estimator looks around to when trying to decide
  # how granular the image is.  500 or so seems to be fine, but pushing up as high as 200 doesn't
  # seem to slow things down much for these file sizes.
  bandwidth = estimate_bandwidth(LocalMaxCoords, quantile=0.4, n_samples=500)
  print "\tThe bandwidth for image number " + str(imageNumber) + " was estemated at: " + str(bandwidth) + "."
  # Construct the MeanShift object.
  ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
  # This actually does the clustering.
  ms.fit(LocalMaxCoords)
  # Get the labels for each point.
  labels = ms.labels_
  # And pull up the cluster centers.
  cluster_centers = ms.cluster_centers_
  # Figure out how many unique clusters we tracked down.
  labels_unique = np.unique(labels)
  n_clusters_ = len(labels_unique)
  print"\tFound " + str(n_clusters_), "clusters in image", imageNumber
  RealClusters = []
  colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
  for k, col in zip(range(n_clusters_), colors):
    cluster_center = cluster_centers[k]
    my_members = labels == k
    #plt.plot(LocalMaxCoords[my_members, 0], LocalMaxCoords[my_members, 1], col + '.')
    ClusterIntegral = dcimage[LocalMaxCoords[my_members, 0].astype(int), LocalMaxCoords[my_members, 1].astype(int)].sum()
    nPixelsInCluster = dcimage[LocalMaxCoords[my_members, 0].astype(int), LocalMaxCoords[my_members, 1].astype(int)].size
    GoodCluster = False
    if((ClusterIntegral > ClusterIntegralThreshold) and (nPixelsInCluster > ClusterSizeThreshold)): GoodCluster = True
    print "\t\tIntegral of pixel values in cluster", k, "is", ClusterIntegral, "in", nPixelsInCluster, "pixels."
    if GoodCluster:
      print "\t\t\tWe're keeping this one."
      ClusterDisplayX = cluster_center[1] * ((lSensorX) / dcimage.shape[1]) + (0.5 * lPixelX)
      ClusterDisplayY = cluster_center[0] * ((lSensorY) / dcimage.shape[0]) + (0.5 * lPixelY)
      RealClusters.append([ClusterDisplayX, ClusterDisplayY])
      plt.plot(ClusterDisplayX, ClusterDisplayY, 'o', markerfacecolor=col, markeredgecolor='k', markersize=14)
    else:
      print "\t\t\tAnd we're discarding it."
  print"\tAnd", len(RealClusters), "of those", str(n_clusters_), "passed our cuts..."
  plt.title("Found " + str(len(RealClusters)) + " Clusters")
  ImagePlotFilePath = DCOutputDir + "/" + thisDCImageName + ".Clusters.png"
  plt.savefig(ImagePlotFilePath)
  plt.clf()
  #################################################################################################
  # Done screwing around with clustering...                                                       #
  #################################################################################################
  # Get the values of the local maxima and the surrounding pixels so that we can build up the Sum(N) spectra.
  thisSum1Spectrum = []
  thisSum9Spectrum = []
  thisSum25Spectrum = []
  for coord in LocalMaxCoords:
    if(dcimage[coord[0]][coord[1]] > LocalMaxThreshold):
      thisSum1Spectrum.append(dcimage[coord[0]][coord[1]])
      # Build up the Sum(9) spectrum by collecting the values of the neighboring pixels.
      thisSum9Val = 0.
      GoodCluster = False
      # Step over the 3x3 cluster and make some quality checks: make sure that at least one of the neighboring pixels is above the SumNThreshold
      for i in range(-1, 2):
        for j in range(-1, 2):
          if(dcimage[coord[0] + i][coord[1] + j] > SumNThreshold): GoodCluster = True
      # If we still have a good cluster, then add up the neighborhood and append it to the histogram queue.
      if(GoodCluster):
        for i in range(-1, 2):
          for j in range(-1, 2):
            thisSum9Val += dcimage[coord[0] + i][coord[1] + j]
      thisSum9Spectrum.append(thisSum9Val)
      # And now the Sum(25) spectrum and the next-to-neighboring pixels.
      thisSum25Val = 0.
      GoodCluster = False
      # Step over the 5x5 cluster and make some quality checks: make sure that at least one of the neighboring pixels is above the SumNThreshold
      for i in range(-2, 3):
        for j in range(-2, 3):
          if(dcimage[coord[0] + i][coord[1] + j] > SumNThreshold): GoodCluster = True
      # If we still have a good cluster, then add up the neighborhood and append it to the histogram queue.
      if(GoodCluster):
        for i in range(-2, 3):
          for j in range(-2, 3):
            thisSum25Val += dcimage[coord[0] + i][coord[1] + j]
      thisSum25Spectrum.append(thisSum25Val)
  Sum01Spectra.append(thisSum1Spectrum)
  Sum09Spectra.append(thisSum9Spectrum)
  Sum25Spectra.append(thisSum25Spectrum)

# Now, make some histograms
xLo, xHi, xStep = -500., 3950., 10.
OutputFile.create_dataset('xAxisParams',data=[xLo, xHi, xStep])
xBins = np.arange(xLo, xHi + xStep, xStep)
xBinCenters = [x + (0.5 * xStep) for x in xBins[:-1]]
xAxisTitle, yAxisTitle = 'Background Corrected ADC Value', 'Counts per ' + str(xStep) + ' ADC Unit Bin'
OutputFile.create_dataset('xHistoAxisTitle', data=[xAxisTitle])
OutputFile.create_dataset('yHistoAxisTitle', data=[yAxisTitle])
for imageNumber in range(nImagesPerFile / 2):
  print "\tCreating sum spectra for image", imageNumber
  # Sum(1)
  thisPlotTitle = 'Sum(1) Spectrum from TEAM Detector at Time Stamp: ' + str(ImagesInThisFile[imageNumber][0][1])
  ImagePlotFilePath = SSOutputDir + "/Sum01Spectrum" + str(imageNumber) + ".pdf"
  if(len(Sum01Spectra[imageNumber]) == 0): Sum01Spectra[imageNumber].append(0.)
  Sum01Vals = PythonTools.PlotHistogram(Sum01Spectra[imageNumber], xBins, 'r', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
  OutputFile.create_dataset('Sum01HistoVals_' + str(imageNumber), data=Sum01Vals)
  # Sum(9)
  thisPlotTitle = 'Sum(9) Spectrum from TEAM Detector at Time Stamp: ' + str(ImagesInThisFile[imageNumber][0][1])
  ImagePlotFilePath = SSOutputDir + "/Sum09Spectrum" + str(imageNumber) + ".pdf"
  if(len(Sum09Spectra[imageNumber]) == 0): Sum09Spectra[imageNumber].append(0.)
  Sum09Vals = PythonTools.PlotHistogram(Sum09Spectra[imageNumber], xBins, 'g', thisPlotTitle, xAxisTitle, yAxisTitle, ImagePlotFilePath)
  OutputFile.create_dataset('Sum09HistoVals_' + str(imageNumber), data=Sum09Vals)
  # Sum(25)
  thisPlotTitle = 'Sum(25) Spectrum from TEAM Detector at Time Stamp: ' + str(ImagesInThisFile[imageNumber][0][1])
  ImagePlotFilePath = SSOutputDir + "/Sum25Spectrum" + str(imageNumber) + ".pdf"
  if(len(Sum25Spectra[imageNumber]) == 0): Sum25Spectra[imageNumber].append(0.)
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
