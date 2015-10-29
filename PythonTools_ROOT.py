#!/usr/bin/python

###################################################################################################
# Support functions for TEAM detector file analysis.                                              #
###################################################################################################

import ROOT
import array
import numpy

def MakeTH2DObject(npixelsx, npixelsy, lpixelx, lpixely):
  xlo = -0.5 * lpixelx
  xhi = (float(npixelsx) - 0.5) * lpixelx
  ylo = -0.5 * lpixely
  yhi = (float(npixelsy) - 0.5) * lpixely
  thisTH2D = ROOT.TH2D("thisTH2D", "2D Histogram for a TEAM Detector Image",
                       npixelsx, xlo, xhi, npixelsy, ylo, yhi)
  thisTH2D.GetXaxis().SetTitle("x Position [mm]")
  thisTH2D.GetYaxis().SetTitle("y Position [mm]")
  thisTH2D.GetYaxis().SetTitleOffset(0.7)
  return thisTH2D

# Create a TH1D object for pixel value histograms of various sorts.
def MakePixValHisto(histoname, histotitle, nbins, xlo, xhi, color):
  AxisTitleSize = 0.05
  AxisTitleOffset = 0.7
  AxisLabelSize = 0.03
  PixValHisto = ROOT.TH1D(histoname, histotitle, nbins, xlo, xhi)
  PixValHisto.SetLineColor(color)
  PixValHisto.GetXaxis().SetTitle("Background Corrected ADC Value")
  PixValHisto.GetXaxis().SetTitleSize(AxisTitleSize)
  PixValHisto.GetXaxis().SetTitleOffset(AxisTitleOffset)
  PixValHisto.GetXaxis().SetLabelSize(AxisLabelSize)
  BinWidth = (xhi - xlo) / float(nbins)
  TitleString = "Counts per " + "{:0.1f}".format(BinWidth) + " ADC Unit Bin"
  PixValHisto.GetYaxis().SetTitle(TitleString)
  PixValHisto.GetYaxis().SetTitleSize(AxisTitleSize)
  PixValHisto.GetYaxis().SetTitleOffset(1.0 * AxisTitleOffset)
  PixValHisto.GetYaxis().SetLabelSize(AxisLabelSize)
  return PixValHisto

# Construct a peak model taken from RadWare, a tool often used in HPGe detector analysis...
def GetRWFitModel(fitmodelname, templatehisto, mean, sigm):
  # Set some basic parameter limits...
  MeanHalfWindow = 100.
  LoFrac =  0.5
  HiFrac =  1.5
  # Build up the peak model
  GausString = "([0] * exp(-0.5 * ((x - [1]) / [2])^2))"
  SkSfString = "([3] * TMath::Erfc((x - [1]) / (sqrt(2.) * [2])))"
  BkGdString = "[4] + ([5] * x) + ([6] * (x^2))"
  PeakModelString = GausString + " + " + SkSfString + " + " + BkGdString
  PeakModel = ROOT.TF1(fitmodelname, PeakModelString, 
                       templatehisto.GetXaxis().GetXmin(), templatehisto.GetXaxis().GetXmax())
  PeakModel.SetLineColor(templatehisto.GetLineColor())
  PeakModel.SetLineStyle(2)
  PeakModel.SetLineWidth(4)
  # Calculate the initial guesses for the model parameters from the histogram
  templatehisto.GetXaxis().SetRangeUser(mean - (5. * sigm), mean + (5. * sigm))
  SpecMax = templatehisto.GetMaximum()
  templatehisto.GetXaxis().UnZoom()
  # Set up the peak model...
  # Gaussian bit:
  PeakModel.SetParName(  0, "Gaus. Nor.")
  PeakModel.SetParLimits(0, 0., 2. * SpecMax)
  PeakModel.SetParameter(0, SpecMax)
  PeakModel.SetParName(  1, "Peak Mean")
  PeakModel.SetParLimits(1, mean - MeanHalfWindow, mean + MeanHalfWindow)
  PeakModel.SetParameter(1, mean)
  PeakModel.SetParName(  2, "Gaus. sig.")
  PeakModel.SetParLimits(2, LoFrac * sigm, HiFrac * sigm)
  PeakModel.SetParameter(2, sigm)
  # Sigmoid function:
  PeakModel.SetParName(  3, "SF Nor.")
  PeakModel.SetParLimits(3, 0., 0.5 * SpecMax)
  PeakModel.SetParameter(3, 0.)
  # Polynomial background:
  PeakModel.SetParName(  4, "BG Cnst.")
  PeakModel.SetParLimits(4, 0., 2. * SpecMax)
  PeakModel.SetParameter(4, templatehisto.GetBinContent(templatehisto.FindBin(mean - 10.)))
  PeakModel.SetParName(  5, "BG Lin.")
  PeakModel.FixParameter(5, 0.)
  PeakModel.SetParName(  6, "BG Quad.")
  PeakModel.FixParameter(6, 0.)
  return PeakModel

# And get the components of the RadWare peak model...
def GetRWFitModelComponents(fitmodel, plotlo, plothi):
  # Peak model components...
  GausString = "([0] * exp(-0.5 * ((x - [1]) / [2])^2))"
  SkSfString = "([3] * TMath::Erfc((x - [1]) / (sqrt(2.) * [2])))"
  BkGdString = "[4] + ([5] * x) + ([6] * (x^2))"
  # Isolate the Gaussian component of the fit model
  GausModel = ROOT.TF1("GausModel", GausString, plotlo, plothi)
  GausModel.SetTitle("Gaussian Peak")
  GausModel.SetLineColor(fitmodel.GetLineColor() - 1)
  GausModel.SetLineWidth(fitmodel.GetLineWidth())
  GausModel.SetLineStyle(fitmodel.GetLineStyle())
  GausModel.FixParameter(0, fitmodel.GetParameter(0))
  GausModel.FixParameter(1, fitmodel.GetParameter(1))
  GausModel.FixParameter(2, fitmodel.GetParameter(2))
  # Isolate the sigmoid component
  SgmdModel = ROOT.TF1("SkGsModel", SkSfString, plotlo, plothi)
  SgmdModel.SetTitle("Sigmoid Fcn.")
  SgmdModel.SetLineColor(fitmodel.GetLineColor() + 1)
  SgmdModel.SetLineWidth(fitmodel.GetLineWidth())
  SgmdModel.SetLineStyle(fitmodel.GetLineStyle())
  SgmdModel.FixParameter(3, fitmodel.GetParameter(3))
  SgmdModel.FixParameter(1, fitmodel.GetParameter(1))
  SgmdModel.FixParameter(2, fitmodel.GetParameter(2))
  # Isolate the sigmoid component
  BkGdModel = ROOT.TF1("BkGdModel", BkGdString, plotlo, plothi)
  BkGdModel.SetTitle("Pol. Bg.")
  BkGdModel.SetLineColor(fitmodel.GetLineColor() + 2)
  BkGdModel.SetLineWidth(fitmodel.GetLineWidth())
  BkGdModel.SetLineStyle(fitmodel.GetLineStyle())
  BkGdModel.FixParameter(4, fitmodel.GetParameter(4))
  BkGdModel.FixParameter(5, fitmodel.GetParameter(5))
  BkGdModel.FixParameter(6, fitmodel.GetParameter(6))
  return [GausModel, SgmdModel, BkGdModel]

# Create an annotation to display the parameters from the Radware fit
def MakeFitAnnotationRW(fitmodel):
  thisChi2 = fitmodel.GetChisquare()
  thisNDF  = fitmodel.GetNDF()
  if(thisNDF == 0): thisNDF = 1
  thisPVal = fitmodel.GetProb()
  Gnor = fitmodel.GetParameter(0)
  GnEr = fitmodel.GetParError(0)
  Mean = fitmodel.GetParameter(1)
  MeEr = fitmodel.GetParError(1)
  Sigm = fitmodel.GetParameter(2)
  SiEr = fitmodel.GetParError(2)
  SmdN = fitmodel.GetParameter(3)
  SNEr = fitmodel.GetParError(3)
  AnnotationLeft  = 0.672
  AnnotationRight = 0.972
  AnnotationTop   = 0.915
  AnnotationBottom = 0.515
  thisAnnotation = ROOT.TPaveText(AnnotationLeft,AnnotationBottom,AnnotationRight,AnnotationTop,"blNDC")
  thisAnnotation.SetTextFont(42)
  thisAnnotation.SetName(fitmodel.GetName() + "_AnnotationText")
  thisAnnotation.SetBorderSize(1)
  thisAnnotation.SetFillColor(ROOT.kWhite)
  ThisLine = "#chi^{2} per DoF = " + "{:6.1f}".format(thisChi2) + " / " + str(thisNDF) + " = " + "{:6.2f}".format(thisChi2 / float(thisNDF))
  thisAnnotation.AddText(ThisLine)
  ThisLine = "(Probability = " + "{:2.6f}".format(thisPVal) + ")"
  thisAnnotation.AddText(ThisLine)
  ThisLine = "Mean = "    + "{:6.2f}".format(Mean) + " #pm " + "{:1.2f}".format(MeEr)
  thisAnnotation.AddText(ThisLine)
  ThisLine = "Sigma =  "  + "{:6.2f}".format(Sigm) + " #pm " + "{:1.2f}".format(SiEr)
  thisAnnotation.AddText(ThisLine)
  ThisLine = "Gaus. Norm. = "   + "{:6.2f}".format(Gnor) + " #pm " + "{:1.2f}".format(GnEr)
  thisAnnotation.AddText(ThisLine)
  ThisLine = "Sigmoid Norm. = "   + "{:6.2f}".format(SmdN) + " #pm " + "{:1.2f}".format(SNEr)
  thisAnnotation.AddText(ThisLine)
  return thisAnnotation
