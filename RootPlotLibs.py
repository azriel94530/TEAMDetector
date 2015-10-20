####################################################################################################
# Some commonly used functions to make my root plots beautiful                                     #
####################################################################################################

# Header, import statements etc.
import ROOT

# Create a simple TCanvas object
def ASimpleCanvas():
  ACanvas = ROOT.TCanvas("ACanvas","Look!  It's a Canvas!",0,0,1024,768)
  ACanvas.SetHighLightColor(ROOT.kWhite)
  ACanvas.SetBorderSize(0)
  ACanvas.SetGridx(0)
  ACanvas.SetGridy(0)
  ACanvas.SetFrameFillColor(ROOT.kWhite)
  ACanvas.SetFillColor(ROOT.kWhite)
  return ACanvas

# Create a TCanvas object for two-panel plots
def ATwoPanelCanvas():
  ACanvas = ROOT.TCanvas("ACanvas","Look!  It's a Canvas!",0,0,1024,371)
  ACanvas.SetHighLightColor(ROOT.kWhite)
  ACanvas.SetBorderSize(0)
  ACanvas.SetGridx(0)
  ACanvas.SetGridy(0)
  ACanvas.SetFrameFillColor(ROOT.kWhite)
  ACanvas.SetFillColor(ROOT.kWhite)
  return ACanvas

# Create a TCanvas object for 720p HD plots
def A720pCanvas():
  ACanvas = ROOT.TCanvas("ACanvas","Look!  It's a Canvas!",0,0,1280,720)
  ACanvas.SetHighLightColor(ROOT.kWhite)
  ACanvas.SetBorderSize(0)
  ACanvas.SetGridx(0)
  ACanvas.SetGridy(0)
  ACanvas.SetFrameFillColor(ROOT.kWhite)
  ACanvas.SetFillColor(ROOT.kWhite)
  return ACanvas

# Create a TCanvas object for 720p HD plots
def A1080pCanvas():
  ACanvas = ROOT.TCanvas("ACanvas","Look!  It's a Canvas!",0,0,1920,1080)
  ACanvas.SetHighLightColor(ROOT.kWhite)
  ACanvas.SetBorderSize(0)
  ACanvas.SetGridx(0)
  ACanvas.SetGridy(0)
  ACanvas.SetFrameFillColor(ROOT.kWhite)
  ACanvas.SetFillColor(ROOT.kWhite)
  return ACanvas

# Create a simple TPad object
def ASimplePad():
  APad = ROOT.TPad("APad", "It's a Pad", 0.,0.,1.,1., ROOT.kWhite, 0,0)
  APad.SetTopMargin(0.01)
  APad.SetLeftMargin(0.06)
  APad.SetBottomMargin(0.08)
  APad.SetRightMargin(0.03)
  APad.SetLogx(0)
  APad.SetLogy(0)
  APad.SetLogz(0)
  APad.SetGridx(1)
  APad.SetGridy(1)

  return APad

# Create the left hand TPad object for a two panel plot
def TwoPanelPad_Left():
  LeftPad = ROOT.TPad("LeftPad", "It's a Pad", 0.,0.,0.5,1., ROOT.kWhite, 0,0)
  LeftPad.SetTopMargin(0.08)
  LeftPad.SetLeftMargin(0.12)
  LeftPad.SetBottomMargin(0.10)
  LeftPad.SetRightMargin(0.01)
  LeftPad.SetLogx(0)
  LeftPad.SetLogy(0)
  LeftPad.SetLogz(0)
  LeftPad.SetGridx(1)
  LeftPad.SetGridy(1)

  return LeftPad

# Create the right hand TPad object for a two panel plot
def TwoPanelPad_Right():
  RightPad = ROOT.TPad("RightPad", "It's another Pad", 0.5,0.,1.,1., ROOT.kWhite, 0,0)
  RightPad.SetTopMargin(0.08)
  RightPad.SetLeftMargin(0.12)
  RightPad.SetBottomMargin(0.10)
  RightPad.SetRightMargin(0.01)
  RightPad.SetLogx(0)
  RightPad.SetLogy(0)
  RightPad.SetLogz(0)
  RightPad.SetGridx(1)
  RightPad.SetGridy(1)

  return RightPad

# Adjust the binning of a TH1D so that it matches the spacing in a TGraph object
def AdjustTH1DBinningToTGraph(Histo, Graph, XLo, XHi):
  # Histo: TH1D object to be adjusted
  # Graph: TGraph object used as a template
  # XLo, XHi: Range of Histo over which we will preserve the adjusted histogram's integral

  # An array to hold the low edges of the bins in the adjusted histogram
  AdjustedHistoLoEdges = []
  
  # Set up the bin low edges for the new histogram
  PointSpacing = 0.
  for i in range(Graph.GetN()):
    if (i == 0):
      PointSpacing = float(Graph.GetX()[1] - Graph.GetX()[0])
    else:
      PointSpacing = float(Graph.GetX()[i] - Graph.GetX()[i - 1])
    AdjustedHistoLoEdges.append(Graph.GetX()[i] - (0.5 * PointSpacing))
  AdjustedHistoLoEdges = array.array("f", AdjustedHistoLoEdges)
  
  # An new name for the adjusted histogram so that ROOT doesn't get scared about memory leaks
  AdjustedHistoName = Histo.GetName() + "_Adjusted"

  # Create the TH1D object
  AdjustedHisto = ROOT.TH1D(AdjustedHistoName, Histo.GetTitle(), Graph.GetN() - 1, AdjustedHistoLoEdges)
  
  # Step over the points in Graph and set the value matching bin in AdjustedHisto to the corresponding bin in Histo
  for i in range(Graph.GetN()):
    AdjustedHisto.SetBinContent(AdjustedHisto.FindBin(Graph.GetX()[i]), Histo.GetBinContent(Histo.FindBin(Graph.GetX()[i])))
  
  # Scale AdjustedHisto to preserve Histo's area
  #print "Scaling", Histo.GetName(), "by", Histo.Integral(Histo.FindBin(XLo), Histo.FindBin(XHi)) / AdjustedHisto.Integral(AdjustedHisto.FindBin(XLo), AdjustedHisto.FindBin(XHi)), "after rebinning"
  AdjustedHisto.Scale(Histo.Integral(Histo.FindBin(XLo), Histo.FindBin(XHi)) / AdjustedHisto.Integral(AdjustedHisto.FindBin(XLo), AdjustedHisto.FindBin(XHi)))
  #AdjustedHisto.Draw()
  #Histo.Draw("same")
  #Graph.Draw("samelp")
  #raw_input()
  
  return AdjustedHisto

# Set up some stylistic flourishes for nicer ROOT plots, and return ROOT TCanvas and TPad objects
# suitable for simple, single-panal plots.
def GetReadyToPlot():
  ROOT.gStyle.SetOptTitle(1)
  ROOT.gStyle.SetOptFit(0)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetPalette(1)
  # and create a canvas on which to plot things...
  canvas = A720pCanvas()
  # and a pad to set up the canvas...
  TopMargin    = 0.08
  LeftMargin   = 0.04
  BottomMargin = 0.09
  RightMargin  = 0.11
  pad = ASimplePad()
  pad.SetTopMargin(TopMargin)
  pad.SetLeftMargin(LeftMargin)
  pad.SetBottomMargin(BottomMargin)
  pad.SetRightMargin(RightMargin)
  pad.SetGridx(1)
  pad.SetGridy(1)
  return canvas, pad
