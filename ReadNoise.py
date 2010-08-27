"""-------------------------------------------------------------------------------------
Use this script to read in however many noise waveforms desirable and average over them.
Most useful to print power spectrum to file for quick reading later, since it takes
a long time to run through them otherwise.

Author: Greg Dooley
Contact: gdooley@princeton.edu
Date: August 2010
"""


from RawNoise import *
from NoiseAnalysis import *
from SimulateNoise import *
from ROOT import *
#import NoiseAnalysis

## number of waveforms to average over.
n = 1
rFile = "struck_run346.root"
rTree = "SIS3302Decoder"

a = RawNoise(n,rFile,rTree)

freq = a.GetXValues() ## extract frequency values from power spectrum
power = a.GetYValues() ##extract power from power spectrum

## uncomment to write to file
"""
fout = open("frequency2.dat", "w")
for i in range(freq.size):
	fout.write(str(freq[i]))
	fout.write("\n")
fout.close()

fout2 = open("power2.dat", "w")
for i in range(power.size):
	fout2.write(str(power[i]))
	fout2.write("\n")
fout2.close()
"""

# Create object for analyzing data.
b = NoiseAnalysis(power, freq, 0)
sTimes = array([4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8, 7e-8, 5e-8, 3e-8, 2e-8, 1.3e-8, 1e-8]) #good with analog filter
sTimes2 = array([2.5e-6,2e-6,1.4e-6,1e-6,9e-7,7e-7,4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8, 7e-8, 5e-8, 3e-8]) # good with trap filter
sTimes3 = array([4e-6,2.5e-6,2e-6,1.4e-6,1e-6,9e-7,7e-7,4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8]) # Good for hand interpolated data
sTimes4 = array([5e-6, 1e-6, 9e-7, 4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8, 7e-8, 5e-8, 3e-8, 2e-8, 1.3e-8, 1e-8, 9e-9, 7e-9, 3e-9, 1e-9, 7e-10])
guessA = 1.8e12
guessB = 3.21e11
guessC = 1.49e11


"""
hist = a.GimmeCumFTHist()
c1 = TCanvas("c1", "Signal Analysis", 200,10,600,480)
#pad = TPad("pad2", "Frequency Domain", .05, .02,.95,.98)
gPad.SetFillColor(0)
gPad.SetLogx()
gPad.SetLogy()
hist.Draw()
c1.Update()
raw_input("Press Enter to Close")
#b = NoiseAnalysis(hist2)
#b.PlotNoise()
"""


"""
arr = hist2.GetArray()
print arr[10]
print hist2.GetBinContent(10)
print hist2.GetBinWidth(0)
print hist2.GetNbinsX()
print hist2.GetNbinsX()*hist2.GetBinWidth(0)

c1 = TCanvas("c1", "Signal Analysis", 200,10,600,480)
pad1 = TPad("pad1", "Time Domain", .05,.5,.95,.97)
pad2 = TPad("pad2", "Frequency Domain", .05, .02,.95,.47)
pad1.SetFillColor(0)
pad2.SetFillColor(0)
pad1.Draw()
pad2.Draw()
pad1.cd()
hist.Draw()
pad2.cd()
hist2.Draw()
"""