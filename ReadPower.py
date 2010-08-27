"""-------------------------------------------------------------------------------------
Read power spectrum data from file and create NoiseAnalysis object for analyzing.
Run within python shell to access all methods dynamically.

Author: Greg Dooley
Contact: gdooley@princeton.edu
Date: August 2010
"""


from RawNoise import *
from NoiseAnalysis import *
from SimulateNoise import *
from ROOT import *

## Read data
power = zeros(17233)
fin = open("power.dat", "r")
for i in range(17233):
	x = fin.readline()
	x = float(x)
	power[i] = x
fin.close()

freq = zeros(17233)
fin2 = open("frequency.dat", "r")
for i in range(17233):
	x = fin2.readline()
	x = float(x)
	freq[i] = x
fin2.close()
print power
print freq

## Create object
a = RawNoise(1)
b = NoiseAnalysis(power, freq, 1)
print power
print freq
sTimes = array([4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8, 7e-8, 5e-8, 3e-8, 2e-8, 1.3e-8, 1e-8])
sTimes2 = array([2.5e-6,2e-6,1.4e-6,1e-6,9e-7,7e-7,4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8, 7e-8, 5e-8, 3e-8]) # good with trap filter
sTimes3 = array([4e-6,2.5e-6,2e-6,1.4e-6,1e-6,9e-7,7e-7,4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8]) # Good for hand interpolated data
sTimes4 = array([5e-6, 1e-6, 9e-7, 4e-7, 3e-7,2e-7, 1.3e-7, 1e-7, 9e-8, 7e-8, 5e-8, 3e-8, 2e-8, 1.3e-8, 1e-8, 9e-9, 7e-9, 3e-9, 1e-9, 7e-10])

## must guess fit parameters for FWHM^2 vs Peaking time curve
guessA = 1.8e12
guessB = 3.21e11
guessC = 1.49e11

##Some examples of using code
#b.PlotNoise()
#b.PlotFWHMvsShapingTime(sTimes2,guessA,guessB,guessC)
#samp = a.GetMaxFreq()*2.0*1e-9 #1e-9 converts to GHz
#a.PlotSimulatedNoise()

#  5e-7 and 1e-8 are good boundaries for having the analog shaper intersect the power spectrum.
#b.PlotFWHMvsShapingTime(sTimes)
#b.PlotNoise()
#b.PlotSteps(1e-7)
#b.PlotShaper(10e-9)
