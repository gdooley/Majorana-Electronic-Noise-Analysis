"""--------------------------------------------------------------------------------------
  Given a sample power spectrum from freq_input and noise_input, this simulates a noise
  pulse of a desired duration using a method described in Wan's thesis. Implements it
  both by a direct means, and using MGDO.

   Date:   August-2010
   Author: Greg Dooley
   Contact: gdooley@princeton.edu
---------------------------------------------------------------------------------------"""

import ROOT
from numpy import *
ROOT.gApplication.ExecuteFile("$MGDODIR/Majorana/LoadMGDOMJClasses.C")

## Constants to Set
Fc =  50000 #hertz - cut off frequency
To = 1 #seconds - length of signal to simulate
nbins = 2*Fc*To # number of interpolated data points to take
nbins = 25000
## Noise in V/sqrt(Hz) - power spectrum to sample to produce noise
freq_input = array([1.1, 30.0, 400., 3311., 39810., 707945., 19054607.17])
noise_input = array([pow(10.,-5.82), pow(10., -5.95), pow(10., -6.4), pow(10, -7.2), pow(10, -7.8), pow(10,-8.15), pow(10,-8.22)])
power_spectral_input = noise_input*noise_input

npts = freq_input.size

## interpolate to get more data points (needs to sample frequency linearly)
freq = linspace(1, Fc, nbins)
noise = zeros(nbins)
grLog = ROOT.TGraph(npts, log(freq_input), log(noise_input))
for i in range(nbins):
   noise[i] = grLog.Eval(log(freq[i]),0,"S") ## Gives a terrible interpolation sans logscale - need to take log first, then raise back in power.
## Reset back to non-log scale
noise = exp(noise)
power_spectral = noise*noise


## My implementation of method to simulate noise described in Wan's thesis.
sigma = noise/sqrt(2) #power_spectral_density =  2*sigma^2

# Genearte Xi and Yi - Re and Im fourier components for each frequency.
freqAmp = empty(nbins,cfloat)
rand = ROOT.TRandom3(0)
for i in range(nbins):
	freqAmp[i] = complex(rand.Gaus(0, sigma[i]), rand.Gaus(0, sigma[i]))

## Create Waveform. Magnitude of FT components is V/sqrt(Hz)
waveform = ROOT.MGTWaveform()
waveformFT = ROOT.MGTWaveformFT()
fft = ROOT.MGWFFastFourierTransformFFTW()
waveformFT.SetDataLength(nbins)
waveformFT.SetSamplingFrequency(Fc*2*1e-9) #1e-9 converts to GHz

for i in range(nbins):
	waveformFT[i] = ROOT.complex('double')(freqAmp[i].real, freqAmp[i].imag)

#for i in range(1000):
#	print freqAmp[i*8], "Mag: ", abs(freqAmp[i*8]), "Sigma:", sigma[i*8]

## Plot Power Spectrum stored in MGTWaveformFT
c1 = ROOT.TCanvas("c1", "Signal Analysis", 200,10,600,480)
pad1 = ROOT.TPad("pad1", "Frequency Domain", .05,.65,.95,.97)
pad2 = ROOT.TPad("pad2", "Time Domain", .05, .35,.95,.63)
pad3 = ROOT.TPad("pad3", "Time Domain", .05, .02, .95, .33)
pad1.SetFillColor(0)
pad2.SetFillColor(0)
pad3.SetFillColor(0)
pad1.Draw()
pad2.Draw()
pad3.Draw()
pad1.cd()
waveformFT.GimmeHist(ROOT.MGTWaveformFT.kPower).Draw()

### Compute Simulated Noise using MGDO method

## Inverse fourier transform (Destroys data in MGTWaveformFT)
fft.PerformInverseFFT(waveform, waveformFT)
waveform /= waveform.GetLength(); #Normalize

## Generate noise from waveformFT
waveform2 = ROOT.MGTWaveform()
waveform2FT = ROOT.MGTWaveformFT()
waveform2FT.SetDataLength(nbins)
waveform2FT.SetSamplingFrequency(Fc*2*1e-9)
waveform2.SetLength(waveform2FT.GetTDomainLength())
waveform2.SetSamplingFrequency(Fc*2*1e-9)

#waveform2.MakeSimilarTo(waveform2FT)

for i in range(nbins):
	waveform2FT[i] = noise[i]
#waveform2FT.SetData(noise)

noise_transform = ROOT.MGWFAddNoiseFromFT()
noise_transform.SetNoiseWaveform(waveform2FT)
noise_transform.Transform(waveform2) # Add noise to waveform

# Draw
pad2.cd()
waveform.GimmeHist().Draw()
pad3.cd()
waveform2.GimmeUniqueHist().Draw()
c1.Update()
"""
c3 = ROOT.TCanvas("c3", "Test",800, 500, 600, 480)
hist = ROOT.TH1D("name", "title", waveformFT.GetDataLength(), 0., .5*waveformFT.GetSamplingFrequency()*1e3)
print "waveformFT data length: ", waveformFT.GetDataLength()
binWidth = .5*waveformFT.GetSamplingFrequency()*1e3/waveformFT.GetDataLength()
hist.SetBins( waveformFT.GetDataLength(), 0.-binWidth/2, 0.5*waveformFT.GetSamplingFrequency()*1e3);
normalization = 1e6 * binWidth # why the normalization?

hist.SetBinContent( 1, pow(abs(freqAmp[0]),2))
hist.SetBinContent(waveformFT.GetDataLength(), pow(abs(freqAmp[waveformFT.GetDataLength()-1]),2));
for i in range(waveformFT.GetDataLength()-1):
	hist.SetBinContent( i+1, 2*pow(abs(freqAmp[i]),2) ) #Why the 2? (2 sigma^2 ?)

hist.Draw()
"""
"""
## Graph of Power Spectral Density
gr1 = ROOT.TGraph(nbins, freq, noise)
gr1.SetLineColor(1)
gr1.SetLineWidth(1)
gr1.SetMarkerStyle(3)

c2 = ROOT.TCanvas("c2", "Power Spectral Density", 300,100,600,480)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
zone = ROOT.TH2F("zone", "Noise Spectral Density; Frequency; Noise [V/sqrt(Hz)]", 1, 1, 100000000, 1, min(noise)/10, max(noise)*10)
zone.DrawCopy()
gr1.Draw("LP")
"""
raw_input("Press Enter to Close")