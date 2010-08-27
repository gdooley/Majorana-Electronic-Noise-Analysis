"""-------------------------------------------------------------------------------------
Plot the power, series and 1/f noise components together on a V/sqrt(Hz) scale.
Choose the constants of the components to see their behavior.
Plots an example spectrum on the same plot in blue.

Author: Greg Dooley
Contact: gdooley@princeton.edu
Date: August 2010
--------------------------------------------------------------------------------------"""


from ROOT import *
from numpy import *
import sys
gApplication.ExecuteFile("$MGDODIR/Majorana/LoadMGDOMJClasses.C")

##Constants
V_to_eV = 1 # Conversion factor = Capacitance / (Coulombs of charge collected per eV of incident energy)
nbins = 50 # number of interpolated data points to take

## Noise in eV/sqrt(Hz)
freq_input = array([1.1, 30.0, 400., 3311., 39810., 707945., 19054607.17])
noise_input = array([pow(10.,-5.82), pow(10., -5.95), pow(10., -6.4), pow(10, -7.2), pow(10, -7.8), pow(10,-8.15), pow(10,-8.22)]) # in V/sqrt(Hz)
noise_input *= V_to_eV # in eV/sqrt(Hz)
npts = freq_input.size

## interpolate to get more data points
freq = linspace(0, 16.76, nbins)
noise = zeros(nbins)
grLog = TGraph(npts, log(freq_input), log(noise_input))
for i in range(nbins):
   noise[i] = grLog.Eval(freq[i],0,"S") # Gives a terrible interpolation on logscale - need to take log first, then raise back in power.

series = zeros(nbins)
parallel = zeros(nbins)
over = zeros(nbins)
### Make Adjustments to constants here to see how the graph of the components change ###
# Equations done in log-log space. See my REU paper or Spieler for further explanation on the equation forms.
for i in range(nbins):
	parallel[i] = sqrt(1e-12 / (1+pow(1e-2*exp(freq[i]) ,2)))
	series[i] = log(1*pow(10,-8.3))
	over[i] = (-.5*freq[i]+log(pow(10,-6.1)) )

# Reset back to non-log scale
freq = exp(freq)
noise = exp(noise)
series = exp(series)
#parallel = exp(parallel)
over = exp(over)

gr2 = TGraph(nbins, freq, series)
gr2.SetLineWidth(1)
gr2.SetLineColor(1)
gr2.SetMarkerColor(1)

gr3 = TGraph(nbins, freq, parallel)
gr3.SetLineColor(1)
gr3.SetLineWidth(1)
gr3.SetMarkerStyle(3)

gr4 = TGraph(nbins, freq, over)
gr4.SetLineColor(1)
gr4.SetLineWidth(1)
gr4.SetMarkerStyle(3)

gr5 = TGraph(nbins, freq, over+series+parallel)
gr5.SetLineColor(2)
gr5.SetLineWidth(1)
gr5.SetMarkerStyle(3)


# Graph of Noise
gr1 = TGraph(npts, freq_input, noise_input)
gr1.SetLineColor(1)
gr1.SetLineWidth(1)
gr1.SetMarkerStyle(3)

# Graph of Interpolated Noise
grInter = TGraph(nbins, freq, noise)
grInter.SetMarkerStyle(24)
grInter.SetLineColor(4)
grInter.SetMarkerColor(4)

canvas = TCanvas("canvas", "Noise Analysis", 200,10,600,480)
canvas.SetFillColor(0)
gStyle.SetOptStat(0);

gPad.SetLogx()
gPad.SetLogy()

zone1 = TH2F("zone1", "Pre and Post amp Noise Spectral Density; Frequency; Noise [eV/sqrt(Hz)]", 1, 1, 100000000, 1, min(noise_input)/10, max(noise_input)*10)
zone1.DrawCopy()
gr1.Draw("P")
grInter.Draw("PL")
gr2.Draw("L")
gr4.Draw("L")
gr3.Draw("L")
gr5.Draw("L")

raw_input("enter to close")