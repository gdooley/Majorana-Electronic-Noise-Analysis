"""--------------------------------------------------------------------------------------
   This script reads a noise power spectrum and provides functions to analyze the data
   namely to plot FWHM^2 vs Peaking Time, and to try to pick off the parallel, series
   and 1/f noise components of the power spectrum.

   Usage: name = NoiseAnalysis(noise data (in charge^2), frequency, option to interpolate)

   Date:   August-2010
   Author: Greg Dooley
   Contact: gdooley@princeton.edu
---------------------------------------------------------------------------------------"""

from ROOT import *
from numpy import *
import sys
gApplication.ExecuteFile("$MGDODIR/Majorana/LoadMGDOMJClasses.C")

## From Trapezoidal Filter transfer function
## Gap time set here
def shaper(f, P=1e-6, G=1e-8):
   response = abs(1/(pi*P*f) * (cos(pi*f*G)-cos(pi*f*(2*P+G))))
   return response

## Analog Filter transfer function. Only the one name shaper will be used in code. Rename to swap.
def shaper2(f, time=1e-7, n=4):
   w = 2*pi*f
   x = time*w
   shapeResponse = abs(x*1j / pow(1+x*1j, n+1))
   return shapeResponse
	
##Constants
ADC_to_eV = 224  ## Conversion value not guaranteed - measured by Greg Dooley on August 17, 2010
nbins = 100 # number of interpolated data points to take when using interpolation option
hand_freq = array([1e4,3e4,1e5,3e5,1e6,2e6,4e6,1e7,2e7,4e7]) ## hand chosen frequencies to match power spectrum
hand_noise = array([3,2.6,1.3,.45,.3,.3,.28,.36,.25,.21]) ## hand chosen noise values in eV/sqrt(Hz)


class NoiseAnalysis:
	## n should be array of the power of the FT of raw noise data. (X^2 + Y^2). Retrieved from RawNoise.
	## f is the corresponding array of frequency values
	def __init__(self, n, f, interpolate=0):
		## Attempt to normalize histogram in units of eV^2/Hz
		self.noise = n
		self.noise *= 1/((f[2]-f[1])*pow(10,6.0)) * 1/sqrt((2*(f.size-1.0))) ## Very important normalization by 1/sqrt(n) and divide by bin width - may not be right!!
		self.noise = sqrt(self.noise) #ADC/sqrt(Hz)
		self.noise *= ADC_to_eV # convert to eV
		
		self.npts = len(n)
		self.freq = f
		# FWHM^2(t) = at + b/t + c
		self.a = 0 # Parallel constant
		self.b = 0 # Series Constant
		self.c = 0 # 1/f Constant 
		self.ParallelConst = 0 # N(f) = u / f.    a = ParallelConst * u^2  See final REU paper for better explanation if unclear.
		self.SeriesConst = 0 # N(f) = v.          b = SeriesConst * v^2
		self.OverConst = 0 # N(f) = w / sqrt(f).  c = OverConst *  w^2
		
		self.u = 0 # for parallel noise	
		self.v = 0 # for series noise
		self.w = 0 # for 1/f noise
		
		## First data point is 0 (undershoot bin of histogram), second is small - probably not physical

		## interpolate hand input values to get more data points
		self.freqI = linspace(9.0, 17.727, nbins)
		self.noiseI = zeros(nbins)
		
		gr2 = TGraph(hand_freq.size, log(hand_freq), log(hand_noise))
		gr2.SetMarkerColor(4)
		gr2.SetMarkerStyle(5)
		
		for i in range(nbins):
			self.noiseI[i] = gr2.Eval(self.freqI[i],0,"S") # Gives a terrible interpolation sans logscale - need to take log first, then raise back in power.
		# Reset back to non-log scale
		self.freqI = exp(self.freqI)
		self.noiseI = exp(self.noiseI)
		
		self.freq_input = self.freq
		self.noise_input = self.noise
		self.npts_input = self.npts
		self.interpolate = interpolate
		if(interpolate):
			self.freq = self.freqI
			self.noise = self.noiseI
			self.npts = nbins
	

	## Plot input and interpolated noise (if option chosen) (eV / sqrt(Hz)) on log-log scale
	def PlotNoise(self):
		# Graph of input noise
		gr1 = TGraph(self.npts_input, self.freq_input, self.noise_input)
		gr1.SetLineColor(1)
		gr1.SetLineWidth(1)
		gr1.SetMarkerStyle(3)
		gr2 = TGraph(hand_freq.size, hand_freq, hand_noise)
		gr2.SetMarkerColor(2)
		gr2.SetMarkerStyle(5)
		
		# Graph of Interpolated Noise
		grInter = TGraph(self.freqI.size, self.freqI, self.noiseI)
		grInter.SetMarkerStyle(24)
		grInter.SetMarkerColor(4)
		grInter.SetMarkerStyle(2)
		grInter.SetLineColor(4)
		grInter.SetLineWidth(2)
		
		canvas = TCanvas("canvas", "Noise Analysis", 200,10,800,400)
		canvas.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		gPad.SetLogy()	
		zone1 = TH2F("zone1", "Electronic Noise Spectral Density; Frequency; Noise [eV/sqrt(Hz)]", 1, min(self.freq)/5, max(self.freq)*10, 1, min(self.noise[2:self.npts-1])/10, max(self.noise)*10)
		zone1.DrawCopy()
		gr1.Draw("L")
		if(self.interpolate):
			gr2.Draw("P")
			grInter.Draw("L")   ## creating a second NoiseAnalysis object in a python script with interpolate = true prevents gr1 from drawing. Bizarre problem.
		raw_input("Enter to Close Noise Plot")
	
	## Fit three noise components directly to the power spectrum	
	def PlotFitNoise(self):
		# Graph of input noise
		gr1 = TGraph(self.npts, self.freq, self.noise)
		gr1.SetLineColor(1)
		gr1.SetLineWidth(1)
		gr1.SetMarkerStyle(3)
				
		xmin = min(self.freq)/10.0
		xmax = max(self.freq)*10.0		
		
		canvas = TCanvas("canvas", "Noise Analysis", 200,10,800,400)
		canvas.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		gPad.SetLogy()
		zone1 = TH2F("zone1", "Electronic Noise Spectral Density; Frequency; Noise [eV/sqrt(Hz)]", 1, min(self.freq)/10.0, max(self.freq)*10, 1, min(self.noise[2:self.npts-1])/10, max(self.noise)*10)
		zone1.DrawCopy()
		gr1.Draw("L")
		
		##Fit to parallel, series and parallel noise.
		func = TF1("func","sqrt([0]/(1+[1]*x*x)) + [2] + [3]/sqrt(x)",xmin,xmax)
		func.SetParameters(self.GetK1(), self.GetK2(), self.GetV(), self.GetU()) #Need to guess parameters here
		func.SetParNames("parallel_1", "parallel_2", "series", "1/f");
		func.SetLineWidth(2);
		func.SetLineStyle(2);
		gr1.Fit("func");
		print "ChiSquare: " + str(func.GetChisquare()) + "  DOF: " + str(func.GetNDF());

		parallel_1 = func.GetParameter(0);
		parallel_2 = func.GetParameter(1)
		series = func.GetParameter(2);
		frequency = func.GetParameter(3);
		
		parFit = TF1("parFit", "sqrt([0]/(1+[1]*x*x))", xmin,xmax);
		parFit.SetParameter(0,parallel_1);
		parFit.SetParameter(1,parallel_2);
		parFit.SetLineWidth(2);
		parFit.SetLineStyle(1);
		parFit.SetLineColor(3);
		parFit.SetParName(0,"Parallel Noise");

		seriesFit = TF1("seriesFit",  "[0]", xmin,xmax);
		seriesFit.SetParameter(0,series);
		seriesFit.SetLineWidth(2);
		seriesFit.SetLineStyle(3);
		seriesFit.SetLineColor(2);
		seriesFit.SetParName(0,"Series Noise");

		freqFit = TF1("freqFit",  "[0]/sqrt(x)", xmin,xmax);
		freqFit.SetParameter(0,frequency);
		freqFit.SetLineWidth(2);
		freqFit.SetLineStyle(4);
		freqFit.SetLineColor(4);
		freqFit.SetParName(0,"1/f Noise");

		parFit.Draw("same");
		seriesFit.Draw("same");
		freqFit.Draw("same");

		#Add legend
		leg = TLegend(.5,.63,.90,.90);
		title1 = "Parallel Noise Consts = %3.3g, %3.3g, Error = %3.3g, %3.3g" %(parallel_1, parallel_2, func.GetParError(0),func.GetParError(1))
		title2 = "Series Noise = %3.3g, Error = %3.3g" %(series, func.GetParError(1))
		title3 = "1/f Noise = %3.3g, Error = %3.3g" %(frequency, func.GetParError(2))
		leg.AddEntry(parFit, title1, "l");
		leg.AddEntry(seriesFit, title2, "l");
		leg.AddEntry(freqFit, title3, "l");
		leg.SetFillColor(0);  
		leg.Draw("same");
		raw_input("Enter to Close Noise Plot")	
	
	# Graph of shaper curve. Uses input frequency as the x-axis range
	def PlotShaper(self, time, gap):
		shapeResponse = shaper(self.freq, time, gap)
		grShaper = TGraph(self.npts, self.freq, shapeResponse)
		grShaper.SetLineWidth(2)
		c = TCanvas("c", "Shaper", 200,10,800,400)
		c.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		zone = TH2F("zone", "Filter Frequency Response; Frequency; Fraction of signal passing", 1,1,max(self.freq)*10,1,0,2)
		zone.DrawCopy()
		grShaper.Draw("C")
		raw_input("Enter to Exit Shaper Plot")
	
	## Noise after Shaper for a given shaper time/peaking time
	def PostShaper(self, time):
		shapeResponse = shaper(self.freq, time)
		return shapeResponse*self.noise
	
	## Plot noise spectral density after shaper/filter has been applied for a given peaking time.	
	def PlotPostShaper(self, time):
		postShaper = self.PostShaper(time)
		grShaper = TGraph(self.npts, self.freq, postShaper)
		grShaper.SetLineWidth(2)
		c = TCanvas("c", "Post Shaper", 200,10,800,400)
		c.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		gPad.SetLogy()
		zone = TH2F("zone", "Noise Spectral Density After Filter; Frequency; Noise [eV/sqrt(Hz)]", 1,1,max(self.freq)*10, 1, min(postShaper[2:self.npts-1])/10, max(postShaper)*10)
		zone.DrawCopy()
		grShaper.Draw("C")
		raw_input("Enter to Exit PostShaper Plot")
	
	## RMS Noise - integrate power spectral noise density and take sqrt
	def RMSNoise(self, time):
		noise_shaper = self.PostShaper(time)
		power_spectral_shaped = pow(noise_shaper, 2) # eV^2 / Hz
		cumInt = zeros(self.npts)
		for i in range(self.npts):
   			cumInt[i] = sqrt(trapz(power_spectral_shaped[0:i+1], self.freq[0:i+1])) # eV
		return cumInt
	

	# Plots Cumulative integral over frequency of RMS noise
	def PlotRMSNoise(self, time):
		cumInt = self.RMSNoise(time)
		grInt = TGraph(self.npts, self.freq, cumInt)
		grInt.SetLineWidth(2)
		grInt.SetMarkerStyle(3)
		c = TCanvas("c", "RMS Noise Cumulative Integral", 200,10,800,400)
		c.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		zone = TH2F("zone4", "Total noise as cumulative integral over frequency; Frequency; RMS Noise (eV)", 1,1,max(self.freq)*10 ,1, cumInt[1]/2,max(cumInt)*1.5)
		zone.DrawCopy()
		grInt.Draw("L") 
		raw_input("Enter to Exit RMSNoise Plot")
	

	## Generate total noise for a variety of shaper times.
	## sTimes is array of shaper times in seconds	
	def ShaperCurve(self, sTimes, n = array([0])):
		if n.size == 1:
			n = self.noise
		total_noise = zeros(sTimes.size)
		for i in range(sTimes.size):
   			noise_shaper2 = n * shaper(self.freq, sTimes[i])
   			power2 = pow(noise_shaper2, 2)
   			total_noise[i] = sqrt(trapz(power2, self.freq))
		return total_noise
	
	
	# Returns array of FWHM for given shaping times from sTimes
	def GetFWHM(self, sTimes):
		total_noise = self.ShaperCurve(sTimes)
		return 2.35*total_noise
	
		
	# Returns array of FWHM^2 for given shaping times from sTimes	
	def GetFWHM2(self, sTimes):
		total_noise = self.ShaperCurve(sTimes)
		return pow(2.35*total_noise,2)
	
		
	## Draw Plot of FWHM^2 vs Shaping Time. Must guess the paralle, series and 1/f constants to reasonable
	## accuracy for fit to work. Must also ensure sTimes is within accptable range of the frequency used.
	def PlotFWHMvsPeakingTime(self, sTimes, guessA, guessB, guessC):
		total_noise = self.ShaperCurve(sTimes)
		FWHM2 = pow(2.35*total_noise,2) # Multiply by 2.35 to get FWHM from RMS
		times = sTimes*1e6 # Multiply by 1e6 to get uS time.
		#print FWHM2
		grFWHM = TGraph(sTimes.size, times, FWHM2) 
		grFWHM.SetLineColor(1)
		grFWHM.SetLineWidth(2)
		grFWHM.SetMarkerStyle(3)
		### TEST CODE HERE##
		"""
		n2 = zeros(self.npts)
		w = 430158
		for i in range(self.npts):
			n2[i] = w / self.freq[i]
		t2 = self.ShaperCurve(sTimes, n2)
		FWHM2_2 = pow(2.35*t2,2) 
		grTEST = TGraph(sTimes.size, times, FWHM2_2)
		grTEST.SetLineColor(8)
		grTEST.SetLineWidth(8)
		grTEST.SetMarkerStyle(2)
		print "w = ", w
		print FWHM2_2
		a = FWHM2_2[0]
		print "a =", a
		print "const = ", a/(w*w)
		"""
		###END TEST CODE###
		c1 = TCanvas("c1", "Noise", 200,10,1100,800)
		c1.SetFillColor(0)
		gStyle.SetOptStat(0)
		gPad.SetLogy()
		gPad.SetLogx()
		
		# Set scale for plots
		xmin = min(times)/10.0
		xmax = max(times)*10.0
		
		zone = TH2F("zone","FWHM^2 vs Peak Shaping Time; Shaping Peak Time [us]; FWHM^2 (eV^2)", 1, xmin, xmax, 1, min(FWHM2)/10, max(FWHM2)*10);
		zone.DrawCopy();
		grFWHM.Draw("PL")
		#grTEST.Draw("PL")  ## TEST CODE
		
		##Fit to parallel, series and parallel noise. 

		func = TF1("func","[0]*x+[1]/x+[2]",xmin,xmax)
		func.SetParameters(guessA, guessB, guessC) #Need to guess parameters here
		func.SetParNames("parallel", "series", "1/f");
		func.SetLineWidth(2);
		func.SetLineStyle(2);
		#func.SetParLimits(2, 0, 1e20);  ## Newly added untested code!!
		#func->Draw("same");
		grFWHM.Fit("func");
		print "ChiSquare: " + str(func.GetChisquare()) + "  DOF: " + str(func.GetNDF());

		parallel = func.GetParameter(0);
		series = func.GetParameter(1);
		frequency = func.GetParameter(2);
		
		self.a = parallel *1e6 # [eV^2 / sec] Need 1e6 to convert from 1/us to 1/S
		self.b = series *1e-6 # [eV^2 * sec] Need 1e-6 to convert from uS to S.
		self.c = frequency # This is just a constant independent of time [eV^2]

		parFit = TF1("parFit",  "[0]*x", xmin,xmax);
		parFit.SetParameter(0,parallel);
		parFit.SetLineWidth(2);
		parFit.SetLineStyle(1);
		parFit.SetLineColor(3);
		parFit.SetParName(0,"Parallel Noise");

		seriesFit = TF1("seriesFit",  "[0]/x", xmin,xmax);
		seriesFit.SetParameter(0,series);
		seriesFit.SetLineWidth(2);
		seriesFit.SetLineStyle(3);
		seriesFit.SetLineColor(2);
		seriesFit.SetParName(0,"Series Noise");

		freqFit = TF1("freqFit",  "[0]", xmin,xmax);
		freqFit.SetParameter(0,frequency);
		freqFit.SetLineWidth(2);
		freqFit.SetLineStyle(4);
		freqFit.SetLineColor(4);
		freqFit.SetParName(0,"1/f Noise");

		parFit.Draw("same");
		seriesFit.Draw("same");
		freqFit.Draw("same");

		#Add legend
		leg = TLegend(.3,.63,.7,.87);
		title1 = "Parallel Noise = %3.3g, Error = %3.3g" %(parallel, func.GetParError(0))
		title2 = "Series Noise = %3.3g, Error = %3.3g" %(series, func.GetParError(1))
		title3 = "1/f Noise = %3.3g, Error = %3.3g" %(frequency, func.GetParError(2))
		leg.AddEntry(parFit, title1, "l");
		leg.AddEntry(seriesFit, title2, "l");
		leg.AddEntry(freqFit, title3, "l");
		leg.SetFillColor(0);  
		leg.Draw("same");
		raw_input("enter to close")
	

	## Draw 4 plots in process of computing FWHM for 1 shaper time
	def PlotSteps(self, time):
		# Graph of Noise
		gr1 = TGraph(self.npts, self.freq, self.noise)
		gr1.SetLineColor(1)
		gr1.SetLineWidth(1)
		gr1.SetMarkerStyle(3)
		"""
		# Graph of Interpolated Noise
		grInter = TGraph(nbins, freq, noise)
		grInter.SetMarkerStyle(24)
		grInter.SetLineColor(4)
		grInter.SetMarkerColor(4)
		"""
		# Graph of shaper curve
		shapeResponse = shaper(self.freq, time)
		grShaper = TGraph(self.npts, self.freq, shapeResponse)
		grShaper.SetLineWidth(2)
		
		# Graph of post shaper Noise
		noise_shaper = self.PostShaper(time)
		grPostShaper = TGraph(self.npts, self.freq, noise_shaper)
		grPostShaper.SetLineColor(1)
		grPostShaper.SetLineWidth(2)
		grPostShaper.SetMarkerStyle(20)
		
		# Cumulative integral curve
		cumInt = self.RMSNoise(time)
		grInt = TGraph(self.npts, self.freq, cumInt)
		grInt.SetLineWidth(2)
		grInt.SetMarkerStyle(3)
		canvas = TCanvas("canvas", "Noise Analysis", 200,10,1600,800)
		canvas.SetFillColor(0)
		gStyle.SetOptStat(0);
		pad1 = TPad("pad1", "a", .04,.5,.48,.97)
		pad2 = TPad("pad2", "b", .04,.02,.48,.47)
		pad3 = TPad("pad3", "c", .52, .5, .95, .97)
		pad4 = TPad("pad4", "d", .52, .02, .95, .47)
		pad1.SetLogx()
		pad1.SetLogy()
		pad2.SetLogx()
		pad2.SetLogy()
		pad3.SetLogx()
		pad4.SetLogx()
		pad1.SetFillColor(0)
		pad2.SetFillColor(0)
		pad3.SetFillColor(0)
		pad4.SetFillColor(0)
		pad1.Draw()
		pad2.Draw()
		pad3.Draw()
		pad4.Draw()
		pad1.cd()
		
		zone1 = TH2F("zone1", "Electronic Noise Spectral Density; Frequency; Noise [eV/sqrt(Hz)]", 1, min(self.freq)/5, max(self.freq)*10, 1, min(self.noise[2:self.npts-1])/10, max(self.noise)*10) 
		zone1.DrawCopy()
		gr1.Draw("L")
		#grInter.Draw("PL")
		pad2.cd()
		zone2 = TH2F("zone2", "Noise Spectral Density After Filter; Frequency; Noise [eV/sqrt(Hz)]", 1,min(self.freq)/5,max(self.freq)*10, 1, min(noise_shaper[2:self.npts-1])/10, max(noise_shaper)*10)
		zone2.DrawCopy()
		grPostShaper.Draw("C")
		
		pad3.cd()
		zone3 = TH2F("zone3", "Filter Frequency Response; Frequency; % signal passing", 1,1,max(self.freq)*10,1,0,2)
		zone3.DrawCopy()
		grShaper.Draw("C")
		pad4.cd()
		
		zone4 = TH2F("zone4", "Total noise as cumulative integral over frequency; Frequency; RMS Noise (eV)", 1,min(self.freq)/2,max(self.freq)*10 ,1, cumInt[1]/2,max(cumInt)*1.5)
		zone4.DrawCopy()
		grInt.Draw("L") 
		raw_input("enter to close")
	

	# Get Parallel fit constant a [eV^2 / sec] Need 1e6 to convert from 1/us to 1/S
	def GetParallel(self):
		if self.a == 0:
			print "Error: must run PlotFWHMvsShapingTime(sTimes, guessA, guessB, guessC) first"
		return self.a
	
	# Get Series fit constant b [eV^2 * sec] Need 1e-6 to convert from uS to S.
	def GetSeries(self):
		if self.b == 0:
			print "Error: must run PlotFWHMvsShapingTime(sTimes, guessA, guessB, guessC) first"
		return self.b
	
	# Get 1/f fit constant c
	def GetOneOverF(self):
		if self.c == 0:
			print "Error: must run PlotFWHMvsShapingTime(sTimes, guessA, guessB, guessC) first"
		return self.c
	
	## Uses sample power spectrum of the parallel noise given a constant u and finds the corresponding value of a
	## in the FWHM^2 space in order to find the constant relation between u and a.
	def GetParallelConst(self):
		if self.ParallelConst !=0:
			return self.ParallelConst
		nbins2 = 50
		freq = linspace(0, 16.8, nbins2)
		freq = exp(freq)
		parallel = zeros(nbins2)
		sTimes = array([100e-6, 31.6e-6, 10e-6, 7e-6, 3.16e-6,2e-6, 1e-6, 316e-9, 100e-9])
		u = pow(10, 5) # pick an arbitrary value for u, then find a to find const. empirically does not depend on choice of u.
		for i in range(nbins2):
			parallel[i] = u/(freq[i])
		
		## Get FWHM for shaping times
		total_noise = zeros(sTimes.size)
		for i in range(sTimes.size):
			noise_shaper2 = parallel * shaper(freq, sTimes[i])
			power2 = pow(noise_shaper2, 2)
			total_noise[i] = sqrt(trapz(power2, freq))
				
		FWHM2 = pow(2.35*total_noise,2) # Multiply by 2.35 to get FWHM from RMS
		guess = exp(log(FWHM2[5]) - log(sTimes[5]))
		
		# Use FWHM^2 vs Shaping Time to find a
		grFWHM = TGraph(sTimes.size, sTimes, FWHM2) # keep in units of seconds
		
		## Uncomment sections to plot curves for debugging purposes
		"""
		grFWHM.SetLineColor(1)
		grFWHM.SetLineWidth(2)
		grFWHM.SetMarkerStyle(3)
		
		c1 = TCanvas("c1", "Noise", 200,10,1100,800)
		c1.SetFillColor(0)
		gStyle.SetOptStat(0)
		gPad.SetLogy()
		gPad.SetLogx()

		zone = TH2F("zone","FWHM^2 vs Peak Shaping Time; Shaping Peak Time [us]; FWHM^2 (eV^2)", 1, min(sTimes), max(sTimes), 1, min(FWHM2)/10, max(FWHM2)*10);
		zone.DrawCopy();
		grFWHM.Draw("PL");
		"""
		##Fit to parallel, series and parallel noise. 
		func = TF1("func","[0]*x",min(sTimes),max(sTimes))
		func.SetParameter(0,guess); #Need to guess parameters here
		func.SetParName(0,"parallel");
		func.SetLineWidth(2);
		func.SetLineStyle(2);
		grFWHM.Fit("func");
		a = func.GetParameter(0);    
		self.ParallelConst = a/(u*u)
		"""		
		parFit = TF1("parFit",  "[0]*x",min(sTimes),max(sTimes));
		parFit.SetParameter(0,a);
		parFit.SetLineWidth(2);
		parFit.SetLineStyle(1);
		parFit.SetLineColor(3);
		parFit.SetParName(0,"Parallel Noise");
		parFit.Draw("same");
		"""
		return a/(u*u)
	

	# N(f) = u / f. a = ParallelConst *  u^2 for parallel noise (V/sqrt(Hz))
	def GetU(self):
		const = self.GetParallelConst()
		a = self.GetParallel()
		self.u = sqrt(a/const)
		return self.u
	
	# Uses sample power spectrum of the series noise given a constant v and finds the corresponding value of b
	## in the FWHM^2 space in order to find the constant relation between v and b.
	def GetSeriesConst(self):
		if self.SeriesConst !=0:
			return self.SeriesConst
		nbins2 = 50
		freq = linspace(0, 16.8, nbins2)
		freq = exp(freq)
		series = zeros(nbins2)
		sTimes = array([100e-6, 31.6e-6, 10e-6, 7e-6, 3.16e-6,2e-6, 1e-6, 316e-9, 100e-9])
		v = pow(10, 1) # pick an arbitrary value for v, then find b to find const. empirically does not depend on choice of v.
		for i in range(nbins2):
			series[i] = v
		
		## Get FWHM for shaping times
		total_noise = zeros(sTimes.size)
		for i in range(sTimes.size):
			noise_shaper2 = series * shaper(freq, sTimes[i])
			power2 = pow(noise_shaper2, 2)
			total_noise[i] = sqrt(trapz(power2, freq))
				
		FWHM2 = pow(2.35*total_noise,2) # Multiply by 2.35 to get FWHM from RMS
		guess = exp(log(FWHM2[5]) + log(sTimes[5]))
		
		# Use FWHM^2 vs Shaping Time to find b
		grFWHM = TGraph(sTimes.size, sTimes, FWHM2) # keep in units of seconds

        ## uncomment for plotting purposes
		"""
		grFWHM.SetLineColor(1)
		grFWHM.SetLineWidth(2)
		grFWHM.SetMarkerStyle(3)

		c1 = TCanvas("c1", "Noise", 200,10,1100,800)
		c1.SetFillColor(0)
		gStyle.SetOptStat(0)
		gPad.SetLogy()
		gPad.SetLogx()

		zone = TH2F("zone","FWHM^2 vs Peak Shaping Time; Shaping Peak Time [us]; FWHM^2 (eV^2)", 1, min(sTimes), max(sTimes), 1, min(FWHM2)/10, max(FWHM2)*10);
		zone.DrawCopy();
		grFWHM.Draw("PL");
		"""
		##Fit to parallel, series and parallel noise. 
		func = TF1("func","[0]/x",min(sTimes),max(sTimes))
		func.SetParameter(0,guess); #Need to guess parameters here
		func.SetParName(0,"series");
		func.SetLineWidth(2);
		func.SetLineStyle(2);
		grFWHM.Fit("func");
		b = func.GetParameter(0);    
		self.SeriesConst = b/(v*v)
		"""
		parFit = TF1("parFit",  "[0]/x",min(sTimes),max(sTimes));
		parFit.SetParameter(0,b);
		parFit.SetLineWidth(2);
		parFit.SetLineStyle(1);
		parFit.SetLineColor(3);
		parFit.SetParName(0,"Series Noise");
		parFit.Draw("same");
		raw_input("wait!")
		"""
		return b/(v*v)
	
	# N(f) = v for Series Noise (V/ sqrt(Hz)) b = SeriesConst * v^2
	def GetV(self):
		const = self.GetSeriesConst()
		b = self.GetSeries()
		self.v = sqrt(b/const)
		return self.v
	
	
	# Uses sample power spectrum of the 1/f noise given a constant w and finds the corresponding value of c
	## in the FWHM^2 space in order to find the constant relation between w and c.	
	def GetOneOverFConst(self):
		if self.OverConst !=0:
			return self.OverConst
		nbins2 = 50
		freq = linspace(0, 16.8, nbins2)
		freq = exp(freq)
		over = zeros(nbins2)
		sTimes = array([100e-6, 31.6e-6, 10e-6, 7e-6, 3.16e-6,2e-6, 1e-6, 316e-9, 100e-9])
		w = pow(10, 2) # pick an arbitrary value for w, then find c to find const. empirically does not depend on choice of w.
		for i in range(nbins2):
			over[i] = w/sqrt(freq[i])
		
		## Get FWHM for shaping times
		total_noise = zeros(sTimes.size)
		for i in range(sTimes.size):
			noise_shaper2 = over * shaper(freq, sTimes[i])
			power2 = pow(noise_shaper2, 2)
			total_noise[i] = sqrt(trapz(power2, freq))
				
		FWHM2 = pow(2.35*total_noise,2) # Multiply by 2.35 to get FWHM from RMS
		guess = FWHM2[5]
		
		# Use FWHM^2 vs Shaping Time to find c
		grFWHM = TGraph(sTimes.size, sTimes, FWHM2) # keep in units of seconds
		"""
		grFWHM.SetLineColor(1)
		grFWHM.SetLineWidth(2)
		grFWHM.SetMarkerStyle(3)

		c1 = TCanvas("c1", "Noise", 200,10,1100,800)
		c1.SetFillColor(0)
		gStyle.SetOptStat(0)
		gPad.SetLogy()
		gPad.SetLogx()

		zone = TH2F("zone","FWHM^2 vs Peak Shaping Time; Shaping Peak Time [us]; FWHM^2 (eV^2)", 1, min(sTimes), max(sTimes), 1, min(FWHM2)/10, max(FWHM2)*10);
		zone.DrawCopy();
		grFWHM.Draw("PL");
		"""
		##Fit to parallel, series and parallel noise. 
		func = TF1("func","[0]",min(sTimes),max(sTimes))
		func.SetParameter(0,guess); #Need to guess parameters here
		func.SetParName(0,"1/f");
		func.SetLineWidth(2);
		func.SetLineStyle(2);
		grFWHM.Fit("func");
		c = func.GetParameter(0);
		self.OverConst = c/(w*w)
		"""
		parFit = TF1("parFit", "[0]",min(sTimes),max(sTimes));
		parFit.SetParameter(0,c);
		parFit.SetLineWidth(2);
		parFit.SetLineStyle(1);
		parFit.SetLineColor(3);
		parFit.SetParName(0,"1/f Noise");
		parFit.Draw("same");
		raw_input("wait!")
		"""
		return c/(w*w)
	
		
	# N(f) = w / sqrt(f). c = OverConst *  w^2 for 1/f noise (V/sqrt(Hz))
	def GetW(self):
		const = self.GetOneOverFConst()
		c = self.GetOneOverF()
		self.w = sqrt(c/const)
		return self.w
	
	## Parallel noise requires 2 consants - not 1. Find them from w and noise curve calibration
	def GetK2(self):
		u = self.GetU()
		v = self.GetV()
		w = self.GetW()
		P = self.noise[3]-v-w/sqrt(self.freq[3]) ## Find remaining power needed to calibrate curves to fit exactly at freq[3] 
		k2 = P*P/(u*u-P*P*self.freq[3]*self.freq[3]) ##Not a good method to do this! not robust
		return k2
	
	## Parallel noise requires 2 consants - not 1. Find them from w and noise curve calibration	
	def GetK1(self):
		u = self.GetU()
		return self.GetK2()*u*u
	
	# Plot Power specturm with its components	
	def PlotPowerComponents(self):
		u = self.GetU()
		v = self.GetV()
		w = self.GetW()
		series = zeros(self.npts)
		over = zeros(self.npts)
		parallel = zeros(self.npts)
		## Parallel noise requires 2 consants - not 1. Find them from w and noise curve calibration
		k1 = self.GetK1()
		k2 = self.GetK2()
		
		for i in range(self.npts):
			parallel[i] = sqrt(k1/(1+k2*pow(self.freq[i],2)))
			series[i] = v
			over[i] = w/sqrt(self.freq[i])
		
		# Graph of Noise
		gr1 = TGraph(self.npts, self.freq, self.noise)
		gr1.SetLineColor(1)
		gr1.SetLineWidth(1)
		gr1.SetMarkerStyle(3)
		
		# Graph of components
		gr2 = TGraph(self.npts, self.freq, series)
		gr2.SetLineWidth(2)
		gr2.SetLineColor(2)
		gr2.SetMarkerColor(1)

		gr3 = TGraph(self.npts, self.freq, over)
		gr3.SetLineColor(3)
		gr3.SetLineWidth(2)
		gr3.SetMarkerStyle(3)

		gr4 = TGraph(self.npts, self.freq, parallel)
		gr4.SetLineColor(4)
		gr4.SetLineWidth(2)
		gr4.SetMarkerStyle(3)
		
		gr5 = TGraph(self.npts, self.freq, parallel+series+over)
		gr5.SetLineColor(6)
		gr5.SetLineWidth(2)
		gr5.SetMarkerStyle(3)

		canvas = TCanvas("canvas", "Noise Analysis", 200,10,600,480)
		canvas.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		gPad.SetLogy()

		zone1 = TH2F("zone1", "Reverse fit of Electronic Noise Components to Power Spectrum; Frequency; Noise [eV/sqrt(Hz)]", 1, min(self.freq)/10.0, max(self.freq)*10, 1, min(self.noise[2:self.npts-1])/10, max(self.noise)*10)
		zone1.DrawCopy()
		gr1.Draw("L")
		gr2.Draw("L")
		gr3.Draw("L")
		gr4.Draw("L")
		gr5.Draw("L")

		#Add legend
		leg = TLegend(.5,.63,.9,.87);
		title1 = "Parallel Noise"
		title2 = "Series Noise"
		title3 = "1/f Noise"
		leg.AddEntry(gr4, title1, "l");
		leg.AddEntry(gr2, title2, "l");
		leg.AddEntry(gr3, title3, "l");
		leg.AddEntry(gr5, "Sum of Components", "l")
		leg.AddEntry(gr1, "Actual Power Spectrum", "l")
		leg.SetFillColor(0);  
		leg.Draw("same");
		raw_input("Enter to Close Power Decomposition")	
	

