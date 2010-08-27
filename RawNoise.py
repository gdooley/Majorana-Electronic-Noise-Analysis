"""--------------------------------------------------------------------------------------
   This script reads waveforms from a root file. It subtracts the baseline from signals to 
   center them at 0. Computes fourier transform of each signal, combines all
   in a histogram averaged over the input data. Vetoes out all waveforms that have a
   |value| > MaxNoiseADC as they are likely a physics pulse or pulse reset. Does not rule
   out the possibility of very small events getting through in the data. Averaging many 
   histograms is very time consuming - so run it once and write GetXValues, GetYValues 
   to file.

   This file makes use of MGDO classes.

   Usage: name = RawNoise(number of samples to read, name of .root file, name of Tree in root file)

   Date:   August-2010
   Author: Greg Dooley
   Contact: gdooley@princeton.edu
----------------------------------------------------------------------------------"""

from ROOT import *
from numpy import *
gApplication.ExecuteFile("$MGDODIR/Majorana/LoadMGDOMJClasses.C")

class RawNoise:
	## Open Root file rFile and read Tree rTree. Read samples number of waveforms in average
	## power spectrum
	def __init__(self, samples = 0, rFile = "struck_run346.root", rTree = "SIS3302Decoder"):
		self.rootFile = rFile
		self.rootTree = rTree
		self.samples = samples
		## Read in waveforms
		self.f = TFile(self.rootFile)
		self.tree = self.f.Get(self.rootTree)
		self.wf = MGTWaveform()
		self.tree.SetBranchAddress("MGTWaveformBranch", self.wf)
		self.baseline = MGWFBaselineRemover()
		self.fft = MGWFFastFourierTransformFFTW()
		self.hist = TH1D() #Empty histogram to store cumulative histogram
		self.simWave = MGTWaveform() #store simulated wave
		self.MaxNoiseADC = 30  ## All noise events fall within +- 30 units.

	## Return MGTWaveform corresponding to id with baseline removed
	def GetMGTWaveform(self, id):
		if id < 0 or id >= self.tree.GetEntries():
			print "Error: id out of bounds max is %g" %(self.tree.GetEntries()-1)
			return
		else:
			self.tree.GetEntry(id)
			self.baseline.SetBaselineTime( self.wf.GetLength()/self.wf.GetSamplingFrequency() )
			self.baseline.Transform(self.wf)
			return self.wf
			
	## Gives Histogram of waveform corresponding to id
	def GimmeWFHist(self, id):
		wave = self.GetMGTWaveform(id)
		return wave.GimmeHist(str(id))
	
	## Plot Waveform corresponding to id	
	def PlotWF(self, id):
		h = self.GimmeWFHist(id)
		c = TCanvas("c1", "Raw Waveform", 200,10,600,480)
		c.SetFillColor(0)
		gStyle.SetOptStat(0);
		h.Draw()
		raw_input("Enter to close waveform")
		
	## Return number of waveforms in data
	def GetNumWaveforms(self):
		return self.tree.GetEntries()
	
	## Give FT of waveform corresponding to ID
	def GetMGTWaveformFT(self, id):
		wave = self.GetMGTWaveform(id)
		waveformFT = MGTWaveformFT()
		self.fft.PerformFFT(wave, waveformFT)
		return waveformFT

	## Histogram of FT of id
	def GimmeFTHist(self, id):
		waveFT1 = self.GetMGTWaveformFT(id)
		return waveFT1.GimmeHist(MGTWaveformFT.kPower, str(id))

	## Plot FT of waveform corresponding to ID
	def PlotFT(self, id):
		h = self.GimmeFTHist(id)
		c = TCanvas("c1", "FT of Raw Waveform", 200,10,600,480)
		c.SetFillColor(0)
		gStyle.SetOptStat(0);
		gPad.SetLogx()
		gPad.SetLogy()
		h.Draw()
		raw_input("Enter to close FT")
		
	## Add power spectrum histograms for all pulses and take the average
	def GimmeCumFTHist(self):
		#empty hist has 1 bin. Return existing histogram if already computed.
		if self.hist.GetNbinsX() != 1:
			return self.hist
		count = 0 # for first loop to find first good noise
		notFound = 1
		vetoes = 0 # count number of physics pulses or reset pulses
		
		# Get first acceptable waveform to make histogram
		while notFound:
			waveFT = self.GetMGTWaveformFT(count)
			m = max(fabs(self.wf.GetVectorData()))
			if m < self.MaxNoiseADC:  #May want to comment out later
				self.hist = waveFT.GimmeHist(MGTWaveformFT.kPower,"-1") # -1 is unique (doesn't interfere with PlotFT ids)
				print "Histogram %g added" %count
				notFound = 0
			else:
				print "Histogram %g vetoed. Max = %g" %(count, m)
				vetoes+=1
			count+=1

		if self.samples == 0:
			self.samples = self.tree.GetEntries()
		for i in range(count,self.samples):
			waveFT = self.GetMGTWaveformFT(i)
			m = max(fabs(self.wf.GetVectorData()))
			if m < self.MaxNoiseADC: # May want to comment out later
				self.hist.Add(waveFT.GimmeHist(MGTWaveformFT.kPower, str(i)))
				print "Histogram %g added" %i
			else:
				print "Histogram %g vetoed. Max = %g" %(i,m)
				vetoes+=1
		self.hist.Scale(1.0/(self.samples-vetoes))
		return self.hist
	
	## Marino suggests using AddNormsSquared from his MGDO library, but I felt safer with
	## my original method. - This function may not be the way to use it correctly, but its
	## a start.
	def GimmeCumFTHist2(self):
		waveFT2 = self.GetMGTWaveformFT(0)
		print "Histogram 0 added"
		if self.samples == 0:
			self.samples = self.tree.GetEntries()
		vetoes2 = 0
		for i in range(1,self.samples):
			waveFTtemp = self.GetMGTWaveformFT(i)
			m2 = max(fabs(self.wf.GetVectorData()))
			if m2 < 30:
				waveFT2.AddNormsSquared(waveFTtemp)
				print "Histogram %g added" %i
			else:
				print "Histogram %g vetoed. Max = %g" %(i,m2)
				vetoes2+=1
			waveFT2 /= sqrt(self.samples - vetoes2)
		return waveFT2.GimmeHist(MGTWaveformFT.kPower, "-2")
		
	## Plot the averaged cumulative FT
	def PlotCumFT(self):
		h = self.GimmeCumFTHist()
		c1 = TCanvas("c1", "Signal Analysis", 200,10,600,480)
		c1.SetFillColor(0)
		gPad.SetFillColor(0)
		gPad.SetLogx()
		gPad.SetLogy()
		h.Draw()
		raw_input("Enter to close Cumulative FT")
		
	## Maximum Frequency available in FT
	def GetMaxFreq(self):
		self.tree.GetEntry(0)
		return self.wf.GetSamplingFrequency()*1e9/2.0  # convert from GHz to Hz
	
	## Minimum frequency possible in FT
	def GetMinFreq(self):
		self.tree.GetEntry(0)
		return 1.0/self.GetTime()
	
	## Number of samples in the waveform
	def GetLength(self):
		self.tree.GetEntry(0)
		return self.wf.GetLength()
	
	## Time interval of measurement in seconds
	def GetTime(self):
		return self.GetLength()/2.0/self.GetMaxFreq()
			
	## Return array of Y values (noise power)
	def GetYValues(self):
		l = self.hist.GetNbinsX()
		if l == 1:
			self.GimmeCumFTHist()
		arr = self.hist.GetArray()
		power = zeros(l)
		for i in range(l):
			power[i] = arr[i]
		return power
	
	## Return array of X values in hertz. (frequency)
	def GetXValues(self):
		if self.hist.GetNbinsX() == 1:
			self.GimmeCumFTHist()
		width = self.hist.GetBinWidth(0) * pow(10,6) #convert from MH to Hz
		npts = self.hist.GetNbinsX()
		freq = zeros(npts)
		for i in range(npts):
			freq[i] = i*width+width/2.0
		return freq
	
	## Returns MGTWaveform of simulated noise.
	## If already computed, returns the most recent waveform in memory
	def GetSimulatedNoise(self):
		if self.simWave.GetLength() != 0:
			return self.simWave
		else:
			h = self.GimmeCumFTHist()
			waveFT = self.GetMGTWaveformFT(0)
			waveTemp = self.GetMGTWaveform(0)
			self.simWave.SetLength(waveFT.GetTDomainLength())
			self.simWave.SetSamplingFrequency(waveTemp.GetSamplingFrequency())
			for i in range(waveFT.GetDataLength()):
				waveFT[i] = sqrt(h[i]) # Take sqrt to get V/sqrt(Hz)
			noise_transform = MGWFAddNoiseFromFT()
			noise_transform.SetNoiseWaveform(waveFT)
			noise_transform.Transform(self.simWave) # Add noise to waveform
			return self.simWave
	
	## Create new simulated waveform
	def GetNewSimulatedNoise(self):
		h = self.GimmeCumFTHist()
		waveFT = self.GetMGTWaveformFT(0)
		waveTemp = self.GetMGTWaveform(0)
		self.simWave = MGTWaveform()
		self.simWave.SetLength(waveFT.GetTDomainLength())
		self.simWave.SetSamplingFrequency(waveTemp.GetSamplingFrequency())
		for i in range(waveFT.GetDataLength()):
			waveFT[i] = sqrt(h[i]) # Take sqrt to get V/sqrt(Hz)
		noise_transform = MGWFAddNoiseFromFT()
		noise_transform.SetNoiseWaveform(waveFT)
		noise_transform.Transform(self.simWave) # Add noise to waveform
		return self.simWave
	
	## Plot most recent simulated noise
	def PlotSimulatedNoise(self):
		wave = self.GetSimulatedNoise()
		c = TCanvas("c1", "Simulated Noise", 200,10,600,480)
		c.SetFillColor(0)
		gStyle.SetOptStat(0);
		wave.GimmeHist().Draw()
		raw_input("Enter to close Simulated Noise")

	## Plot Imaginary vs Real Components of Fourier Transform of waveform id
	def PlotYvsX(self, id):
		wave = self.GetMGTWaveformFT(id)
		hX = wave.GimmeHist(MGTWaveformFT.kReal,str(id))
		hY = wave.GimmeHist(MGTWaveformFT.kImaginary,str(id))
		X = hX.GetArray()
		Y = hY.GetArray()
		gr = TGraph(wave.GetDataLength(),X,Y)
		gr.SetMarkerStyle(1)
		c = TCanvas("c", "Distribution of fourier components", 200, 10, 600, 480)
		c.SetFillColor(0)
		gStyle.SetOptStat(0)
		zone = TH2F("zone","Power Spectrum Fourier Components; Real Xi; Imaginary Yi", 1,-2000,2000,1,-2000,2000)
		zone.DrawCopy()
		gr.Draw("P")
		raw_input("Enter to close")

	## Create histogram of phase distribution of phase of waveformFT id
	def PlotPhaseDistribution(self, id):
		wave = self.GetMGTWaveformFT(id)
		hX = wave.GimmeHist(MGTWaveformFT.kReal,str(id))
		hY = wave.GimmeHist(MGTWaveformFT.kImaginary,str(id))
		X = hX.GetArray()
		Y = hY.GetArray()
		h = TH1F("h", "Phase Distribution of Power Spectrum; Phase (degrees); Counts", 90, -90, 90)
		h.SetMinimum(0) #Force Y axis to draw 0
		for i in range(wave.GetDataLength()):
			h.Fill(arctan2(Y[i],X[i])*180/pi) 	
		c = TCanvas("c", "Distribution of fourier components", 200, 10, 600, 480)
		c.SetFillColor(0)
		gStyle.SetOptStat(0)
		h.Draw()
		raw_input("Enter to close")

	## Create histogram of distribution of Xi or Yi Fourier components of waveformFT id
	## Default is imaginary components. type=0 for real components. Used to show they
	## are gaussian distributed, justifying the simulation method.
	def PlotFTDistribution(self, id, type):
		wave = self.GetMGTWaveformFT(id)
		if type == 0:
			h = wave.GimmeHist(MGTWaveformFT.kReal,str(id))
		else:
			h = wave.GimmeHist(MGTWaveformFT.kImaginary,str(id))
		X = h.GetArray()
		h = TH1F("h", "Power Spectrum Fourier Component Distribution; Xi or Yi; Counts", 300, -3000, 3000)
		h.SetMinimum(0) #Force Y axis to draw 0
		for i in range(wave.GetDataLength()):
			h.Fill(X[i]) 	
		c = TCanvas("c", "Distribution of Xi or Yi", 200, 10, 600, 480)
		c.SetFillColor(0)
		gStyle.SetOptStat(0)
		h.Draw()
		h.Fit("gaus", "Q")
		raw_input("Enter to close")