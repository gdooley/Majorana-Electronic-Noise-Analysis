## Plots 1000 waveforms to file.ps 4 waves per page.
#  Use to quickly study the pulses
# Author: Greg Dooley
# Contact: gdooley@princeton.edu
# Date: August 2010

from RawNoise import *
from ROOT import *

a = RawNoise(1)
c1 = TCanvas("c1");
ps = TPostScript("file.ps",112);
c1.Divide(2,2);

for i in range(250):
	ps.NewPage();
	for j in range(4):
		w = a.GetMGTWaveform(i*4+j)
		p = w.GimmeHist(str(i*4+j))
		c1.cd(j+1)
		c1.Update()        
		p.Draw()
	c1.Update()
ps.Close();

# Alternate method
"""
c1 = TCanvas("c1");
c1.Divide(2,2);
w = a.GetMGTWaveform(0)
p = w.GimmeHist("0")
p.Draw();
c1.Print("c1.ps(");  # write canvas and keep the ps file open

for i in range(1,40):
	w = a.GetMGTWaveform(i)
	p = w.GimmeHist(str(i))
	p.Draw();
	c1.Print("c1.ps");

w = a.GetMGTWaveform(41)
p = w.GimmeHist("41")
p.Draw();
c1.Print("c1.ps)");  # canvas is added to "c1.ps" and ps file is closed
"""