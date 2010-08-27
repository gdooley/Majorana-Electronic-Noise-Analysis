/************************************************************************
* This script is used to measure the FWHM of a number of energy
* spectra directly and plot it as a function of peaking time
* It requires a degree of manual labor in opening the spectra
* and picking off the points on the pulser peak and two Cobalt peaks
* as to where to apply gaussian fits. Uses the 1.33 MeV cobalt peak and
* 0 point to calibrate the pulser and 1.17 MeV peak. Read through
* comments before using.
*
* Requires use of the file read.C
*
* Author: Greg Dooley
* Contact: gdooley@princeton.edu
* Date: August 2010
************************************************************************/

#include "read.C"
#include <sstream>
#include "math.h"
#include "TGraphErrors.h"
void FWHM2()
{
  int npts = 6; // number of energy spectra taken
                // tried to use this variable throughout, but doesn't seem to work with root interpreter
               // in allocating arrays 
  double peakTimes[6] = {1,10,50,100,500,1000}; //Peak times used for the different spectra as taken with Orca
  for (int i=0;i<6;i++){
    peakTimes[i] *= 80*1e-3;  // Calibrate peak times: *decimation * 10 nanoseconds per decimation and convert to uS
  }

  double FWHM2_pulser[6]; // Store FWHM^2
  double Error_pulser[6]; // array of errors on pulser FWHM
  double Error_peak1[6]; // 1.17 MeV peak
  double Error_peak2[6]; //1.33 MeV peak
  double FWHM2_peak1[6];
  double FWHM2_peak2[6];
  double Energy_pulser[6];

  double* data;
  // This is the tedios part. For every energy spectra must specify the value in ADC of the min and max location
  // where the gaussian is in order to fit a curve to that portion of the graph alone. Must specify min,max for
  // pulser, 1.17, 1.33 peaks, thus 6 numbers per spectrum.
  int bounds[6][6] = {{2500,2530,10185,10240, 11565,11625},
		      {13430,13500,50940,51160,57840,58100},
		      {90850,91120,255200,256500,290000,291200},
		      {239450,240050,511500,513500,581000,583250},
		      {2073000,2087000, 2557500,2565000,2905000,2914000},
		      {1,100,200,322,455,659}};

  // Choose the number of bins to use in the histogram of the data for each spectrum
  int bins[6] = {14000,25000,25000,25000,25000,25000};
  // Choose the maximum ADC data range to use in each energy spectra. - again, a manual process.
  int range[6] = {14000,70000,400000,700000,4000000,40000000};

  char fname[56];
  int run;
  int dummy;
  string runString;
  stringstream out;
  for (int i = 0; i < npts; i++){
    run = 500 + i; // Choose the starting point run number. .root files must be named consecutively
    dummy = sprintf(fname, "~/Documents/noise/struckDataGregRun3/struck_run%d.root", run); // put location of your files here
    //cout <<fname<<endl;                                                                                                                    
    data = read(fname,bounds[i],bins[i],range[i]);
    // store data                                                                                                                            
    FWHM2_pulser[i] = data[0]*data[0];
    Error_pulser[i] = sqrt( pow((data[2]/data[0]),2) *2) *data[0]*data[0]; // Error propagation in taking the square
    cout <<"Pulser FWHM^2 = "<<FWHM2_pulser[i]<<endl;                                                                             
    FWHM2_peak1[i] = data[3]*data[3];
    Error_peak1[i] = sqrt( pow((data[5]/data[3]),2) *2) *data[3]*data[3];    
    FWHM2_peak2[i] = data[6]*data[6];
    Error_peak2[i] = sqrt( pow((data[8]/data[6]),2) *2) *data[6]*data[6];
    Energy_pulser[i] = data[1];
    cout <<"Pulser Energy ="<<Energy_pulser[i]<<endl; 
  }

  for (int i=0;i<npts;i++){
    cout <<"error = "<<Error_pulser[i]<<endl;
  }

  // Plot Data (FWHM^2 vs peaking time) with Error bars
  TGraphErrors* gr1 = new TGraphErrors(npts,peakTimes,FWHM2_pulser, 0,Error_pulser);
  gr1->SetLineWidth(1);
  gr1->SetMarkerStyle(0); 

  TGraphErrors* gr2 = new TGraphErrors(npts,peakTimes,FWHM2_peak1, 0, Error_peak1);
  gr2->SetLineWidth(1);
  gr2->SetMarkerStyle(0);

  TGraphErrors* gr3 = new TGraphErrors(npts,peakTimes,FWHM2_peak2,0,Error_peak2);
  gr3->SetLineWidth(1);
  gr3->SetMarkerStyle(0);

  TGraph* gr4 = new TGraph(npts,peakTimes,Energy_pulser);
  gr4->SetLineWidth(1);
  gr4->SetMarkerStyle(6);

  TCanvas *c2 = new TCanvas("c1","Energy Plots",200,10,700,500);
  c2->SetFillColor(0);
  c2->Divide(2,2);
  gStyle->SetOptStat(0);
  c2->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  TH2F* zone1 = new TH2F("zone1", "Pulser FWHM; Peaking Time(uS); FWHM^2 (eV^2)", 1, .01, 100, 1, 300000,30000000);
  zone1->DrawCopy();
  gr1->Draw("PL");

  c2->cd(2);
  gPad->SetLogx();
  gPad->SetLogy();
  TH2F* zone2 = new TH2F("zone2", "1.17 MeV FWHM; Peaking Time(uS); FWHM^2 (eV^2)", 1, .01,1100, 1, 100000, 100000000);
  zone2->DrawCopy();
  gr2->Draw("PL");

  c2->cd(3);
  gPad->SetLogx();
  gPad->SetLogy();
  TH2F* zone3 = new TH2F("zone3", "1.33 MeV FWHM; Peaking Time (uS); FWHM^2 (eV^2)", 1, .01,100, 1, 100000, 100000000);
  zone3->DrawCopy();
  gr3->Draw("PL");

  c2->cd(4);
  gPad->SetLogx();
  gPad->SetLogy();
  TH2F* zone4 = new TH2F("zone4", "Pulser Energy; Peaking Time (uS); Energy (MeV)", 1, .01, 100, 1, .1, 1.2);
  zone4->DrawCopy();
  gr4->Draw("PL");
  c2->Update();
}
