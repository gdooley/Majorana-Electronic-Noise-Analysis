/************************************************************************
* Reads a .root file and plots the energy spectrum with range and number
* of bins specified by user. Also fits gaussians to 3 regions of interest
* designated by bounds. 
* Requires use of the file read.C
* Return: array of length 9 containing FWHM, mean, FWHM error in that
* order for each peak.
*
* Author: Greg Dooley
* Contact: gdooley@princeton.edu
* Date: August 2010
************************************************************************/

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH2F.h"
#include <iostream>

double* read(const char* fname = "~/Documents/noise/struckDataGregRun3/struck_run500.root", int bounds[], int bins=25000, int range=400000)
{
	//Uncomment to use this file by itself to view energy spectra
    //int bounds[] = {2500,2530,10185,10240, 11570,11625};

	// Read From ROOT
	TFile* f = new TFile(fname,"READ");
	TTree* t = (TTree*)f->Get("SIS3302Decoder");
	assert(t);
	//cout <<t->GetEntries()<<endl;
	
	UInt_t energy;
	t->SetBranchStatus("*",0);
	t->SetBranchStatus("energy",1);
	t->SetBranchAddress("energy",&energy);
	
	// t->Scan("energy");
	int n = t->GetEntries();
	int nbins = bins;
	
	TH1F* h = new TH1F("h", "Energy Spectrum", nbins, 0, range);
	for (int i=0; i<n; ++i) {
		t->GetEntry(i);
		h->Fill((Double_t)energy);
		//cout <<i<<"\t"<<energy<<endl;
	}
	
	// Draw Energy Spectrum
	TCanvas *c1 = new TCanvas("c1","Energy Spectrum",200,10,700,500);
	c1->SetFillColor(0);
	TH2F* zone = new TH2F("zone", "Gaussian; ADU; Energy", 1, 0, 100000000, 1, 0, 150);
	zone->DrawCopy();
	h->Draw();
	
	// Integrate over gaussian to find # of counts
	//cout <<"counts = "<<h->Integral(7943, 7973)<<endl;	
	
	// Fit Gaussians
	double *width = new double[9];
	h->Fit("gaus", "RQ", "", bounds[0], bounds[1]);
	double constant = gaus->GetParameter(0);
	double mean = gaus->GetParameter(1);
	double sigma = gaus->GetParameter(2);
	double error = gaus->GetParError(2);
	double FWHM = 2.35*sigma;
	width[0] = FWHM;
	width[1] = mean;
	width[2] = error*2.35;
	//cout <<"error = "<<error<<endl;

	h->Fit("gaus", "RQ", "", bounds[2],bounds[3]);
	constant = gaus->GetParameter(0);
	mean = gaus->GetParameter(1);
	sigma = gaus->GetParameter(2);
	error = gaus->GetParError(2);
	FWHM = 2.35*sigma;
	width[3] = FWHM;
	width[4] = mean;
	width[5] = error*2.35;

//histogram->GetFunction  how to define gaus

	h->Fit("gaus", "RQ", "", bounds[4],bounds[5]);
	constant = gaus->GetParameter(0);
	mean = gaus->GetParameter(1);
	sigma = gaus->GetParameter(2);
	error = gaus->GetParError(2);
 	FWHM = 2.35*sigma;
	width[6] = FWHM;
	width[7] = mean;
	width[8] = error*2.35;
	
	double deltaADC = width[7];
	double deltaeV = 1330000;
	double ADC_to_eV = deltaeV/deltaADC;
	cout <<"means - "<<width[1]<<", "<<width[4]<<", "<<width[7]<<endl;
	// Choose 1.33 MeV line as calibration point for other peaks
	width[1] = ADC_to_eV*(width[1]) *1e-6;
	width[4] = ADC_to_eV*(width[4]) *1e-6;
	width[7] = 1.33;
	width[0] *= ADC_to_eV;
	width[2] *= ADC_to_eV;
	width[3] *= ADC_to_eV;
	width[5] *= ADC_to_eV;
	width[6] *= ADC_to_eV;
	width[8] *= ADC_to_eV;

	cout <<"means - "<<width[1]<<", "<<width[4]<<", "<<width[7]<<endl;
	return width;
	char dummy = cin.get();
}
