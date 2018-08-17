import math
import ROOT
from ROOT import *
from predict import *
from array import array

def plot_iso(isolation_list1, isolation_list2, bins, low, high, filename):
	"""
	docstring
	"""
	hNum1 = TH1F("hNum1", "hNum1", bins, low, high) # correct tau
	hNum2 = TH1F("hNum2", "hNum2", bins, low, high) # incorrect tau

	for (jetID, iso) in isolation_list1:
		if (abs(jetID) == 15):
			hNum1.Fill(iso)
	for (jetID, iso) in isolation_list2:
		if (abs(jetID) != 15):
			hNum2.Fill(iso)

	hNum1.Sumw2() # handle errors
	hNum2.Sumw2()

	hNum1.Scale(1./hNum1.Integral()) # scale histograms to unity	
	hNum2.Scale(1./hNum2.Integral())

	max_ = [hNum1.GetMaximum(), hNum2.GetMaximum()]
	hNum1.SetMaximum(max(max_)*1.4)
	hNum1.SetStats(False)

	cst = TCanvas("cst","cst", 10, 10, 700, 700) # define canvas	
	
	hNum1.SetLineWidth(2) # set width + colors
	hNum1.SetLineColor(1)
	hNum2.SetLineWidth(2)
	hNum2.SetLineColor(2)

	hNum1.Draw() # draw histogram
	hNum2.Draw('same')
	hNum1.SetTitle("Isolation Plot") # set title and axis labels
	hNum1.GetXaxis().SetTitle("Isolation")
	hNum1.GetYaxis().SetTitle("Number")
	
	leg = TLegend(0.6, 0.7, 0.9, 0.9) # create legend
	leg.SetNColumns(1)
	
	entry1 = leg.AddEntry("hNum1","Tau","l")
	entry1.SetLineWidth(2)
	entry1.SetLineColor(1)
	
	entry2 = leg.AddEntry("hNum2","Background","l")
	entry2.SetLineWidth(2)
	entry2.SetLineColor(2)
	
	leg.Draw()

	cst.Update()
	cst.SaveAs(filename + ".png") # save as png

	return None

def plot_roc(isolation_list1, isolation_list2, bins, low, high, filename, decaymode, logy = False):
	"""
	docstring
	"""
	hNum1 = TH1F("hNum1", "hNum1", bins, low, high) # correct tau
	hNum2 = TH1F("hNum2", "hNum2", bins, low, high) # incorrect tau

	for (jetID, iso) in isolation_list1:
		if (abs(jetID) == 15):
			hNum1.Fill(iso)
	for (jetID, iso) in isolation_list2:  
		if (abs(jetID) != 15): 
			hNum2.Fill(iso)
		
	hNum1_total = hNum1.Integral() # get totals
	hNum2_total = hNum2.Integral(1, bins + 1)
	
	x = array('d'); y = array('d') # ROC curve arrays
	x.append(0.); y.append(0.)
	
	lx = array('d'); ly = array('d') # Reference line arrays
	lx.append(0.); lx.append(1.); ly.append(0.); ly.append(1.)

	actual_x = array('d'); actual_y = array('d')
	if (decaymode == 'GluGluHToTauTau'):
		actual_x.append(49.0/100); actual_x.append(40.8/100); actual_x.append(38.1/100)
		actual_y.append(0.00386); actual_y.append(0.00206); actual_y.append(0.00175)

	if (decaymode == 'ZPrimeToTauTau'):
		actual_x.append(58.9/100); actual_x.append(50.8/100); actual_x.append(48.1/100)
		actual_y.append(0.00386); actual_y.append(0.00206); actual_y.append(0.00175)

	for i in range(bins):
		x.append(hNum1.Integral(1, i + 1)/hNum1_total) # signal efficiency
		y.append(hNum2.Integral(1, i + 1)/hNum2_total) # background rejection

	t1 = TGraph(bins, x, y)
	t2 = TGraph(2, lx, ly)
	t3 = TGraph(3, actual_x, actual_y)	

	cst = TCanvas("cst","cst", 10, 10, 1000, 1000) # define canvas	
	
	if (logy): 
		cst.SetLogy() # set y-axis to log-scale

	t1.SetLineColor(2)
	t1.SetLineWidth(2)
	t1.SetMarkerColor(4)
	t1.SetMarkerStyle(22)
	t2.SetLineColor(1)
	t2.SetLineWidth(2)
	if (logy):
		t2.SetLineColor(0)
		t2.SetLineWidth(0)
	t3.SetLineColor(0)
	t3.SetLineWidth(0)
	t3.SetMarkerColor(3)
	t3.SetMarkerStyle(29)

	mg = TMultiGraph() # create multigraph
	mg.Add(t1)
	mg.Add(t2)
	mg.Add(t3)
	mg.Draw('alp')
	mg.SetTitle("ROC Curve") # set title and axis labels
	mg.GetXaxis().SetTitle("Signal Efficiency")
	mg.GetYaxis().SetTitle("Background Efficiency")

	cst.Update()
	cst.SaveAs(filename + ".png") # save as png

#----------Plotting-----------
filename1 = "/home/drankin/TauID/GenNtuple_GluGluHToTauTau_M125_13TeV_powheg_pythia8.root"
filename2 = "/home/drankin/TauID/GenNtuple_ZprimeToTauTau_M-3000_TuneCP5_13TeV-pythia8-tauola.root"
filename3 = "/home/drankin/TauID/GenNtuple_QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root"
filename4 = "/home/drankin/TauID/GenNtuple_QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root"
treeNum = "GenNtupler/gentree"
iterations = 25000

isolation_list1 = predict(filename1, treeNum, iterations, 0.5, 0.5)
isolation_list2 = predict(filename2, treeNum, iterations, 0.5, 0.5)
isolation_list3 = predict(filename3, treeNum, iterations, 0.5, 0.5)
isolation_list4 = predict(filename4, treeNum, iterations, 0.5, 0.5)

plot_roc(isolation_list1, isolation_list3, 200, 0, 2, "roc_GluGluHToTauTau", "GluGluHToTauTau")
plot_roc(isolation_list2, isolation_list4, 200, 0, 2, "roc_ZPrimeToTauTau", "ZPrimeToTauTau")
plot_roc(isolation_list1, isolation_list3, 200, 0, 2, "roc_GluGluHToTauTau_log", "GluGluHToTauTau", True)
plot_roc(isolation_list2, isolation_list4, 200, 0, 2, "roc_ZPrimeToTauTau_log", "ZPrimeToTauTau", True)

plot_iso(isolation_list1, isolation_list3, 40, 0, 2, "iso_GluGluHToTauTau")
plot_iso(isolation_list1, isolation_list3, 40, 0, 0.5, "iso_GluGluHToTauTau_zoom")
plot_iso(isolation_list2, isolation_list4, 40, 0, 2, "iso_ZPrimeToTauTau")
plot_iso(isolation_list2, isolation_list4, 40, 0, 0.5, "iso_ZPrimeToTauTau_zoom")


