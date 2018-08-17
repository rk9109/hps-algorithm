import math
import ROOT
from ROOT import *
from predict import predict

def parameter_stack(pred_tau, bins, low, high, parameter, xlabel, ylabel, title, filename):
	"""
	Return: None

	Input: pred_tau   | Predicted tau candidates
		   bins       | Number of bins
		   low        | Lower limit
		   high       | Upper limit
		   parameter  | Lambda function for parameter
		   xlabel     | X-axis title
		   ylabel     | Y-axis title
		   title      | Plot title
		   filename   | Name to save file 
	"""
	stack = THStack("hs", "hs")
	hNum1 = TH1F("hNum1", "hNum1", bins, low, high) # three hadrons
	hNum2 = TH1F("hNum2", "hNum2", bins, low, high) # hadrons + photons
	hNum3 = TH1F("hNum3", "hNum3", bins, low, high) # one hadron

	for (jetID, guessID, vec) in pred_tau:
		if (abs(jetID) == 15):
			if (guessID == '1'): hNum1.Fill(parameter(vec))
			elif (guessID in ['2', '3']): hNum2.Fill(parameter(vec))
			elif (guessID == '4'): hNum3.Fill(parameter(vec))
	
	hNum1.SetFillColor(15) # set colors
	hNum2.SetFillColor(46)
	hNum3.SetFillColor(9)

	stack.Add(hNum2) # add histograms to stack
	stack.Add(hNum1)
	stack.Add(hNum3)

	cst = TCanvas("cst","cst", 10, 10, 700, 700) # define canvas	
	
	stack.Draw("hist")
	stack.SetTitle(title) # set title and axis labels
	stack.GetXaxis().SetTitle(xlabel)
	stack.GetYaxis().SetTitle(ylabel)

	leg = TLegend(0.6, 0.7, 0.9, 0.9) # create legend
	leg.SetNColumns(1)
	
	entry1 = leg.AddEntry("hNum1","Three hadrons","f")
	entry2 = leg.AddEntry("hNum2","Hadrons + strips","f")
	entry3 = leg.AddEntry("hNum3","One hadron","f")

	entry1.SetFillColor(15)
	entry2.SetFillColor(46)
	entry3.SetFillColor(9)

	leg.Draw()

	cst.Update()
	cst.SaveAs(filename + ".png") # save as png
	
	return None

def parameter_hist(pred_tau, bins, low, high, parameter, xlabel, ylabel, title, filename, pred_other = None, logy = False):
	"""
	Return: None

	Input: pred_tau   | Predicted tau candidates
		   bins       | Number of bins
		   low        | Lower limit
		   high       | Upper limit
		   parameter  | Lambda function for parameter
		   xlabel     | X-axis title
		   ylabel     | Y-axis title
		   title      | Plot title
		   filename   | Name to save file 
		   pred_other | Predicted other candidates
		   logy       | Include logy scaling (True/False)
	"""
	hNum1 = TH1F("hNum1", "hNum1", bins, low, high) # three hadrons
	hNum2 = TH1F("hNum2", "hNum2", bins, low, high) # hadrons + photons
	hNum3 = TH1F("hNum3", "hNum3", bins, low, high) # one hadron
	hNum4 = TH1F("hNum4", "hNum4", bins, low, high) # background

	for (jetID, guessID, vec) in pred_tau:
		if (abs(jetID) == 15):
			if (guessID == '1'): hNum1.Fill(parameter(vec))
			elif (guessID in ['2', '3']): hNum2.Fill(parameter(vec))
			elif (guessID == '4'): hNum3.Fill(parameter(vec))

		elif pred_other: continue

		else: hNum4.Fill(parameter(vec))
	
	if pred_other:
		for (jetID, vec) in pred_other:
			if (abs(jetID) != 15):
				hNum4.Fill(parameter(vec))

	hNum1.Sumw2() # handle errors
	hNum2.Sumw2()
	hNum3.Sumw2()
	hNum4.Sumw2()

	hNum1.Scale(1./hNum1.Integral()) # scale histograms to unity	
	hNum2.Scale(1./hNum2.Integral()) 
	hNum3.Scale(1./hNum3.Integral())
	hNum4.Scale(1./hNum4.Integral())
	
	max_ = [hNum1.GetMaximum(), hNum2.GetMaximum(), hNum3.GetMaximum(), hNum4.GetMaximum()]
	hNum1.SetMaximum(max(max_)*1.4)
	hNum1.SetStats(False)

	cst = TCanvas("cst","cst", 10, 10, 700, 700) # define canvas	
	
	if (logy): 
		cst.SetLogy() # set y-axis to log-scale
	
	hNum1.SetLineWidth(2) # set width + colors
	hNum1.SetLineColor(1)
	hNum2.SetLineWidth(2)
	hNum2.SetLineColor(2)
	hNum3.SetLineWidth(2)
	hNum3.SetLineColor(3)
	hNum4.SetLineWidth(2)
	hNum4.SetLineColor(4)

	hNum1.Draw() # draw histogram
	hNum2.Draw('same')
	hNum3.Draw('same')
	hNum4.Draw('same')
	hNum1.SetTitle(title) # set title and axis labels
	hNum1.GetXaxis().SetTitle(xlabel)
	hNum1.GetYaxis().SetTitle(ylabel)
	
	leg = TLegend(0.6, 0.7, 0.9, 0.9) # create legend
	leg.SetNColumns(1)
	
	entry1 = leg.AddEntry("hNum1","Three hadrons","l")
	entry1.SetLineWidth(2)
	entry1.SetLineColor(1)
	
	entry2 = leg.AddEntry("hNum2","Hadrons + strips","l")
	entry2.SetLineWidth(2)
	entry2.SetLineColor(2)
	
	entry3 = leg.AddEntry("hNum3","One hadron","l")
	entry3.SetLineWidth(2)
	entry3.SetLineColor(3)
	
	entry4 = leg.AddEntry("hNum4","Background","l")
	entry4.SetLineWidth(2)
	entry4.SetLineColor(4)
	
	leg.Draw()

	cst.Update()
	cst.SaveAs(filename + ".png") # save as png
    
	return None

def eff_hist(pred_tau1, pred_other1, pred_tau2, pred_other2, bins, low, high, parameter, xlabel, ylabel, title, filename):
	"""
	Return: None

	Input: pred_tau   | Predicted tau candidates
		   pred_other | Predicted other candidates
		   bins       | Number of bins
		   low        | Lower limit
		   high       | Upper limit
		   parameter  | Lambda function for parameter
		   xlabel     | X-axis title
		   ylabel     | Y-axis title
		   title      | Plot title
		   filename   | Name to save file 
	"""
	hNum1 = TH1F("hNum1", "hNum1", bins, low, high) # correct tau
	hNum2 = TH1F("hNum2", "hNum2", bins, low, high) # all tau
	hNum3 = TH1F("hNum3", "hNum3", bins, low, high) # incorrect other
	hNum4 = TH1F("hNum4", "hNum4", bins, low, high) # all other
	
	for (jetID, _, vec) in pred_tau1:
		if (abs(jetID) == 15): 
			hNum1.Fill(parameter(vec))
			hNum2.Fill(parameter(vec))
	
	for (jetID, _, vec) in pred_tau2:
		if (abs(jetID) != 15):
			hNum3.Fill(parameter(vec))
			hNum4.Fill(parameter(vec))

	for (jetID, vec) in pred_other1:
		if (abs(jetID) == 15):
			hNum2.Fill(parameter(vec))
		
	for (jetID, vec) in pred_other2:
		if (abs(jetID) != 15):
			hNum4.Fill(parameter(vec))
					
	cst = TCanvas("cst","cst", 10, 10, 1000, 1000) # define canvas	

	T1 = TGraphAsymmErrors() # create graphs
	T2 = TGraphAsymmErrors()
	T1.Divide(hNum1, hNum2, "cl = 0.683")
	T2.Divide(hNum3, hNum4, "cl = 0.683")

	T1.SetMarkerColor(4)
	T2.SetMarkerColor(2)
	T1.SetMarkerStyle(8)
	T2.SetMarkerStyle(22)
	T1.SetLineColor(4)
	T2.SetLineColor(2)
					
	mg = TMultiGraph() # create multigraph
	mg.Add(T1)
	mg.Add(T2)
	mg.Draw('ap')
	mg.SetTitle(title) # set title and axis labels
	mg.GetXaxis().SetTitle(xlabel)
	mg.GetYaxis().SetTitle(ylabel)
	mg.GetYaxis().SetRangeUser(0., 1.2)
	
	leg = TLegend(0.6, 0.8, 0.9, 0.9) # create legend
	leg.SetNColumns(1)
	
	entry1 = leg.AddEntry("T1","Signal Efficiency","p")
	entry1.SetMarkerColor(4)
	entry1.SetMarkerStyle(8)
	
	entry2 = leg.AddEntry("T2","Background Efficiency","p")
	entry2.SetMarkerColor(2)
	entry2.SetMarkerStyle(22)
	
	leg.Draw()

	cst.Update()
	cst.SaveAs(filename + ".png") # save as png
    
	return None

def eff_hist_detail(pred_tau1, pred_other1, pred_tau2, pred_other2, bins, low, high, parameter, xlabel, ylabel, title, filename):
	"""
	Return: None

	Input: pred_tau   | Predicted tau candidates
		   pred_other | Predicted other candidates
		   bins       | Number of bins
		   low        | Lower limit
		   high       | Upper limit
		   parameter  | Lambda function for parameter
		   xlabel     | X-axis title
		   ylabel     | Y-axis title
		   title      | Plot title
		   filename   | Name to save file 
	"""
	hNum1 = TH1F("hNum1", "hNum1", bins, low, high) # correct tau 1
	hNum2 = TH1F("hNum2", "hNum2", bins, low, high) # correct tau 2/3
	hNum3 = TH1F("hNum3", "hNum3", bins, low, high) # correct tau 4
	hNum4 = TH1F("hNum4", "hNum4", bins, low, high) # all tau
	hNum5 = TH1F("hNum1", "hNum1", bins, low, high) # incorrect other 1
	hNum6 = TH1F("hNum2", "hNum2", bins, low, high) # incorrect other 2/3
	hNum7 = TH1F("hNum3", "hNum3", bins, low, high) # incorrect other 4
	hNum8 = TH1F("hNum4", "hNum4", bins, low, high) # all other
	
	for (jetID, guessID, vec) in pred_tau1:
		if (abs(jetID) == 15): 
			if (guessID == '1'): hNum1.Fill(parameter(vec))
			elif (guessID in ['2', '3']): hNum2.Fill(parameter(vec))
			elif (guessID == '4'): hNum3.Fill(parameter(vec))
			hNum4.Fill(parameter(vec))

	for (jetID, vec) in pred_other1:
		if (abs(jetID) == 15):
			hNum4.Fill(parameter(vec))
	
	for (jetID, guessID, vec) in pred_tau2:
		if (abs(jetID) != 15):	
			if (guessID == '1'): hNum5.Fill(parameter(vec))
			elif (guessID in ['2', '3']): hNum6.Fill(parameter(vec))
			elif (guessID == '4'): hNum7.Fill(parameter(vec))
			hNum8.Fill(parameter(vec))
	
	for (jetID, vec) in pred_other2:
		if (abs(jetID) != 15):
			hNum8.Fill(parameter(vec))

	cst = TCanvas("cst","cst", 10, 10, 1000, 1000) # define canvas	

	T1 = TGraphAsymmErrors() # create graphs
	T2 = TGraphAsymmErrors()
	T3 = TGraphAsymmErrors()
	T4 = TGraphAsymmErrors()
	T5 = TGraphAsymmErrors()
	T6 = TGraphAsymmErrors()
	T1.Divide(hNum1, hNum4, "cl = 0.683")
	T2.Divide(hNum2, hNum4, "cl = 0.683")
	T3.Divide(hNum3, hNum4, "cl = 0.683")
	T4.Divide(hNum5, hNum8, "cl = 0.683")
	T5.Divide(hNum6, hNum8, "cl = 0.683")
	T6.Divide(hNum7, hNum8, "cl = 0.683")
	
	T1.SetMarkerColor(4)
	T2.SetMarkerColor(2)
	T3.SetMarkerColor(3)
	T4.SetMarkerColor(38)
	T5.SetMarkerColor(46)
	T6.SetMarkerColor(30)

	T1.SetMarkerStyle(8)
	T2.SetMarkerStyle(8)
	T3.SetMarkerStyle(8)
	T4.SetMarkerStyle(8)
	T5.SetMarkerStyle(8)
	T6.SetMarkerStyle(8)

	T1.SetLineColor(4)
	T2.SetLineColor(2)
	T3.SetLineColor(3)
	T4.SetLineColor(38)
	T5.SetLineColor(46)
	T6.SetLineColor(30)
					
	mg = TMultiGraph() # create multigraph
	mg.Add(T1)
	mg.Add(T2)
	mg.Add(T3)
	mg.Add(T4)
	mg.Add(T5)
	mg.Add(T6)
	mg.Draw('ap')
	mg.SetTitle(title) # set title and axis labels
	mg.GetXaxis().SetTitle(xlabel)
	mg.GetYaxis().SetTitle(ylabel)
	mg.GetYaxis().SetRangeUser(0., 1.2)
	
	leg = TLegend(0.1, 0.8, 0.9, 0.9) # create legend
	leg.SetNColumns(6)
	
	entry1 = leg.AddEntry("T1","Eff 1","p")
	entry1.SetMarkerColor(4)
	entry1.SetMarkerStyle(8)
	
	entry2 = leg.AddEntry("T2","Eff 2/3","p")
	entry2.SetMarkerColor(2)
	entry2.SetMarkerStyle(8)
	
	entry3 = leg.AddEntry("T3","Eff 4","p")
	entry3.SetMarkerColor(3)
	entry3.SetMarkerStyle(8)

	entry4 = leg.AddEntry("T4","FP 1","p")
	entry4.SetMarkerColor(38)
	entry4.SetMarkerStyle(8)

	entry5 = leg.AddEntry("T5","FP 2/3","p")
	entry5.SetMarkerColor(46)
	entry5.SetMarkerStyle(8)

	entry6 = leg.AddEntry("T6","FP 4","p")
	entry6.SetMarkerColor(30)
	entry6.SetMarkerStyle(8)
	
	leg.Draw()

	cst.Update()
	cst.SaveAs(filename + ".png") # save as png
    
	return None

#------------PLOTTING------------

filename1 = "/home/drankin/TauID/GenNtuple_GluGluHToTauTau_M125_13TeV_powheg_pythia8.root"
filename2 = "/home/drankin/TauID/GenNtuple_ZprimeToTauTau_M-3000_TuneCP5_13TeV-pythia8-tauola.root"
filename3 = "/home/drankin/TauID/GenNtuple_QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root"
filename4 = "/home/drankin/TauID/GenNtuple_QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root"
treeNum = "GenNtupler/gentree"
iterations = 25000

# working points
working_points1 = [('low', 0.05), ('medium', 0.03), ('tight', 0.02)]
working_points2 = [('low', 0.20), ('medium', 0.12), ('tight', 0.08)]

# get predictions
for name, wp in working_points2:
	pred_tau1, pred_other1 = predict(filename1, treeNum, iterations, 0.5, 0.5, wp)
	pred_tau3, pred_other3 = predict(filename3, treeNum, iterations, 0.5, 0.5, wp)

	# efficiency plots
	eff_hist(pred_tau1, pred_other1, pred_tau3, pred_other3, 15, 0., 200., lambda x: x.Pt(), 
	         "Pt", "Efficiency", "GluGluHToTauTau WP " + name, name + "eff_GluGluHToTauTau_Pt")
	
	eff_hist(pred_tau1, pred_other1, pred_tau3, pred_other3, 10, -3., 3., lambda x: x.Eta(), 
			 "Eta", "Efficiency", "GluGluHToTauTau WP " + name, name + "eff_GluGluHToTauTau_Eta")

	eff_hist_detail(pred_tau1, pred_other1, pred_tau3, pred_other3, 15, 0., 200., lambda x: x.Pt(), 
					"Pt", "Efficiency", "GluGluHToTauTau WP " + name, name + "effdet_GluGluHToTauTau_Pt")
	
	eff_hist_detail(pred_tau1, pred_other1, pred_tau3, pred_other3, 10, -3., 3., lambda x: x.Eta(), 
					"Eta", "Efficiency", "GluGluHToTauTau WP " + name, name + "effdet_GluGluHToTauTau_Eta")

"""
# stack plots
parameter_stack(pred_tau1, 100, 0., 2., lambda x: x.M(), "Invariant Mass", "Events", "GluGluHToTauTau", "stack_GluGluHToTauTau")
parameter_stack(pred_tau2, 100, 0., 2., lambda x: x.M(), "Invariant Mass", "Events", "GluGluHToTauTau", "stack_ZPrimeToTauTau")

# hist plots
parameter_hist(pred_tau1, 75, 0., 150., lambda x: x.Pt(), "Pt", "", "GluGluHToTauTau_Pt", "hist_GluGluHToTauTau_Pt")
parameter_hist(pred_tau1, 75, 0., 150., lambda x: x.Pt(), "Pt", "", "GluGluHToTauTau_Pt", "hist_GluGluHToTauTau_Pt_log", None, True)
parameter_hist(pred_tau1, 50, -3., 3., lambda x: x.Eta(), "Eta", "", "GluGluHToTauTau_Eta", "hist_GluGluHToTauTau_Eta")
parameter_hist(pred_tau1, 50, -3., 3., lambda x: x.Eta(), "Eta", "", "GluGluHToTauTau_Eta", "hist_GluGluHToTauTau_Eta_log", None, True)
parameter_hist(pred_tau1, 50, -3., 3., lambda x: x.Phi(), "Phi", "", "GluGluHToTauTau_Phi", "hist_GluGluHToTauTau_Phi")
parameter_hist(pred_tau1, 50, -3., 3., lambda x: x.Phi(), "Phi", "", "GluGluHToTauTau_Phi", "hist_GluGluHToTauTau_Phi_log", None, True)
parameter_hist(pred_tau1, 100, 0., 300., lambda x: x.E(), "Energy", "", "GluGluHToTauTau_Energy", "hist_GluGluHToTauTau_Energy")
parameter_hist(pred_tau1, 100, 0., 300., lambda x: x.E(), "Energy", "", "GluGluHToTauTau_Energy", "hist_GluGluHToTauTau_Energy_log", None, True)
parameter_hist(pred_tau1, 100, 0., 2., lambda y: y.M(), "Mass", "", "GluGluHToTauTau_Mass", "hist_GluGluHToTauTau_Mass", pred_other1)
parameter_hist(pred_tau1, 100, 0., 2., lambda y: y.M(), "Mass", "", "GluGluHToTauTau_Mass", "hist_GluGluHToTauTau_Mass_log", pred_other1, True)

parameter_hist(pred_tau2, 75, 0., 1750., lambda x: x.Pt(), "Pt", "", "ZPrimeToTauTau_Pt", "hist_ZPrimeToTauTau_Pt")
parameter_hist(pred_tau2, 75, 0., 1750., lambda x: x.Pt(), "Pt", "", "ZPrimeToTauTau_Pt", "hist_ZPrimeToTauTau_Pti_log", None, True)
parameter_hist(pred_tau2, 50, -3., 3., lambda x: x.Eta(), "Eta", "", "ZPrimeToTauTau_Eta", "hist_ZPrimeToTauTau_Eta")
parameter_hist(pred_tau2, 50, -3., 3., lambda x: x.Eta(), "Eta", "", "ZPrimeToTauTau_Eta", "hist_ZPrimeToTauTau_Eta_log", None, True)
parameter_hist(pred_tau2, 50, -3., 3., lambda x: x.Phi(), "Phi", "", "ZPrimeToTauTau_Phi", "hist_ZPrimeToTauTau_Phi")
parameter_hist(pred_tau2, 50, -3., 3., lambda x: x.Phi(), "Phi", "", "ZPrimeToTauTau_Phi", "hist_ZPrimeToTauTau_Phi_log", None, True)
parameter_hist(pred_tau2, 100, 0., 2500., lambda x: x.E(), "Energy", "", "ZPrimeToTauTau_Energy", "hist_ZPrimeToTauTau_Energy")
parameter_hist(pred_tau2, 100, 0., 2500., lambda x: x.E(), "Energy", "", "ZPrimeToTauTau_Energy", "hist_ZPrimeToTauTau_Energy_log", None, True)
parameter_hist(pred_tau2, 100, 0., 2., lambda y: y.M(), "Mass", "", "ZPrimeToTauTau_Mass", "hist_ZPrimeToTauTau_Mass", pred_other2)
parameter_hist(pred_tau2, 100, 0., 2., lambda y: y.M(), "Mass", "", "ZPrimeToTauTau_Mass", "hist_ZPrimeToTauTau_Mass_log", pred_other2, True)
"""


