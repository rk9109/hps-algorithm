import json, math
import ROOT, sys, os, re, string
from ROOT import *
import numpy as n
# load libraries

ROOT.gROOT.SetBatch(ROOT.kTRUE) # do not print outputs of draw or load graphics
ROOT.gStyle.SetOptStat(0)       # do not draw stat box
#ROOT.gStyle.SetPalette(55)

def compplot(filename1, treeNum, varname, nbins, low, high, logy, xlabel, outname, title, cuts1, label1, cuts2, label2, cuts3="", label3="", cuts4="", label4="", cuts5="", label5="", cuts6="", label6=""): #plot same variable with different selection
	
	f1 = TFile(filename1) #open file
	tNum1 = f1.Get(treeNum) #get TTree from file
	hNum1 = TH1F("hNum1","hNum1",nbins,low,high) #define histograms
	hNum2 = TH1F("hNum2","hNum2",nbins,low,high)
	hNum3 = TH1F("hNum3","hNum3",nbins,low,high)
	hNum4 = TH1F("hNum4","hNum4",nbins,low,high)
	hNum5 = TH1F("hNum5","hNum5",nbins,low,high)
	hNum6 = TH1F("hNum6","hNum6",nbins,low,high)
	if (cuts3 is ""): cuts3 = cuts1 #only worry about cuts that have been defined
	if (cuts4 is ""): cuts4 = cuts1
	if (cuts5 is ""): cuts5 = cuts1
	if (cuts6 is ""): cuts6 = cuts1
	tNum1.Draw(varname+" >> hNum1",cuts1,"goff") #draw TTree into histograms
	tNum1.Draw(varname+" >> hNum2",cuts2,"goff")
	tNum1.Draw(varname+" >> hNum3",cuts3,"goff")
	tNum1.Draw(varname+" >> hNum4",cuts4,"goff")
	tNum1.Draw(varname+" >> hNum5",cuts5,"goff")
	tNum1.Draw(varname+" >> hNum6",cuts6,"goff")

	hNum1.Sumw2() #handle errors correctly
	hNum2.Sumw2()
	hNum3.Sumw2()
	hNum4.Sumw2()
	hNum5.Sumw2()
	hNum6.Sumw2()

	hNum1.Scale(1./hNum1.Integral()) #scale all histograms unity to make comparisons easier
	hNum2.Scale(1./hNum2.Integral())
	hNum3.Scale(1./hNum3.Integral())
	hNum4.Scale(1./hNum4.Integral())
	hNum5.Scale(1./hNum5.Integral())
	hNum6.Scale(1./hNum6.Integral())

	#set colors and line thickness
	hNum1.SetLineWidth(2)
	hNum1.SetLineColor(1)

	hNum2.SetLineWidth(2)
	hNum2.SetLineColor(2)

	hNum3.SetLineWidth(2)
	hNum3.SetLineColor(3)

	hNum4.SetLineWidth(2)
	hNum4.SetLineColor(4)

	hNum5.SetLineWidth(2)
	hNum5.SetLineColor(6)

	hNum6.SetLineWidth(2)
	hNum6.SetLineColor(7)

	if (logy): #set y scale properly
		hNum1.SetMaximum(max([hNum1.GetMaximum(),hNum2.GetMaximum(),hNum3.GetMaximum(),hNum4.GetMaximum(),
							  hNum5.GetMaximum(),hNum6.GetMaximum()])*80.)
	else:
		hNum1.SetMaximum(max([hNum1.GetMaximum(),hNum2.GetMaximum(),hNum3.GetMaximum(),hNum4.GetMaximum(),
							  hNum5.GetMaximum(),hNum6.GetMaximum()])*1.3)
	
	cst = TCanvas("cst","cst",10,10,700,700) #define canvas so it can be saved later
	
	if (logy): cst.SetLogy() #set y-axis to log scale if needed

	hNum1.Draw() #draw first histogram
	hNum1.SetTitle(title) #set title and axis labels
	hNum1.GetXaxis().SetTitle(xlabel)
	hNum1.GetYaxis().SetTitle("arb.")

	hNum2.Draw("same") #draw other histograms without overwriting
	if (label3 is not ""): hNum3.Draw("same")
	if (label4 is not ""): hNum4.Draw("same")
	if (label5 is not ""): hNum5.Draw("same")
	if (label6 is not ""): hNum6.Draw("same")

	leg = TLegend(0.1,0.75,0.9,0.9) #create legend
	leg.SetNColumns(2)
	entry1 = leg.AddEntry("hNum1",label1,"l")
	entry1.SetLineWidth(2)
	entry1.SetLineColor(1)
	entry2 = leg.AddEntry("hNum2",label2,"l")
	entry2.SetLineWidth(2)
	entry2.SetLineColor(2)

	if (label3 is not ""):
		entry3 = leg.AddEntry("hNum3",label3,"l")
		entry3.SetLineWidth(2)
		entry3.SetLineColor(3)

	if (label4 is not ""):
		entry4 = leg.AddEntry("hNum4",label4,"l")
		entry4.SetLineWidth(2)
		entry4.SetLineColor(4)

	if (label5 is not ""):
		entry5 = leg.AddEntry("hNum5",label5,"l")
		entry5.SetLineWidth(2)
		entry5.SetLineColor(6)

	if (label6 is not ""):
		entry6 = leg.AddEntry("hNum6",label6,"l")
		entry6.SetLineWidth(2)
		entry6.SetLineColor(7)

	leg.Draw()

	cst.Update()
	cst.SaveAs("EGMComp_"+outname+".png") #save as png, can add other lines for different file types (.pdf, .C, .root, ...)

def distplot(filename1, treeNum, varname, nbins, low, high, cuts, logy, xlabel, outname, title, markcut = 0.):
    
	f1 = TFile(filename1) #open file
	tNum1 = f1.Get(treeNum) #get TTree
	hNum1 = TH1F("hNum1","hNum1",nbins,low,high) #define histogram
	tNum1.Draw(varname+" >> hNum1",cuts,"goff") #draw TTree into histogram

	hNumM = TH1F("hNumM","hNumM",20,-10.*markcut,10.*markcut) #for drawing vertical line somewhere, doesnt get used in markcut is 0.

	hNum1.Sumw2() #for proper error handling
	hNum1.Scale(1./hNum1.Integral()) #scale to unity

	hNumM.SetBinContent(10,1.)
	hNumM.SetBinContent(11,1.)

	hNum1.SetLineWidth(2) #set colors/thicknes
	hNum1.SetLineColor(2)
	hNumM.SetLineWidth(2)
	hNumM.SetLineColor(1)
	hNumM.SetLineStyle(2)
	hNumM.SetFillStyle(3004)
	hNumM.SetFillColor(1)

#   if (logy):
#       hNum1.SetMaximum(hNum1.GetMaximum()*80.)
#   else:
#       hNum1.SetMaximum(hNum1.GetMaximum()*1.3)

	cst = TCanvas("cst","cst",10,10,700,700) #define canvas
	if (logy): cst.SetLogy() #set y-axis to log scale

	hNum1.Draw() #draw histogram
	hNum1.SetTitle(title) #set title and axis labels
	hNum1.GetXaxis().SetTitle(xlabel)
	hNum1.GetYaxis().SetTitle("arb.")
	if (markcut!=0.): hNumM.Draw("same hist") #draw line for cut if asked for

	latex = TLatex() #add text for cut
	latex.SetTextSize(0.05)
	latex.SetTextAlign(21)
	if (markcut!=0.): latex.DrawLatex(0.,hNum1.GetMinimum(),"98% Cut")

	cst.Update()
	cst.SaveAs("EGMDist_"+outname+".png") #save as png, can save as other types too (.pdf, .C, .root, ...)

# ---- Plotting ----

filename1 = "/home/drankin/TauID/GenNtuple_GluGluHToTauTau_M125_13TeV_powheg_pythia8.root"
filename2 = "/home/drankin/TauID/GenNtuple_ZprimeToTauTau_M-3000_TuneCP5_13TeV-pythia8-tauola.root"
filename3 = "/home/drankin/TauID/GenNtuple_QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root"
filename4 = "/home/drankin/TauID/GenNtuple_QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root"
filename5 = "/home/pharris/GenNtuple.root"

"""
prefix1 = "_GluGluHToTauTau"
prename1 = ", GluGlu to H to TauTau"
prefix2 = "_ZPrimeToTauTau"
prename2 = ", ZPrime to TauTau"
basecut = ""
barcut = basecut+" && abs(genjeteta)<1.5"
endcut = basecut+" && abs(genjeteta)>1.5"

varlist = [
	 ['genjetpt','gen jet p_{T}','genjetPt',50,0.,100.,False]
	,['abs(genjeteta)','gen jet |#eta|','genjetEta',50,0.,3.,False]
	,['genjetphi','gen jet #phi','genjetPhi',64,-3.2,3.2,False]
	,['abs(genjetid)','gen jet ID','genjetID',30,0,30,True]]


for var in varlist:
    distplot(filename2,"GenNtupler/gentree",var[0],var[3],var[4],var[5],"1"+barcut,var[6],var[1],var[2]+prefix2+"_Barrel",
			 var[1]+" Distribution Barrel"+prename2)

    distplot(filename2,"GenNtupler/gentree",var[0],var[3],var[4],var[5],"1"+endcut,var[6],var[1],var[2]+prefix2+"_Endcap",
			 var[1]+" Distribution Endcap"+prename2)

    compplot(filename2,"GenNtupler/gentree",var[0],var[3],var[4],var[5],var[6],var[1],var[2]+prefix2+"_Barrel",
			 var[1]+" Comparison Barrel"+prename2,"(abs(genjetid)<6 || abs(genjetid)==21)"+barcut,"Q/G Jet",
			 "abs(genjetid)==15"+barcut,"Tau","(abs(genjetid)>=6 && abs(genjetid)!=21 && abs(genjetid)!=15)"+barcut,"Other")

    compplot(filename2,"GenNtupler/gentree",var[0],var[3],var[4],var[5],var[6],var[1],var[2]+prefix2+"_Endcap",
			 var[1]+" Comparison Endcap"+prename2,"(abs(genjetid)<6 || abs(genjetid)==21)"+endcut,"Q/G Jet",
			 "abs(genjetid)==15"+endcut,"Tau","(abs(genjetid)>=6 && abs(genjetid)!=21 && abs(genjetid)!=15)"+endcut,"Other")
"""

# ---- Dummy Functions ----

def list_particles(filename, treeNum, eventNum):

    f1 = TFile(filename) #open file
    tNum = f1.Get(treeNum) #get TTree

    i = 0
    for event in tNum:
        if (i >= eventNum): break
        for j in range(len(event.genjetid)):
            print('Event '+str(i)+' jet '+str(j)+' contains')
            for k in range(len(event.genindex)):
                if (event.genindex[k] is j): print '\t ID '+str(event.genid[k])
        i += 1

def list_pt(filename, treeNum, eventNum, cutoff=0):
	
	f1 = TFile(filename) 
	tNum = f1.Get(treeNum)

	i = 0
	for event in tNum:
		if (i >= eventNum): break
		for j in range(len(event.genjetid)):
			if (abs(event.genjetid[j]) == 15): #only examine tau jets
				print('Event '+str(i)+' jet '+str(j)+' PT '+str(event.genjetpt[j]))
				for k in range(len(event.genindex)):
					if (event.genindex[k] == j and event.genpt[k] > cutoff):
						print('\t PT '+str(event.genpt[k])+' ID '+str(event.genid[k]))
		i += 1

def list_parent(filename, treeNum, eventNum):
	
	f1 = TFile(filename)
	tNum = f1.Get(treeNum)

	i = 0
	for event in tNum:
		if (i >= eventNum): break
		for j in range(len(event.genjetid)):
			if (abs(event.genjetid[j]) == 15):
				print('Event '+str(i)+' jet '+str(j)+' contains')
				for k in range(len(event.genindex)):
					if (event.genindex[k] == j):
						print('\t PartID '+str(event.genpartid[k])+' Parent '+str(event.genparent[k])+' Status '+str(event.genstatus[k]))
						print('\t Particle '+str(event.genid[k])+'\n')
		i += 1

def find_parent(filename, treeNum, eventNum):
	
	f1 = TFile(filename)
	tNum = f1.Get(treeNum)

	i = 0
	for event in tNum:
		if (i >= eventNum): break
		for j in range(len(event.genjetid)): 
			if (abs(event.genjetid[j]) == 15): #only examine tau jets
				print('Event '+str(i)+' jet '+str(j)+' contains') 
				for k in range(len(event.genindex)): 
					if (event.genindex[k] == j): 
						if (event.genstatus[k] == 1):
							index = k
							while ((event.genparent[index] != -2) and (abs(event.genid[index]) != 15)):
								index = event.genparent[index]
							if (abs(event.genid[index]) == 15):
								print('\t'+str(event.genid[index])+' decays to '+str(event.genid[k]))
		i += 1

#list_particles(filename5, "GenNtupler/gentree", 5)
#list_pt(filename5, "GenNtupler/gentree", 5)
#list_parent(filename5, "GenNtupler/gentree", 100)
find_parent(filename5, "GenNtupler/gentree", 5)



