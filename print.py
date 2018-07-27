import sys, os, math
import ROOT, sys, os, re, string
from ROOT import *

filename1 = "/data/t3home000/rinik/GenNtuple_GluGluHToTauTau_M125_13TeV_powheg_pythia8.root"
filename2 = "/data/t3home000/rinik/GenNtuple_ZprimeToTauTau_M-3000_TuneCP5_13TeV-pythia8-tauola.root"
filename3 = "/data/t3home000/rinik/GenNtuple_QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root"
filename4 = "/data/t3home000/rinik/GenNtuple_QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root"
filename5 = "/mnt/hadoop/scratch/drankin/TauID/GenNtuple_parentage_GluGluHToTauTau_M125_13TeV_powheg_pythia8.root"
filename6 = "/mnt/hadoop/scratch/drankin/TauID/GenNtuple_parentage_ZprimeToTauTau_M-3000_TuneCP5_13TeV-pythia8-tauola.root"

# ---- Print Functions ----

def list_particles(filename, treeNum, totalNum):

    f1 = TFile(filename)   #open file
    tNum = f1.Get(treeNum) #get TTree

    eventNum = 0
    for event in tNum:
        if (eventNum >= totalNum): break
        for j in range(len(event.genjetid)):
            print('Event '+str(eventNum)+' jet '+str(j)+' contains')
            for k in range(len(event.genindex)):
                if (event.genindex[k] is j): print('\t ID '+str(event.genid[k]))
        eventNum += 1

def list_taus(filename, treeNum, totalNum, cutoff=0):
	
	f1 = TFile(filename) 
	tNum = f1.Get(treeNum)

	eventNum = 0
	for event in tNum:
		if (eventNum >= totalNum): break
		for j in range(len(event.genjetid)):
			if (abs(event.genjetid[j]) == 15): #only examine tau jets
				print('Event '+str(eventNum)+' jet '+str(j)+' contains')
				for k in range(len(event.genindex)):
					if (event.genindex[k] == j and event.genpt[k] > cutoff):
						print('\t ID ' +str(event.genid[k]))
						print('\t Pt ' +str(event.genpt[k]))
						print('\t Eta '+str(event.geneta[k]))
						print('\t Phi '+str(event.genphi[k]))
						print('\t E '  +str(event.genenergy[k])+'\n')
		eventNum += 1

def list_parent(filename, treeNum, totalNum):
	
	f1 = TFile(filename)
	tNum = f1.Get(treeNum)

	eventNum = 0
	for event in tNum:
		if (eventNum >= totalNum): break
		for j in range(len(event.genjetid)):
			if (abs(event.genjetid[j]) == 15): #only examine tau jets
				print('Event '+str(eventNum)+' jet '+str(j)+' contains')
				for k in range(len(event.genindex)):
					if (event.genindex[k] == j):
						print('\t PartID '+str(event.genpartid[k])+' Parent '+str(event.genparent[k])+
							  ' Status '+str(event.genstatus[k]))
						print('\t Particle '+str(event.genid[k])+'\n')
		eventNum += 1

def find_parent(filename, treeNum, totalNum):
	
	f1 = TFile(filename)
	tNum = f1.Get(treeNum)

	eventNum = 0
	jetNum = 0
	counter = 0
	
	for event in tNum:
		if (eventNum >= totalNum): break
		for j in range(len(event.genjetid)): 
			claim_tau = False; actual_tau = False
			
			if (abs(event.genjetid[j]) == 15): #only examine tau jets
				claim_tau = True	
			#print('Event '+str(eventNum)+' jet '+str(j)+' contains') 	
			
			for k in range(len(event.genindex)): 
				if (event.genindex[k] == j): 
					if (event.genstatus[k] == 1):
						index = k
						while ((event.genparent[index] != -2) and (abs(event.genid[index]) != 15)):
							index = event.genparent[index]
						if (abs(event.genid[index]) == 15):
							#print('\t'+str(event.genid[index])+' decays to '+str(event.genid[k]))
							actual_tau = True
			
			if claim_tau != actual_tau:
				counter += 1
				print(claim_tau, actual_tau)
				print('Event '+str(eventNum)+' jet '+str(j)+' misidentified')
			
			jetNum += 1
		eventNum += 1
	
	print('Number of misidentified events: ' + str(counter))
	print('Number of jets: ' + str(jetNum))

# ---- Run Functions ----

#list_particles(filename1, "GenNtupler/gentree", 5)
#list_taus(filename1, "GenNtupler/gentree", 5)
#list_parent(filename5, "GenNtupler/gentree", 10)
find_parent(filename5, "GenNtupler/gentree", 1000)


