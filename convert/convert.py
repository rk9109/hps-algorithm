import argparse
import numpy as np
import rootpy
from ROOT import *
from utilities import progress

def create_jets(event):
	"""
	Create jets from particle flow inputs
	Return: List of jet candidates (seed, [jet_particles])
	"""
	part_candidates = []
	for part in event.pf:
		part_candidates.append(part)
	
	# Sort particle candidates by Pt
	part_candidates_sorted = sorted(part_candidates, key=lambda x: x[0].Pt())[::-1]

	used_candidates = []
	jet_candidates = []
	for part in part_candidates_sorted:
		if len(jet_candidates) >= 5: break # Maximum 5 jets/event
		# Seed critera: Charged hadron
		#               Pt > 20
		#               Eta < 2.5
		if (part[1] == 211) and (abs(part[0].Eta() < 2.5)) and (part not in used_candidates):
			jet = []; seed = part[0]
			jet_vecSum = TLorentzVector()
			jet_vecSum.SetPtEtaPhiE(0., 0., 0., 0.)

			for cand in part_candidates:
				if (seed.DeltaR(cand[0]) < 0.4) and (cand not in used_candidates):
					jet_vecSum += cand[0]
					jet.append(cand)
					used_candidates.append(cand)

			if jet_vecSum.Pt() > 20:
				jet_candidates.append((seed, jet))
	
	return jet_candidates

def create_taus(event):
	"""
	Create taus from gen inputs
	Return: List of tau candidate 4-vectors 
	"""
	indices = []           # List of indices
	gen_candidates = []    # List of gen candidates: (cand, index)
	for idx, cand in enumerate(event.gen):
		index = event.tauIndex[idx]
		if index not in indices: indices.append(index)
		gen_candidates.append((cand, index))

	tau_candidates = []    # List of tau 4-vectors	
	for index in indices:
		tau_vecSum = TLorentzVector()
		tau_vecSum.SetPtEtaPhiE(0., 0., 0., 0.)

		hadron_decay = False
		for gen_cand in gen_candidates:
			if gen_cand[1] == index:
				tau_vecSum += gen_cand[0][0]
				if abs(gen_cand[0][1]) == 211:
					hadron_decay = True

		# Tau criteria: Hadronic decay
		#               Pt < 20
		#               Eta < 2.5
		if (hadron_decay) and (tau_vecSum.Pt() > 20) and (abs(tau_vecSum.Eta()) < 2.5):
			tau_candidates.append(tau_vecSum)
	
	return tau_candidates

def match_taus(jet_candidates, tau_candidates):
	"""
	Match reconstructed taus to jets
	Return: List of jets ([jet_particles], tau)
	"""
	jets = []              # List of jets
	used_candidates = []
	for seed, jet in jet_candidates:
		tau = None
		for vec in tau_candidates:
			if (seed.DeltaR(vec) < 0.4) and (vec not in used_candidates):
				tau = vec
				used_candidates.append(vec)
				break
		jets.append((jet, tau))

	return jets

def convert_tree(tree, number=None):
	"""
	Particle flow inputs => HPS inputs
	"""
	event_num = 0.         # Event counter 
	if number: total_num = number
	else: total_num = int(tree.GetEntries())

	jets_array = []

	for event in tree:
		if event_num == total_num: break

		jet_candidates = create_jets(event)
		tau_candidates = create_taus(event)
		jets = match_taus(jet_candidates, tau_candidates)
		jets_array.append(jets)

		event_num += 1
	
	#cleanup
	jets_nparray = np.array(jets_array)
	return jets_nparray
			
if __name__ == "__main__":	
	parser = argparse.ArgumentParser()
	parser.add_argument('filename', help='ROOT filename')
	parser.add_argument('-n', '--number', dest='number', default=0, help='number of events')
	parser.add_argument('-t', '--tree', dest='tree', default='dumpP4/objects', help='tree name')
 	options = parser.parse_args()

	# Get ROOT TTree	
	filename = options.filename
	rf = TFile(filename)              # open file
	tree = rf.Get(options.tree)       # get TTree

	# Generate and save files
