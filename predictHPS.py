import sys, os
import math
import itertools
import ROOT
from ROOT import *

def lepton_decay(jet):
	"""
	Return: Boolean
	Input: jet  | List of jet particles (TLorentzVector, ID)
	"""
	electrons = [] # electron candidates
	muons = []     # muon candidates

	for vec, ID in jet:
		if abs(ID) == 11: electrons.append(vec)
		if abs(ID) == 13: muons.append(vec)
	
	if len(electrons) == 0:
		if len(muons) != 0: return True
		return False

	if len(electrons) % 2 != 0: return True
	
	return False

def gen_strips(ep_4v):
	"""
	Return: List of strip four vectors
	Input: ep_4v | List of electron/photon candidate four vectors
	"""
	total_num = len(ep_4v) # total number of candidates
	strip_4v = []         
	
	for i in range(min(6, total_num)): # consider 
		indices = [i]

		# initialize parameters
		strip_vec = ep_4v[i]
		for j in range(i + 1, 6):
			eta = strip_vec.Eta()

			# if particle in strip
			if ((eta-0.025) < ep_4v[j].Eta() < (eta+0.025)) and (ep_4v[j].DeltaPhi(strip_vec) < 0.10):
				indices.append(j)

				# initialize new parameters
				new_eta = 0; new_phi = 0; total_Pt = 0; total_E = 0

				for idx in indices: # update parameters
					new_eta += ep_4v[idx].Eta()*ep_4v[idx].Pt()
					new_phi += ep_4v[idx].Phi()*ep_4v[idx].Pt()
					total_Pt += ep_4v[idx].Pt()
					total_E += ep_4v[idx].E()
				
				strip_vec.SetPtEtaPhiE(total_Pt, new_eta/total_Pt, new_phi/total_Pt, total_E)
		
		if strip_vec.Pt() > 2.5:
			strip_4v.append(strip_vec)
	
	return strip_4v

def gen_candidates(jet, hadron_cut, ep_cut):
	"""
	Return: (List of hadron candidates, List of strips)
	Input: jet        | List of jet particles (TLorentzVector, ID)
		   hadron_cut | Pt cutoff for hadrons
		   ep_cut     | Pt cutoff for electrons/photons
	"""
	hadron_4v = [] # hadron candidate 4vec
	ep_4v = []     # electron/photon candidate 4vec

	for vec, ID in jet:
		if (ID == 211) and (vec.Pt() > hadron_cut): 
			hadron_4v.append(vec)
		if (abs(ID) in [11, 22]) and (vec.Pt() > ep_cut):
			ep_4v.append(vec)

	hadron_4v = sorted(hadron_4v, key = lambda x: x.Pt())[-5:]  # sort hadrons
	ep_4v = sorted(ep_4v, key = lambda x: x.Pt())[::-1]         # sort strips
	
	return (hadron_4v, gen_strips(ep_4v))

def hypothesis1(event, hadron_4v):
	"""
	Return: Guesses for h+/-, h-/+, h-/+ hypothesis
	Input: event     | ROOT event
		   hadron_4v | List of hadron candidate (indexes, four vectors)
	"""
	guesses = []

	for triple in itertools.combinations(hadron_4v, 3):
		if (any(event.genid[e[0]] < 0 for e in triple) and any(event.genid[e[0]] > 0 for e in triple)):
			if (not any(event.gencharge[t[0]] == 0 for t in triple)): # charge check
				candidates = (triple[0][1], triple[1][1], triple[2][1])
				vec_sum = candidates[0] + candidates[1] + candidates[2]

				if 0.8 < vec_sum.M() < 1.5: # mass check
					DeltaR_cut = 3.0/vec_sum.Pt()
					if (not any(vec_sum.DeltaR(c2) > DeltaR_cut for c2 in candidates)): # deltaR check
						pt_sum = candidates[0].Pt() + candidates[1].Pt() + candidates[2].Pt()
						guesses.append((pt_sum, vec_sum))
	return guesses

def hypothesis2(event, hadron_4v, strip_4v):
	"""
	Return: Guesses for h+/-, pi0, pi0 hypothesis
	Input: event     | ROOT event
	       hadron_4v | List of hadron candidate (indexes, four vectors)
		   strip_4v  | List of strip four vectors
	"""
	guesses = []

	for h in hadron_4v
		if (abs(event.gencharge[h[0]]) == 1): # charge check
			for pair in itertools.combinations(strip_4v, 2):
				candidates = (h[1], pair[0], pair[1])
				vec_sum = candidates[0] + candidates[1] + candidates[2]
				cutoff = min(max(1.2*math.sqrt(vec_sum.Pt()/100), 1.2), 4.0)
	
				if 0.4 < vec_sum.M() < cutoff: # mass check
					DeltaR_cut = 3.0/vec_sum.Pt()
					if (not any(vec_sum.DeltaR(c2) > DeltaR_cut for c2 in candidates)): # deltaR check 
						pt_sum = candidates[0].Pt() + candidates[1].Pt() + candidates[2].Pt()
						guesses.append((pt_sum, vec_sum))
	return guesses

def hypothesis3(event, hadron_4v, strip_4v):
	"""
	Return: Guesses for h+/-, pi0 hypothesis
	Input: event     | ROOT event
	       hadron_4v | List of hadron candidate (indexes, four vectors)
		   strip_4v  | List of strip four vectors
	"""
	guesses = []

	for h in hadron_4v:
		if (abs(event.gencharge[h[0]]) == 1): # charge check
			for s in strip_4v:
				candidates = (h[1], s)
				vec_sum = candidates[0] + candidates[1]
				cutoff = min(max(1.3*math.sqrt(vec_sum.Pt()/100), 1.3), 4.2)

				if 0.3 < vec_sum.M() < cutoff: # mass check
					DeltaR_cut = 3.0/vec_sum.Pt()
					if (not any(vec_sum.DeltaR(c2) > DeltaR_cut for c2 in candidates)): 
						pt_sum = candidates[0].Pt() + candidates[1].Pt()
						guesses.append((pt_sum, vec_sum))
	return guesses

def hypothesis4(event, hadron_4v, strip_4v):
	"""
	Return: Guesses for h+/-, pi0 hypothesis
	Input: event     | ROOT event
	       hadron_4v | List of hadron candidate (indexes, four vectors)
		   strip_4v  | List of strip four vectors
	"""
	guesses = []
	
	if (len(strip_4v) == 0 and len(hadron_4v) == 1):
		if (abs(event.gencharge[hadron_4v[0][0]]) == 1): # charge check
			guesses.append((hadron_4v[0][1].Pt(), hadron_4v[0][1]))
	return guesses

def isolation(event, pair, cutoff):
	"""
	Return: Boolean
	Input: event  | ROOT event
		   pair   | Tuple: (pt_sum, 4vec)
		   cutoff | Isolation cutoff
	"""	
	vec_ptsum, vec = pair
	ptsum = 0

	for k, _ in enumerate(event.genisoid):
		if ((event.genisoid[k] > 40) and (abs(event.genisocharge[k]) == 1)) or (event.genisoid[k] == 22):
			vec_test = TLorentzVector()
			vec_test.SetPtEtaPhiE(event.genisopt[k], event.genisoeta[k], event.genisophi[k], event.genisoenergy[k])
			if ((vec.DeltaR(vec_test) < 0.4) and (vec.Pt() > 0.5)):
				ptsum += vec_test.Pt()
	
	iso = (ptsum - vec_ptsum)/vec.Pt()

	if (iso > cutoff): return (iso, False)

	return (iso, True)

def predict(array, #TODO):	
	
	for event in range(array):
		for jet, tau in jet_candidates:

		if lepton_decay(jet):
			continue

		# generate list of hadronic tau candidates
		# generate list of photons/electrons



		
