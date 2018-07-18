import math, itertools, time
import ROOT, sys, os, re, string
from ROOT import *

def update_progress(progress):
	"""
	Return: Progress Bar
	Input: progress | Fraction completed from 0->1
	"""
	barLength = 20
	status = ''

	if (progress >= 1):
		progress = 1
		status = "Done...\r\n"

	block = int(round(barLength*progress))
	text = '\rPercent: [{0}] {1}% {2}'.format('#'*block + '-'*(barLength - block), progress*100, status)
	sys.stdout.write(text)
	sys.stdout.flush()

def remove_leptons(event, jet_num):
	"""
	Return: Boolean
	Input: event   | ROOT event
	       jet_num | Jet number
	"""
	ep = [] # electron/positron candidates
	other = [] # other lepton candidates

	for k, _ in enumerate(event.genindex):
		if (event.genindex[k] == jet_num) and (10 < abs(event.genid[k]) < 20): 
			if (abs(event.genid[k]) != 11):
				other.append(event.genid[k])
			else:
				ep.append(event.genid[k])
	
	if (len(ep) == 0):
		if (len(other) != 0): return True
		return False

	if (len(ep) % 2 != 0): return True
	
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
		for j in range(i + 1, total_num):
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

def gen_candidates(event, jet_num, hadron_cut, ep_cut):
	"""
	Return: (List of hadron candidates, List of strips)
	Input: event      | ROOT event
	       jet_num    | Jet number
		   hadron_cut | Pt cutoff for hadrons
		   ep_cut     | Pt cutoff for electrons/photons
	"""

	hadron_4v = [] # hadron candidate 4vec
	ep_4v = []     # electron/photon candidate 4vec

	for k, _ in enumerate(event.genindex):
		if (event.genindex[k] == jet_num):
			vec = TLorentzVector()
			vec.SetPtEtaPhiE(event.genpt[k], event.geneta[k], event.genphi[k], event.genenergy[k])

			if ((abs(event.genid[k]) > 40) and (event.genpt[k] > hadron_cut)):
				hadron_4v.append((k, vec))
			if ((event.genid[k] in ([22, 11, -11])) and (event.genpt[k] > ep_cut)):
				ep_4v.append(vec)

			hadron_4v = sorted(hadron_4v, key = lambda x: x[1].Pt())[-5:] # sort hadrons
			ep_4v = sorted(ep_4v, key = lambda x: x.Pt())[::-1]           # sort strips

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

	for h in hadron_4v:
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

def predict(filename, treeNum, iterations, hadron_cut = 0, ep_cut = 0, iso_cutoff = float('inf')):
	"""
	Return: List of predictions
	Input: filename   | ROOT file
	       treeNum    | ROOT tree number
	       iterations | Number of events to consider
		   hadron_cut | Pt cutoff for hadrons
		   ep_cut     | Pt cutoff for electrons/photons
	"""
	rf = TFile(filename)      # open file
	tree = rf.Get(treeNum)    # get TTree
	
	predictions_tau = [] # list of (ID, 4vec)
	predictions_other = [] # list of (ID, 4vec)
	isolation_list = [] # list of (ID, iso)

	accuracy = [0., 0., 0., 0., 0., 0.]
	event_num = 0.
	
	t0 = time.time()
		
	for event in tree:
		if (event_num >= iterations): break
		for jet_num, _  in enumerate(event.genjetid):

			if remove_leptons(event, jet_num): # disregard lepton decays
				continue

			tau_present = False; max_pt = []

			hadron_4v, strip_4v = gen_candidates(event, jet_num, hadron_cut, ep_cut)
			
			guesses1 = hypothesis1(event, hadron_4v)
			guesses2 = hypothesis2(event, hadron_4v, strip_4v)
			guesses3 = hypothesis3(event, hadron_4v, strip_4v)
			guesses4 = hypothesis4(event, hadron_4v, strip_4v)

			if guesses1: max_pt.append(('1', max(guesses1, key = lambda x: x[1].Pt())))
			if guesses2: max_pt.append(('2', max(guesses2, key = lambda x: x[1].Pt())))
			if guesses3: max_pt.append(('3', max(guesses3, key = lambda x: x[1].Pt())))	
			if guesses4: max_pt.append(('4', max(guesses4, key = lambda x: x[1].Pt())))
		
			# isolist contains (all tau reconstructions (tau + not tau) +  

			if max_pt:
				max_ = max(max_pt, key = lambda x: x[1][1].Pt())
				if (max_[1][1].Pt() > 20): # Only consider candidates with Pt > 20
					iso, truth = isolation(event, max_[1], iso_cutoff)
					isolation_list.append((event.genjetid[jet_num], iso))

					if truth:
						tau_present = True
						predictions_tau.append((event.genjetid[jet_num], max_[0], max_[1][1]))
				
					else:
						predictions_other.append((event.genjetid[jet_num], max_[1][1]))
			
			elif (event.genjetpt[jet_num] > 20): # Pt > 20 (Adjust value?)
				vec = TLorentzVector()
				vec.SetPtEtaPhiE(event.genjetpt[jet_num], event.genjeteta[jet_num], event.genjetphi[jet_num], event.genjetenergy[jet_num])
				predictions_other.append((event.genjetid[jet_num], vec))
				isolation_list.append((event.genjetid[jet_num], 1000)) # arbitrary large iso
			
			# accuracy
			if (abs(event.genjetid[jet_num]) == 15) and (event.genjetpt[jet_num] > 20):
				if tau_present:
					accuracy[0] += 1 # total correct
					accuracy[2] += 1 # total taus identified
				accuracy[3] += 1 # total taus
			elif (event.genjetpt[jet_num] > 20):
				if (not tau_present):
					accuracy[0] += 1 # total correct
					accuracy[4] += 1 # total not taus identified
				accuracy[5] += 1 # total not taus
			accuracy[1] += 1 # total
			
			"""
			# debug code here
			"""

		event_num += 1
		update_progress(event_num/iterations)

	t1 = time.time()

	total_accuracy = accuracy[0]/accuracy[1]*100
	if accuracy[3]:
		efficiency = accuracy[2]/accuracy[3]*100
		false_positive_rate = (1 - accuracy[4]/accuracy[5])*100

	print 'Runtime: ', (t1 - t0), 'seconds\n'
	print 'Total Accuracy: ', total_accuracy, '%'

	if accuracy[3]:
		print 'Efficency: ', efficiency, '%'
		print 'False Positive Rate: ', false_positive_rate, '%\n'
	
	#return isolation_list
	return predictions_tau, predictions_other



