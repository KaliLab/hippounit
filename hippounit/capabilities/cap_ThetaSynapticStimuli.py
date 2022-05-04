import sciunit
from sciunit import Capability

class ThetaSynapticStimuli(sciunit.Capability):
	"""Indicates that the model receives multiple synaptic inputs mimicing synaptic inputs (as pathway stimuli) during exploratory theta. The model must be initialised befor applying this method."""

	def activate_theta_stimuli(self, dend_loc, AMPA_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train):
		""" Doesn't return anything, places, sets the parameters and activates synapses on the appropriate dendritic locations ."""
		raise NotImplementedError()

