import sciunit
from sciunit import Capability

class ReceivesSynapse(sciunit.Capability):
	"""Indicates that the model receives synapse"""

	def run_syn(self, dend_loc, interval, number, AMPA_weight):
		""" Must return numpy arrays containing the time and voltage values (at the soma and at the synaptic location )"""
		raise NotImplementedError()

	def run_synapse_get_vm(self, dend_loc, interval, number, AMPA_weight):

		t, v, v_dend = self.run_syn(dend_loc, interval, number, AMPA_weight)

		return t, v, v_dend
