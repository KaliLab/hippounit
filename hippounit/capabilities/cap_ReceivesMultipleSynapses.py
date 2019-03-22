import sciunit
from sciunit import Capability

class ReceivesMultipleSynapses(sciunit.Capability):
	"""Indicates that the model receives one or multiple synapses"""

	def run_multiple_syn(self, dend_loc, interval, number, weight):
		""" Must return numpy arrays containing the time and voltage values (at the soma and at the synaptic location )"""
		raise NotImplementedError()

	def run_multiple_synapse_get_vm(self, dend_loc, interval, number, weight):

		t, v, v_dend = self.run_multiple_syn(dend_loc, interval, number, weight)

		return t, v, v_dend
