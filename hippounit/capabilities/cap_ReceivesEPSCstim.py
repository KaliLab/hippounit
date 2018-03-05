import sciunit
from sciunit import Capability

class ReceivesEPSCstim(sciunit.Capability):
	"""Indicates that the model receives synapse"""

	def run_EPSCstim(self, dend_loc, weight, tau1, tau2):
		""" Must return numpy arrays containing the time and voltage values (at the soma and at the synaptic location )"""
		raise NotImplementedError()

	def run_EPSC_stim_get_vm(self, dend_loc, weight, tau1, tau2):

		t, v, v_dend = self.run_EPSCstim(dend_loc, weight, tau1, tau2)

		return t, v, v_dend
