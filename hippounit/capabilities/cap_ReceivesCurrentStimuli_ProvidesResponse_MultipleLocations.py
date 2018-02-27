import sciunit
from sciunit import Capability


class ReceivesSquareCurrent_ProvidesResponse_MultipleLocations(sciunit.Capability):
	"""Indicates that current can be injected into the model as
    a square pulse. And records at multiple locations."""

	def inject_current_record_respons_multiple_loc(self, amp, delay, dur, section_stim, loc_stim, dend_locations):
		""" Must return numpy arrays containing the time and voltage values"""
		raise NotImplementedError()

	def get_multiple_vm(self, amp, delay, dur, section_stim, loc_stim, dend_locations):
		# v : dictionary -  keys: dendritic location, values: the voltage trace for each recording locations
		t, v_stim, v = self.inject_current_record_respons_multiple_loc(amp, delay, dur, section_stim, loc_stim, dend_locations)
		return t, v_stim, v
