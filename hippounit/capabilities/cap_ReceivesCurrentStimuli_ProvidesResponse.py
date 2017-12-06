import sciunit
from sciunit import Capability


class ReceivesSquareCurrent_ProvidesResponse(sciunit.Capability):
	"""Indicates that current can be injected into the model as
    a square pulse. """

	def inject_current(self, amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec):
		""" Must return numpy arrays containing the time and voltage values"""
		raise NotImplementedError()

	def get_vm(self, amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec):

		t, v = self.inject_current(amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec)
		return t, v
