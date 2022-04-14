import sciunit
from sciunit import Capability

class RunSimulation_ReturnTraces(sciunit.Capability):
	"""Runs the simulation on an initialised model and records the voltage trace on the some and one selected dendritic location. Any stimulus must by set and activted previously """

	def run_simulation(self, dend_loc, recording_loc, tstop):
		""" Returns time vector and voltage traces of the soma and a selected location (t, v, v_dend - as numpy arrays)."""

		raise NotImplementedError()

