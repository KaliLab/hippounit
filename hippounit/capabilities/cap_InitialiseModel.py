import sciunit
from sciunit import Capability

class InitialiseModel(sciunit.Capability):
	"""Initialises the model. In the case of NEURON models: compiles (if needed) and loads the mod files, and loads the hoc file"""

	def initialise(self):
		""" Doesn't return anything, but loads the model to be ready for simulations."""
		raise NotImplementedError()

