
import sciunit
from sciunit import Capability
import multiprocessing


class ProvidesRecordingLocationsOnTrunk(sciunit.Capability):
	""" Indicates that the model provides a list of locations on the trunk (primary apical dendrite) to be tested"""

	def find_trunk_locations(self,distances, tolerance):
		""" Must provide a dictionary - keys: distances, values: corresponding locations on the trunk (primary apical dendrite) in list
		at 50, 105, 250, 350 um distances from the soma
		The form must be: dend_loc = (dist1, ['trunk_segment1_1',location],['trunk_segment1_2',location]), (dist2, ['trunk_segment2',location]),(dist3, ['trunk_segment3',location]), (dist4, ['trunk_segment4',location])
		E.g. : OrderedDict([(50, ['dendrite[0]', 0.6956994222486329]), (150, ['dendrite[81]', 0.5557523508251703]), (250, ['dendrite[109]', 0.33250043844278565])]) """

		raise NotImplementedError()

	def find_trunk_locations_multiproc(self, distances, tolerance):
		""" Used to keep all NEURON related tasks in independent processes, to avoid errors like 'template can not be redefined'"""
		pool_trunk = multiprocessing.Pool(1, maxtasksperchild = 1)
		self.dend_locations, actual_distances = pool_trunk.apply(self.find_trunk_locations, (distances, tolerance,))  # this way model.dend_loc gets the values
		pool_trunk.terminate()
		pool_trunk.join()
		del pool_trunk
		return self.dend_locations, actual_distances
