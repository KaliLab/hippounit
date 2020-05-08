
import sciunit
from sciunit import Capability
import multiprocessing


class ProvidesRecordingLocationsOnTrunk(sciunit.Capability):
    """ Indicates that the model provides a list of locations on the trunk (primary apical dendrite) to be tested"""

    def find_trunk_locations(self,distances, tolerance, trunk_origin):
        """
        This function must be implemented by the model.

        Must return two dictionaries

        (1) keys: distances, values: corresponding locations on the trunk (primary apical dendrite) in list
        at 50, 105, 250, 350 um distances from the soma
        The form must be: dend_locations = (dist1, ['trunk_segment1_1',location],['trunk_segment1_2',location]), (dist2, ['trunk_segment2',location]),     (dist3, ['trunk_segment3',location]), (dist4, ['trunk_segment4',location])
        E.g. : OrderedDict([(50, ['dendrite[0]', 0.6956994222486329]), (150, ['dendrite[81]', 0.5557523508251703]), (250, ['dendrite[109]', 0.33250043844278565])])

        (2) keys: locations on the trunk, values: its actual distance from the soma.
        The form must be: actial distances = (('trunk_segment1_1',location), distance), (('trunk_segment2',location), distance),    (['trunk_segment3',location], distance)
        E.g. : {('dendrite[95]', 0.5): 191.4537639215934, ('dendrite[91]', 0.5): 186.10161451767556}

        Parameters
        ----------
        distances : list
            list of dinstances from the soma to be considered. Eg.: [50,150,250,350]
        tolerance : int
            number indicating the range around the distance values in the distances list to be considered. If tolerance = 20, first distance range will be 50+-20 um.
        trunk_origin : list
            first element : name of the section from which the trunk originates, second element : position on section (E.g. ['soma[5]', 1]). If not set by the user, the end of the default soma section is used.      
        """

        raise NotImplementedError()

    def find_trunk_locations_multiproc(self, distances, tolerance, trunk_origin):
        """
        This function is called by the test and calls the find_trunk_locations() function.
        Used to keep all NEURON related tasks in independent processes, to avoid errors like 'template can not be redefined'
        """
        pool_trunk = multiprocessing.Pool(1, maxtasksperchild = 1)
        self.dend_locations, actual_distances = pool_trunk.apply(self.find_trunk_locations, (distances, tolerance, trunk_origin,))  # this way model.dend_loc gets the values
        pool_trunk.terminate()
        pool_trunk.join()
        del pool_trunk
        return self.dend_locations, actual_distances
