
import sciunit
from sciunit import Capability
import multiprocessing


class ProvidesRandomDendriticLocations(sciunit.Capability):
    """ Indicates that the model provides a list of randomly selected locations on the trunk (primary apical dendrite) to be tested"""

    def get_random_locations(self, num, seed, dist_range, trunk_origin):
        """
        This function must be implemented by the model.

        Must return a list of lists [dendrite, seg.x]. Eg. : [['dendrite[31]', 0.5], ['dendrite[117]', 0.8333333333333333], ['dendrite[117]', 0.16666666666666666], ['dendrite[77]', 0.5], ['dendrite[99]', 0.5]],
        and a dictionary where the keys are the locations, the value is the actual distance of the location from the soma. Eg.: {('dendrite[95]', 0.5): 191.4537639215934, ('dendrite[91]', 0.5): 186.10161451767556}

        Parameters
        ----------
        num : int
            the number of dendritic locations to be selected
        seed : float
            the random seed
        dist_range : list
                containing the mimnimum and maximum distance from the soma. Eg.: [50,150]
        trunk_origin : list
            first element : name of the section from which the trunk originates, second element : position on section (E.g. ['soma[5]', 1]). If not  set by the user, the end of the default soma section is used.
        """

        raise NotImplementedError()

    def get_random_locations_multiproc(self, num, seed, dist_range, trunk_origin):
        """
        This function is called by the test and calls the get_random_locations() function.
        Used to keep all NEURON related tasks in independent processes, to avoid errors like 'template can not be redefined'
        """
        pool = multiprocessing.Pool(1, maxtasksperchild = 1)
        self.dend_locations, actual_distances = pool.apply(self.get_random_locations, (num, seed, dist_range, trunk_origin,))  # this way model.dend_loc gets the values
        pool.terminate()
        pool.join()
        del pool
        return self.dend_locations, actual_distances
