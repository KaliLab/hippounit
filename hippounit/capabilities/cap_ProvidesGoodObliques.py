
import sciunit
from sciunit import Capability
import multiprocessing


class ProvidesGoodObliques(sciunit.Capability):
    """ Indicates that the model provides a list of oblique dendrites and locations to be tested"""

    def find_good_obliques(self, trunk_origin):
        """
        This function must be implemented by the model.

        Must provide a list of oblique dendrites
        that meet the criteria of the experimental protocol (Losonczy, Magee 2006),
        and also proximal and distal locations on them.
        Criteria: originate from the trunk, have no child, close to the soma (at most 120 microns)

        The form must be: dend_loc = [['name_of_dend1',prox_location, "prox"],['name_of_dend1',dist_location, "dist"],['name_of_dend2',prox_location, "prox"] ['name_of_dend2',dist_location, "dist"]]
        E.g. : [['CCell[0].apic[47]', 0.5, "prox"], ['CCell[0].apic[47]', 0.8333333333333333, "dist"]]

        Parameters
        ----------
        trunk_origin : list
            first element : name of the section from which the trunk originates, second element : position on section (E.g. ['soma[5]', 1]). If not set by the user, the end of the default soma section is used.
        """

        raise NotImplementedError()

    def find_obliques_multiproc(self, trunk_origin):
        """
        This function is called by the test and calls the find_good_obliques() function.
        Used to keep all NEURON related tasks in independent processes, to avoid errors like 'template can not be redefined'
        """
        pool_obl = multiprocessing.Pool(1, maxtasksperchild = 1)
        self.dend_loc = pool_obl.apply(self.find_good_obliques, (trunk_origin,))  # this way model.dend_loc gets the values
        pool_obl.terminate()
        pool_obl.join()
        del pool_obl
