import sciunit
from sciunit import Capability

class ReceivesSynapse(sciunit.Capability):
    """Indicates that the model receives synapse"""

    def run_syn(self, dend_loc, interval, number, AMPA_weight):
        """
        This function must be implemented by the model.

        Must return numpy arrays containing the time and voltage values (at the soma and at the synaptic location )

        Parameters
        ----------
        dend_loc : list
            containing the name of the section (string) and the location (float) where input is received. Eg.: ['dendrite[3]', 0.5]
        number: int
            number of synaptic input
        interval : float
            time interval between the synaptic inputs
        weight : float
            weight of the synaptic input
        """

        raise NotImplementedError()

    def run_synapse_get_vm(self, dend_loc, interval, number, AMPA_weight):

        """
        This function is called by the test and calls the run_syn() function.
        """
        t, v, v_dend = self.run_syn(dend_loc, interval, number, AMPA_weight)

        return t, v, v_dend
