import sciunit
from sciunit import Capability

class ReceivesEPSCstim(sciunit.Capability):
    """
    Indicates that the model receives an excitatory post-synaptic current (EPSC) shaped input
    """

    def run_EPSCstim(self, dend_loc, weight, tau1, tau2):
        """
        This function must be implemented by the model.

        Must return numpy arrays containing the time and voltage values (at the soma and at the synaptic location )

        Parameters
        ----------
        dend_loc : list
            containing the name of the section (string) and the location (float) where input is received. Eg.: ['dendrite[3]', 0.5]
        weight : float
            weight of the synaptic input
        tau1 : float
            rising time constant of the synaptic input
        tau2 : float
            decay time constant of the synaptic input
        """

        raise NotImplementedError()

    def run_EPSC_stim_get_vm(self, dend_loc, weight, tau1, tau2):
        """
        This function is called by the test and calls the run_EPSCstim() function.
        """

        t, v, v_dend = self.run_EPSCstim(dend_loc, weight, tau1, tau2)

        return t, v, v_dend
