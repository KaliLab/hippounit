import sciunit
from sciunit import Capability

class ReceivesSynapseGivenPathway(sciunit.Capability):
    """Indicates that the model receives synapse on via a given pathway (eg. PP or SC)"""

    def run_syn_pathway(self, dend_loc, AMPA_weight, pathway):
        """ Must return numpy arrays containing the time and voltage values (at the soma and at the synaptic location )"""
        raise NotImplementedError()

    def run_synapse_pathway_get_vm(self, dend_loc, AMPA_weight, pathway):

        t, v, v_dend = self.run_syn_pathway(dend_loc, AMPA_weight, pathway)

        return t, v, v_dend
