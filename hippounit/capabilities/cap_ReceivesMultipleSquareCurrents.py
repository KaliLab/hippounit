import sciunit
from sciunit import Capability


class ReceivesMultipleSquareCurrents(sciunit.Capability):
    """Indicates that current can be injected into the model as
    a square pulse. """

    def activate_current_stimuli(self, amp, delay, dur, number, interval_bw_stimuli, section_stim, loc_stim):
        """ Doesn't return anything, places, sets the parameters and activates step current stimulus at the given dendritic location ."""
        raise NotImplementedError()

