import sciunit
from sciunit import Capability


class ReceivesSquareCurrent_ProvidesResponse_MultipleLocations(sciunit.Capability):
    """Indicates that current can be injected into the model as
    a square pulse. And records at multiple locations."""

    def inject_current_record_respons_multiple_loc(self, amp, delay, dur, section_stim, loc_stim, dend_locations):
        """
        This function must be implemented by the model.

        Must return numpy arrays containing the time vector and the voltage values recorded on the stimulus location (soma), and a nested dictionary containing voltage vectors of the recorded dendritic locations at the examined distances in this test.

        eg.: {dist1: { ('trunk_segment1',location1): numpy.array(voltage trace),
                       ('trunk_segment2',location2): numpy.array(voltage trace)
                      },
                      { ('trunk_segment3',location3): numpy.array(voltage trace),
                       ('trunk_segment4',location4): numpy.array(voltage trace)
                      },

        Parameters
        ----------
        amp : float
            amplitude of the current injection (mV)
        delay : float
            delay before the current  injection (ms)
        duration : float
            duration of the current pulse
        section_stim : str
            the name of the stimulated section (eg. "soma")
        loc_stim : float
            location on the stimulated section (eg. 0.5)
        dend_locations : list
            list of recording locations in the form:  dend_loc = (dist1, ['trunk_segment1_1',location],['trunk_segment1_2',location]), (dist2, ['trunk_segment2',location]),(dist3, ['trunk_segment3',location]), (dist4, ['trunk_segment4',location])
        """

        raise NotImplementedError()

    def get_multiple_vm(self, amp, delay, dur, section_stim, loc_stim, dend_locations):
        """
        This function is called by the test and calls the inject_current_record_respons_multiple_loc() function.
        """
        # v : dictionary -  keys: dendritic location, values: the voltage trace for each recording locations
        t, v_stim, v = self.inject_current_record_respons_multiple_loc(amp, delay, dur, section_stim, loc_stim, dend_locations)
        return t, v_stim, v
