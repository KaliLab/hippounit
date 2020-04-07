import sciunit
from sciunit import Capability


class ReceivesSquareCurrent_ProvidesResponse(sciunit.Capability):
    """Indicates that current can be injected into the model as
    a square pulse. """

    def inject_current(self, amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec):
        """
        This function must be implemented by the model.
        Must return numpy arrays containing the time and voltage values recorded at the location described by section_rec and dend_loc_rec

        Parameters
        ----------
        amp : float
            amplitude of the current injection (mV)
        delay : float
            delay before the current  injection (ms)
        duration : float
            duration of the current pulse
        section_stim : string
            the name of the stimulated section (eg. "soma")
        loc_stim : float
            location on the stimulated section (eg. 0.5)
        section_rec : string
            the name of the section whose response is recorded (eg. "soma")
        loc_rec : float
            location on the section from where the response is recorded (eg. 0.5)
        """

        raise NotImplementedError()

    def get_vm(self, amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec):

        """
        This function is called by the test and calls the inject_current() function.
        """

        t, v = self.inject_current(amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec)
        return t, v
