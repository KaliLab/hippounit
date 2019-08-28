from __future__ import division
from builtins import str
from builtins import range
from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless
import collections

class ZScore_PSPAttenuation(Score):
    """
    Average of Z scores. A float indicating the average of standardized difference
    from reference means for PSP attenuation values.
    """

    def __init__(self, score, related_data={}):

        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_PSPAttenuation,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction, distances):
        """Computes average of z-scores from observation and prediction for PSP attenuation values"""

        errors = collections.OrderedDict()

        for i in range (0, len(distances)):
            p_value = prediction['mean_attenuation_soma/dend_'+str(distances[i])+'_um']['mean']
            o_mean = observation['mean_attenuation_soma/dend_'+str(distances[i])+'_um']
            o_std = observation['std_attenuation_soma/dend_'+str(distances[i])+'_um']

            try:
                error = abs(p_value - o_mean)/o_std
                error = assert_dimensionless(error)
            except (TypeError,AssertionError) as e:
                error = e

            errors['error_attenuation_soma/dend_'+str(distances[i])+'_um'] = error

        error_list = []
        for key, value in errors.items():
            error_list.append(value)

        score_avg = numpy.nanmean(error_list)

        return score_avg, errors

    def __str__(self):

        return 'ZScore_avg = %.2f' % self.score
