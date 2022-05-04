from __future__ import division
from builtins import str
from builtins import range

from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless
import collections

class ZScore_PathwayInteraction(Score):
    """
    Average of Z scores. A float indicating the average of standardized difference
    from reference means for Pathway Interaction features.
    """

    def __init__(self, score, related_data={}):

        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_PathwayInteraction,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        """Computes average of z-scores from observation and prediction for Pathway Interaction features"""

        errors_dict = collections.OrderedDict()
        errors=[]

        for pathway, value in prediction.items():
            if pathway not in list(errors_dict.keys()):
                errors_dict[pathway] = {}
            for feat_name, feat_value in prediction[pathway].items():
                p_mean = prediction[pathway][feat_name]['mean']
                o_mean = observation[pathway][feat_name]['mean']
                o_std = observation[pathway][feat_name]['std']

                try:
                    feature_error = abs(p_mean - o_mean)/o_std
                    feature_error = assert_dimensionless(feature_error)

                except (TypeError,AssertionError) as e:
                    feature_error = e
                errors.append(feature_error)
                errors_dict[pathway].update({feat_name : feature_error})
                
        penalty_PP_depol = 0    
        penalty_SC_PP = 0      
        if numpy.isnan(errors_dict['PP+depol']['plateau duration']):
            penalty_PP_depol += 100
        if numpy.isnan(errors_dict['SC+PP']['plateau duration']):
            penalty_SC_PP += 100

        score_avg=numpy.nanmean(errors) + penalty_PP_depol + penalty_SC_PP

        return score_avg, errors_dict, penalty_PP_depol, penalty_SC_PP



    def __str__(self):

        return 'ZScore_avg = %.2f' % self.score


