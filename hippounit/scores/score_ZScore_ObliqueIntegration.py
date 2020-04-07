from __future__ import division
from sciunit import Score
import numpy
import collections
from sciunit.utils import assert_dimensionless

class ZScore_ObliqueIntegration(Score):
    """
    Average of Z scores. A float indicating the average of standardized difference
    from reference means for oblique integration features.
    """

    def __init__(self, score, related_data={}):

        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_ObliqueIntegration,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        """Computes average of z-scores from observation and prediction for features of dendritic integration in oblique dendrites"""
        #print observation
        #print prediction


        errors_dict=collections.OrderedDict()
        errors=[]

        for feat_name, value in observation.items():
            if 'mean' in feat_name:
                p_mean = prediction['model_' + feat_name]
                o_mean = observation[feat_name]
                o_std = observation[feat_name[5:] + '_std']

                try:
                    feature_error = abs(p_mean - o_mean)/o_std
                    feature_error = assert_dimensionless(feature_error)

                except (TypeError,AssertionError) as e:
                    feature_error = e
                errors.append(feature_error)
                errors_dict[feat_name[5:] + '_error'] = feature_error

        #print errors_dict

        score_avg=numpy.nanmean(errors)

        return score_avg, errors_dict

    def __str__(self):

        return 'ZScore_avg = %.2f' % self.score
