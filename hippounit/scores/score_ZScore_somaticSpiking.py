from __future__ import division
from builtins import range
from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless

class ZScore_somaticSpiking(Score):
    """
    Mean of Z scores. A float indicating the sum of standardized difference
    from reference means for somatic spiking features.
    """

    def __init__(self, score, related_data={}):

        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_somaticSpiking,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        """Computes average of z-scores from observation and prediction for somatic spiking features"""

        feature_errors=numpy.array([])
        features_names=(list(observation.keys()))
        feature_results_dict={}
        bad_features = []

        for i in range (0, len(features_names)):
            p_value = prediction[features_names[i]]['feature mean']
            o_mean = float(observation[features_names[i]]['Mean'])
            o_std = float(observation[features_names[i]]['Std'])

            p_std = prediction[features_names[i]]['feature sd']


            try:
                feature_error = abs(p_value - o_mean)/o_std
                feature_error = assert_dimensionless(feature_error)
            except ZeroDivisionError:
                feature_error = float("inf")
                feature_error = float("inf")
            except (TypeError,AssertionError) as e:
                feature_error = e
            #feature_errors=numpy.append(feature_errors,feature_error)
            feature_result={features_names[i]: feature_error}

            feature_results_dict.update(feature_result)

            if numpy.isnan(feature_error) or numpy.isinf(feature_error):
                bad_features.append(features_names[i])
            else:
                feature_errors=numpy.append(feature_errors,feature_error)
        score_avg=numpy.nanmean(feature_errors)

        return score_avg, feature_results_dict, features_names, bad_features

    def __str__(self):

        return 'ZScore_avg = %.2f' % self.score
