from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless

class ZScore_somaticSpiking(Score):
    """
    Sum of Z scores. A float indicating the sum of standardized difference
    from reference means for somatic spiking features.
    """

    def __init__(self, score, related_data={}):

	    if not isinstance(score, Exception) and not isinstance(score, float):
	        raise InvalidScoreError("Score must be a float.")
	    else:
	        super(ZScore_somaticSpiking,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        """Computes sum of z-scores from observation and prediction for somatic spiking features"""

        feature_error_means=numpy.array([])
        feature_error_stds=numpy.array([])
        features_names=(observation.keys())
        feature_results_dict={}

        for i in range (0, len(features_names)):
            p_value = prediction[features_names[i]]['feature mean']
            o_mean = float(observation[features_names[i]]['Mean'])
            o_std = float(observation[features_names[i]]['Std'])

            p_std = prediction[features_names[i]]['feature sd']


            try:
                feature_error = abs(p_value - o_mean)/o_std
                feature_error = assert_dimensionless(feature_error)
                feature_error_mean=numpy.mean(feature_error)
                feature_error_sd=numpy.std(feature_error)

            except (TypeError,AssertionError) as e:
                feature_error = e
            feature_error_means=numpy.append(feature_error_means,feature_error_mean)
            feature_result={features_names[i]:{'mean feature error':feature_error_mean,
                                            'feature error sd':feature_error_sd}}

            feature_results_dict.update(feature_result)


        score_sum=numpy.sum(feature_error_means)

        return score_sum, feature_results_dict, features_names

    def __str__(self):

		return 'Z_sum = %.2f' % self.score
