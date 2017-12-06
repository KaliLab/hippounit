from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless
import collections

class ZScore_backpropagatingAP(Score):
    """
    Sum of Z scores. A float indicating the sum of standardized difference
    from reference means for back-propagating AP amplitudes.
    """

    def __init__(self, score, related_data={}):

        self.score_l=[]
        for i in range(0, len(score)):
            if not isinstance(score[i], Exception) and not isinstance(score[i], float):
                raise InvalidScoreError("Score must be a float.")
            else:
                super(ZScore_backpropagatingAP,self).__init__(score[i], related_data=related_data)
                self.score_l.append(score[i])

    @classmethod
    def compute(cls, observation, prediction, distances):
        """Computes sum of z-scores from observation and prediction for back-propagating AP amplitudes"""

        errors = collections.OrderedDict()

        for i in range (0, len(distances)):
            if 'mean_AP1_amp_strong_propagating_at_'+str(distances[i])+'um' in observation.keys() or 'mean_AP1_amp_weak_propagating_at_'+str(distances[i])+'um' in observation.keys():
                p_value = prediction['model_AP1_amp_at_'+str(distances[i])+'um']['mean']
                o_mean = observation['mean_AP1_amp_strong_propagating_at_'+str(distances[i])+'um']
                o_std = observation['std_AP1_amp_strong_propagating_at_'+str(distances[i])+'um']

                try:
                    error = abs(p_value - o_mean)/o_std
                    error = assert_dimensionless(error)
                except (TypeError,AssertionError) as e:
                    error = e
                errors['AP1_amp_strong_propagating_at_'+str(distances[i])] = error


                o_mean = observation['mean_AP1_amp_weak_propagating_at_'+str(distances[i])+'um']
                o_std = observation['std_AP1_amp_weak_propagating_at_'+str(distances[i])+'um']

                try:
                    error = abs(p_value - o_mean)/o_std
                    error = assert_dimensionless(error)
                except (TypeError,AssertionError) as e:
                    error = e
                errors['AP1_amp_weak_propagating_at_'+str(distances[i])] = error

            else:
                p_value = prediction['model_AP1_amp_at_'+str(distances[i])+'um']['mean']
                o_mean = observation['mean_AP1_amp_at_'+str(distances[i])+'um']
                o_std = observation['std_AP1_amp_at_'+str(distances[i])+'um']

                try:
                    error = abs(p_value - o_mean)/o_std
                    error = assert_dimensionless(error)
                except (TypeError,AssertionError) as e:
                    error = e
                errors['AP1_amp_at_'+str(distances[i])] = error

        for i in range (0, len(distances)): # to keep better order: first all AP1, then all APlast
            p_value_l = prediction['model_APlast_amp_at_'+str(distances[i])+'um']['mean']
            o_mean_l = observation['mean_APlast_amp_at_'+str(distances[i])+'um']
            o_std_l = observation['std_APlast_amp_at_'+str(distances[i])+'um']

            try:
                error_l = abs(p_value_l - o_mean_l)/o_std_l
                error_l = assert_dimensionless(error_l)
            except (TypeError,AssertionError) as e:
                error_l = e
            errors['APlast_amp_at_'+str(distances[i])] = error_l

        score_sum_strong_propagating = 0.0
        score_sum_weak_propagating = 0.0

        for key, value in errors.iteritems():
            if 'strong' not in key:
                score_sum_weak_propagating += value
        for key, value in errors.iteritems():
            if 'weak' not in key:
                score_sum_strong_propagating += value
        return [score_sum_strong_propagating, score_sum_weak_propagating], errors

    def __str__(self):

		return 'Z_strong_propagating = %.2f, Z_weak_propagating = %.2f' % (self.score_l[0], self.score_l[1])
