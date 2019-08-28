from __future__ import division
from builtins import str
from builtins import range

from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless
import collections

class ZScore_backpropagatingAP(Score):
    """
    Average of Z scores. A float indicating the average of standardized difference
    from reference means for back-propagating AP amplitudes.
    """

    def __init__(self, score, related_data={}):

        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_backpropagatingAP,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction, distances):
        """Computes average of z-scores from observation and prediction for back-propagating AP amplitudes"""

        errors = collections.OrderedDict()

        for i in range (0, len(distances)):
            if 'mean_AP1_amp_strong_propagating_at_'+str(distances[i])+'um' in list(observation.keys()) or 'mean_AP1_amp_weak_propagating_at_'+str(distances[i])+'um' in list(observation.keys()):
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

        score_strong_propagating = []
        score_weak_propagating = []

        for key, value in errors.items():
            if 'strong' not in key:             # everything except 'strong'
                score_weak_propagating.append(value)
        for key, value in errors.items():
            if 'weak' not in key:
                score_strong_propagating.append(value)


        score_avg_weak_propagating = numpy.nanmean(score_weak_propagating)
        score_avg_strong_propagating = numpy.nanmean(score_strong_propagating)


        if score_avg_weak_propagating < score_avg_strong_propagating:
            cls.strong = False
        elif score_avg_weak_propagating > score_avg_strong_propagating:
            cls.strong = True
        elif score_avg_weak_propagating == score_avg_strong_propagating:
            cls.strong = None

        return [score_avg_strong_propagating, score_avg_weak_propagating], errors

    def __str__(self):

        if ZScore_backpropagatingAP.strong:
            return 'Z_score_avg_STRONG_propagating = %.2f' % self.score
        elif ZScore_backpropagatingAP.strong is False:
            return 'Z_score_avg_WEAK_propagating = %.2f' % self.score
        elif ZScore_backpropagatingAP.strong is None:
            return 'Z_score_avg = %.2f' % self.score

        '''
        if self.score_l[0] < self.score_l[1]:
            self.score = self.score_l[0]
            return 'Z_score_avg_STRONG_propagating = %.2f' % self.score_l[0]
        elif self.score_l[1] < self.score_l[0]:
            self.score = self.score_l[1]
            return 'Z_score_avg_WEAK_propagating = %.2f' % self.score_l[1]
        elif self.score_l[1] == self.score_l[0]:
            self.score = self.score_l[0]
            return 'Z_score_avg = %.2f' % self.score_l[0]
        '''
