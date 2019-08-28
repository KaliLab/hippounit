from __future__ import division
from builtins import range
from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless
from scipy import stats


class P_Value_ObliqueIntegration(Score):
    """
    P valuee from t-test. Used in ObliqueIntegrationTest.
    """

    def __init__(self, score, related_data={}):

        self.score_l=[]

        for i in range(0, len(score)):
            if not isinstance(score[i], Exception) and not isinstance(score[i], float):
                raise InvalidScoreError("Score must be a float.")
            else:
                super(P_Value_ObliqueIntegration,self).__init__(score[i], related_data=related_data)
                self.score_l.append(score[i])

    @classmethod
    def ttest(cls, exp_mean, model_mean, exp_sd, model_sd, exp_n, model_n):
        m1 = exp_mean
        m2 = model_mean
        v1 = exp_sd**2
        v2 = model_sd**2
        n1 = exp_n
        n2 = model_n

        if n2 != 0 and v2 != 0:
            vn1 = v1 / n1
            vn2 = v2 / n2

            df = ((vn1 + vn2)**2) / ((vn1**2) / (n1 - 1) + (vn2**2) / (n2 - 1))

            denom = numpy.sqrt(vn1 + vn2)
            d = m1 - m2
            t = numpy.divide(d, denom)

            prob = stats.t.sf(numpy.abs(t), df) * 2  # use np.abs to get upper tail
        else:
            prob = float('NaN')

        return prob

    @classmethod
    def ttest_calc(cls, observation, prediction):

        exp_means=[observation['mean_threshold'], observation['mean_prox_threshold'], observation['mean_dist_threshold'], observation['mean_peak_deriv'], observation['mean_nonlin_at_th'], observation['mean_nonlin_suprath'],  observation['mean_amp_at_th'], observation['mean_time_to_peak'], observation['mean_async_nonlin']]
        exp_SDs=[observation['threshold_std'], observation['prox_threshold_std'], observation['dist_threshold_std'], observation['peak_deriv_std'], observation['nonlin_at_th_std'], observation['nonlin_suprath_std'],  observation['amp_at_th_std'], observation['time_to_peak_std'], observation['async_nonlin_std']]
        exp_Ns=[observation['exp_n'], observation['prox_n'], observation['dist_n'], observation['exp_n'], observation['exp_n'], observation['exp_n'], observation['exp_n'], observation['exp_n'], observation['async_n']]

        model_means = [prediction['model_mean_threshold'], prediction['model_mean_prox_threshold'], prediction['model_mean_dist_threshold'], prediction['model_mean_peak_deriv'], prediction['model_mean_nonlin_at_th'], prediction['model_mean_nonlin_suprath'],  prediction['model_mean_amp_at_th'], prediction['model_mean_time_to_peak'], prediction['model_mean_async_nonlin']]
        model_SDs = [prediction['model_threshold_std'], prediction['model_prox_threshold_std'], prediction['model_dist_threshold_std'], prediction['model_peak_deriv_std'], prediction['model_nonlin_at_th_std'], prediction['model_nonlin_suprath_std'], prediction['model_amp_at_th_std'], prediction['model_time_to_peak_std'], prediction['model_async_nonlin_std']]
        model_N= [prediction['model_n'], prediction['model_prox_n'], prediction['model_dist_n'], prediction['model_n'], prediction['model_n'], prediction['model_n'], prediction['model_n'], prediction['model_n'], prediction['model_n']]

        p_values=[]

        for i in range (0, len(exp_means)):

            try:
                ttest_result = cls.ttest(exp_means[i], model_means[i], exp_SDs[i], model_SDs[i], exp_Ns[i], model_N[i])
                ttest_result = assert_dimensionless(ttest_result)
                p_values.append(ttest_result)

            except (TypeError,AssertionError) as e:
                ttest_result = e

        return p_values

    def __str__(self):

        return '\n p_value_threshold = %.2f,\n p_value_prox_threshold  = %.2f,\n p_value_dist_threshold = %.2f,\n p_value_peak_dV/dt_at_threshold = %.2f,\n p_value_nonlin_at_th = %.2f,\n p_value_suprath_nonlin = %.2f,\n p_value_amplitude_at_th = %.2f,\n p_value_time_to_peak_at = %.2f,\n p_value_nonlin_at_th_asynch = %.2f\n' % (self.score_l[0], self.score_l[1],self.score_l[2], self.score_l[3],self.score_l[4], self.score_l[5],self.score_l[6], self.score_l[7],self.score_l[8])
