from sciunit import Score
import numpy
from sciunit.utils import assert_dimensionless


class ZScore_depolblock(Score):
    """
    Z scores. Float indicating standardized difference
    from reference means for current threshold (Ith), and equilibrium potential (Veq) during depolarization block.
    """

    def __init__(self, score, related_data={}):

        self.score_l=[]
        for i in range(0, len(score)):
	        if not isinstance(score[i], Exception) and not isinstance(score[i], float):
	            raise InvalidScoreError("Score must be a float.")
	        else:
	            super(ZScore_depolblock,self).__init__(score[i], related_data=related_data)
	            self.score_l.append(score[i])

    @classmethod
    def compute(cls, observation, prediction):
        """Computes a z-scores from an observation and a prediction."""

        p_value_Ith = prediction['model_Ith']
        o_mean_Ith = observation['mean_Ith']
        o_std_Ith = observation['Ith_std']
        p_value_Veq = prediction['model_Veq']
        o_mean_Veq = observation['mean_Veq']
        o_std_Veq = observation['Veq_std']

        try:
            result_Ith = (p_value_Ith - o_mean_Ith)/o_std_Ith
            result_Ith = assert_dimensionless(result_Ith)
            result_Veq = (p_value_Veq - o_mean_Veq)/o_std_Veq
            result_Veq = assert_dimensionless(result_Veq)

        except (TypeError,AssertionError) as e:
            result_Ith = e
            result_Veq = e

        return [result_Ith, result_Veq]

    def __str__(self):

		return 'Z_Ith = %.2f, Z_Veq = %.2f' % (self.score_l[0], self.score_l[1])
