from __future__ import division
from sciunit import Score
import numpy
import math
from sciunit.utils import assert_dimensionless
from quantities import nA

class ZScore_depolblock(Score):
    """
    Z scores. Float indicating standardized difference
    from reference means for current threshold (Ith), and equilibrium potential (Veq) during depolarization block.
    """

    def __init__(self, score, related_data={}):

        #self.score_l=[]
        if not isinstance(score, Exception) and not isinstance(score, float):
            raise InvalidScoreError("Score must be a float.")
        else:
            super(ZScore_depolblock,self).__init__(score, related_data=related_data)

    @classmethod
    def compute(cls, observation, prediction):
        """Computes a z-scores from an observation and a prediction."""

        p_value_I_maxNumAP = prediction['model_I_maxNumAP']
        p_value_I_below_depol_block = prediction['model_I_below_depol_block']
        o_mean_Ith = observation['mean_Ith']
        o_std_Ith = observation['Ith_std']
        p_value_Veq = prediction['model_Veq']
        o_mean_Veq = observation['mean_Veq']
        o_std_Veq = observation['Veq_std']

        try:
            result_I_maxNumAP = abs(p_value_I_maxNumAP - o_mean_Ith)/o_std_Ith
            result_I_maxNumAP = assert_dimensionless(result_I_maxNumAP)
            if not math.isnan(p_value_I_below_depol_block):
                result_I_below_depol_block = abs(p_value_I_below_depol_block - o_mean_Ith)/o_std_Ith
            else:
                result_I_below_depol_block = float('NaN')
            result_I_below_depol_block = assert_dimensionless(result_I_below_depol_block)
            if not math.isnan(p_value_Veq):
                result_Veq = abs(p_value_Veq - o_mean_Veq)/o_std_Veq
            else:
                result_Veq = float('NaN')
            result_Veq = assert_dimensionless(result_Veq)

        except (TypeError,AssertionError) as e:
            result_I_maxNumAP = e
            result_I_below_depol_block = e
            result_Veq = e

        if p_value_I_maxNumAP != p_value_I_below_depol_block and not math.isnan(result_I_below_depol_block): # according to the experiment thesetwo should be equal
            I_diff_penalty = 200.0*(abs(p_value_I_maxNumAP - p_value_I_below_depol_block)/(1*nA))   # divided be (1*nA) to make it dimensionless
            I_diff_penalty = assert_dimensionless(I_diff_penalty)
        else:
            I_diff_penalty = 0

        if math.isnan(result_I_below_depol_block) or math.isnan(result_Veq) :
            final_score = 100.0
        else:
            final_score = numpy.nanmean([result_I_maxNumAP, result_I_below_depol_block, result_Veq]) + I_diff_penalty

        return final_score, result_I_maxNumAP, result_I_below_depol_block, result_Veq, I_diff_penalty

    def __str__(self):

        #return 'Z_Ith = %.2f, Z_Veq = %.2f' % (self.score_l[0], self.score_l[1])
        return 'ZScore_avg+penalty = %.2f' % self.score
