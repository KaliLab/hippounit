from __future__ import division
from sciunit import Score
import math
from sciunit.utils import assert_dimensionless

class ZScore_APfreq(Score):
    """
    Z scores. Float indicating standardized difference
    from reference means for current threshold (Ith), and equilibrium potential (Veq) during depolarization block.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dont_hide = ["related_data"]
    
    @classmethod
    def compute(cls, observation, prediction):
        """Computes average of scores from an observation and a prediction."""

        stim_spikecounts = []
        sum_scores = 0
        for obs, pred in zip(observation, prediction):
            # evaluate individual Z-scores
            # Note: if mean, std zero, then taking absolute difference (see if better option)
            score = None
            if obs["std"] == 0:
                score = float(abs(pred["freq"] - obs["mean"]))
            else:
                score = float((pred["freq"] - obs["mean"]) / obs["std"])
            sum_scores = sum_scores + abs(score)

            stim_spikecounts.append({
                "i_inj": obs["i_inj"],
                "mean": obs["mean"],
                "std": obs["std"],
                "freq": pred["freq"],
                "score": score
            })
        avg_score = sum_scores/len(observation)
        return ZScore_APfreq(avg_score), stim_spikecounts

    def __str__(self):
        return '%.2f' % self.score