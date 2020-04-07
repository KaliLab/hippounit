######
Scores
######

Feature error scores are computed as the absolute difference between the feature value of the model (extracted from its response to the stimuli) and the experimental mean feature value, divided by the experimental standard deviation (*Z-score*).

.. automodapi:: hippounit.scores
    :nosignatures:
    :no-main-docstr:
    :skip: {{ score_classes|join(', ') }}