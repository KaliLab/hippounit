Introduction
============

This package contains validation tests for models of hippocampus, 
based on the SciUnit framework and the NeuronUnit package.

Current tests cover somatic behavior and signal propagation and integration in apical dendrites of
hippocampal CA1 pyramidal cell models.

Feature errors: Z-score (the difference of model feature value from the exerimental mean feature value in the unit of the experimental SD.)

SomaticFeaturesTest
-------------------
 
The Somatic Features Test can be used for both pyramidal cells and interneurons. It evaluates the model against various eFEL features under somatic current injection of varying amplitudes.

Score type: average of Z-Scores.


DepolarizationBlockTest
-----------------------

Aims to determine whether the model enters depolarization block to prolonged, high intensity somatic current stimulus.

Features tested:

* Ith: the threshold current to reach the depolarization block, the amplitude of the current injection at which the cell exhibits
  the maximum number of spikes.
  In the test two separate features are evaluated: the current intensity to which the model fires the maximum number of action potentials, and 	 the current intensity before the model enters depolarization block. If these two are not equal a penalty is added to the score.
* Veq: the average equilibrium value during the depolarization block, 
  average of the membrane potential over the last 100 ms of a 
  current pulse 50 pA above Ith.

(Bianchi et al. (2012) J Comput Neurosci, 33: 207-225)

Score type: Sum of Z-Scores. If the model did not enter depolarization block the Score is a 100 penalty.

BackpropagatingAPTest
----------------------

The test evaluates the mode and efficacy of back-propagating action potentials on the apical trunk in locations of different distances from the soma. The amplitude of the first and last AP of around 15 Hz train is compared to experimental data from Golding et al. 2001 (J Neurophysiol 86:2998-3010).

Score type: average of Z-scores.

PSPAttenuationTest 
----------------------

Evaluates how much the post synaptic potential (using EPSC stimulus) attenuates from the dendrite (different distances) to the soma. The soma/dendrite attenuation is compared to data from Magee & Cook 2000 (Nat Neurosci 3:895-903).

Score type: average of Z-scores.

ObliqueIntegrationTest
----------------------

Tests signal integration in oblique dendrites for increasing number of synchronous and asynchronous inputs.

Features tested:

* voltage threshold: the upper limit of the threshold for dendritic spike initiation (an expected somatic depolarization at wich step like dV/dt increase occur)
* proximal threshold
* distal threshold
* degree of nonlinearity at threshold
* suprathreshold degree of nonlinearity
* peak derivative of (somatic voltage) at threshold
* peak amplitude of somatic EPSP
* time to peak of somatic EPSP
* degree of nonlinearity in the case of asynchronous inputs

(A. Losonczy and J. C. Magee (2006) Neuron, 50: 291-307)

Score used: Sum of Z-scores. (Also available: p value from T-test. If the p value < 0.05, the model mean differs from the experimental mean)


