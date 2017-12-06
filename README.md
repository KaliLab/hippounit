Introduction
============

This package contains validation tests for models of hippocampus, 
based on the SciUnit framework and the NeuronUnit package.

It was derived from the `hippounit` module of https://github.com/BlueBrain/neuronunit 
by separating out the validation tests from the models.

Current tests cover somatic behavior and signal integration in radial oblique dendrites of
hippocampal CA1 pyramidal cell models.

DepolarizationBlockTest
-----------------------

Tests if the model enters depolarization block under current injection of
increasing amplitudes to the soma.

Features tested:

* Ith: the threshold current to reach the depolarization block, 
  the amplitude of the current injection at which the cell exhibits
   the maximum number of spikes
* Veq: the average equilibrium value during the depolarization block, 
  average of the membrane potential over the last 100 ms of a 
  current pulse 50 pA above Ith.

(Bianchi et al. (2012) J Comput Neurosci, 33: 207-225)

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

Score used: p value from T-test. If the p value < 0.05, the model mean differs from the experimental mean

SomaticFeaturesTest
-------------------
 
Tests some somatic features using the BBP eFel
