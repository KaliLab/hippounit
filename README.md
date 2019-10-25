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

Score type: Average of Z-Scores plus penalty if the two htreshold features (described above) differ. If the model did not enter depolarization block the Score is a 100 penalty.

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

Tests signal integration in oblique dendrites for increasing number of synchronous and asynchronous synaptic inputs.

A defalt synapse model is provided in this tets which consists of the Exp2Syn built in synapse of NEURON (as the AMPA component) and an additional NMDA receptor model with Jahr & Stevens voltage dependance and rise and decay time constants of 3.3 and 102.38 ms respectively. The time constant values used here are Q10 corrected values from McDermott et al. 2006. Q10 values for the rise and decay time constants are 2.2 (Hestrin et al., 1990) and 1.7 (Korinek et al., 2010) respectively. 

It is also possible to use the model's own AMPA and NMDA mechanisms.

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

Score used: Average of Z-scores. (Also available: p value from T-test. If the p value < 0.05, the model mean differs from the experimental mean)


Install HippoUnit
------------------

Install `git` and type:

    git clone https://github.com/KaliLab/hippounit.git

After cloning the repository you can install it by the standard installation method for Python packages:

    sudo python setup.py install

or as a user

    python setup.py install --user


Run HippoUnit
-------------------

See the example jupyter notebooks at https://github.com/KaliLab/HippoUnit_demo


Test Platforms
--------------

    1. Ubuntu 16.04 LTS
      - python 2.7.12
      - sciunit 0.2.0.2
      - efel 2.13.1
      - numpy 1.14.2
      - quantities 0.12.1
      - scipy 1.0.1
      - matplotlib 2.2.2
      - neuron 7.4

    2. Ubuntu 14.04.4 LTS 
      - python 2.7.6 
      - sciunit 0.2.0.2
      - efel 3.0.16
      - numpy 1.8.2'
      - quantities 0.12.1
      - scipy 0.18.1
      - matplotlib 2.0.2
      - neuron 7.4

    3. Ubuntu 16.04.6 LTS
      - python 2.7.12
      - sciunit 
      - efel 3.0.22
      - numpy 1.15.1
      - quantities 0.12.2
      - scipy 1.17.0
      - matplotlib 2.0.2, 2.2.4
      - neuron 7.4


    4. Ubuntu 16.04.6 LTS
      - python 3.5.2
      - sciunit 0.2.1.1
      - efel 3.0.58
      - numpy 1.16.4
      - quantities 0.12.1
      - scipy 1.3.0
      - matplotlib 3.0.3
      - neuron 7.6.2

Acknowledgments
-----------------

This open source software code was developed in part in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1 and SGA2).

SUPPORTED BY THE ÚNKP-19-3-III NEW NATIONAL EXCELLENCE PROGRAM OF THE MINISTRY FOR INNOVATION AND TECHNOLOGY.  

This research has been partially supported by the European Union, co-financed by the European Social Fund (EFOP-3.6.3-VEKOP- 16-2017-00002 ).
