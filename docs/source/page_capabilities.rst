############
Capabilities
############

The capabilities are the interface between the tests and the models. The ``ModelLoader()`` class (in ``utils.py``) inherits from the capabilities and must implement the methods of the capability. The test can only be run on a model if the necessary capability methods are implemented in the ``ModelLoader()``. All communication between the test and the model happens through the capabilities, therefore all functions that are needed to run simulations on the models and record their response are implemented as capabilities.

**SomaticFeaturesTest:** 
 - ``ReceivesSquareCurrent_ProvidesResponse()``  

**DepolarizationBlockTest:**
 - ``ReceivesSquareCurrent_ProvidesResponse()``

**BackPropagatingAPTest:** 
 - ``ReceivesSquareCurrent_ProvidesResponse()``  
 - ``ProvidesRecordingLocationsOnTrunk()``
 - ``ReceivesSquareCurrent_ProvidesResponse_MultipleLocations()``  

**PSPAttenuationTest:** 
 - ``ProvidesRandomDendriticLocation()``
 - ``ReceivesEPSCstim()``

**ObliqueIntegrationTest:** 
 - ``ProvidesGoodObliques()``
 - ``ReceivesMultipleSynapse()``

.. automodapi:: hippounit.capabilities
    :nosignatures:
    :no-main-docstr:
    :skip: {{ capability_classes|join(', ') }}