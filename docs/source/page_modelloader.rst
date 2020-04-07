####################
ModelLoader class
####################

HippoUnit contains a general ``ModelLoader()`` class in ``utils.py`` which is capable of inheriting from the ``SciUnit.Model`` class. This class is implemented in a way that it is able to load and deal with most types of models defined in the HOC language of the NEURON simulator (either as standalone HOC models or as HOC templates). It implements all model-related (capability) methods that are needed to simulate these kinds of neural models in order to generate the prediction without any further coding required from the user. For neural models developed using other software and methods, the user needs to implement the capabilities through which the tests of HippoUnit perform the simulations and recordings on the model.

The path to the MOD (mechanism) files is argument to the ``ModelLoader()`` class, that hase to be given when instantiating the class.

General instance variables of the ``ModelLoader()`` class
############################################

These instance variables are used by each tests of HippoUnit and can/needed to be set after instantiating the *ModelLoader()* class.

* ``self.hocpath``:  String variable. The path to the model (HOC) file. This must be set by the user.

* ``self.cvode_active``: Boolean variable. If set to ``True`` simulation is run using variable time step integration. Default value is ``True``.

* ``self.template_name`` : String variable. Must be set if the model is implemented as a template, otherwise its value should be ``None``. Default value is ``None``.

* ``self.SomaSecList_name``: String variable. Must be set if the soma of the model is not a single section , otherwise its value should be ``None``. Default value is ``None``.

* ``self.soma``:  String variable. Name of the somatic section. Must be set by the user, except when the self.SomaSecList_name is set.

* ``self.v_init``: Voltage value from which the simulation is initiated. Default value is -70 mV.

* ``self.celsius``: Temperature of the simulation. Default value is 34 celsius degrees.

* ``self.name``: String variable. Name of the model. It is important to set, as outputs will be saved into a folder named like this.

* ``self.base_directory``: String variable.  Results will be saved in this directory. Default value is './validation_results/'.

* ``self.find_section_lists``: Boolean variable. If set to ``True`` , the different dendritic types of the apical tree (main apical trunk, apical tuft dendrites, radial oblique dendrites) are automatically classified. Can only be used for HOC models that load their morphologies from a separate morphology file (e.g., ASC, SWC). In this case the name of the section list should be set to None. Default value is ``False``


Instance variables of the ``Model Loader()`` class that are specific to a given test
###########################################################################################

* **Oblique Integration Test**

    * ``self.max_dist_from_soma``: The maximum distance from the soma on the apical trunk where oblique dendrites can be selected. Default value is 150 um.
    * ``self.c_minmax``: numpy array. Range in which the binary search algorithm searches for the proper synaptic weight.  Default value is ``numpy.array([0.00004, 0.04])``.
    * ``self.c_step_start``: Float variable.  Initial step size to create an array between the values of ``self.c_minmax``. Default value is 0.00004.
    * ``self.c_step_stop``: Float variable.  Edge condition  for step size to create an array between the values of *self.c_minmax*. Default value is 0.000004.
    * ``self.threshold``: Threshold for spike detection by eFEL.  Default value is -20.
    * ``self.ObliqueSecList_name``: String variable. The name of the section list containing the oblique dendrites implemented in the HOC file. Must be set by the user.
    * ``self.NMDA_name``:  String variable. Name of the NMDA receptor model included in the model. If its value is ``None`` the default NMDA model of the test is used.
    * ``self.AMPA_name`` :  String variable. Name of the AMPA receptor model included in the model. If its value is ``None`` the default *Exp2Syn* of NEURON is used.
    * ``self.AMPA_NMDA_ratio``: Float variable. The ration of the NMDA component of the synapse relative to AMPA.  Default value is AMPA/NMDA = 2.0.
    * ``self.AMPA_tau1``: Float variable. The rising time constant of the AMPA. Default value is 0.1 ms.
    * ``self.AMPA_tau2`` : Float variable. The decay time constant of the AMPA. Default value is 2.0 ms.
    * ``self.start``: time of the start of the synaptic stimulation. Default value is 150 ms.


* **PSP Attenuation Test**

    * ``self.TrunkSecList_name``: String variable. The name of the section lis containing the trunk dendrites implemented in the HOC file. Must be set by the user.
    * ``self.AMPA_tau1``: Float variable. The rising time constant of the AMPA. Default value is 0.1 ms.
    * ``self.AMPA_tau2``: Float variable. The decay time constant of the AMPA. Default value is 2.0 ms.

* **Backpropagating AP Test**

    * ``self.TrunkSecList_name``: String variable. The name of the section list containing the trunk dendrites implemented in the HOC file. Must be set by the user.

       