from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
#from builtins import str
from builtins import range
from quantities.quantity import Quantity
import sciunit
from sciunit import Test,Score
try:
    from sciunit import ObservationError
except:
    from sciunit.errors import ObservationError
import hippounit.capabilities as cap
from sciunit.utils import assert_dimensionless# Converters.
from sciunit.scores import BooleanScore,ZScore # Scores.

try:
    import numpy
except:
    print("NumPy not loaded.")

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#from neuron import h
import collections
import efel
import os
import multiprocessing
import multiprocessing.pool
import functools
import math
from scipy import stats

from scipy.optimize import fsolve, curve_fit
import scipy.interpolate as interpolate
from scipy.signal import find_peaks

import json
from hippounit import plottools
import collections


try:
    import pickle as pickle
except:
    import pickle
import gzip

try:
    import copy_reg
except:
    import copyreg

from types import MethodType

from quantities import mV, nA, ms, V, s

from hippounit import scores
import copy


def _pickle_method(method):
    func_name = method.__func__.__name__
    obj = method.__self__
    cls = method.__self__.__class__
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)


class NonDaemonPool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwds):
        proc = super(NonDaemonPool, self).Process(*args, **kwds)

        class NonDaemonProcess(proc.__class__):
            """Monkey-patch process to ensure it is never daemonized"""

            @property
            def daemon(self):
                return False

            @daemon.setter
            def daemon(self, val):
                pass

        proc.__class__ = NonDaemonProcess

        return proc

try:
    copy_reg.pickle(MethodType, _pickle_method, _unpickle_method)
except:
    copyreg.pickle(MethodType, _pickle_method, _unpickle_method)


class PathwayInteraction(Test):
    """ """

    def __init__(self, config = {},
                 observation = {},
                 name="Pathway Interaction test" ,
                 force_run=False,
                 force_run_adjust_syn_weight=False,
                 base_directory= None,
                 num_of_dend_locations = 15,
                 random_seed = 1,
                 show_plot=True,
                 save_all = True,
                 AMPA_weight_init = 0.000748,
                 trunk_origin = None):

        self.num_of_dend_locations = num_of_dend_locations
        self.random_seed = random_seed

        observation = self.format_data(observation)
        observation = self.add_std_to_observation(observation)

        Test.__init__(self, observation, name)

        self.required_capabilities = (cap.ProvidesRandomDendriticLocations, cap.ProvidesRecordingLocationsOnTrunk, cap.ReceivesSynapse, cap.InitialiseModel, cap.ThetaSynapticStimuli, cap.RunSimulation_ReturnTraces, cap.NumOfPossibleLocations, cap.ReceivesSynapseGivenPathway, cap.ReceivesMultipleSquareCurrents) # +=

        self.force_run_adjust_syn_weight = force_run_adjust_syn_weight
        self.force_run = force_run
        self.show_plot = show_plot
        self.save_all = save_all

        self.base_directory = base_directory

        self.path_figs = None    #added later, because model name is needed
        self.path_results = None
        self.trunk_origin = trunk_origin

        self.logFile = None
        self.test_log_filename = 'test_log.txt'
        self.message_to_logFile = ''

        self.config = config

        self.npool = multiprocessing.cpu_count() - 1

        self.AMPA_weight_init = AMPA_weight_init

        description = ""

    score_type = scores.ZScore_PathwayInteraction

    def format_data(self, observation):

        for key, val in list(observation.items()):
            for ke, va in list(val.items()):
                for k, v in list(va.items()):
                    try:
                        assert type(observation[key][ke][k]) is Quantity
                    except Exception as e:
                        try:
                            observation[key][ke][k] = float(val)
                        except Exception as e:
                            quantity_parts = v.split(" ")
                            number = float(quantity_parts[0])
                            units = " ".join(quantity_parts[1:])
                            observation[key][ke][k] = Quantity(number, units)
        return observation

    def add_std_to_observation (self, observation):

        for key, val in list(observation.items()):
            for ke, va in list(val.items()):
                observation[key][ke]['std'] = float(observation[key][ke]['sem'] * numpy.sqrt(observation[key][ke]['mean'])) * observation[key][ke]['mean'].units
        #print(observation)
        return observation

        
    def analyse_syn_traces(self, model, t, v, t_no_input, v_no_input):

        if not numpy.array_equal(t, t_no_input):    #if the  time vectors are not equal, the traces are resampled with fixed time step
            dt = 0.025
            time_vector = numpy.arange(t[0], t[-1], dt)  #from the first to the last element of the original time vector

            interp_trace = numpy.interp(time_vector, t, v)
            interp_trace_no_input = numpy.interp(time_vector, t, v_no_input)

            depol = interp_trace - interp_trace_no_input


            print("Voltage traces are resampled using linear interpolation")

        else:
            depol = v - v_no_input
            time_vector = t

        max_depol = max(depol)

        return max_depol  


    def synapse(self, model, t_no_input, v_no_input, weight, path_adjust_syn_weight, pathway, dend_loc0):
        file_name = path_adjust_syn_weight + 'Trace_' + str(dend_loc0[0]) + '(' + str(dend_loc0[1]) + ')_' + 'weight_' + str(weight) + '.p'

        if self.force_run_adjust_syn_weight or (os.path.isfile(file_name) is False):

            t, v, v_dend = model.run_synapse_pathway_get_vm(dend_loc0, weight, pathway)
            if self.save_all:
                pickle.dump([t, v, v_dend], gzip.GzipFile(file_name, "wb"))

        else:
            [t, v, v_dend] = pickle.load(gzip.GzipFile(file_name, "rb"))

        max_soma_depol = self.analyse_syn_traces(model, t, v, t_no_input, v_no_input)

        return max_soma_depol 

    def adjust_syn_weight(self, model, dend_loc, pathway):

        if self.base_directory:
            path_adjust_syn_weight = self.base_directory + 'temp_data/' + 'pathway_interaction/' + model.name + '/adjust_syn_weight/'
        else:
            path_adjust_syn_weight = model.base_directory + 'temp_data/' + 'pathway_interaction/adjust_syn_weight/'

        try:
            if not os.path.exists(path_adjust_syn_weight) and self.save_all:
                os.makedirs(path_adjust_syn_weight)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        if pathway == 'SC':
            file_name = path_adjust_syn_weight + 'SC_weight.p'
            desired_somatic_depol = 0.2
        if pathway == 'PP':
            file_name = path_adjust_syn_weight + 'PP_weight.p'
            desired_somatic_depol = 0.2


        # SC_desired_somatic_depol = 0.2
        # PP_desired_somatic_depol = 0.2

        if self.force_run_adjust_syn_weight or (os.path.isfile(file_name) is False):

            file_name_no_input = path_adjust_syn_weight + 'Traces_no_input.p'

            if os.path.isfile(file_name_no_input) is False:

                pool_syn_ = multiprocessing.Pool(1, maxtasksperchild = 1)    # I use multiprocessing to keep every NEURON related task in independent processes
                t_no_input, v_no_input, v_dend_no_input = pool_syn_.apply(model.run_synapse_pathway_get_vm, args = (dend_loc[0], 0.0, pathway))
                # plt.plot(t_no_input, v_no_input)
                # plt.show()
                pool_syn_.terminate()
                pool_syn_.join()
                del pool_syn_
                if self.save_all:
                    pickle.dump([t_no_input, v_no_input, v_dend_no_input], gzip.GzipFile(file_name_no_input, "wb"))

            else:
                [t_no_input, v_no_input, v_dend_no_input] = pickle.load(gzip.GzipFile(file_name_no_input, "rb"))


            synapse_ = functools.partial(self.synapse, model, t_no_input, v_no_input, self.AMPA_weight_init, path_adjust_syn_weight, pathway)

            pool_syn = multiprocessing.Pool(self.npool, maxtasksperchild = 1)    # I use multiprocessing to keep every NEURON related task in independent processes
            max_soma_depols = pool_syn.map(synapse_, dend_loc, chunksize=1)
            pool_syn.terminate()
            pool_syn.join()
            del pool_syn

            print("before:" , max_soma_depols)
            avg_max_soma_depols = numpy.mean(max_soma_depols)
            print('avg before', avg_max_soma_depols)

            scale_factor = desired_somatic_depol / avg_max_soma_depols
            print('scale_factor', scale_factor)

            synapse_ = functools.partial(self.synapse, model, t_no_input, v_no_input, self.AMPA_weight_init * scale_factor, path_adjust_syn_weight, pathway)

            pool_syn = multiprocessing.Pool(self.npool, maxtasksperchild = 1)    # I use multiprocessing to keep every NEURON related task in independent processes
            max_soma_depols = pool_syn.map(synapse_, dend_loc, chunksize=1)
            pool_syn.terminate()
            pool_syn.join()
            del pool_syn

            print("after:" , max_soma_depols)
            avg_max_soma_depols = numpy.mean(max_soma_depols)
            print('avg after', avg_max_soma_depols)

            AMPA_weight_final = self.AMPA_weight_init * scale_factor


            pickle.dump(AMPA_weight_final, gzip.GzipFile(file_name, "wb"))


        else:
            AMPA_weight_final = pickle.load(gzip.GzipFile(file_name, "rb"))


        return AMPA_weight_final

    def adjust_num_syn(self, model, SC_weight, PP_weight, recording_loc, stimuli_params, t_no_input_rec_dend, v_no_input_rec_dend, pathway):
        interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train = stimuli_params

        new_stimuli_params = [interval_bw_trains, interval_bw_stimuli_in_train, 1, num_stimuli_in_train] 

        dist_range = [0,9999999999]
        random_seed = self.random_seed


        if self.base_directory:
            path = self.base_directory + 'temp_data/' + 'pathway_interaction/' + model.name + '/'
        else:
            path = model.base_directory + 'temp_data/' + 'pathway_interaction/'

        file_name = path + pathway + '_dendritic_locations.p'

        if self.force_run_adjust_syn_weight or (os.path.isfile(file_name) is False):

            if pathway == 'SC':

                model.SecList = model.ObliqueSecList_name
                dend_loc, locations_distances = model.get_random_locations_multiproc(10, self.random_seed, dist_range, self.trunk_origin) # number of random locations , seed
                PP_dend_loc =[] 
                num_of_loc = model.get_num_of_possible_locations()
                
                exp_depol = 16.0
                exp_depol_sd = 1.6

                # traces = self.theta_pathway_stimulus(model, SC_weight, PP_weight, dend_loc, PP_dend_loc, recording_loc, new_stimuli_params, 600, pathway, save_traces=False))


            elif pathway == 'PP':
                model.SecList = model.TuftSecList_name
                dend_loc, locations_distances = model.get_random_locations_multiproc(10, self.random_seed, dist_range, self.trunk_origin) # number of random locations , seed
                
                SC_dend_loc =[] 
                num_of_loc = model.get_num_of_possible_locations()

                exp_depol = 10.2
                exp_depol_sd = 1.0

                # traces = self.theta_pathway_stimulus(model, SC_weight, PP_weight, SC_dend_loc, dend_loc, recording_loc, new_stimuli_params, 600, pathway, save_traces=False))


            # max_depol = self.analyse_syn_traces(model, traces[pathway]['t'], traces[pathway]['v_dend'], t_no_input_rec_dend, v_no_input_rec_dend)
            # plt.figure()
            # plt.plot(traces[pathway]['t'], traces[pathway]['v_dend'])
            # plt.show()
            # print(max_depol)

            found = False
            prev_max_depol = None

            print(pathway, 'num_of_loc', num_of_loc)
            
            #"""
            if pathway == 'PP' or pathway == 'SC':      #changing this I can play with which pathway to tune automatically and which not
                while not found and len(dend_loc) > 1 and len(dend_loc) <= num_of_loc:

                    random_seed += 1

                    if prev_max_depol:
                        prev_max_depol = max_depol    # if it already has a value (we are not in the first iteration), it gets the value of the previous iteration 

                    pool = multiprocessing.Pool(1, maxtasksperchild = 1)   # multiprocessing pool is used so that the model can be killed after the simulation, avoiding pickle errors
                
                    if pathway == 'SC':
                        traces = pool.apply(self.theta_pathway_stimulus, args = (model, SC_weight, PP_weight, dend_loc, PP_dend_loc, recording_loc, new_stimuli_params, 1600, 0, pathway, False))   # , save_traces=False because we don't want to save all the traces during adjustment

                    elif pathway == 'PP':
                        traces = pool.apply(self.theta_pathway_stimulus, args = (model, SC_weight, PP_weight, SC_dend_loc, dend_loc, recording_loc, new_stimuli_params, 1600, 0, pathway, False))

                    pool.terminate()
                    pool.join()
                    del pool


                    max_depol = self.analyse_syn_traces(model, traces[pathway]['t'], traces[pathway]['v_dend'], t_no_input_rec_dend, v_no_input_rec_dend)
                    print(pathway, ': ', max_depol)


                    if not prev_max_depol:         # if it has a value of None (we are in the first iteration), it gets the same value as the max_depol 
                        prev_max_depol = max_depol


                    if max_depol < exp_depol - exp_depol_sd and prev_max_depol < exp_depol - exp_depol_sd:
                        if pathway == 'SC':
                            model.SecList = model.ObliqueSecList_name
                        elif pathway == 'PP':
                            model.SecList = model.TuftSecList_name
                        
                        prev_dend_loc = list(dend_loc)
             
                        dend_loc_, locations_distances_ = model.get_random_locations_multiproc(1, random_seed, dist_range, self.trunk_origin) # select one more location

                        while dend_loc_[0] in dend_loc and len(dend_loc) <= num_of_loc: 
                            random_seed += 1
                            dend_loc_, locations_distances_ = model.get_random_locations_multiproc(1, random_seed, dist_range, self.trunk_origin) # select one more location
                        dend_loc.append(dend_loc_[0]) 
                        print(pathway, ': ', dend_loc)

                    elif max_depol < exp_depol - exp_depol_sd and prev_max_depol > exp_depol + exp_depol_sd:
                        #print(pathway, ' koztes1')
                        
                        #print('depols', max_depol, prev_max_depol)
                        #print('dend_locs')
                        #print(dend_loc)
                        #print(prev_dend_loc)
                        accepted_depol_diff= min(abs(max_depol-exp_depol), abs(prev_max_depol-exp_depol))
                        if accepted_depol_diff == abs(prev_max_depol-exp_depol):
                            dend_loc = prev_dend_loc
                            found = True
                        else:
                            # dend_loc remains
                            found = True 
                        #print('chosen dend_loc', dend_loc)

                    elif max_depol > exp_depol + exp_depol_sd and prev_max_depol > exp_depol + exp_depol_sd:
                        prev_dend_loc = list(dend_loc)

                        dend_loc.pop()  #removing last element
                        print(pathway, ': ', dend_loc)

                    elif max_depol > exp_depol + exp_depol_sd and prev_max_depol < exp_depol - exp_depol_sd:
                        #print(pathway, ' koztes2')

                        #print('depols', max_depol, prev_max_depol)
                        #print('dend_locs')
                        #print(dend_loc)
                        #print(prev_dend_loc)

                        accepted_depol_diff= min(abs(max_depol-exp_depol), abs(prev_max_depol-exp_depol))
                        if accepted_depol_diff == abs(prev_max_depol-exp_depol):
                            dend_loc = list(prev_dend_loc)
                            found = True
                        else:
                            # dend_loc remains
                            found = True 
                        #print('chosen dend_loc', dend_loc)

                    elif exp_depol - exp_depol_sd < max_depol < exp_depol + exp_depol_sd:

                        found = True
                        #print(pathway, ': ', dend_loc)
                        
                    


                if not found:
                    print("The number of activated synapses could not be adjusted properly on pathway:", pathway)
                    print("Maximum depolarization achieved:", max_depol, "mV")
                    print("Stimulated dendritic locations:", dend_loc)
                    self.message_to_logFile += "The number of activated synapses could not be adjusted properly on pathway: " + pathway + "\n" + "Maximum depolarization achieved: " +  str(max_depol) + " mV \n" + "Stimulated dendritic locations :" + str(dend_loc) + "\n"
                
            #"""

            pathway_dend_locs = {pathway: dend_loc}
            print("final dend_loc " + pathway + " : ", dend_loc)
            pickle.dump(pathway_dend_locs, gzip.GzipFile(file_name, "wb"))

        else:
            
            pathway_dend_locs = pickle.load(gzip.GzipFile(file_name, "rb"))
            print("final dend_loc ", pathway_dend_locs)

        return pathway_dend_locs    


    def generate_no_input_traces(self, model, recording_loc):


        if self.base_directory:
            path = self.base_directory + 'temp_data/' + 'pathway_interaction/' + model.name + '/'
        else:
            path = model.base_directory + 'temp_data/' + 'pathway_interaction/'

        file_name_no_input = path + 'traces_no_input.p'

        if self.force_run or (os.path.isfile(file_name_no_input) is False):
            model.initialise()
            t, v, v_dend,v_stim_locs = model.run_simulation([], recording_loc, 1600)

            if self.save_all:
                pickle.dump([t, v, v_dend], gzip.GzipFile(file_name_no_input, "wb"))

        else:
            [t, v, v_dend] = pickle.load(gzip.GzipFile(file_name_no_input, "rb"))

        return t, v, v_dend

    def spikecount(self, delay, duration, v_trace):

        efel.setThreshold(-25)

        trace = {}
        traces=[]
        trace['T'] = v_trace[0]
        trace['V'] = v_trace[1]
        trace['stim_start'] = [delay]
        trace['stim_end'] = [delay + duration]
        traces.append(trace)

        traces_results = efel.getFeatureValues(traces, ['Spikecount'])

        spikecount = traces_results[0]['Spikecount'][0]

        return spikecount

    def extract_efel_features(self, delay, duration, v_trace, features):

        efel.setThreshold(-25)

        trace = {}
        traces=[]
        trace['T'] = v_trace[0]
        trace['V'] = v_trace[1]
        trace['stim_start'] = [delay]
        trace['stim_end'] = [delay + duration]
        traces.append(trace)

        traces_results = efel.getFeatureValues(traces, features)

        return traces_results[0]

    def run_current_stim(self, model, path_adjust_current_amplitude, amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x):

        file_name = path_adjust_current_amplitude + 'dend_trace_' + str(amplitude) + '_nA.p'

        if os.path.isfile(file_name) is False:

            pool = multiprocessing.Pool(1, maxtasksperchild = 1)   # multiprocessing pool is used so that the model can be killed after the simulation, avoiding pickle errors
            t, v = pool.apply(model.get_vm, args = (amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x))
            pool.terminate()
            pool.join()
            del pool
        
            if self.save_all:
                pickle.dump([t, v], gzip.GzipFile(file_name, "wb"))

        else:
            [t, v] = pickle.load(gzip.GzipFile(file_name, "rb"))

        return t, v

    def binsearch(self, model, path_adjust_current_amplitude, stim_range, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x, desired_depol,v_no_input):
        c_minmax = stim_range
        c_step_start = 0.01
        c_step_stop= 0.002

        found = False
        spikecounts = []
        amplitudes = []
        max_depols = []

        tolerance = 2.0

        print("DOING BINARY SEARCH")

        while c_step_start >= c_step_stop and not found:

            c_stim = numpy.arange(c_minmax[0], c_minmax[1], c_step_start)
            print('c_stim: ', c_stim)

            first = 0
            last = numpy.size(c_stim, axis=0)-1

            while first <= last and not found:

                midpoint = (first + last)//2
                amplitude = c_stim[midpoint]
                #print('INFO: ', first, c_stim[first], last, c_stim[last])

                result=[]

                t, v = self.run_current_stim(model, path_adjust_current_amplitude, amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x)

                depol_dend = v - v_no_input

                max_depol = numpy.max(depol_dend)
                spike_count = self.spikecount(delay, duration, [t, v])
                print("amp: ", amplitude, "depol: ", max_depol, "spike count: ", spike_count)

                amplitudes.append(amplitude)
                spikecounts.append(spike_count)
                max_depols.append(max_depol)

                if spike_count == 0 and max_depol <= desired_depol + tolerance and max_depol >= desired_depol - tolerance:
                    found = True
                else:
                    if spike_count == 1 or (spike_count == 0 and max_depol > desired_depol + tolerance):
                        last = midpoint-1
                    elif spike_count == 0 and max_depol < desired_depol - tolerance:
                        first = midpoint+1
            c_step_start=c_step_start/2

        if not found:
            amp_index = min((p for p in range(len(spikecounts)) if spikecounts[p] == 0), key=lambda i: abs(max_depols[i]-desired_depol)) # we choose the one that is nearest to the desired depol, but no APs
            amplitude = amplitudes[amp_index]
            spike_count = spikecounts[amp_index]
            max_depol = max_depols[amp_index]


        binsearch_result=[found, amplitude, max_depol, spike_count]
        print("binsearch result: ", binsearch_result)

        return binsearch_result


    def adjust_current_amplitude(self, model, stimuli_list):

        amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x = stimuli_list

        if self.base_directory:
            path_adjust_current_amplitude = self.base_directory + 'temp_data/' + 'pathway_interaction/' + model.name + '/adjust_current_amplitude/'
        else:
            path_adjust_current_amplitude = model.base_directory + 'temp_data/' + 'pathway_interaction/adjust_current_amplitude/'


        try:
            if not os.path.exists(path_adjust_current_amplitude) and self.save_all:
                os.makedirs(path_adjust_current_amplitude)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        desired_depol = 15.8 # mV

        file_name_current_amp = path_adjust_current_amplitude + 'current_amp.p'

        if self.force_run_adjust_syn_weight or (os.path.isfile(file_name_current_amp) is False):

            t_no_input, v_no_input = self.run_current_stim(model, path_adjust_current_amplitude, 0.0, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x)

            spike_count_no_input = self.spikecount(delay, duration, [t_no_input, v_no_input])

            if spike_count_no_input > 0:
                print("Cell fires spontaneously")
                current_amp_final = floeat('nan')

            else:

                amplitude = 0.25

                t, v = self.run_current_stim(model, path_adjust_current_amplitude, amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x)

                depol_dend = v - v_no_input

                max_depol = numpy.max(depol_dend)
                spike_count = self.spikecount(delay, duration, [t, v])
                print("amp: ", amplitude, "depol: ", max_depol, "spike count: ", spike_count)


                while spike_count > 0:

                    amplitude = amplitude/2.0

                    t, v = self.run_current_stim(model, path_adjust_current_amplitude, amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x)

                    depol_dend = v - v_no_input

                    max_depol = numpy.max(depol_dend)
                    spike_count = self.spikecount(delay, duration, [t, v])
                    print("amp: ", amplitude, "depol: ", max_depol, "spike count: ", spike_count)


                scale_factor = desired_depol / max_depol
                print("scale_factor: ", scale_factor)


                amplitude = amplitude * scale_factor

                t, v = self.run_current_stim(model, path_adjust_current_amplitude, amplitude, delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x)

                depol_dend = v - v_no_input

                max_depol = numpy.max(depol_dend)
                spike_count = self.spikecount(delay, duration, [t, v])
                print("amp: ", amplitude, "depol: ", max_depol, "spike count: ", spike_count)

                
                if spike_count == 0:
                    current_amp_final = amplitude

                    pickle.dump(current_amp_final, gzip.GzipFile(file_name_current_amp, "wb"))
                else:
                    binsearc_result = self.binsearch(model, path_adjust_current_amplitude, [0, amplitude], delay, duration, stim_section_name, stim_location_x, rec_section_name, rec_location_x, desired_depol, v_no_input)

                    current_amp_final = binsearc_result[1]

                    pickle.dump(current_amp_final, gzip.GzipFile(file_name_current_amp, "wb"))


        else:
            current_amp_final = pickle.load(gzip.GzipFile(file_name_current_amp, "rb"))

        print("current_amp_final: ", current_amp_final)

        self.message_to_logFile += "current_amp_final: " + str(current_amp_final) + "\n"

        return current_amp_final
        
    def theta_pathway_stimulus(self, model, SC_weight, PP_weight, SC_dend_loc, PP_dend_loc, recording_loc, stimuli_params, tstop, depol_amp, pathway, save_traces):

        """Simulates pathway stimulation of the Schaffer-collateral or the Perforant Path, or both at the same time. The simultaneous activation of the 2 pathways is solved by calling the same Capability function but with different arguments (section list, synaptic parameters). For this to be feasible, the model must be loaded first, and therefore separate capability methods are needed to (1) load the model, (2) define the synaptic stimulus and (3) to run the simulation and make the recordings. In other tests all of these were done through a single capability method."""

        '''
        interval_bw_trains = 1/ self.config["frequency of stimulus sequence"] * 1000
        interval_bw_stimuli_in_train = 1/ self.config["frequency of trains"] * 1000
        num_trains = self.config["number of trains"]
        num_stimuli_in_train = self.config["number of stimuli in a train"]
        '''

        interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train = stimuli_params

        if self.base_directory:
            path = self.base_directory + 'temp_data/' + 'pathway_interaction/' + model.name + '/'
        else:
            path = model.base_directory + 'temp_data/' + 'pathway_interaction/'

        file_name = path + pathway + '_traces.p'

        if (self.force_run or (os.path.isfile(file_name) is False)) and save_traces:


            model.initialise()  # should be solved more like using Capabilities (problem: to add synapses, model should be loaded, but can not be reloaded. We want to be able to add PP and SC stimulation separately. (Later different synaptic parameters, different delay etc) 

            if pathway == 'SC':
                model.activate_theta_stimuli(SC_dend_loc, SC_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                t, v, v_dend, v_stim_locs = model.run_simulation(SC_dend_loc, recording_loc, tstop)
                '''
                plt.figure()
                plt.plot(t,v)
                plt.plot(t,v_dend)
                plt.title('SC stimulus')
                '''
            elif pathway == 'PP':
                model.activate_theta_stimuli(PP_dend_loc, PP_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                t, v, v_dend, v_stim_locs = model.run_simulation(PP_dend_loc, recording_loc, tstop)
                '''
                plt.figure()
                plt.plot(t,v)
                plt.plot(t,v_dend)
                plt.title('PP stimulus')
                '''
            elif pathway == 'SC+PP':
                model.activate_theta_stimuli(PP_dend_loc, PP_weight, 'PP', interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                # model.activate_theta_stimuli(SC_dend_loc + PP_dend_loc, PP_weight, 'PP')            
                model.activate_theta_stimuli(SC_dend_loc, SC_weight, 'SC', interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                t, v, v_dend, v_stim_locs = model.run_simulation(SC_dend_loc + PP_dend_loc, recording_loc, tstop)
                '''
                plt.figure()
                plt.plot(t,v)
                plt.plot(t,v_dend)
                plt.title('SC+PP stimulus')
                '''
            # plt.show()

            elif pathway == 'depol':
                (rec_ndend, xloc), distance = recording_loc
                model.activate_current_stimuli(depol_amp, model.start, num_stimuli_in_train * interval_bw_stimuli_in_train, num_trains, interval_bw_trains, rec_ndend, xloc)
                t, v, v_dend, v_stim_locs = model.run_simulation([], recording_loc, tstop)

            elif pathway == 'PP+depol':
                (rec_ndend, xloc), distance = recording_loc
                model.activate_theta_stimuli(PP_dend_loc, PP_weight, 'PP', interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                model.activate_current_stimuli(depol_amp, model.start, num_stimuli_in_train * interval_bw_stimuli_in_train, num_trains, interval_bw_trains, rec_ndend, xloc)
                t, v, v_dend, v_stim_locs = model.run_simulation(PP_dend_loc, recording_loc, tstop)

            elif pathway == 'SC+depol':
                (rec_ndend, xloc), distance = recording_loc
                model.activate_theta_stimuli(SC_dend_loc, SC_weight, 'SC', interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                model.activate_current_stimuli(depol_amp, model.start, num_stimuli_in_train * interval_bw_stimuli_in_train, num_trains, interval_bw_trains, rec_ndend, xloc)
                t, v, v_dend, v_stim_locs = model.run_simulation(SC_dend_loc, recording_loc, tstop)

            traces = {pathway: {'t' : t, 'v_soma' : v, 'v_dend' : v_dend, 'v_stim_locs' : v_stim_locs}} 

            if self.save_all:
                pickle.dump(traces, gzip.GzipFile(file_name, "wb"))


        elif save_traces is False:

            model.initialise()  # should be solved more like using Capabilities (problem: to add synapses, model should be loaded, but can not be reloaded. We want to be able to add PP and SC stimulation separately. (Later different synaptic parameters, different delay etc) 

            if pathway == 'SC':
                model.activate_theta_stimuli(SC_dend_loc, SC_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                t, v, v_dend, v_stim_locs = model.run_simulation(SC_dend_loc, recording_loc, tstop)

            elif pathway == 'PP':
                model.activate_theta_stimuli(PP_dend_loc, PP_weight, pathway, interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train)
                t, v, v_dend, v_stim_locs = model.run_simulation(PP_dend_loc, recording_loc, tstop)

            traces = {pathway: {'t' : t, 'v_soma' : v, 'v_dend' : v_dend, 'v_stim_locs' : v_stim_locs}} 


        elif (self.force_run is False and (os.path.isfile(file_name))) and save_traces:

            traces = pickle.load(gzip.GzipFile(file_name, "rb"))

        return traces

    def plot_traces(self, model, traces_dict):
    
    
        if self.base_directory:
            self.path_figs = self.base_directory + 'figs/' + 'pathway_interaction/' + model.name + '/'
        else:
            self.path_figs = model.base_directory + 'figs/' + 'pathway_interaction/'

        try:
            if not os.path.exists(self.path_figs) and self.save_all:
                os.makedirs(self.path_figs)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        print("The figures are saved in the directory: ", self.path_figs)
        
        fig= plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)

        ax1.plot(traces_dict['SC']['t'], traces_dict['SC']['v_soma'], label = 'soma')
        ax1.plot(traces_dict['SC']['t'], traces_dict['SC']['v_dend'], label = 'distal dendrite')
        ax1.set_xlabel('Time (ms)')
        ax1.set_ylabel('Voltage (mV)')
        ax1.title.set_text('SC stimulus')

        ax2.plot(traces_dict['PP']['t'], traces_dict['PP']['v_soma'])
        ax2.plot(traces_dict['PP']['t'], traces_dict['PP']['v_dend'])
        ax2.set_xlabel('Time (ms)')
        ax2.set_ylabel('Voltage (mV)')
        ax2.title.set_text('PP stimulus')

        ax3.plot(traces_dict['SC+PP']['t'], traces_dict['SC+PP']['v_soma'])
        ax3.plot(traces_dict['SC+PP']['t'], traces_dict['SC+PP']['v_dend'])
        ax3.set_xlabel('Time (ms)')
        ax3.set_ylabel('Voltage (mV)')
        ax3.title.set_text('SC+PP stimulus')
  
        fig.subplots_adjust(wspace = 0.5, hspace = 0.6)
        handles, labels = ax1.get_legend_handles_labels()
        lgd=fig.legend(handles, labels, bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_subplots', bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure()
        plt.plot(traces_dict['SC']['t'], traces_dict['SC']['v_soma'], label = 'soma')
        plt.plot(traces_dict['SC']['t'], traces_dict['SC']['v_dend'], label = 'distal dendrite')
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('SC stimulus')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_SC_stimulus', bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure()
        plt.plot(traces_dict['PP']['t'], traces_dict['PP']['v_soma'], label = 'soma')
        plt.plot(traces_dict['PP']['t'], traces_dict['PP']['v_dend'], label = 'distal dendrite')
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('PP stimulus')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_PP_stimulus', bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure()
        plt.plot(traces_dict['SC+PP']['t'], traces_dict['SC+PP']['v_soma'], label = 'soma')
        plt.plot(traces_dict['SC+PP']['t'], traces_dict['SC+PP']['v_dend'], label = 'distal dendrite')
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('SC+PP stimulus')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_SC_PP_stimulus', bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure()
        plt.plot(traces_dict['SC+depol']['t'], traces_dict['SC+depol']['v_soma'], label = 'soma')
        plt.plot(traces_dict['SC+depol']['t'], traces_dict['SC+depol']['v_dend'], label = 'distal dendrite')
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('SC stimulus + depolarization')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_SC_stimulus+depol', bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure()
        plt.plot(traces_dict['PP+depol']['t'], traces_dict['PP+depol']['v_soma'], label = 'soma')
        plt.plot(traces_dict['PP+depol']['t'], traces_dict['PP+depol']['v_dend'], label = 'distal dendrite')
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('PP stimulus + depolarization')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_PP_stimulus+depol', bbox_extra_artists=(lgd,), bbox_inches='tight')

        plt.figure()
        plt.plot(traces_dict['depol']['t'], traces_dict['depol']['v_soma'], label = 'soma')
        plt.plot(traces_dict['depol']['t'], traces_dict['depol']['v_dend'], label = 'distal dendrite')
        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')
        plt.title('depolarization')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        plt.savefig(self.path_figs + 'trace_only_depol', bbox_extra_artists=(lgd,), bbox_inches='tight')
        
        ncols = 3
        nrows = int(numpy.ceil(len(list(traces_dict['SC']['v_stim_locs'].keys()))/float(ncols)))
        fig2, axs2 = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
        fig2.subplots_adjust(wspace = 0.4, hspace = 0.5)
        axs2=axs2.flatten()
        for i, (key, value) in enumerate(traces_dict['SC']['v_stim_locs'].items()):
            axs2[i].plot(traces_dict['SC']['t'], value)
            axs2[i].set_title(str(key))
            axs2[i].set_xlabel('Time (ms)')
            axs2[i].set_ylabel('Voltage (mV)')
        fig2.suptitle('SC stimulus')
        plt.savefig(self.path_figs + 'local_traces_SC_stimulus', bbox_inches='tight')
	
        ncols = 3
        nrows = int(numpy.ceil(len(list(traces_dict['PP']['v_stim_locs'].keys()))/float(ncols)))
        fig3, axs3 = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
        fig3.subplots_adjust(wspace = 0.4, hspace = 0.5)
        axs3=axs3.flatten()
        for i, (key, value) in enumerate(traces_dict['PP']['v_stim_locs'].items()):
            axs3[i].plot(traces_dict['PP']['t'], value)
            axs3[i].set_title(str(key))
            axs3[i].set_xlabel('Time (ms)')
            axs3[i].set_ylabel('Voltage (mV)')
        fig3.suptitle('PP stimulus')
        plt.savefig(self.path_figs + 'local_traces_PP_stimulus', bbox_inches='tight')
	
        ncols = 3
        nrows = int(numpy.ceil(len(list(traces_dict['SC+PP']['v_stim_locs'].keys()))/float(ncols)))
        fig4, axs4 = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
        fig4.subplots_adjust(wspace = 0.4, hspace = 0.5)
        axs4=axs4.flatten()
        for i, (key, value) in enumerate(traces_dict['SC+PP']['v_stim_locs'].items()):
            axs4[i].plot(traces_dict['SC+PP']['t'], value)
            axs4[i].set_title(str(key))
            axs4[i].set_xlabel('Time (ms)')
            axs4[i].set_ylabel('Voltage (mV)')
        fig4.suptitle('SC+PP stimulus')
        plt.savefig(self.path_figs + 'local_traces_SC_PP_stimulus', bbox_inches='tight')

        ncols = 3
        nrows = int(numpy.ceil(len(list(traces_dict['SC+depol']['v_stim_locs'].keys()))/float(ncols)))
        fig5, axs5 = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
        fig5.subplots_adjust(wspace = 0.4, hspace = 0.5)
        axs5=axs5.flatten()
        for i, (key, value) in enumerate(traces_dict['SC+depol']['v_stim_locs'].items()):
            axs5[i].plot(traces_dict['SC+depol']['t'], value)
            axs5[i].set_title(str(key))
            axs5[i].set_xlabel('Time (ms)')
            axs5[i].set_ylabel('Voltage (mV)')
        fig5.suptitle('SC stimulus + depolarization')
        plt.savefig(self.path_figs + 'local_traces_SC_stimulus+depol', bbox_inches='tight')

        ncols = 3
        nrows = int(numpy.ceil(len(list(traces_dict['PP+depol']['v_stim_locs'].keys()))/float(ncols)))
        fig6, axs6 = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
        fig6.subplots_adjust(wspace = 0.4, hspace = 0.5)
        axs6=axs6.flatten()
        for i, (key, value) in enumerate(traces_dict['PP+depol']['v_stim_locs'].items()):
            axs6[i].plot(traces_dict['PP+depol']['t'], value)
            axs6[i].set_title(str(key))
            axs6[i].set_xlabel('Time (ms)')
            axs6[i].set_ylabel('Voltage (mV)')
        fig6.suptitle('PP stimulus + depolarization')
        plt.savefig(self.path_figs + 'local_traces_PP_stimulus+depol', bbox_inches='tight')


    def extract_plateau_features(self, model, traces, traces_no_input, stimuli_params, pathway): 

        interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train = stimuli_params


        time = traces['t']
        v_dend = traces['v_dend']
        v_soma = traces['v_soma']

        v_dend_no_input = traces_no_input[2]
        v_soma_no_input = traces_no_input[1]     #[t, v_soma, v_dend]

        depol_dend = v_dend-v_dend_no_input

        ''' Remove points with high derivative, and around the peak, interpolate remaining points to get the plateau without the bAPs'''

        dt=numpy.diff(time)
        dV=numpy.diff(depol_dend)

        deriv_dend = dV/dt

        indices_to_keep = numpy.where((abs(deriv_dend) < 1))[0] 

        peaks_ind, _ = find_peaks(depol_dend)

        peaks_v = depol_dend

        #print(peaks_ind)


        for i in peaks_ind:

            indices_to_keep = numpy.setdiff1d(indices_to_keep, numpy.where((time >= time[i] - 0.2) & (time <= time[i] + 0.2))[0])  # remove peaks and points 0.2 ms around the peaks


        depol_dend_plateau = depol_dend[indices_to_keep]

        time_plateau = time[indices_to_keep]

        interp_depol_dend_plateau = numpy.interp(time, time_plateau, depol_dend_plateau)


        dV_interp_depol_dend_plateau=numpy.diff(interp_depol_dend_plateau)
        deriv_interp_depol_dend_plateau = dV_interp_depol_dend_plateau/dt


        ''' start figure'''
        plt.figure()

        stim_start = model.start

        start_indices = []
        stop_indices = []
        plateau_amplitudes = []
        plateau_durations = []
        for i in range(num_trains):
            start = stim_start + i * interval_bw_trains
            start_index = numpy.where(time >= start)[0][0]
            stop = start + interval_bw_trains
            stop_index = numpy.where(time >= stop)[0][0]
            start_indices.append(start_index)
            stop_indices.append(stop_index)

            ''' extract plateau duration and plateau amplitude'''

            amplitude = numpy.max(interp_depol_dend_plateau[start_index:stop_index])
            plateau_amplitudes.append(amplitude)

            half_amplitude = amplitude/2.0 


            cross_idx = numpy.argwhere(numpy.diff(numpy.sign(interp_depol_dend_plateau[start_index:stop_index] - half_amplitude))).flatten()   # First the interp_depol_dend_plateau[start_index:stop_index] - half_amplitude and the corresponding signs  are calculated using numpy.sign;  numpy.diff gives the positions, where the sign changes (e.g. the lines cross); numpy.argwhere gives the indices.

            #print(cross_idx)

            plt.plot(time[start_index:stop_index][cross_idx], interp_depol_dend_plateau[start_index:stop_index][cross_idx], 'oc')


            durs = []
            if len(cross_idx) >1:
                for j in range(0, len(cross_idx), 2):
                    #print(j)
                    if j+1 < len(cross_idx):
                        if numpy.sign(deriv_interp_depol_dend_plateau[start_index:stop_index][cross_idx])[j] > 0 and numpy.sign(deriv_interp_depol_dend_plateau[start_index:stop_index][cross_idx])[j+1] < 0:
                            duration = time[start_index:stop_index][cross_idx][j+1] - time[start_index:stop_index][cross_idx][j]
                            durs.append(duration)
                #print(durs)
                if len(durs) > 0:
                    plt.plot([time[start_index:stop_index][cross_idx][numpy.argmax(durs)*2], time[start_index:stop_index][cross_idx][numpy.argmax(durs)*2 +1]], [interp_depol_dend_plateau[start_index:stop_index][cross_idx][numpy.argmax(durs)*2], interp_depol_dend_plateau[start_index:stop_index][cross_idx][numpy.argmax(durs)*2+1]], 'o-g')
                    plateau_durations.append(numpy.max(durs))
                else:
                    plateau_durations.append(float('nan'))
            else:
                plateau_durations.append(float('nan'))
                 

        #print('plateau amplitudes: ', plateau_amplitudes)
        #print('plateau durations: ', plateau_durations)

        ''' plot info'''
        plt.plot(time, depol_dend, color='orange')
        #plt.plot(time_plateau, depol_dend_plateau, 'm*')
        #plt.plot(time[:-1], deriv_dend, color='black')
        #plt.plot(time[peaks_ind], depol_dend[peaks_ind], 'go')
        plt.plot(time, interp_depol_dend_plateau, color='blue')
        plt.plot(time[start_indices], interp_depol_dend_plateau[start_indices], 'or')
        plt.plot(time[stop_indices], interp_depol_dend_plateau[stop_indices], 'ok')
        plt.title(pathway + ' - Distal dendrite')
        plt.savefig(self.path_figs + pathway + '_plateau_half_dur', bbox_inches='tight')

        for i, plateau_dur in enumerate(plateau_durations):         # if there is no interpretable plateau, we neither interpret its amplitude
            if numpy.isnan(plateau_dur):
                plateau_amplitudes[i] = float('nan')


        return plateau_amplitudes, plateau_durations


    def extract_features(self, model, traces, traces_no_input, stimuli_params, pathway):

        #print(pathway)

        interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train = stimuli_params

        time = traces['t']
        v_dend = traces['v_dend']
        v_soma = traces['v_soma']

        start_indices = []
        stop_indices = []
        plateau_amplitudes = []
        plateau_durations = []

        num_APs_list = []
        ISIs_list = []
        bAP_amp_list = []

        stim_start = model.start

        for i in range(num_trains):
            start = stim_start + i * interval_bw_trains
            start_index = numpy.where(time >= start)[0][0]
            stop = start + interval_bw_trains
            stop_index = numpy.where(time >= stop)[0][0]
            start_indices.append(start_index)
            stop_indices.append(stop_index)

            ''' extract number of APs, bAP amp, soma ISI'''

            efel_results = self.extract_efel_features(0, time[stop_index] - time[start_index], [time[start_index:stop_index], v_dend[start_index:stop_index]], ['Spikecount', 'AP_amplitude', 'AP_begin_voltage', 'peak_voltage'])
            #num_APs = spikecount(0, time[stop_index] - time[start_index], [time[start_index:stop_index], v_dend[start_index:stop_index]])
            #print(efel_results)
            if pathway == 'SC' or pathway == 'PP' or pathway == 'SC+depol' or pathway == 'PP+depol':
                num_APs = efel_results['Spikecount'][0]
                num_APs_list.append(num_APs)
                #print('num APs: ', num_APs)
            
            if  pathway == 'SC+depol' or pathway == 'PP+depol':
                try:    # sometimes eFEL finds less AP_begin_voltage than peak_voltage
                    bAP_amp = efel_results['peak_voltage'] - efel_results['AP_begin_voltage']   #efel_results['AP_amplitude'] somehow this often gives empty array
                except:
                    bAP_amp = float('nan')  
                bAP_amp_list.append(numpy.mean(bAP_amp))
                #print('bAP amp', bAP_amp)
            

            if pathway == 'SC+PP':
                efel_results_soma = self.extract_efel_features(0, time[stop_index] - time[start_index], [time[start_index:stop_index], v_soma[start_index:stop_index]], ['ISI_values'])
                ISI = efel_results_soma['ISI_values']
                if ISI is None:
                    ISI = []
                #print('ISI: ', ISI)
                ISIs_list.append(numpy.mean(ISI))
                #print('ISIs_list: ', ISIs_list)


        if pathway == 'SC' or pathway == 'PP':
            features = {pathway: {'num AP' : {'mean' : numpy.mean(num_APs_list), 'std' : numpy.std(num_APs_list)}}}

        if pathway == 'SC+depol':
            #later this should get a bAP_amp feature
            plateau_amplitudes, plateau_durations = self.extract_plateau_features(model, traces, traces_no_input, stimuli_params, pathway)
            features = {pathway: {'num AP' : {'mean' : numpy.mean(num_APs_list), 'std' : numpy.std(num_APs_list)}, 'bAP amp' : {'mean' : numpy.mean(bAP_amp_list) * mV, 'std' : numpy.std(bAP_amp_list) * mV}, '1st plateau duration' : {'mean' : numpy.nanmean(plateau_durations[0]) * ms, 'std' : numpy.nanstd(plateau_durations[0]) * ms}, '3-5th plateau duration' : {'mean' : numpy.nanmean(plateau_durations[2:]) * ms, 'std' : numpy.nanstd(plateau_durations[2:]) * ms}}}

        if pathway == 'PP+depol':
            #later this should get a bAP_amp feature
            plateau_amplitudes, plateau_durations = self.extract_plateau_features(model, traces, traces_no_input, stimuli_params, pathway)
            features = {pathway: {'num AP' : {'mean' : numpy.mean(num_APs_list), 'std' : numpy.std(num_APs_list)}, 'bAP amp' : {'mean' : numpy.mean(bAP_amp_list) * mV, 'std' : numpy.std(bAP_amp_list) * mV}, 'plateau duration' : {'mean' : numpy.nanmean(plateau_durations) * ms, 'std' : numpy.nanstd(plateau_durations) * ms}, '1st plateau duration' : {'mean' : numpy.nanmean(plateau_durations[0]) * ms, 'std' : numpy.nanstd(plateau_durations[0]) * ms}, '3-5th plateau duration' : {'mean' : numpy.nanmean(plateau_durations[2:]) * ms, 'std' : numpy.nanstd(plateau_durations[2:]) * ms}}}

        if pathway == 'SC+PP':
            plateau_amplitudes, plateau_durations = self.extract_plateau_features(model, traces, traces_no_input, stimuli_params, pathway)
            features = {pathway: {'plateau amplitude' : {'mean' : numpy.nanmean(plateau_amplitudes) * mV, 'std' : numpy.nanstd(plateau_amplitudes) * mV}, 'plateau duration' : {'mean' : numpy.nanmean(plateau_durations) * ms, 'std' : numpy.nanstd(plateau_durations) * ms}, '1st plateau duration' : {'mean' : numpy.nanmean(plateau_durations[0]) * ms, 'std' : numpy.nanstd(plateau_durations[0]) * ms}, '3-5th plateau duration' : {'mean' : numpy.nanmean(plateau_durations[2:]) * ms, 'std' : numpy.nanstd(plateau_durations[2:]) * ms}, 'somatic AP ISI' : {'mean' : numpy.mean(ISIs_list) * ms, 'std' : numpy.std(ISIs_list)* ms}}}
            #print('plateau features: ', plateau_amplitudes, plateau_durations)

        #print(features)

        return features


    def plot_features(self, prediction):
        
        observation = self.observation
        '''
        feat_means = []
        feat_stds =  []
        labels = []
        plt.figure()
        for pathway, feats in prediction.items():
            for feat, v in feats.items():
                feat_means.append(v['mean'])
                feat_stds.append(v['std'])
                labels.append(pathway + ' - ' + feat)

        y = range(len(feat_means))

        plt.errorbar(feat_means, y, xerr=feat_stds, linestyle='none', marker='o', color='blue')
        plt.yticks(y, labels)
        plt.savefig(self.path_figs + 'feature_values', bbox_inches='tight')
        '''

        model_num_AP_means = []
        model_plateau_duration_means = []
        model_plateau_amplitude_means = []
        model_ISI_means = []
        model_bAP_amp_means = []
        model_1st_plateau_duration_means = []
        model_3_5th_plateau_duration_means = []
        
        model_num_AP_stds = []
        model_plateau_duration_stds = []
        model_plateau_amplitude_stds = []
        model_ISI_stds = []
        model_bAP_amp_stds = []
        model_1st_plateau_duration_stds = []
        model_3_5th_plateau_duration_stds = []

        exp_num_AP_means = []
        exp_plateau_duration_means = []
        exp_plateau_amplitude_means = []
        exp_ISI_means = []
        exp_bAP_amp_means = []
        exp_1st_plateau_duration_means = []
        exp_3_5th_plateau_duration_means = []        

        exp_num_AP_stds = []
        exp_plateau_duration_stds = []
        exp_plateau_amplitude_stds = []
        exp_ISI_stds = []
        exp_bAP_amp_stds = []
        exp_1st_plateau_duration_stds = []
        exp_3_5th_plateau_duration_stds = []

        labels_num_AP = []
        labels_plateau_duration = []
        labels_plateau_amplitude = []
        labels_ISI = []
        labels_bAP_amp = []
        labels_1st_plateau_duration = []
        labels_3_5th_plateau_duration = []
        

        for pathway, feats in prediction.items():
            if 'num AP' in list(feats.keys()):
                model_num_AP_means.append(prediction[pathway]['num AP']['mean'])
                model_num_AP_stds.append(prediction[pathway]['num AP']['std'])

                exp_num_AP_means.append(observation[pathway]['num AP']['mean'])
                exp_num_AP_stds.append(observation[pathway]['num AP']['std'])

                labels_num_AP.append(pathway + ' - ' + 'num AP')
                
            if 'bAP amp' in list(feats.keys()):
                model_bAP_amp_means.append(prediction[pathway]['bAP amp']['mean'])
                model_bAP_amp_stds.append(prediction[pathway]['bAP amp']['std'])

                exp_bAP_amp_means.append(observation[pathway]['bAP amp']['mean'])
                exp_bAP_amp_stds.append(observation[pathway]['bAP amp']['std'])

                labels_bAP_amp.append(pathway + ' - ' + 'bAP amp')                

            if 'plateau duration' in list(feats.keys()):
                model_plateau_duration_means.append(prediction[pathway]['plateau duration']['mean'])
                model_plateau_duration_stds.append(prediction[pathway]['plateau duration']['std'])

                exp_plateau_duration_means.append(observation[pathway]['plateau duration']['mean'])
                exp_plateau_duration_stds.append(observation[pathway]['plateau duration']['std'])

                labels_plateau_duration.append(pathway + ' - ' + 'plateau duration')
                
            if '1st plateau duration' in list(feats.keys()):
                model_1st_plateau_duration_means.append(prediction[pathway]['1st plateau duration']['mean'])
                model_1st_plateau_duration_stds.append(prediction[pathway]['1st plateau duration']['std'])

                exp_1st_plateau_duration_means.append(observation[pathway]['1st plateau duration']['mean'])
                exp_1st_plateau_duration_stds.append(observation[pathway]['1st plateau duration']['std'])

                labels_1st_plateau_duration.append(pathway + ' - ' + '1st plateau duration')
                
            if '3-5th plateau duration' in list(feats.keys()):
                model_3_5th_plateau_duration_means.append(prediction[pathway]['3-5th plateau duration']['mean'])
                model_3_5th_plateau_duration_stds.append(prediction[pathway]['3-5th plateau duration']['std'])

                exp_3_5th_plateau_duration_means.append(observation[pathway]['3-5th plateau duration']['mean'])
                exp_3_5th_plateau_duration_stds.append(observation[pathway]['3-5th plateau duration']['std'])

                labels_3_5th_plateau_duration.append(pathway + ' - ' + '3-5th plateau duration')

            if 'plateau amplitude' in list(feats.keys()):
                model_plateau_amplitude_means.append(prediction[pathway]['plateau amplitude']['mean'])
                model_plateau_amplitude_stds.append(prediction[pathway]['plateau amplitude']['std'])

                exp_plateau_amplitude_means.append(observation[pathway]['plateau amplitude']['mean'])
                exp_plateau_amplitude_stds.append(observation[pathway]['plateau amplitude']['std'])

                labels_plateau_amplitude.append(pathway + ' - ' + 'plateau amplitude')

            if 'somatic AP ISI' in list(feats.keys()):
                model_ISI_means.append(prediction[pathway]['somatic AP ISI']['mean'])
                model_ISI_stds.append(prediction[pathway]['somatic AP ISI']['std'])

                exp_ISI_means.append(observation[pathway]['somatic AP ISI']['mean'])
                exp_ISI_stds.append(observation[pathway]['somatic AP ISI']['std'])

                labels_ISI.append(pathway + ' - ' + 'somatic AP ISI')
                
        model_all_plateau_duration_means = model_plateau_duration_means + model_1st_plateau_duration_means + model_3_5th_plateau_duration_means
        model_all_plateau_duration_stds = model_plateau_duration_stds + model_1st_plateau_duration_stds + model_3_5th_plateau_duration_stds
        exp_all_plateau_duration_means = exp_plateau_duration_means + exp_1st_plateau_duration_means + exp_3_5th_plateau_duration_means
        exp_all_plateau_duration_stds = exp_plateau_duration_stds + exp_1st_plateau_duration_stds + exp_3_5th_plateau_duration_stds
        labels_all_plateau_duration = labels_plateau_duration + labels_1st_plateau_duration + labels_3_5th_plateau_duration
        
        fig, axs = plt.subplots(3,2, figsize=(2*4, 2*4))
        plt.subplots_adjust(wspace = 0.5, hspace = 0.8)

        axs[0,0].errorbar(range(len(labels_num_AP)), model_num_AP_means, yerr=model_num_AP_stds, linestyle='none', marker='o', color='blue')
        axs[0,0].errorbar(range(len(labels_num_AP)), exp_num_AP_means, yerr=exp_num_AP_stds, linestyle='none', marker='o', color='red')
        axs[0,0].set_xticks(range(len(labels_num_AP)))
        axs[0,0].set_xticklabels(labels_num_AP, rotation = 20)
        axs[0,0].set_ylabel('# APs')

        axs[0,1].errorbar(range(len(labels_ISI)), model_ISI_means, yerr=model_ISI_stds, linestyle='none', marker='o', color='blue', label = 'model')
        axs[0,1].errorbar(range(len(labels_ISI)), exp_ISI_means, yerr=exp_ISI_stds, linestyle='none', marker='o', color='red', label = 'experiment')
        axs[0,1].set_xticks(range(len(labels_ISI)))
        axs[0,1].set_xticklabels(labels_ISI, rotation = 20)
        axs[0,1].set_ylabel('Somatic AP ISI (ms)')
        
        axs[1,0].errorbar(range(len(labels_plateau_amplitude)), model_plateau_amplitude_means, yerr=model_plateau_amplitude_stds, linestyle='none', marker='o', color='blue')
        axs[1,0].errorbar(range(len(labels_plateau_amplitude)), exp_plateau_amplitude_means, yerr=exp_plateau_amplitude_stds, linestyle='none', marker='o', color='red')
        axs[1,0].set_xticks(range(len(labels_plateau_amplitude)))
        axs[1,0].set_xticklabels(labels_plateau_amplitude, rotation = 20)
        axs[1,0].set_ylabel('Plateau amplitude (mV)')
        
        '''
        axs[1,1].errorbar(range(len(labels_plateau_duration)), model_plateau_duration_means, yerr=model_plateau_duration_stds, linestyle='none', marker='o', color='blue')
        axs[1,1].errorbar(range(len(labels_plateau_duration)), exp_plateau_duration_means, yerr=exp_plateau_duration_stds, linestyle='none', marker='o', color='red')
        axs[1,1].set_xticks(range(len(labels_plateau_duration)))
        axs[1,1].set_xticklabels(labels_plateau_duration, rotation = 20)
        axs[1,1].set_ylabel('Plateau duration (ms)')
        '''

        axs[1,1].errorbar(range(len(labels_all_plateau_duration)), model_all_plateau_duration_means, yerr=model_all_plateau_duration_stds, linestyle='none', marker='o', color='blue')
        axs[1,1].errorbar(range(len(labels_all_plateau_duration)), exp_all_plateau_duration_means, yerr=exp_all_plateau_duration_stds, linestyle='none', marker='o', color='red')
        axs[1,1].set_xticks(range(len(labels_all_plateau_duration)))
        axs[1,1].set_xticklabels(labels_all_plateau_duration, rotation = 90)
        axs[1,1].set_ylabel('Plateau duration (ms)')
        
        axs[2,0].errorbar(range(len(labels_bAP_amp)), model_bAP_amp_means, yerr=model_bAP_amp_stds, linestyle='none', marker='o', color='blue')
        axs[2,0].errorbar(range(len(labels_bAP_amp)), exp_bAP_amp_means, yerr=exp_bAP_amp_stds, linestyle='none', marker='o', color='red')
        axs[2,0].set_xticks(range(len(labels_bAP_amp)))
        axs[2,0].set_xticklabels(labels_bAP_amp, rotation = 20)
        axs[2,0].set_ylabel('bAP amplitude (mV)')
        
        axs[2,1].set_axis_off()
        

        lgd=axs[0,1].legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')

        fig.suptitle('Feature values')

        if self.save_all:
            plt.savefig(self.path_figs + 'feature_values', bbox_extra_artists=(lgd,), bbox_inches='tight')


    def plot_errors(self, errors):

        feat_errors = []
        labels = []
        plt.figure()
        for pathway, feats in errors.items():
            for feat, v in feats.items():
                feat_errors.append(v)
                labels.append(pathway + ' - ' + feat)

        y = range(len(feat_errors))

        plt.plot(feat_errors, y, linestyle='none', marker='o', color='blue')
        plt.yticks(y, labels)
        plt.title('Feature errors')
        plt.xlabel('# SDs')
        plt.savefig(self.path_figs + 'feature_errors', bbox_inches='tight')


    def validate_observation(self, observation):
        pass


    def generate_prediction(self, model, verbose=False):
        """Implementation of sciunit.Test.generate_prediction."""
        
        model.start = 400

        efel.reset()
        plt.close('all') #needed to avoid overlapping of saved images when the test is run on multiple models in a for loop

        if self.base_directory:
            self.path_results = self.base_directory + 'results/' + 'pathway_interaction/' + model.name + '/'
        else:
            self.path_results = model.base_directory + 'results/' + 'pathway_interaction/'

        try:
            if not os.path.exists(self.path_results):
                os.makedirs(self.path_results)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        dist_range = [0,9999999999]

        model.SecList = model.ObliqueSecList_name
        SC_dend_loc, SC_locations_distances = model.get_random_locations_multiproc(self.num_of_dend_locations, self.random_seed, dist_range, self.trunk_origin) # number of random locations , seed

        model.SecList = model.TuftSecList_name
        PP_dend_loc, PP_locations_distances = model.get_random_locations_multiproc(self.num_of_dend_locations, self.random_seed, dist_range, self.trunk_origin) # number of random locations , seed

        """Finding recording location on Trunk whose distance is closest to 300 um"""
        distances = [self.config["distance of recording location"]]
        tolerance = self.config["distance tolerance"]

        rec_locs, rec_locs_actual_distances = model.find_trunk_locations_multiproc(distances, tolerance, self.trunk_origin)
        print("recording locs", rec_locs, rec_locs_actual_distances)

        # recording_loc = min(rec_locs_actual_distances, key=abs(distances[0] - rec_locs_actual_distances.get))
        recording_loc = min(rec_locs_actual_distances.items(), key=lambda kv : abs(kv[1] - distances[0]))
        #print(recording_loc, type(recording_loc))


        if not model.AMPA_name:
            print('')
            print('The built in Exp2Syn is used as the AMPA component. Tau1 =', model.AMPA_tau1, ',Tau2 =', model.AMPA_tau2 , '.')
            print('')
        if not model.NMDA_name: 
            print('')
            print('The default NMDA model of HippoUnit is used with Jahr, Stevens voltage dependence.')
            print('')

        print("Adjusting synaptic weights ...")

        SC_weight = self.adjust_syn_weight(model, SC_dend_loc, pathway = 'SC') #0.000748
        print('SC AMPA weight', SC_weight)
        self.message_to_logFile += "SC AMPA weight: " + str(SC_weight) + "\n"

        PP_weight = self.adjust_syn_weight(model, PP_dend_loc, pathway = 'PP') #0.000748
        print('PP AMPA weight', PP_weight)
        self.message_to_logFile += "PP AMPA weight: " + str(PP_weight) + "\n"
        
        pool = multiprocessing.Pool(1, maxtasksperchild = 1)
        t_no_input_rec_dend, v_soma_no_input, v_no_input_rec_dend  = pool.apply(self.generate_no_input_traces, (model, recording_loc,)) # this is run in multiprocessing pool so that the model can be completely killed after done 
        pool.terminate()
        pool.join()
        del pool

        interval_bw_trains = 1/ self.config["frequency of stimulus sequence"] * 1000
        interval_bw_stimuli_in_train = 1/ self.config["frequency of trains"] * 1000
        num_trains = self.config["number of trains"]
        num_stimuli_in_train = self.config["number of stimuli in a train"]


        stimuli_params =[interval_bw_trains, interval_bw_stimuli_in_train, num_trains, num_stimuli_in_train] 

        # self.adjust_num_syn(model, SC_weight, PP_weight, recording_loc, stimuli_params, t_no_input_rec_dend, v_no_input_rec_dend, 'SC')
       
        pool = NonDaemonPool(self.npool, maxtasksperchild=1)  # NoDeamonPool is needed because Random locations are needed to be chosen several times, for which the model is loaded in a multiprocessing pool 
        adjust_num_syn_= functools.partial(self.adjust_num_syn, model, SC_weight, PP_weight, recording_loc, stimuli_params, t_no_input_rec_dend, v_no_input_rec_dend)
        dend_locs = pool.map(adjust_num_syn_, ['SC', 'PP'], chunksize=1)

        pool.terminate()
        pool.join()
        del pool

        dend_locs_dict = {} 
        for locs in dend_locs:
            dend_locs_dict.update(locs)
        
        (rec_ndend, xloc), distance = recording_loc
        stimuli_list = [0.25, 400, num_stimuli_in_train * interval_bw_stimuli_in_train, rec_ndend, xloc, rec_ndend, xloc]

        current_amp  = self.adjust_current_amplitude(model, stimuli_list)
        
        SC_dend_loc = dend_locs_dict['SC'] #[['dendrite[52]', 0.1], ['dendrite[112]', 0.3], ['dendrite[107]', 0.5], ['dendrite[54]', 0.7], ['dendrite[90]', 0.07142857142857142], ['dendrite[84]', 0.9], ['dendrite[98]', 0.7], ['dendrite[107]', 0.9285714285714286]] #dend_locs_dict['SC'] 
        PP_dend_loc = dend_locs_dict['PP'] #[['dendrite[119]', 0.5], ['dendrite[155]', 0.5], ['dendrite[152]', 0.7], ['dendrite[130]', 0.5]] #dend_locs_dict['PP']

        self.message_to_logFile += "SC dencd_loc: " + str(SC_dend_loc) + "\n"
        self.message_to_logFile += "PP dencd_loc: " + str(PP_dend_loc) + "\n"

        tstop = 1600

        pool = multiprocessing.Pool(self.npool, maxtasksperchild=1)
        theta_pathway_stimulus_= functools.partial(self.theta_pathway_stimulus, model, SC_weight, PP_weight, SC_dend_loc, PP_dend_loc, recording_loc, stimuli_params, tstop, current_amp, save_traces=True)  # save_traces=True because we want to save traces into pickle files for later use
        traces = pool.map(theta_pathway_stimulus_, ['SC', 'PP', 'SC+PP', 'depol', 'SC+depol', 'PP+depol'], chunksize=1)

        pool.terminate()
        pool.join()
        del pool

        traces_dict = {} 
        for trace in traces:
            traces_dict.update(trace)
        # print(traces_dict)


        self.plot_traces(model, traces_dict)

        print('Extracting features')

        prediction = {}
        for pathway, traces in traces_dict.items():
            if pathway != 'depol':
                features = self.extract_features(model, traces, [t_no_input_rec_dend, v_soma_no_input, v_no_input_rec_dend], stimuli_params, pathway)
                prediction.update(features)
        print(prediction)


        ''' printing to logFile'''

        filepath = self.path_results + self.test_log_filename
        self.logFile = open(filepath, 'w') # if it is opened before multiprocessing, the multiporeccing won't work under python3


        if not model.AMPA_name:
            self.logFile.write('The built in Exp2Syn is used as the AMPA component. Tau1 = ' + str(model.AMPA_tau1) + ', Tau2 = ' + str(model.AMPA_tau2) + '.\n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")
        if not model.NMDA_name:
            self.logFile.write('The default NMDA model of HippoUnit is used with Jahr, Stevens voltage dependence.\n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        self.logFile.write(self.message_to_logFile)
        self.message_to_logFile = ''


        prediction_json = copy.deepcopy(prediction)

        for key, val in list(prediction.items()):
            for ke, va in list(val.items()):
                for k, v in list(va.items()):
                    try:
                        v = str(v)
                        quantity_parts = v.split("*")
                        prediction_json[key][ke][k] = " ".join(quantity_parts)
                    except:
                        prediction_json[key][ke][k] = str(v)



        file_name_json = self.path_results + 'pathway_interaction_model_features.json'

        json.dump(prediction_json, open(file_name_json, "w"), indent=4)

        self.plot_features(prediction)



        print("Results are saved in the directory: ", self.path_results)

        efel.reset()

        return prediction

    def compute_score(self, observation, prediction, verbose=False):
        """Implementation of sciunit.Test.score_prediction."""


        score, errors, penalty_PP_depol, penalty_SC_PP = scores.ZScore_PathwayInteraction.compute(observation,prediction)

        score=scores.ZScore_PathwayInteraction(score)



        file_name_errors = self.path_results + 'pathway_interaction_errors.json'
        json.dump(errors, open(file_name_errors, "w"), indent=4)


        if self.show_plot:
            plt.show()

        final_score={'score' : str(score)}
        file_name_score= self.path_results + 'final_score.json'
        json.dump(final_score, open(file_name_score, "w"), indent=4)

        if penalty_PP_depol > 0:
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")
            self.logFile.write('PP + depolarization stimulus didn\'t generate interpretable plateau potential. A penalty (100) is added to the final score. Please have a look at the traces. \n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        if penalty_SC_PP > 0:
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")
            self.logFile.write('SC+PP stimulus didn\'t generate interpretable plateau potential. A penalty (100) is added to the final score. Please have a look at the traces. \n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        self.logFile.write("---------------------------------------------------------------------------------------------------\n")
        self.logFile.write(str(score)+'\n')
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        self.logFile.close()

        self.logFile = self.path_results + self.test_log_filename
        
        self.plot_errors(errors)

        return score

    def bind_score(self, score, model, observation, prediction):

        if self.path_figs is not None:
             score.related_data["figures"] = [self.path_figs + 'feature_values.png', self.path_figs + 'feature_errors.png',
                                        self.path_figs + 'SC+PP_plateau_half_dur.png', self.path_figs + 'PP+depol_plateau_half_dur.png',
                                        self.path_figs + 'local_traces_SC_stimulus+depol.png', self.path_figs + 'local_traces_PP_stimulus+depol.png',
                                        self.path_figs + 'local_traces_SC_PP_stimulus.png', self.path_figs + 'local_traces_PP_stimulus.png',
                                        self.path_figs + 'local_traces_SC_stimulus.png', self.path_figs + 'trace_SC_stimulus+depol.png',
                                        self.path_figs + 'trace_SC_stimulus.png', self.path_figs + 'trace_SC_PP_stimulus.png',
                                        self.path_figs + 'trace_PP_stimulus+depol.png', self.path_figs + 'trace_PP_stimulus.png',
                                        self.path_figs + 'trace_only_depol.png', self.path_figs + 'trace_subplots.png']
        score.related_data["results"] = [self.path_results + 'pathway_interaction_model_features.json', self.path_results + 'pathway_interaction_errors.json',
                                        self.path_results + 'final_score.json', self.path_results + 'test_log.txt']
        return score
