from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
#from builtins import str
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

import matplotlib
matplotlib.use('Agg')
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


try:
    copy_reg.pickle(MethodType, _pickle_method, _unpickle_method)
except:
    copyreg.pickle(MethodType, _pickle_method, _unpickle_method)


class DepolarizationBlockTest(Test):
    """
    Tests if the model enters depolarization block under current injection of increasing amplitudes.

    Parameters
    ----------
    observation : dict
        dictionary loaded from a JSON file, containing the experimental mean and std values for the features to be tested
    force_run : boolean
        If True and the pickle files containing the model's response to the simulation exists, the simulation won't be run again, traces are loaded from the pickle file
    base_directory : str
        Results will be saved here
    show_plot : boolean
        If False, plots are not displayed but still saved
    save_all : boolean
        If False, only the JSON files containing the absolute feature values, the feature error scores and the final scores, and a log file are saved, but the figures and pickle files are not.
    """

    def __init__(self,
                 observation = {'mean_Ith':None, 'Ith_std':None, 'mean_Veq': None, 'Veq_std': None}  ,
                 name="Depolarization block test" ,
                 force_run=False,
                 base_directory= None,
                show_plot=True,
                save_all=True):

        observation = self.format_data(observation)

        Test.__init__(self,observation,name)

        self.required_capabilities += (cap.ReceivesSquareCurrent_ProvidesResponse,)

        self.force_run = force_run
        self.base_directory = base_directory
        self.save_all = save_all

        self.show_plot = show_plot

        self.path_temp_data = None #added later, because model name is needed
        self.path_figs = None
        self.path_results = None

        self.figures = []

        self.npool = multiprocessing.cpu_count() - 1

        self.logFile = None
        self.test_log_filename = 'test_log.txt'

        description = "Tests if the model enters depolarization block under current injection of increasing amplitudes."

    score_type = scores.ZScore_depolblock

    def format_data(self, observation):

        for key, val in list(observation.items()):
            try:
                assert type(observation[key]) is Quantity
            except Exception as e:
                quantity_parts = val.split(" ")
                number = float(quantity_parts[0])
                units = " ".join(quantity_parts[1:])
                observation[key] = Quantity(number, units)
        return observation

    def cclamp(self, model, amp, delay, dur):

        if self.base_directory:
            self.path_temp_data = self.base_directory + 'temp_data/' + 'depol_block/' + model.name + '/'
        else:
            self.path_temp_data = model.base_directory + 'temp_data/' + 'depol_block/'

        try:
            if not os.path.exists(self.path_temp_data) and self.save_all:
                os.makedirs(self.path_temp_data)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        file_name = self.path_temp_data + 'cclamp_' + str(amp) + '.p'

        if self.force_run or (os.path.isfile(file_name) is False):

            trace = {}
            traces=[]

            t, v = model.get_vm(amp, delay, dur, 'soma', 0.5, 'soma', 0.5)

            #print "- running amplitude: " + str(amp)  + " on model: " + model.name + " at: " + str(model.soma) + "(" + str(0.5) + ")"


            #t, v = model.run_cclamp()



            trace['T'] = t
            trace['V'] = v
            trace['stim_start'] = [delay]
            trace['stim_end'] = [delay + dur]
            traces.append(trace)

            traces_results = efel.getFeatureValues(traces,
                                        ['Spikecount'])

            traces.append(traces_results)
            if self.save_all:
                pickle.dump(traces, gzip.GzipFile(file_name, "wb"))

        else:
            traces = pickle.load(gzip.GzipFile(file_name, "rb"))

        return traces

    def analyse_trace_end(self, results, Veq_index):

        Veq_trace=results[Veq_index][0]['V']
        time=numpy.array(results[Veq_index][0]['T'])
        indices1 = numpy.where(1400<=time)[0]           # we need the last 100ms of the current pulse
        indices2 = numpy.where(1500>=time)[0]
        trace_end_index_beginning = numpy.amin(indices1)
        trace_end_index_end=numpy.amax(indices2)
        trace_end=Veq_trace[trace_end_index_beginning:trace_end_index_end]      # THIS CONTAINS the last 100ms of the current pulse
        time_end=time[trace_end_index_beginning:trace_end_index_end]

        trace = {}
        traces = []
        trace['T'] = time_end
        trace['V'] = trace_end
        trace['stim_start'] = [0]
        trace['stim_end'] = [100]
        traces.append(trace)

        result = efel.getFeatureValues(traces,
                                        ['Spikecount'])
        spikecount = result[0]['Spikecount']

        return trace_end, spikecount

    def find_Ith_Veq(self, model, results, amps):

        if self.base_directory:
            self.path_figs = self.base_directory + 'figs/' + 'depol_block/' + model.name + '/'
        else:
            self.path_figs = model.base_directory + 'figs/' + 'depol_block/'

        try:
            if not os.path.exists(self.path_figs) and self.save_all:
                os.makedirs(self.path_figs)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        print("The figures are saved in the directory: ", self.path_figs)

        if self.base_directory:
            self.path_results = self.base_directory + 'results/' + 'depol_block/' + model.name + '/'
        else:
            self.path_results = model.base_directory + 'results/' + 'depol_block/'

        try:
            if not os.path.exists(self.path_results):
                os.makedirs(self.path_results)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        filepath = self.path_results + self.test_log_filename
        self.logFile = open(filepath, 'w')

        spikecount_array=numpy.array([])

        for i, amp in enumerate(amps):

            spikecount_array=numpy.append(spikecount_array, results[i][1][0]['Spikecount'])


        maximum=numpy.amax(spikecount_array)


        I_maxNumAP_index = numpy.where(spikecount_array==maximum)[0]        # this is an array if there are a lot of same values in spike_count_array

        if I_maxNumAP_index.size > 1 or I_maxNumAP_index == (spikecount_array.size) -1:     # If the max num AP is the last element, it didn`t enter depol. block
            I_maxNumAP=amps[I_maxNumAP_index][0]                                 # If Ith == None, it means it didn`t enter depol block!!!
            I_above_maxNumAP_index = None #int(I_maxNumAP_index[0])
            I_above_maxNumAP = None #I_maxNumAP
            I_below_depol_block = float('NaN')
            I_below_depol_block_index = None
            Veq=float('NaN')
            Veq_index= None#I_maxNumAP_index
            #Veq_index = int(Veq_index[0])
            '''
            plt.figure()
            plt.plot(results[spikecount_array.size-1][0]['T'],results[spikecount_array.size-1][0]['V'])
            plt.title("somatic response to the highest current intensity ("+str(amps[-1])+ " nA) \n (The model did not enter depol. block.)")
            print " the model did not enter depolarization block"
            plt.savefig(self.path_figs + 'somatic_response' + '.pdf', dpi=300)
            '''

        else:
            I_maxNumAP=amps[I_maxNumAP_index]
            I_maxNumAP=I_maxNumAP[0]
            I_above_maxNumAP_index=I_maxNumAP_index+1
            I_above_maxNumAP_index = int(I_above_maxNumAP_index[0])
            '''
            plt.figure()
            plt.plot(results[I_above_maxNumAP_index][0]['T'],results[I_above_maxNumAP_index][0]['V'])
            plt.title("somatic response at I_maxNumAP + 0.05 nA - "+str(amps[I_above_maxNumAP_index])+" nA")
            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            #plt.tick_params(labelsize=18)
            plt.savefig(self.path_figs + 'somatic_resp_at_above_I_maxNumAP' + '.pdf', dpi=300)
            '''

            trace_end, spikecount_end_of_trace = self.analyse_trace_end(results, I_above_maxNumAP_index)

            """ If there are APs at the end of the trace - the model didn't really entered depolarization block"""
            if spikecount_end_of_trace != 0:
                Veq = float('NaN')
                index = I_above_maxNumAP_index

                while spikecount_end_of_trace != 0 and index+1 < len(amps):
                    #if index == I_above_maxNumAP_index:
                    index += 1
                    trace_end, spikecount_end_of_trace = self.analyse_trace_end(results, index)
                    '''
                    plt.figure()
                    plt.plot(results[index][0]['T'],results[index][0]['V'])
                    plt.title("somatic response at "+str(amps[index])+" nA")
                    plt.xlabel("time (ms)")
                    plt.ylabel("Somatic voltage (mV)")
                    '''
                    #print spikecount_end_of_trace
                    #index += 1
                if amps[index] == amps[-1] and spikecount_end_of_trace != 0:
                    Veq = float('NaN')
                    Veq_index = None
                    I_below_depol_block_index = None
                    I_below_depol_block = float('NaN')
                elif amps[index] == amps[-1] and spikecount_end_of_trace == 0:
                    depol_block_index = index
                    Veq = numpy.average(trace_end)
                    Veq_index=index
                    I_below_depol_block_index = index-1
                    I_below_depol_block = amps[index-1]
                else:
                    depol_block_index = index
                    Veq = numpy.average(trace_end)
                    Veq_index=index
                    I_below_depol_block_index = index-1
                    I_below_depol_block = amps[index-1]
            else:
                I_below_depol_block = amps[I_maxNumAP_index[0]]
                I_below_depol_block_index = I_maxNumAP_index[0]
                Veq_index=I_above_maxNumAP_index
                Veq = numpy.average(trace_end)


        print("I_maxNumAP (the current intensity for which the model exhibited the maximum number of APs):", I_maxNumAP *nA)
        print("I_below_depol_block (the current intensity before the model enters depolarization block):", I_below_depol_block *nA)
        print("Veq (the equilibrium value during the depolarization block):", Veq * mV)

        self.logFile.write("I_maxNumAP (the current intensity for which the model exhibited the maximum number of APs) :" + str(I_maxNumAP) + " nA\n")
        self.logFile.write("I_below_depol_block (the current intensity before the model enters depolarization block): " + str(I_below_depol_block) + " nA\n")
        self.logFile.write("Veq (the equilibrium value during the depolarization block): " + str(Veq) + " mV\n")
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")


        plt.figure()
        fig = plt.gcf()
        #fig.set_size_inches(14, 12)
        plt.plot(amps,spikecount_array,'o-')
        #plt.tick_params(labelsize=20)
        plt.xlabel("I (nA)")
        #plt.ylabel("number of APs")
        plt.ylabel("num. of APs")
        plt.margins(0.01)
        if self.save_all:
            plt.savefig(self.path_figs + 'number_of_APs' + '.pdf', dpi=600, bbox_inches='tight')
            self.figures.append(self.path_figs + 'number_of_APs' + '.pdf')
        #plt.savefig(self.path_figs + 'num. of Aps' + '.pdf')


        spikecount_dict = {"current amplitudes" : list(amps), "spikecounts" : list(spikecount_array)}

        file_name_sc = self.path_results + 'current_amps_spikecounts.json'
        json.dump(spikecount_dict, open(file_name_sc, "w"))


        if I_maxNumAP_index.size > 1:
            I_maxNumAP_index = int(I_maxNumAP_index[-1])
            #Veq_index = int(Veq_index[-1])
        elif I_maxNumAP_index.size == 1:
            I_maxNumAP_index = int(I_maxNumAP_index[0])
            #Veq_index = int(Veq_index[0])


        plt.figure()
        plt.plot(results[I_maxNumAP_index][0]['T'],results[I_maxNumAP_index][0]['V'])
        if I_maxNumAP == amps[-1]:
            plt.title("somatic response to the highest current intensity ("+str(amps[-1])+ " nA) \n (The model did not enter depol. block)")
        else:
            plt.title("somatic response at I, to which the model elicited the maximum number of spikes  ("+str(amps[I_maxNumAP_index])+ " nA)")
        plt.xlabel("time (ms)")
        plt.ylabel("Somatic voltage (mV)")
        #plt.tick_params(labelsize=18)
        if self.save_all:
            plt.savefig(self.path_figs + 'somatic_resp_at_I_maxNumAP' + '.pdf', dpi=600, bbox_inches='tight')
            self.figures.append(self.path_figs + 'somatic_resp_at_I_maxNumAP' + '.pdf')

        if I_above_maxNumAP_index is not None and I_above_maxNumAP_index != Veq_index: # this plot is not needed if to the next step the model enters depol. block
            plt.figure()
            plt.plot(results[I_above_maxNumAP_index][0]['T'],results[I_above_maxNumAP_index][0]['V'])
            plt.title("somatic response at I_maxNumAP + 0.05 nA ("+str(amps[I_above_maxNumAP_index])+" nA)")
            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            #plt.tick_params(labelsize=18)
            if self.save_all:
                plt.savefig(self.path_figs + 'somatic_resp_at_above_I_maxNumAP' + '.pdf', dpi=300)
                self.figures.append(self.path_figs + 'somatic_resp_at_above_I_maxNumAP' + '.pdf')

        if I_below_depol_block_index is not None and I_above_maxNumAP_index != I_below_depol_block_index and I_above_maxNumAP_index != Veq_index:
            plt.figure()
            plt.plot(results[I_below_depol_block_index][0]['T'],results[I_below_depol_block_index][0]['V'])
            #plt.title("somatic response before depol. block - at "+str(amps[I_below_depol_block_index])+" nA")
            if I_below_depol_block == amps[-1]:
                plt.title("somatic response to the highest current intensity ("+str(amps[-1])+ " nA) \n (The model did not enter depol. block)")
            else:
                plt.title("somatic response before depol. block - at "+str(amps[I_below_depol_block_index])+" nA")

            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            #plt.tick_params(labelsize=18)
            if self.save_all:
                plt.savefig(self.path_figs + 'somatic_resp_before_depol_block' + '.pdf', dpi=300, bbox_inches='tight')
                self.figures.append(self.path_figs + 'somatic_resp_before_depol_block' + '.pdf')

        if Veq_index is not None:
            plt.figure()
            plt.plot(results[Veq_index][0]['T'],results[Veq_index][0]['V'])
            plt.title("somatic response - depol. block - at "+str(amps[Veq_index])+" nA")
            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            #plt.tick_params(labelsize=18)
            if self.save_all:
                plt.savefig(self.path_figs + 'somatic_resp_depol_block' + '.pdf', dpi=300, bbox_inches='tight')
                self.figures.append(self.path_figs + 'somatic_resp_depol_block' + '.pdf')



        if I_maxNumAP != amps[-1]:   # plot the voltage response to the highest current intensity, if not plotted yet
            plt.figure()
            plt.plot(results[-1][0]['T'],results[-1][0]['V'])
            plt.title("somatic response to the highest current intensity ("+str(amps[-1])+ " nA)")
            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            #plt.tick_params(labelsize=18)
            if self.save_all:
                plt.savefig(self.path_figs + 'somatic_resp_at_highest_amp' + '.pdf', dpi=600, bbox_inches='tight')
                self.figures.append(self.path_figs + 'somatic_resp_at_highest_amp' + '.pdf')



        if I_maxNumAP != amps[-1]:    # not meaningful if the model did not enter depol. block
            x =numpy.array([1, 2, 3])
            Ith_array = numpy.array([self.observation['mean_Ith'], I_maxNumAP, I_below_depol_block])
            labels = ['Target Ith with SD', 'model_I_maxNumAP', 'model_I_below_depol_block']
            #e = numpy.array([self.observation['Ith_std'], 0.0, 0.0])

            x2 =numpy.array([1])
            y2 = numpy.array([self.observation['mean_Ith']])
            e = numpy.array([self.observation['Ith_std']])

            plt.figure()
            plt.plot(x, Ith_array, 'o')
            plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
            plt.xticks(x, labels, rotation=10)
            plt.margins(0.2)
            plt.ylabel("Ith (nA)")
            if self.save_all:
                plt.savefig(self.path_figs + 'Ith' + '.pdf', dpi=600, bbox_inches='tight')
                self.figures.append(self.path_figs + 'Ith' + '.pdf')



            x =numpy.array([1, 2])
            Veq_array = numpy.array([self.observation['mean_Veq'], Veq])
            labels = ['Target Veq with SD',model.name]
            e = numpy.array([self.observation['Veq_std'], 0.0])

            x2 =numpy.array([1])
            y2 = numpy.array([self.observation['mean_Veq']])
            e = numpy.array([self.observation['Veq_std']])

            plt.figure()
            plt.plot(x, Veq_array, 'o')
            plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
            plt.xticks(x, labels, rotation=10)
            plt.margins(0.2)
            plt.ylabel("Veq (mV)")
            if self.save_all:
                plt.savefig(self.path_figs + 'Veq' + '.pdf', dpi=600, bbox_inches='tight')
                self.figures.append(self.path_figs + 'Veq' + '.pdf')


        #errors
        '''
        if not math.isnan(Veq):
            Veq_error=abs(Veq*mV-self.observation['mean_Veq'])/self.observation['Veq_std']
            print "The error of Veq in units of the experimental SD: ", Veq_error
        else:
            Veq_error=float('NaN')

        if not math.isnan(Ith):
            Ith_error=abs(Ith*nA-self.observation['mean_Ith'])/self.observation['Ith_std']
            print "The error of Ith in units of the experimental SD: ", Ith_error
        else:
            Ith_error=float('NaN')
        '''


        num_AP_min=13
        num_AP_max=82
        labels2 = [model.name]

        plt.figure()
        plt.axhline(y=num_AP_min, label='Min. num.of AP\n observed\n experimentally', color='green')
        plt.axhline(y=num_AP_max, label='Max. num.of AP\n observed\n experimentally', color='red')
        plt.legend()

        plt.plot([1], spikecount_array[I_maxNumAP_index], 'o')
        plt.title("For the models that doesn't enter depol block,\n the plot shows the num of AP-s for the highest current intensity")
        plt.xticks([1], labels2)
        plt.margins(0.2)
        plt.ylabel("Number of AP at Ith  (" +str(amps[I_maxNumAP_index]) + " nA)")
        if self.save_all:
            plt.savefig(self.path_figs + 'num_of_APs_at_Ith' + '.pdf', dpi=600, bbox_inches='tight')
            self.figures.append(self.path_figs + 'num_of_APs_at_Ith' + '.pdf')

        file_name_f = self.path_results + 'depol_block_features_traces.p'
        file_name_json = self.path_results + 'depol_block_model_features.json'


        features_json={}
        features_json['I_maxNumAP']=str(float(I_maxNumAP) * nA)
        features_json['I_below_depol_block']=str(float(I_below_depol_block) * nA)
        features_json['Veq']=str(float(Veq) * mV)
        json.dump(features_json, open(file_name_json, "w"), indent=4)

        features={}
        features['I_maxNumAP']=[I_maxNumAP]
        features['I_below_depol_block']=[I_below_depol_block]
        features['Veq']=[Veq]
        #features['Ith_error']=[float(Ith_error)]
        #features['Veq_error']=[float(Veq_error)]
        features['Spikecount']=spikecount_array
        features['I_maxNumAP_trace']=results[I_maxNumAP_index][0]['V']  # trace where Ith is measured (max num of APs)
        features['I_maxNumAP_time']=results[I_maxNumAP_index][0]['T']
        if I_below_depol_block_index is not None and Veq_index is not None:
            features['I_below_depol_block_trace']=results[I_below_depol_block_index][0]['V']  # trace at current intensity below the one to which the model enters depol. block
            features['I_below_depol_block_time']=results[I_below_depol_block_index][0]['T']
            features['Veq_trace']=results[Veq_index][0]['V']  # trace where Veq is measured
            features['Veq_time']=results[Veq_index][0]['T']

        if self.save_all:
            pickle.dump(features, gzip.GzipFile(file_name_f, "wb"))

        print("Results are saved in the directory: ", self.path_results)


        return I_maxNumAP, I_below_depol_block, Veq

    def validate_observation(self, observation):

        try:
            assert type(observation['mean_Ith']) is Quantity
            assert type(observation['Ith_std']) is Quantity
            assert type(observation['mean_Veq']) is Quantity
            assert type(observation['Veq_std']) is Quantity
        except Exception as e:
            raise ObservationError(("Observation must be of the form "
                                    "{'mean':float*mV,'std':float*mV}"))

    def generate_prediction(self, model, verbose=False):
        """Implementation of sciunit.Test.generate_prediction."""

        efel.reset()

        pool = multiprocessing.Pool(self.npool, maxtasksperchild=1)
        #amps = numpy.arange(0,3.55,0.05)
        amps = numpy.arange(0,1.65,0.05)

        cclamp_ = functools.partial(self.cclamp, model, delay = 500, dur = 1000)
        results = pool.map(cclamp_, amps, chunksize=1)
        #results = result.get()

        pool.terminate()
        pool.join()
        del pool

        plt.close('all') #needed to avoid overlapping of saved images when the test is run on multiple models in a for loop

        I_maxNumAP, I_below_depol_block, Veq = self.find_Ith_Veq(model, results, amps)

        prediction = {'model_I_maxNumAP':float(I_maxNumAP)*nA, 'model_I_below_depol_block':float(I_below_depol_block)*nA, 'model_Veq': float(Veq)*mV}

        efel.reset()

        return prediction

    def compute_score(self, observation, prediction, verbose=False):
        """Implementation of sciunit.Test.score_prediction."""

        final_score, I_maxNumAP_error, I_below_depol_block_error, Veq_error, I_diff_penalty = scores.ZScore_depolblock.compute(observation,prediction)

        file_name_json_err = self.path_results + 'depol_block_model_errors.json'

        errors={}
        errors['I_maxNumAP_error']=I_maxNumAP_error
        errors['I_below_depol_block_error']=I_below_depol_block_error
        errors['Veq_error']=Veq_error
        errors['I_diff_penalty']=I_diff_penalty

        json.dump(errors, open(file_name_json_err, "w"), indent=4)

        plt.figure()

        if final_score == 100.0:
            plt.plot(errors['I_maxNumAP_error'], 0, marker = 'o')
            plt.plot(final_score, 1, marker = 'o')
            plt.yticks([0,1], ['I_maxNumAP_error', 'no_depol_block_penalty'])
        else:
            plt.plot(errors['I_maxNumAP_error'], 0, marker = 'o')
            plt.plot(errors['I_below_depol_block_error'], 1, marker = 'o')
            plt.plot(errors['Veq_error'], 2, marker = 'o')
            plt.plot(errors['I_diff_penalty'], 3, marker = 'o')
            plt.yticks([0,1,2, 3], ['I_maxNumAP_error', 'I_below_depol_block_error', 'Veq error', 'I_diff_penalty'])
        plt.margins(0.2)
        plt.xlabel("error (# sd)")
        plt.title('Errors')
        if self.save_all:
            plt.savefig(self.path_figs + 'Errors' + '.pdf', dpi=600, bbox_inches='tight')
            self.figures.append(self.path_figs + 'Errors' + '.pdf')

        if errors['I_diff_penalty'] != 0:
            print('According to the experiment I_maxNumAP and I_below_depol_block should be equal. If they are not equal a penalty is applied (200 * the difference between them in nA).')
            self.logFile.write('According to the experiment I_maxNumAP and I_below_depol_block should be equal. If they are not equal a penalty is applied  (200 * the difference between them in nA).\n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        score=scores.ZScore_depolblock(final_score)

        score_json= {'score' : final_score}
        file_name_score = self.path_results + 'depol_block_final_score.json'
        json.dump(score_json, open(file_name_score, "w"), indent=4)

        if self.show_plot:
            plt.show()
        if final_score == 100:
            print('The model did not enter depolarization block.')
            self.logFile.write('The model did not enter depolarization block.\n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        self.logFile.write(str(score)+'\n')
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")
        self.logFile.close()

        self.logFile = self.path_results + self.test_log_filename

        return score

    def bind_score(self, score, model, observation, prediction):

        self.figures.append(self.path_results + 'depol_block_model_errors.json')
        self.figures.append(self.path_results + 'depol_block_model_features.json')
        self.figures.append(self.path_results + 'depol_block_final_score.json')
        self.figures.append(self.path_results + self.test_log_filename)

        score.related_data["figures"] = self.figures
        score.related_data["results"] = [self.path_results + 'depol_block_model_errors.json', self.path_results + 'depol_block_model_features.json', self.path_results + 'depol_block_final_score.json', self.path_results + 'depol_block_features_traces.p']
        return score
