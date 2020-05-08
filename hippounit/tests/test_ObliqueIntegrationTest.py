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


class ObliqueIntegrationTest(Test):
    """
    Tests the signal integration in oblique dendrites for increasing number of synchronous and asynchronous inputs.

    Parameters
    ----------
    observation : dict
        dictionary loaded from a JSON file, containing the experimental mean and std values for the features to be tested
    force_run_bin_search : boolean
        If True and the pickle files containing the synaptic weights for each dendritic location exist, the binary search for the proper synaptic weights won't run again, weights are loaded from the pickle file
    force_run_bin_synapse : boolean
        If True and the pickle files containing the model's response to the simulation exists, the simulation won't be run again, traces are loaded from the pickle file
    base_directory : str
        Results will be saved here
    show_plot : boolean
        If False, plots are not displayed but still saved
    save_all : boolean
        If False, only the JSON files containing the absolute feature values, the feature error scores and the final scores, and a log file are saved, but the figures and pickle files are not.
    trunk_origin : list
        first element : name of the section from which the trunk originates, second element : position on section (E.g. ['soma[5]', 1]). If not set by the user, the end of the default soma section is used.
    """

    def __init__(self,
                 observation = {'mean_threshold':None,'threshold_sem':None, 'threshold_std': None,
                                'mean_prox_threshold':None,'prox_threshold_sem':None, 'prox_threshold_std': None,
                                'mean_dist_threshold':None,'dist_threshold_sem':None, 'dist_threshold_std': None,
                                'mean_nonlin_at_th':None,'nonlin_at_th_sem':None, 'nonlin_at_th_std': None,
                                'mean_nonlin_suprath':None,'nonlin_suprath_sem':None, 'nonlin_suprath_std': None,
                                'mean_peak_deriv':None,'peak_deriv_sem':None, 'peak_deriv_std': None,
                                'mean_amp_at_th':None,'amp_at_th_sem':None, 'amp_at_th_std': None,
                                'mean_time_to_peak':None,'time_to_peak_sem':None, 'time_to_peak_std': None,
                                'mean_async_nonlin':None,'async_nonlin_sem':None, 'async_nonlin_std': None}  ,
                 name="Oblique integration test" ,
                 force_run_synapse=False,
                 force_run_bin_search=False,
                 base_directory= None,
                 show_plot=True,
                 save_all = True,
                 trunk_origin = None):

        observation = self.format_data(observation)
        observation = self.add_std_to_observation(observation)

        Test.__init__(self, observation, name)

        self.required_capabilities = (cap.ProvidesGoodObliques, cap.ReceivesSynapse,) # +=

        self.force_run_synapse = force_run_synapse
        self.force_run_bin_search = force_run_bin_search
        self.show_plot = show_plot
        self.save_all = save_all

        self.base_directory = base_directory

        self.path_figs = None    #added later, because model name is needed
        self.path_results = None
        self.trunk_origin = trunk_origin

        self.logFile = None
        self.test_log_filename = 'test_log.txt'

        self.npool = multiprocessing.cpu_count() - 1

        self.max_num_syn = 10

        self.threshold_index = 5  #threshold input number for dendritic spike generation  - index 0 is 0 input

        description = "Tests the signal integration in oblique dendrites for increasing number of synchronous and asynchronous inputs"

    score_type = scores.ZScore_ObliqueIntegration

    def format_data(self, observation):

        for key, val in list(observation.items()):
            try:
                assert type(observation[key]) is Quantity
            except Exception as e:
                try:
                    observation[key] = float(val)
                except Exception as e:
                    quantity_parts = val.split(" ")
                    number = float(quantity_parts[0])
                    units = " ".join(quantity_parts[1:])
                    observation[key] = Quantity(number, units)
        return observation

    def add_std_to_observation (self, observation):

        observation["threshold_std"]=float(observation["threshold_sem"]*math.sqrt(observation["exp_n"]))*mV
        observation["prox_threshold_std"]=float(observation["prox_threshold_sem"]*math.sqrt(observation["prox_n"]))*mV
        observation["dist_threshold_std"]=float(observation["dist_threshold_sem"]*math.sqrt(observation["dist_n"]))*mV
        observation["nonlin_at_th_std"]=float(observation["nonlin_at_th_sem"]*math.sqrt(observation["exp_n"]))
        observation["nonlin_suprath_std"]=float(observation["nonlin_suprath_sem"]*math.sqrt(observation["exp_n"]))
        observation["peak_deriv_std"]=float(observation["peak_deriv_sem"]*math.sqrt(observation["exp_n"]))*V /s
        observation["amp_at_th_std"]=float(observation["amp_at_th_sem"]*math.sqrt(observation["exp_n"]))*mV
        observation["time_to_peak_std"]=float(observation["time_to_peak_sem"]*math.sqrt(observation["exp_n"]))*ms
        observation["async_nonlin_std"]=float(observation["async_nonlin_sem"]*math.sqrt(observation["async_n"]))
        return observation

    def analyse_syn_traces(self, model, t, v, v_dend, threshold):

        trace = {}
        trace['T'] = t
        trace['V'] = v
        trace['stim_start'] = [model.start]
        trace['stim_end'] = [model.start+500]  # should be changed
        traces = [trace]

        trace_dend = {}
        trace_dend['T'] = t
        trace_dend['V'] = v_dend
        trace_dend['stim_start'] = [model.start]
        trace_dend['stim_end'] = [model.start+500]
        traces_dend = [trace_dend]

        efel.setThreshold(threshold)
        traces_results_dend = efel.getFeatureValues(traces_dend,['Spikecount_stimint'], raise_warnings=True)
        traces_results = efel.getFeatureValues(traces,['Spikecount_stimint'], raise_warnings=True)
        spikecount_dend=traces_results_dend[0]['Spikecount_stimint']
        spikecount=traces_results[0]['Spikecount_stimint']

        result = [traces, traces_dend, spikecount, spikecount_dend]

        return result


    def run_synapse(self,model, dend_loc_num_weight, interval):

        ndend, xloc, loc_type, num, weight = dend_loc_num_weight

        if self.base_directory:
            path = self.base_directory + 'temp_data/' + 'oblique_integration/synapse/' + model.name + '/'
        else:
            path = model.base_directory + 'temp_data/' + 'oblique_integration/synapse/'


        try:
            if not os.path.exists(path) and self.save_all:
                os.makedirs(path)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        if interval>0.1:
            file_name = path + 'synapse_async_' + str(num)+ '_' + str(ndend)+ '_' + str(xloc) + '.p'
        else:
            file_name = path + 'synapse_' + str(num)+ '_' + str(ndend)+ '_' + str(xloc) + '.p'

        if self.force_run_synapse or (os.path.isfile(file_name) is False):

            print("- number of inputs:", num, "dendrite:", ndend, "xloc:", xloc)


            t, v, v_dend = model.run_multiple_synapse_get_vm([ndend, xloc, loc_type], interval, num, weight)

            result = self.analyse_syn_traces(model, t, v, v_dend, model.threshold)

            if self.save_all:
                pickle.dump(result, gzip.GzipFile(file_name, "wb"))

        else:
            result = pickle.load(gzip.GzipFile(file_name, "rb"))

        return result


    def syn_binsearch(self, model, dend_loc, interval, number, weight):


        t, v, v_dend = model.run_multiple_synapse_get_vm(dend_loc, interval, number, weight)

        return t, v, v_dend

    def binsearch(self, model, dend_loc0):


        if self.base_directory:
            path_bin_search = self.base_directory + 'temp_data/' + 'oblique_integration/bin_search/' + model.name + '/'
        else:
            path_bin_search = model.base_directory + 'temp_data/' + 'oblique_integration/bin_search/'

        try:
            if not os.path.exists(path_bin_search) and self.save_all:
                os.makedirs(path_bin_search)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        interval=0.1

        file_name = path_bin_search + 'weight_' +str(dend_loc0[0])+ '_' + str(dend_loc0[1]) + '.p'


        if self.force_run_bin_search or (os.path.isfile(file_name) is False):

            c_minmax = model.c_minmax
            c_step_start = model.c_step_start
            c_step_stop= model.c_step_stop
            #c_stim=numpy.arange()


            found = False


            while c_step_start >= c_step_stop and not found:

                c_stim = numpy.arange(c_minmax[0], c_minmax[1], c_step_start)

                #print c_stim
                first = 0
                last = numpy.size(c_stim, axis=0)-1

                num = [self.threshold_index - 1,self.threshold_index]   # [4,5]


                while first<=last and not found:

                    midpoint = (first + last)//2
                    result=[]

                    for n in num:


                        pool_syn = multiprocessing.Pool(1, maxtasksperchild = 1)    # I use multiprocessing to keep every NEURON related task in independent processes

                        t, v, v_dend = pool_syn.apply(self.syn_binsearch, args = (model, dend_loc0, interval, n, c_stim[midpoint]))
                        pool_syn.terminate()
                        pool_syn.join()
                        del pool_syn

                        result.append(self.analyse_syn_traces(model, t, v, v_dend, model.threshold))
                        #print result

                    if result[0][3]==0 and result[1][3]>=1:
                        found = True
                    else:
                        if (result[0][3]>=1 and result[1][3]>=1) or (result[0][3]>=1 and result[1][3]==0):
                            last = midpoint-1

                        elif result[0][3]==0 and result[1][3]==0:
                            first = midpoint+1

                c_step_start=c_step_start/2

                if found:
                    if result[1][2]>=1 :            # somatic AP is generated
                        found = None
                        break

            if not found:
                if result[0][3]>=1 and result[1][3]>=1:
                    found = 'always spike'
                if result[0][3]==0 and result[1][3]==0:
                    found = 'no spike'


            binsearch_result=[found, c_stim[midpoint]]

            if self.save_all:
                pickle.dump(binsearch_result, gzip.GzipFile(file_name, "wb"))

        else:
            binsearch_result = pickle.load(gzip.GzipFile(file_name, "rb"))


        return binsearch_result

    def calcs_plots(self, model, results, dend_loc000, dend_loc_num_weight):

        experimental_mean_threshold=self.observation['mean_threshold']
        threshold_SEM=self.observation['threshold_sem']
        threshold_SD=self.observation['threshold_std']

        threshold_prox=self.observation['mean_prox_threshold']
        threshold_prox_SEM=self.observation['prox_threshold_sem']
        threshold_prox_SD=self.observation['prox_threshold_std']

        threshold_dist=self.observation['mean_dist_threshold']
        threshold_dist_SEM=self.observation['dist_threshold_sem']
        threshold_dist_SD=self.observation['dist_threshold_std']

        exp_mean_nonlin=self.observation['mean_nonlin_at_th']
        nonlin_SEM=self.observation['nonlin_at_th_sem']
        nonlin_SD=self.observation['nonlin_at_th_std']

        suprath_exp_mean_nonlin=self.observation['mean_nonlin_suprath']
        suprath_nonlin_SEM=self.observation['nonlin_suprath_sem']
        suprath_nonlin_SD=self.observation['nonlin_suprath_std']

        exp_mean_peak_deriv=self.observation['mean_peak_deriv']
        deriv_SEM=self.observation['peak_deriv_sem']
        deriv_SD=self.observation['peak_deriv_std']

        exp_mean_amp=self.observation['mean_amp_at_th']
        amp_SEM= self.observation['amp_at_th_sem']
        amp_SD=self.observation['amp_at_th_std']

        exp_mean_time_to_peak=self.observation['mean_time_to_peak']
        exp_mean_time_to_peak_SEM=self.observation['time_to_peak_sem']
        exp_mean_time_to_peak_SD=self.observation['time_to_peak_std']


        if self.base_directory:
            self.path_figs = self.base_directory + 'figs/' + 'oblique_integration/' + model.name + '/'
        else:
            self.path_figs = model.base_directory + 'figs/' + 'oblique_integration/'

        try:
            if not os.path.exists(self.path_figs) and self.save_all:
                os.makedirs(self.path_figs)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        print("The figures are saved in the directory: ", self.path_figs)

        stop=len(dend_loc_num_weight)+1
        sep=numpy.arange(0,stop,self.max_num_syn+1)
        sep_results=[]


        num = numpy.arange(0,self.max_num_syn+1)

        for i in range (0,len(dend_loc000)):
            sep_results.append(results[sep[i]:sep[i+1]])             # a list that contains the results of the 10 locations seperately (in lists)

        # sep_results[0]-- the first location
        # sep_results[0][5] -- the first location at 5 input
        # sep_results[0][1][0] -- the first location at 1 input, SOMA
        # sep_results[0][1][1] -- the first location at 1 input, dendrite
        # sep_results[0][1][1][0] -- just needed


        # subplots settings

        frames = len(dend_loc000)
        columns = 2
        width_ratios=[1]*columns
        rows = int(numpy.ceil(frames/float(columns)))
        height_ratios=[1]*rows

        fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.92, bottom=0.04, left=0.07, right=0.97, hspace=0.75, wspace=0.3)
        t=fig.suptitle('Synchronous inputs (red: dendritic trace, black: somatic trace)')
        for i in range (0,len(dend_loc000)):
            plt.subplot(gs[i])
            #plt.subplots_adjust(top=0.85, bottom=0.04, left=0.07, right=0.97, hspace=0.75, wspace=0.3)
            for j in num:
                plt.plot(sep_results[i][j][0][0]['T'],sep_results[i][j][0][0]['V'], 'k')       # somatic traces
                plt.plot(sep_results[i][j][1][0]['T'],sep_results[i][j][1][0]['V'], 'r')        # dendritic traces
            plt.title(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1])+')')

            plt.xlabel("time (ms)")
            plt.ylabel("Voltage (mV)")
            plt.xlim(sep_results[i][0][0][0]['stim_start'][0] - 10, sep_results[i][0][0][0]['stim_start'][0]+150)
            #plt.tick_params(labelsize=20)

        #fig0 = plt.gcf()
        #fig0.set_size_inches(16, 24)
        if self.save_all:
            plt.savefig(self.path_figs + 'traces_sync' + '.pdf', dpi=600, bbox_extra_artists=(t,), bbox_inches='tight')

        frames = len(dend_loc000)
        columns = 2
        width_ratios=[1]*columns
        rows = int(numpy.ceil(frames/float(columns)))
        height_ratios=[1]*rows
        fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.92, bottom=0.04, left=0.07, right=0.97, hspace=0.75, wspace=0.3)

        fig.suptitle('Synchronous inputs')
        for i in range (0,len(dend_loc000)):
            plt.subplot(gs[i])
            #plt.subplots_adjust(hspace = 0.5)
            for j, number in enumerate(num):
                plt.plot(sep_results[i][j][0][0]['T'],sep_results[i][j][0][0]['V'], 'k')       # somatic traces
            plt.title(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1]) + ')')

            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            plt.xlim(sep_results[i][0][0][0]['stim_start'][0] - 10, sep_results[i][0][0][0]['stim_start'][0]+150)
            #plt.tick_params(labelsize=20)
        #fig0 = plt.gcf()
        #fig0.set_size_inches(16, 24)
        if self.save_all:
            plt.savefig(self.path_figs + 'somatic_traces_sync' + '.pdf', dpi=600, bbox_inches='tight')

        soma_depol=numpy.array([])
        soma_depols=[]
        sep_soma_depols=[]
        #dV_dt=[]
        #sep_dV_dt=[]
        soma_max_depols=numpy.array([])
        soma_expected=numpy.array([])
        sep_soma_max_depols=[]
        sep_soma_expected=[]
        max_dV_dt=numpy.array([])
        sep_max_dV_dt=[]
        max_dV_dt_index=numpy.array([],dtype=numpy.int64)
        sep_threshold=numpy.array([])
        prox_thresholds=numpy.array([])
        dist_thresholds=numpy.array([])
        peak_dV_dt_at_threshold=numpy.array([])
        nonlin=numpy.array([])
        suprath_nonlin=numpy.array([])
        amp_at_threshold=[]
        sep_time_to_peak=[]
        time_to_peak_at_threshold=numpy.array([])
        time_to_peak=numpy.array([])
        threshold_index_dV_dt_jump=numpy.array([])

        for i in range(0, len(sep_results)):
            for j in range(0,self.max_num_syn+1):
        # calculating somatic depolarization and first derivative

                if not numpy.array_equal(sep_results[i][j][0][0]['T'], sep_results[i][0][0][0]['T']):    #if the  time vectors are not equal, the traces are resampled with fixed time step
                    dt = 0.025

                    time_vector = numpy.arange(sep_results[i][j][0][0]['T'][0], sep_results[i][j][0][0]['T'][-1], dt)  #from the first to the last element of the original time vector

                    interp_trace_input = numpy.interp(time_vector, sep_results[i][j][0][0]['T'], sep_results[i][j][0][0]['V'])
                    interp_trace_no_input = numpy.interp(time_vector, sep_results[i][0][0][0]['T'], sep_results[i][0][0][0]['V'])

                    soma_depol=interp_trace_input - interp_trace_no_input

                    dt=numpy.diff(time_vector)
                    dV=numpy.diff(interp_trace_input)

                    if i== 0 and j==1:
                        print("Voltage traces are resampled using linear interpolation")

                else:

                    soma_depol=sep_results[i][j][0][0]['V'] - sep_results[i][0][0][0]['V']
                    time_vector = sep_results[i][j][0][0]['T']

                    dt=numpy.diff(sep_results[i][j][0][0]['T'] )
                    dV=numpy.diff(sep_results[i][j][0][0]['V'] )

                soma_depols.append(soma_depol)

                soma_max_depols=numpy.append(soma_max_depols,numpy.amax(soma_depol))


                deriv=dV/dt
                #dV_dt.append(deriv)

                max_dV_dt=numpy.append(max_dV_dt, numpy.amax(deriv))

                if j==0:
                    soma_expected=numpy.append(soma_expected,0)
                else:
                    soma_expected=numpy.append(soma_expected,soma_max_depols[1]*j)

                if j!=0:
                    peak=numpy.amax(soma_depol)

                    peak_index=numpy.where(soma_depol==peak)[0]

                    peak_time=time_vector[peak_index]   # sep_results[i][j][0][0]['T'][peak_index]

                    t_to_peak=peak_time - model.start#150
                    time_to_peak = numpy.append(time_to_peak, t_to_peak)
                else:
                    time_to_peak = numpy.append(time_to_peak, 0)

            diff_max_dV_dt=numpy.diff(max_dV_dt)

            # calculate the threshold according to the experiment - highest jump in max dV/dt value
            # but synaptic weight is adjusted to have the threshold of dendritic spike generation at 5 inputs,
            # so this will be used in feature calculations
            threshold_index0=numpy.where(diff_max_dV_dt==numpy.amax(diff_max_dV_dt))[0]
            threshold_index0=numpy.add(threshold_index0,1)
            threshold_index0=threshold_index0[0]

            if sep_results[i][threshold_index0][3] > 1 and sep_results[i][threshold_index0-1][3]==1:    #double spikes can cause bigger jump in dV/dt than the first single spike, to find the threshol, we want to eliminate this, but we also need the previous input level to generate spike
                threshold_index_dV_dt_jump=numpy.where(diff_max_dV_dt==numpy.amax(diff_max_dV_dt[1:threshold_index0-1]))
                threshold_index_dV_dt_jump=numpy.add(threshold_index_dV_dt_jump,1)
                threshold_index_dV_dt_jump=threshold_index_dV_dt_jump[0]
            else:
                threshold_index_dV_dt_jump=threshold_index0

            threshold_dV_dt_jump=soma_expected[threshold_index_dV_dt_jump]

            threshold_index = self.threshold_index #threshold of dendritic spike generation is at 5 inputs - weight adjusted for this


            threshold = soma_expected[threshold_index] #threshold of dendritic spike generation is at 5 inputs - weight adjusted for this

            #lists or numpy arrays (for those that are averaged) of features for all the locations
            sep_soma_depols.append(soma_depols)
            #sep_dV_dt.append(dV_dt)
            sep_soma_max_depols.append(soma_max_depols)
            sep_soma_expected.append(soma_expected)
            sep_max_dV_dt.append(max_dV_dt)
            sep_threshold=numpy.append(sep_threshold, threshold)
            peak_dV_dt_at_threshold=numpy.append(peak_dV_dt_at_threshold,max_dV_dt[threshold_index])
            nonlin=numpy.append(nonlin, soma_max_depols[threshold_index]/ soma_expected[threshold_index]*100)  #degree of nonlinearity
            suprath_nonlin=numpy.append(suprath_nonlin, soma_max_depols[threshold_index+1]/ soma_expected[threshold_index+1]*100)  #degree of of nonlinearity
            amp_at_threshold=numpy.append(amp_at_threshold, soma_max_depols[threshold_index])
            sep_time_to_peak.append(time_to_peak)
            time_to_peak_at_threshold=numpy.append(time_to_peak_at_threshold, time_to_peak[threshold_index])

            soma_depols=[]
            #dV_dt=[]
            soma_max_depols=numpy.array([])
            soma_expected=numpy.array([])
            max_dV_dt=numpy.array([])
            threshold_index=numpy.array([])
            threshold_index_dV_dt_jump=numpy.array([])
            threshold_index0=numpy.array([])
            time_to_peak=numpy.array([])

        prox_thresholds = numpy.array([]) #sep_threshold[0::2]
        dist_thresholds = numpy.array([]) #sep_threshold[1::2]
        for i in range(len(sep_threshold)):
            if dend_loc000[i][2] == 'prox':
                prox_thresholds = numpy.append(prox_thresholds, sep_threshold[i])
            if dend_loc000[i][2] == 'dist':
                dist_thresholds = numpy.append(dist_thresholds, sep_threshold[i])


        threshold_errors = numpy.array([abs(experimental_mean_threshold - threshold_errors*mV)/threshold_SD  for threshold_errors in sep_threshold])     # does the same calculation on every element of a list  #x = [1,3,4,5,6,7,8] t = [ t**2 for t in x ]
        prox_threshold_errors=numpy.array([abs(threshold_prox - prox_threshold_errors*mV)/threshold_prox_SD  for prox_threshold_errors in prox_thresholds])        # and I could easily make it a numpy array : t = numpy.array([ t**2 for t in x ])
        dist_threshold_errors=numpy.array([abs(threshold_dist - dist_threshold_errors*mV)/threshold_dist_SD  for dist_threshold_errors in dist_thresholds])
        peak_deriv_errors=numpy.array([abs(exp_mean_peak_deriv - peak_deriv_errors*mV /ms )/deriv_SD  for peak_deriv_errors in peak_dV_dt_at_threshold])
        nonlin_errors=numpy.array([abs(exp_mean_nonlin- nonlin_errors)/nonlin_SD  for nonlin_errors in nonlin])
        suprath_nonlin_errors=numpy.array([abs(suprath_exp_mean_nonlin- suprath_nonlin_errors)/suprath_nonlin_SD  for suprath_nonlin_errors in suprath_nonlin])
        amplitude_errors=numpy.array([abs(exp_mean_amp- amplitude_errors*mV)/amp_SD  for amplitude_errors in amp_at_threshold])
        time_to_peak_errors=numpy.array([abs(exp_mean_time_to_peak- time_to_peak_errors*ms)/exp_mean_time_to_peak_SD  for time_to_peak_errors in time_to_peak_at_threshold])


        # means and SDs
        """ ZScore:obliqueIntegration is used instead"""
        """
        mean_threshold_errors=numpy.mean(threshold_errors)
        mean_prox_threshold_errors=numpy.mean(prox_threshold_errors)
        mean_dist_threshold_errors=numpy.mean(dist_threshold_errors)
        mean_peak_deriv_errors=numpy.mean(peak_deriv_errors)
        mean_nonlin_errors=numpy.mean(nonlin_errors)
        suprath_mean_nonlin_errors=numpy.mean(suprath_nonlin_errors)
        mean_amplitude_errors=numpy.mean(amplitude_errors)
        mean_time_to_peak_errors=numpy.mean(time_to_peak_errors)

        sd_threshold_errors=numpy.std(threshold_errors)
        sd_prox_threshold_errors=numpy.std(prox_threshold_errors)
        sd_dist_threshold_errors=numpy.std(dist_threshold_errors)
        sd_peak_deriv_errors=numpy.std(peak_deriv_errors)
        sd_nonlin_errors=numpy.std(nonlin_errors)
        suprath_sd_nonlin_errors=numpy.std(suprath_nonlin_errors)
        sd_amplitude_errors=numpy.std(amplitude_errors)
        sd_time_to_peak_errors=numpy.std(time_to_peak_errors)
        """

        #try:
        mean_sep_threshold=float(numpy.mean(sep_threshold)) *mV
        mean_prox_thresholds=float(numpy.mean(prox_thresholds)) *mV
        mean_dist_thresholds=numpy.mean(dist_thresholds) *mV
        mean_peak_dV_dt_at_threshold=numpy.mean(peak_dV_dt_at_threshold) *V /s
        mean_nonlin=numpy.mean(nonlin)
        suprath_mean_nonlin=numpy.mean(suprath_nonlin)
        mean_amp_at_threshold=numpy.mean(amp_at_threshold) *mV
        mean_time_to_peak_at_threshold=numpy.mean(time_to_peak_at_threshold) *ms

        sd_sep_threshold=float(numpy.std(sep_threshold)) *mV
        sd_prox_thresholds=float(numpy.std(prox_thresholds)) *mV
        sd_dist_thresholds=numpy.std(dist_thresholds) *mV
        sd_peak_dV_dt_at_threshold=numpy.std(peak_dV_dt_at_threshold) *V /s
        sd_nonlin=numpy.std(nonlin)
        suprath_sd_nonlin=numpy.std(suprath_nonlin)
        sd_amp_at_threshold=numpy.std(amp_at_threshold) *mV
        sd_time_to_peak_at_threshold=numpy.std(time_to_peak_at_threshold) *ms
        '''
        except RuntimeWarning:
            mean_sep_threshold=float('nan') *mV
            mean_prox_thresholds=float('nan') *mV
            mean_dist_thresholds=float('nan') *mV
            mean_peak_dV_dt_at_threshold=float('nan')*V /s
            mean_nonlin=float('nan')
            suprath_mean_nonlin=float('nan')
            mean_amp_at_threshold=float('nan') *mV
            mean_time_to_peak_at_threshold=float('nan')*ms

            sd_sep_threshold=float('nan') *mV
            sd_prox_thresholds=float('nan') *mV
            sd_dist_thresholds=float('nan') *mV
            sd_peak_dV_dt_at_threshold=float('nan')*V /s
            sd_nonlin=float('nan')
            suprath_sd_nonlin=float('nan')
            sd_amp_at_threshold=float('nan') *mV
            sd_time_to_peak_at_threshold=float('nan')*ms
        '''
        depol_input=numpy.array([])
        mean_depol_input=[]
        SD_depol_input=[]
        SEM_depol_input=[]

        expected_depol_input=numpy.array([])
        expected_mean_depol_input=[]

        prox_depol_input=numpy.array([])
        prox_mean_depol_input=[]
        prox_SD_depol_input=[]
        prox_SEM_depol_input=[]

        prox_expected_depol_input=numpy.array([])
        prox_expected_mean_depol_input=[]

        dist_depol_input=numpy.array([])
        dist_mean_depol_input=[]
        dist_SD_depol_input=[]
        dist_SEM_depol_input=[]

        dist_expected_depol_input=numpy.array([])
        dist_expected_mean_depol_input=[]

        peak_deriv_input=numpy.array([])
        mean_peak_deriv_input=[]
        SD_peak_deriv_input=[]
        SEM_peak_deriv_input=[]
        n=len(sep_soma_max_depols)

        prox_sep_soma_max_depols=[]
        dist_sep_soma_max_depols=[]
        prox_sep_soma_expected=[]
        dist_sep_soma_expected=[]

        for i in range(len(sep_soma_max_depols)):
            if dend_loc000[i][2] == 'prox':
                prox_sep_soma_max_depols.append(sep_soma_max_depols[i])
                prox_sep_soma_expected.append(sep_soma_expected[i])
            if dend_loc000[i][2] == 'dist':
                dist_sep_soma_max_depols.append(sep_soma_max_depols[i])
                dist_sep_soma_expected.append(sep_soma_expected[i])

        prox_n=len(prox_sep_soma_max_depols)
        dist_n=len(dist_sep_soma_max_depols)
        # needed for input-output curves (expected-measured EPSP)
        #mean of depols for each input level
        for j in range (0, self.max_num_syn+1):
            for i in range (0, len(sep_soma_max_depols)):
                depol_input=numpy.append(depol_input,sep_soma_max_depols[i][j])
                expected_depol_input=numpy.append(expected_depol_input,sep_soma_expected[i][j])
                peak_deriv_input=numpy.append(peak_deriv_input,sep_max_dV_dt[i][j])
            mean_depol_input.append(numpy.mean(depol_input))
            expected_mean_depol_input.append(numpy.mean(expected_depol_input))
            mean_peak_deriv_input.append(numpy.mean(peak_deriv_input))
            SD_depol_input.append(numpy.std(depol_input))
            SEM_depol_input.append(numpy.std(depol_input)/math.sqrt(n))
            SD_peak_deriv_input.append(numpy.std(peak_deriv_input))
            SEM_peak_deriv_input.append(numpy.std(peak_deriv_input)/math.sqrt(n))
            depol_input=numpy.array([])
            expected_depol_input=numpy.array([])
            peak_deriv_input=numpy.array([])

        for j in range (0, self.max_num_syn+1):
            for i in range (0, len(prox_sep_soma_max_depols)):
                prox_depol_input=numpy.append(prox_depol_input,prox_sep_soma_max_depols[i][j])
                prox_expected_depol_input=numpy.append(prox_expected_depol_input,prox_sep_soma_expected[i][j])
            prox_mean_depol_input.append(numpy.mean(prox_depol_input))
            prox_expected_mean_depol_input.append(numpy.mean(prox_expected_depol_input))
            prox_SD_depol_input.append(numpy.std(prox_depol_input))
            prox_SEM_depol_input.append(numpy.std(prox_depol_input)/math.sqrt(prox_n))
            prox_depol_input=numpy.array([])
            prox_expected_depol_input=numpy.array([])

        for j in range (0, self.max_num_syn+1):
            for i in range (0, len(dist_sep_soma_max_depols)):
                dist_depol_input=numpy.append(dist_depol_input,dist_sep_soma_max_depols[i][j])
                dist_expected_depol_input=numpy.append(dist_expected_depol_input,dist_sep_soma_expected[i][j])
            dist_mean_depol_input.append(numpy.mean(dist_depol_input))
            dist_expected_mean_depol_input.append(numpy.mean(dist_expected_depol_input))
            dist_SD_depol_input.append(numpy.std(dist_depol_input))
            dist_SEM_depol_input.append(numpy.std(dist_depol_input)/math.sqrt(dist_n))
            dist_depol_input=numpy.array([])
            dist_expected_depol_input=numpy.array([])

        EPSPs ={'proximal_expected_mean' : prox_expected_mean_depol_input,
                'proximal_measured_mean' : prox_mean_depol_input,
                'proximal_sd' : prox_SD_depol_input,
                'proximal_sem' : prox_SEM_depol_input,
                'distal_expected_mean' : dist_expected_mean_depol_input,
                'distal_measured_mean' : dist_mean_depol_input,
                'distal_sd' : dist_SD_depol_input,
                'distal_sem' : dist_SEM_depol_input,
                'all_expected_mean' : expected_mean_depol_input,
                'all_measured_mean' : mean_depol_input,
                'all_sd' : SD_depol_input,
                'all_sem' : SEM_depol_input}


        plt.figure(figsize = (210/25.4, 210/25.4))

        plt.title('Synchronous inputs')

        # Expected EPSP - Measured EPSP plot
        try:
            colormap = plt.cm.spectral      #http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
        except:
            colormap = plt.cm.nipy_spectral
        plt.gca().set_prop_cycle(plt.cycler('color', colormap(numpy.linspace(0, 0.9, len(sep_results)))))
        for i in range (0, len(sep_results)):

            plt.plot(sep_soma_expected[i],sep_soma_max_depols[i], '-o', label='dend '+str(dend_loc000[i][0])+ ' loc ' +str(dend_loc000[i][1]))
            plt.plot(sep_soma_expected[i],sep_soma_expected[i], 'k--')         # this gives the linear line
            plt.xlabel("expected EPSP (mV)")
            plt.ylabel("measured EPSP (mV)")
            #plt.legend(loc=2, prop={'size':10})
        #fig = plt.gcf()
        #fig.set_size_inches(12, 12)
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        if self.save_all:
            plt.savefig(self.path_figs + 'input_output_curves_sync' + '.pdf', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')


        columns = 1
        width_ratios=[1]*columns
        rows = 3
        height_ratios=[1]*rows

        fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.92, bottom=0.04, left=0.07, right=0.97, hspace=0.5, wspace=0.3)

        plt.suptitle('Synchronous inputs')

        plt.subplot(gs[0])
        plt.errorbar(expected_mean_depol_input, mean_depol_input, yerr=SD_depol_input, linestyle='-', marker='o', color='red', label='SD')
        plt.errorbar(expected_mean_depol_input, mean_depol_input, yerr=SEM_depol_input, linestyle='-', marker='o', color='blue', label='SEM')
        plt.plot(expected_mean_depol_input,expected_mean_depol_input, 'k--')         # this gives the linear line
        plt.margins(0.1)
        plt.legend(loc=2)
        plt.title("Summary plot of mean input-output curve for all locations")
        plt.xlabel("expected EPSP (mV)")
        plt.ylabel("measured EPSP (mV)")

        plt.subplot(gs[1])
        plt.errorbar(prox_expected_mean_depol_input, prox_mean_depol_input, yerr=prox_SD_depol_input, linestyle='-', marker='o', color='red', label='SD')
        plt.errorbar(prox_expected_mean_depol_input, prox_mean_depol_input, yerr=prox_SEM_depol_input, linestyle='-', marker='o', color='blue', label='SEM')
        plt.plot(prox_expected_mean_depol_input,prox_expected_mean_depol_input, 'k--')         # this gives the linear line
        plt.margins(0.1)
        plt.legend(loc=2)
        plt.title("Summary plot of mean input-output curve for proximal locations")
        plt.xlabel("expected EPSP (mV)")
        plt.ylabel("measured EPSP (mV)")

        plt.subplot(gs[2])
        plt.errorbar(dist_expected_mean_depol_input, dist_mean_depol_input, yerr=dist_SD_depol_input, linestyle='-', marker='o', color='red', label='SD')
        plt.errorbar(dist_expected_mean_depol_input, dist_mean_depol_input, yerr=dist_SEM_depol_input, linestyle='-', marker='o', color='blue', label='SEM')
        plt.plot(dist_expected_mean_depol_input,dist_expected_mean_depol_input, 'k--')         # this gives the linear line

        plt.margins(0.1)
        plt.legend(loc=2)
        plt.title("Summary plot of mean input-output curve for distal locations")
        plt.xlabel("expected EPSP (mV)")
        plt.ylabel("measured EPSP (mV)")

        #fig = plt.gcf()
        #fig.set_size_inches(12, 15)
        if self.save_all:
            plt.savefig(self.path_figs + 'summary_input_output_curve_sync' + '.pdf', dpi=600, bbox_inches='tight')

        plt.figure()

        plt.subplot(2,1,1)
        plt.title('Synchronous inputs')
        #Derivative plot
        try:
            colormap = plt.cm.spectral      #http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
        except:
            colormap = plt.cm.nipy_spectral
        plt.gca().set_prop_cycle(plt.cycler('color', colormap(numpy.linspace(0, 0.9, len(sep_results)))))
        for i in range (0, len(sep_results)):

            plt.plot(num,sep_max_dV_dt[i], '-o', label='dend '+str(dend_loc000[i][0])+ ' loc ' +str(dend_loc000[i][1]))
            plt.xlabel("# of inputs")
            plt.ylabel("dV/dt (V/s)")
            plt.legend(loc=2, prop={'size':10})

        plt.subplot(2,1,2)

        plt.errorbar(num, mean_peak_deriv_input, yerr=SD_peak_deriv_input, linestyle='-', marker='o', color='red', label='SD')
        plt.errorbar(num, mean_peak_deriv_input, yerr=SEM_peak_deriv_input, linestyle='-', marker='o', color='blue', label='SEM')
        plt.margins(0.1)
        plt.legend(loc=2)
        plt.title("Summary plot of mean peak dV/dt amplitude")
        plt.xlabel("# of inputs")
        plt.ylabel("dV/dt (V/s)")

        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        if self.save_all:
            plt.savefig(self.path_figs + 'peak_derivative_plots_sync' + '.pdf', dpi=600, bbox_inches='tight')

        peak_derivatives = {'mean_peak_deriv': mean_peak_deriv_input,
                            'sd_peak_deriv': SD_peak_deriv_input,
                            'sem_peak_deriv': SEM_peak_deriv_input}

        #VALUES PLOT

        columns = 2
        width_ratios=[1]*columns
        rows = 4
        height_ratios=[1]*rows
        plt.figure()
        #fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.95, bottom=0.04, left=0.07, right=0.97, hspace=0.5, wspace=0.3)
        fig.suptitle('Synchronous inputs')

        plt.subplot(gs[0])
        # plot of thresholds
        x =numpy.array([])
        labels = ['exp mean with SD']
        e = numpy.array([threshold_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([experimental_mean_threshold])
        for i in range (0, len(sep_results)+1):
            x=numpy.append(x, i+1)
        for i in range (0, len(sep_results)):
            labels.append(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1])+')')
            plt.plot(x[i+1], sep_threshold[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x, labels, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Threshold (mV)")

        plt.subplot(gs[1])

        # plot of proximal thresholds
        x_prox =numpy.array([])
        labels_prox = ['exp mean with SD']
        e = numpy.array([threshold_prox_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([threshold_prox])
        for i in range(len(dend_loc000)):
            if dend_loc000[i][2] == 'prox':
                labels_prox.append(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1]) + ')')
        for i in range(0, len(prox_thresholds)+1):
            x_prox=numpy.append(x_prox, i+1)
        for i in range(0, len(prox_thresholds)):
            plt.plot(x_prox[i+1], prox_thresholds[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x_prox, labels_prox, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Proximal threshold (mV)")

        plt.subplot(gs[2])

        # plot of distal thresholds
        x_dist =numpy.array([])
        labels_dist = ['exp mean with SD']
        e = numpy.array([threshold_dist_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([threshold_dist])
        for i in range (len(dend_loc000)):
            if dend_loc000[i][2] == 'dist':
                labels_dist.append(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1])+')')
        for i in range (0, len(dist_thresholds)+1):
            x_dist=numpy.append(x_dist, i+1)
        for i in range (0, len(dist_thresholds)):
            plt.plot(x_dist[i+1], dist_thresholds[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x_dist, labels_dist, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Distal threshold (mV)")

        plt.subplot(gs[3])

        # plot of peak derivateives at threshold
        e = numpy.array([deriv_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([exp_mean_peak_deriv])
        for i in range (0, len(sep_results)):
            plt.plot(x[i+1], peak_dV_dt_at_threshold[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x, labels, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("peak derivative at threshold (V/s)")

        plt.subplot(gs[4])

        # plot of degree of nonlinearity at threshold

        e = numpy.array([nonlin_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([exp_mean_nonlin])
        for i in range (0, len(sep_results)):
            plt.plot(x[i+1], nonlin[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x, labels, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("degree of nonlinearity (%)")

        plt.subplot(gs[5])

        # plot of suprathreshold degree of nonlinearity

        e = numpy.array([suprath_nonlin_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([suprath_exp_mean_nonlin])
        for i in range (0, len(sep_results)):
            plt.plot(x[i+1], suprath_nonlin[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x, labels, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("suprath. degree of nonlinearity (%)")

        plt.subplot(gs[6])

        # plot of amplitude at threshold
        e = numpy.array([amp_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([exp_mean_amp])
        for i in range (0, len(sep_results)):
            plt.plot(x[i+1], amp_at_threshold[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x, labels, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Amplitude at threshold (mV)")


        plt.subplot(gs[7])

        # plot of time to peak at threshold
        e = numpy.array([exp_mean_time_to_peak_SD])
        x2 =numpy.array([1])
        y2 = numpy.array([exp_mean_time_to_peak])
        for i in range (0, len(sep_results)):
            plt.plot(x[i+1], time_to_peak_at_threshold[i], 'o')

        plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
        plt.xticks(x, labels, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("time to peak at threshold (ms)")

        fig = plt.gcf()
        fig.set_size_inches(14, 18)
        if self.save_all:
            plt.savefig(self.path_figs + 'values_sync' + '.pdf', dpi=600, bbox_inches='tight')



        # ERROR PLOTS


        columns = 2
        width_ratios=[1]*columns
        rows = 4
        height_ratios=[1]*rows
        plt.figure()
        #fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.95, bottom=0.04, left=0.07, right=0.97, hspace=0.5, wspace=0.3)
        fig.suptitle(' Errors in units of the experimental SD of the feature (synchronous inputs)')

        plt.subplot(gs[0])

        #threshold error plot
        x_error =numpy.array([])
        labels_error = []
        for i in range (0, len(sep_results)):
            labels_error.append(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1])+')')
            x_error=numpy.append(x_error, i+1)

            plt.plot(x_error[i], threshold_errors[i], 'o')
        plt.xticks(x_error, labels_error, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Threshold error")

        plt.subplot(gs[1])

        # proximal threshold error plot

        x_prox_err =numpy.array([])
        labels_prox_err = []
        for i in range (len(dend_loc000)):
            if dend_loc000[i][2] == 'prox':
                labels_prox_err.append(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1])+')')
        for i in range (0, len(prox_threshold_errors)):
            x_prox_err=numpy.append(x, i+1)

            plt.plot(x_prox_err[i], prox_threshold_errors[i], 'o')
        plt.xticks(x_prox_err, labels_prox_err, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Proximal threshold error")

        plt.subplot(gs[2])

    # distal threshold error plot

        x_dist_err =numpy.array([])
        labels_dist_err = []
        for i in range(len(dend_loc000)):
            if dend_loc000[i][2] == 'dist':
                labels_dist_err.append(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1])+')')
        for i in range(0, len(dist_threshold_errors)):
            x_dist_err=numpy.append(x_dist_err, i+1)

            plt.plot(x_dist_err[i], dist_threshold_errors[i], 'o')
        plt.xticks(x_dist_err, labels_dist_err, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Distal threshold error")

        plt.subplot(gs[3])

    #peak deriv error plot

        for i in range (0, len(sep_results)):

            plt.plot(x_error[i], peak_deriv_errors[i], 'o')
        plt.xticks(x_error, labels_error, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Peak derivative error")

        plt.subplot(gs[4])

    #  degree of nonlin. error plot

        for i in range (0, len(sep_results)):
            plt.plot(x_error[i], nonlin_errors[i], 'o')
        plt.xticks(x_error, labels_error, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Degree of nonlinearity error")


        plt.subplot(gs[5])

    # suprathreshold degree of nonlin. error plot

        for i in range (0, len(sep_results)):
            plt.plot(x_error[i], suprath_nonlin_errors[i], 'o')
        plt.xticks(x_error, labels_error, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Suprath. degree of nonlinearity error")

        plt.subplot(gs[6])

    # amplitude error plot

        for i in range (0, len(sep_results)):
            plt.plot(x_error[i], amplitude_errors[i], 'o')
        plt.xticks(x_error, labels_error, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Amplitude error")

        plt.subplot(gs[7])

    # time to peak error plot

        for i in range (0, len(sep_results)):
            plt.plot(x_error[i], time_to_peak_errors[i], 'o')
        plt.xticks(x_error, labels_error, rotation=20)
        #plt.tick_params(labelsize=10)
        plt.margins(0.1)
        plt.ylabel("Time to peak error")

        fig = plt.gcf()
        fig.set_size_inches(14, 18)
        if self.save_all:
            plt.savefig(self.path_figs + 'errors_sync' + '.pdf', dpi=600, bbox_inches='tight')

    # mean values plot

        import matplotlib.gridspec as gridspec
        gs = gridspec.GridSpec(1, 4,width_ratios=[4,1,2,1])
        fig3, axes3 = plt.subplots(nrows=1, ncols=4)
        fig3.tight_layout()
        fig3.suptitle('Synchronous inputs', fontsize=15)
        #plt.subplot(1,4,1)
        plt.subplot(gs[0])

        e_values = numpy.array([threshold_SD, sd_sep_threshold, threshold_prox_SD, sd_prox_thresholds, threshold_dist_SD, sd_dist_thresholds, amp_SD, sd_amp_at_threshold])
        x_values =numpy.array([1,2,4,5,7,8,10,11])
        y_values = numpy.array([experimental_mean_threshold, mean_sep_threshold, threshold_prox, mean_prox_thresholds, threshold_dist, mean_dist_thresholds, exp_mean_amp, mean_amp_at_threshold])
        labels_values=['','threshold', '', 'proximal threshold', '', 'distal threshold','', 'amplitude at th.']
        colors = ['r', 'b', 'r', 'b', 'r', 'b', 'r', 'b']
        for i in range(len(y_values)):
            if i==0:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='experimental mean + SD')
            elif i == 1:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='model mean + SD')
            else:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i])
        plt.xticks(x_values, labels_values, rotation=20)
        plt.tick_params(labelsize=15)
        plt.margins(0.1)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0, prop={'size':12})
        plt.ylabel("Voltage (mV)", fontsize=15)

        #plt.subplot(1,4,2)
        plt.subplot(gs[1])

        e_values = numpy.array([deriv_SD, sd_peak_dV_dt_at_threshold])
        x_values =numpy.array([1,2])
        y_values = numpy.array([exp_mean_peak_deriv, mean_peak_dV_dt_at_threshold])
        labels_values=['', 'peak dV/dt at th.']
        colors=['r', 'b']
        for i in range(len(y_values)):
            if i==0:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='experimental mean')
            elif i == 1:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='model mean')
        plt.xticks(x_values, labels_values, rotation=10)
        plt.tick_params(labelsize=15)
        plt.margins(0.3)
        #plt.legend(loc=1, prop={'size':10})
        plt.ylabel("peak dV/dt(V/s)", fontsize=15)

        #plt.subplot(1,4,3)
        plt.subplot(gs[2])

        e_values = numpy.array([nonlin_SD, sd_nonlin, suprath_nonlin_SD, suprath_sd_nonlin])
        x_values =numpy.array([1,2,4,5])
        y_values = numpy.array([exp_mean_nonlin, mean_nonlin, suprath_exp_mean_nonlin, suprath_mean_nonlin])
        labels_values=['', 'degree of nonlinearity at th.', '', 'suprath. degree of nonlinearity']
        colors=['r', 'b', 'r', 'b']
        for i in range(len(y_values)):
            if i==0:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='experimental mean')
            elif i == 1:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='model mean')
            else:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i])
        plt.xticks(x_values, labels_values, rotation=10)
        plt.tick_params(labelsize=15)
        plt.margins(0.1)
        #plt.legend(loc=2, prop={'size':10})
        plt.ylabel("degree of nonlinearity(%)", fontsize=15)

        #plt.subplot(1,4,4)
        plt.subplot(gs[3])

        e_values = numpy.array([exp_mean_time_to_peak_SD, sd_time_to_peak_at_threshold])
        x_values =numpy.array([1,2])
        y_values = numpy.array([exp_mean_time_to_peak, mean_time_to_peak_at_threshold])
        labels_values=['', 'time to peak at th.']
        colors=['r', 'b']
        for i in range(len(y_values)):
            if i==0:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='experimental mean')
            elif i == 1:
                plt.errorbar(x_values[i], y_values[i], e_values[i], linestyle='None', marker='o', color=colors[i], label='model mean')
        plt.xticks(x_values, labels_values, rotation=10)
        plt.tick_params(labelsize=15)
        plt.margins(0.3)
        #plt.legend(loc=1, prop={'size':10})
        plt.ylabel("time to peak (ms)", fontsize=15)

        fig = plt.gcf()
        fig.set_size_inches(22, 18)
        if self.save_all:
            plt.savefig(self.path_figs + 'mean_values_sync' + '.pdf', dpi=600, bbox_inches='tight')

        """ ZScore:obliqueIntegration is used instead"""
    # mean errors plot
        '''
        plt.figure(figsize = (210/25.4, 210/25.4))
        plt.title('Synchronous inputs')
        e_errors = numpy.array([sd_threshold_errors, sd_prox_threshold_errors, sd_dist_threshold_errors, sd_peak_deriv_errors, sd_nonlin_errors, suprath_sd_nonlin_errors, sd_amplitude_errors, sd_time_to_peak_errors])
        x_errors =numpy.array([1,2,3,4,5,6,7,8])
        y_errors = numpy.array([mean_threshold_errors,  mean_prox_threshold_errors, mean_dist_threshold_errors, mean_peak_deriv_errors, mean_nonlin_errors, suprath_mean_nonlin_errors, mean_amplitude_errors, mean_time_to_peak_errors ])
        labels_errors=['mean threshold error', 'mean proximal threshold error', 'mean distal threshold error', 'mean peak dV/dt at th. error', 'mean degree of nonlinearity at th.error', 'mean suprath. degree of nonlinearity error','mean amplitude at th. error', 'mean time to peak at th. error']
        plt.errorbar(x_errors, y_errors, e_errors, linestyle='None', marker='o')
        plt.xticks(x_errors, labels_errors, rotation=20)
        #plt.tick_params(labelsize=15)
        plt.margins(0.1)
        plt.ylabel("model mean errors in unit of the experimental SD (with SD)")

        fig = plt.gcf()
        #fig.set_size_inches(16, 18)
        plt.savefig(self.path_figs + 'mean_errors_sync' + '.pdf', dpi=600, bbox_inches='tight')


        errors = {'threshold_error' :{'mean': mean_threshold_errors,
                                        'sd':sd_threshold_errors},
                    'proximal_threshold_error' : {'mean':mean_prox_threshold_errors,
                                                  'sd':sd_prox_threshold_errors},
                    'distal_threshold_error' : {'mean':mean_dist_threshold_errors,
                                                'sd':sd_dist_threshold_errors},
                    'peak_dV/dt_at_th_error' : {'mean':mean_peak_deriv_errors,
                                                'sd':sd_peak_deriv_errors},
                    'degree_of_nonlinearity_at_th_error' : {'mean':mean_nonlin_errors,
                                                            'sd':sd_nonlin_errors},
                     'suprath_degree_of_nonlinearity_error' : {'mean':suprath_mean_nonlin_errors,
                                                               'sd':suprath_sd_nonlin_errors},
                     'amplitude_at_th_error' : {'mean':mean_amplitude_errors,
                                                'sd':sd_amplitude_errors},
                     'time_to_peak_at_th_error' : {'mean':mean_time_to_peak_errors,
                                                   'sd':sd_time_to_peak_errors}
                     }
        '''

        model_means = [mean_sep_threshold, mean_prox_thresholds, mean_dist_thresholds, mean_peak_dV_dt_at_threshold , mean_nonlin, suprath_mean_nonlin,  mean_amp_at_threshold, mean_time_to_peak_at_threshold]
        model_SDs = [sd_sep_threshold, sd_prox_thresholds, sd_dist_thresholds, sd_peak_dV_dt_at_threshold , sd_nonlin, suprath_sd_nonlin, sd_amp_at_threshold, sd_time_to_peak_at_threshold]
        model_N= len(sep_results)
        model_prox_N = prox_n
        model_dist_N = dist_n

        return model_means, model_SDs, model_N, model_prox_N, model_dist_N, EPSPs, peak_derivatives


    def calcs_plots_async(self, model, results, dend_loc000, dend_loc_num_weight):

        async_nonlin=self.observation['mean_async_nonlin']
        async_nonlin_SEM=self.observation['async_nonlin_sem']
        async_nonlin_SD=self.observation['async_nonlin_std']

        stop=len(dend_loc_num_weight)+1
        sep=numpy.arange(0,stop,self.max_num_syn+1)
        sep_results=[]

        num = numpy.arange(0,self.max_num_syn+1)

        for i in range (0,len(dend_loc000)):
            sep_results.append(results[sep[i]:sep[i+1]])             # a list that contains the results of the 10 locations seperately (in lists)

        # sep_results[0]-- the first location
        # sep_results[0][5] -- the first location at 5 input
        # sep_results[0][1][0] -- the first location at 1 input, SOMA
        # sep_results[0][1][1] -- the first location at 1 input, dendrite
        # sep_results[0][1][1][0] -- just needed


        frames = len(dend_loc000)
        columns = 2
        width_ratios=[1]*columns
        rows = int(numpy.ceil(frames/float(columns)))
        height_ratios=[1]*rows

        fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.92, bottom=0.04, left=0.07, right=0.97, hspace=0.75, wspace=0.3)

        fig.suptitle('Asynchronous inputs (red: dendritic trace, black: somatic trace)')
        for i in range (0,len(dend_loc000)):
            plt.subplot(gs[i])
            for j, number in enumerate(num):
                plt.plot(sep_results[i][j][0][0]['T'],sep_results[i][j][0][0]['V'], 'k')       # somatic traces
                plt.plot(sep_results[i][j][1][0]['T'],sep_results[i][j][1][0]['V'], 'r')        # dendritic traces
            plt.title(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1]) + ')')

            plt.xlabel("time (ms)")
            plt.ylabel("Voltage (mV)")
            plt.xlim(sep_results[i][0][0][0]['stim_start'][0] - 10, sep_results[i][0][0][0]['stim_start'][0]+250)

            #plt.tick_params(labelsize=20)
        if self.save_all:
            plt.savefig(self.path_figs + 'traces_async' + '.pdf', dpi=600, bbox_inches='tight')

        frames = len(dend_loc000)
        columns = 2
        width_ratios=[1]*columns
        rows = int(numpy.ceil(frames/float(columns)))
        height_ratios=[1]*rows
        fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.92, bottom=0.04, left=0.07, right=0.97, hspace=0.75, wspace=0.3)

        fig.suptitle('Asynchronous inputs')
        for i in range (0,len(dend_loc000)):
            plt.subplot(gs[i])
            for j, number in enumerate(num):
                plt.plot(sep_results[i][j][0][0]['T'],sep_results[i][j][0][0]['V'], 'k')       # somatic traces
            plt.title(str(dend_loc000[i][0])+ '(' +str(dend_loc000[i][1]) +')')

            plt.xlabel("time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            plt.xlim(sep_results[i][0][0][0]['stim_start'][0] - 10, sep_results[i][0][0][0]['stim_start'][0]+250)

            #plt.tick_params(labelsize=20)
        if self.save_all:
            plt.savefig(self.path_figs + 'somatic_traces_async' + '.pdf', dpi=600, bbox_inches='tight')

        soma_depol=numpy.array([])
        soma_depols=[]
        sep_soma_depols=[]
        #dV_dt=[]
        #sep_dV_dt=[]
        soma_max_depols=numpy.array([])
        soma_expected=numpy.array([])
        sep_soma_max_depols=[]
        sep_soma_expected=[]
        max_dV_dt=numpy.array([])
        sep_max_dV_dt=[]
        nonlin=numpy.array([])
        sep_nonlin=[]

        nonlins_at_th=numpy.array([])


        for i in range (0, len(sep_results)):
            for j in range (0,self.max_num_syn+1):

        # calculating somatic depolarization and first derivative

                if not numpy.array_equal(sep_results[i][j][0][0]['T'], sep_results[i][0][0][0]['T']):    #if the  time vectors are not equal, the traces are resampled with fixed time step
                    dt = 0.025

                    time_vector = numpy.arange(sep_results[i][j][0][0]['T'][0], sep_results[i][j][0][0]['T'][-1], dt)  #from the first to the last element of the original time vector

                    interp_trace_input = numpy.interp(time_vector, sep_results[i][j][0][0]['T'], sep_results[i][j][0][0]['V'])
                    interp_trace_no_input = numpy.interp(time_vector, sep_results[i][0][0][0]['T'], sep_results[i][0][0][0]['V'])

                    soma_depol=interp_trace_input - interp_trace_no_input

                    dt=numpy.diff(time_vector)
                    dV=numpy.diff(interp_trace_input)

                else:

                    soma_depol=sep_results[i][j][0][0]['V'] - sep_results[i][0][0][0]['V']
                    dt=numpy.diff(sep_results[i][j][0][0]['T'] )
                    dV=numpy.diff(sep_results[i][j][0][0]['V'] )

                soma_depols.append(soma_depol)

                soma_max_depols=numpy.append(soma_max_depols,numpy.amax(soma_depol))

                deriv=dV/dt

                max_dV_dt=numpy.append(max_dV_dt, numpy.amax(deriv))

                if j==0:
                    soma_expected=numpy.append(soma_expected,0)
                    nonlin=numpy.append(nonlin,100)
                else:
                    soma_expected=numpy.append(soma_expected,soma_max_depols[1]*j)
                    nonlin=numpy.append(nonlin, soma_max_depols[j]/ soma_expected[j]*100)

            sep_soma_depols.append(soma_depols)
            #sep_dV_dt.append(dV_dt)
            sep_soma_max_depols.append(soma_max_depols)
            sep_soma_expected.append(soma_expected)
            sep_max_dV_dt.append(max_dV_dt)
            sep_nonlin.append(nonlin)
            nonlins_at_th=numpy.append(nonlins_at_th, nonlin[5])      #degree of nonlin at 5 inputs, that is the threshold in the synchronous case

            soma_depols=[]
            #dV_dt=[]
            soma_max_depols=numpy.array([])
            soma_expected=numpy.array([])
            max_dV_dt=numpy.array([])
            nonlin=numpy.array([])

        expected_depol_input=numpy.array([])
        expected_mean_depol_input=[]
        depol_input=numpy.array([])
        mean_depol_input=[]
        SD_depol_input=[]
        SEM_depol_input=[]

        peak_deriv_input=numpy.array([])
        mean_peak_deriv_input=[]
        SD_peak_deriv_input=[]
        SEM_peak_deriv_input=[]
        n=len(sep_soma_max_depols)


        for j in range (0, self.max_num_syn+1):
            for i in range (0, len(sep_soma_max_depols)):
                depol_input=numpy.append(depol_input,sep_soma_max_depols[i][j])
                expected_depol_input=numpy.append(expected_depol_input,sep_soma_expected[i][j])
                peak_deriv_input=numpy.append(peak_deriv_input,sep_max_dV_dt[i][j])
            mean_depol_input.append(numpy.mean(depol_input))
            expected_mean_depol_input.append(numpy.mean(expected_depol_input))
            mean_peak_deriv_input.append(numpy.mean(peak_deriv_input))
            SD_depol_input.append(numpy.std(depol_input))
            SEM_depol_input.append(numpy.std(depol_input)/math.sqrt(n))
            SD_peak_deriv_input.append(numpy.std(peak_deriv_input))
            SEM_peak_deriv_input.append(numpy.std(peak_deriv_input)/math.sqrt(n))
            depol_input=numpy.array([])
            expected_depol_input=numpy.array([])
            peak_deriv_input=numpy.array([])

        mean_nonlin_at_th=numpy.mean(nonlins_at_th)
        SD_nonlin_at_th=numpy.std(nonlins_at_th)
        SEM_nonlin_at_th=SD_nonlin_at_th/math.sqrt(n)

        plt.figure()
        plt.subplot(2,1,1)
        plt.title('Asynchronous inputs')
        # Expected EPSP - Measured EPSP plot
        try:
            colormap = plt.cm.spectral      #http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
        except:
            colormap = plt.cm.nipy_spectral
        plt.gca().set_prop_cycle(plt.cycler('color', colormap(numpy.linspace(0, 0.9, len(sep_results)))))
        for i in range (0, len(sep_results)):

            plt.plot(sep_soma_expected[i],sep_soma_max_depols[i], '-o', label='dend '+str(dend_loc000[i][0])+ ' loc ' +str(dend_loc000[i][1]))
            plt.plot(sep_soma_expected[i],sep_soma_expected[i], 'k--')         # this gives the linear line
        plt.xlabel("expected EPSP (mV)")
        plt.ylabel("measured EPSP (mV)")
        plt.legend(loc=2, prop={'size':10})

        plt.subplot(2,1,2)

        plt.errorbar(expected_mean_depol_input, mean_depol_input, yerr=SD_depol_input, linestyle='-', marker='o', color='red', label='SD')
        plt.errorbar(expected_mean_depol_input, mean_depol_input, yerr=SEM_depol_input, linestyle='-', marker='o', color='blue', label='SEM')
        plt.plot(expected_mean_depol_input,expected_mean_depol_input, 'k--')         # this gives the linear line
        plt.margins(0.1)
        plt.legend(loc=2)
        plt.title("Summary plot of mean input-output curve")
        plt.xlabel("expected EPSP (mV)")
        plt.ylabel("measured EPSP (mV)")

        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        if self.save_all:
            plt.savefig(self.path_figs + 'input_output_curves_async' + '.pdf', dpi=600, bbox_inches='tight')

        EPSPs ={'expected_mean' : expected_mean_depol_input,
                'measured_mean' : mean_depol_input,
                'sd' : SD_depol_input,
                'sem' : SEM_depol_input}


        plt.figure()

        plt.subplot(2,1,1)
        plt.title('Asynchronous inputs')

        #Derivative plot
        try:
            colormap = plt.cm.spectral      #http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
        except:
            colormap = plt.cm.nipy_spectral
        plt.gca().set_prop_cycle(plt.cycler('color', colormap(numpy.linspace(0, 0.9, len(sep_results)))))
        for i in range (0, len(sep_results)):

            plt.plot(num,sep_max_dV_dt[i], '-o', label='dend '+str(dend_loc000[i][0])+ ' loc ' +str(dend_loc000[i][1]))
            plt.xlabel("# of inputs")
            plt.ylabel("dV/dt (V/s)")
            plt.legend(loc=2, prop={'size':10})

        plt.subplot(2,1,2)

        plt.errorbar(num, mean_peak_deriv_input, yerr=SD_peak_deriv_input, linestyle='-', marker='o', color='red', label='SD')
        plt.errorbar(num, mean_peak_deriv_input, yerr=SEM_peak_deriv_input, linestyle='-', marker='o', color='blue', label='SEM')
        plt.margins(0.01)
        plt.legend(loc=2)
        plt.title("Summary plot of mean peak dV/dt amplitude")
        plt.xlabel("# of inputs")
        plt.ylabel("dV/dt (V/s)")

        fig = plt.gcf()
        fig.set_size_inches(10, 10)
        if self.save_all:
            plt.savefig(self.path_figs + 'peak_derivative_plots_async' + '.pdf', dpi=600, bbox_inches='tight')

        peak_derivatives = {'mean_peak_deriv': mean_peak_deriv_input,
                            'sd_peak_deriv': SD_peak_deriv_input,
                            'sem_peak_deriv': SEM_peak_deriv_input}

        fig0, axes0 = plt.subplots(nrows=1, ncols=1)
        fig0.tight_layout()
        fig0.suptitle('Asynchronous inputs')
        for j in range (0,len(dend_loc000)):
            plt.subplot(numpy.ceil(len(dend_loc000)/2.0),2,j+1)
            x =numpy.array([])
            labels = ['exp. mean\n with SD']
            e = numpy.array([async_nonlin_SD])
            x2 =numpy.array([1])
            y2 = numpy.array([async_nonlin])
            for i in range (0, len(sep_nonlin[0])+1):
                x=numpy.append(x, i+1)
                labels.append(str(i)+ ' inputs')
            for i in range (0, len(sep_nonlin[j])):
                plt.plot(x[i+1], sep_nonlin[j][i], 'o')

            plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
            plt.xticks(x, labels, rotation=40)
            #plt.tick_params(labelsize=15)
            plt.margins(0.1)
            plt.ylabel("Degree of nonlinearity (%)")
            plt.title(str(dend_loc000[j][0])+ '(' +str(dend_loc000[j][1])+')')

        fig = plt.gcf()
        fig.set_size_inches(12, 12)
        if self.save_all:
            plt.savefig(self.path_figs + 'nonlin_values_async' + '.pdf', dpi=600, bbox_inches='tight')

        async_nonlin_errors=[]
        asynch_nonlin_error_at_th=numpy.array([])

        for i in range (0, len(sep_nonlin)):
            async_nonlin_err = numpy.array([abs(async_nonlin- async_nonlin_err)/async_nonlin_SD  for async_nonlin_err in sep_nonlin[i]])     # does the same calculation on every element of a list  #x = [1,3,4,5,6,7,8] t = [ t**2 for t in x ]
            async_nonlin_errors.append(async_nonlin_err)
            asynch_nonlin_error_at_th=numpy.append(asynch_nonlin_error_at_th, async_nonlin_err[4])

        mean_nonlin_error_at_th=numpy.mean(asynch_nonlin_error_at_th)
        SD_nonlin_error_at_th=numpy.std(asynch_nonlin_error_at_th)
        SEM_nonlin_error_at_th=SD_nonlin_error_at_th/math.sqrt(n)

        error ={'asynch_degree_of_nonlin_at_th_error': {'mean':mean_nonlin_error_at_th,
                                                        'sd':SD_nonlin_error_at_th}
                }

        fig0, axes0 = plt.subplots(nrows=1, ncols=1)
        fig0.tight_layout()
        fig0.suptitle('Asynchronous inputs')
        for j in range (0,len(dend_loc000)):
            plt.subplot(numpy.ceil(len(dend_loc000)/2.0),2,j+1)
            for i in range (0, len(async_nonlin_errors[j])):
                plt.plot(x[i], async_nonlin_errors[j][i], 'o')

            plt.xticks(x, labels[1:-1], rotation=20)
            #plt.tick_params(labelsize=15)
            plt.margins(0.1)
            plt.ylabel("Degree of nonlin. error (%)")
            plt.title(str(dend_loc000[j][0])+ '(' +str(dend_loc000[j][1])+')')
        fig = plt.gcf()
        fig.set_size_inches(12, 12)
        if self.save_all:
            plt.savefig(self.path_figs + 'nonlin_errors_async' + '.pdf', dpi=600, bbox_inches='tight')

        return mean_nonlin_at_th, SD_nonlin_at_th, EPSPs, peak_derivatives


    def validate_observation(self, observation):

        #self.add_std_to_observation(observation)

        try:
            assert type(observation['mean_threshold']) is Quantity
            assert type(observation['threshold_std']) is Quantity
            assert type(observation['mean_prox_threshold']) is Quantity
            assert type(observation['prox_threshold_std']) is Quantity
            assert type(observation['mean_dist_threshold']) is Quantity
            assert type(observation['dist_threshold_std']) is Quantity
            assert type(observation['mean_peak_deriv']) is Quantity
            assert type(observation['peak_deriv_std']) is Quantity
            assert type(observation['mean_amp_at_th']) is Quantity
            assert type(observation['amp_at_th_std']) is Quantity
            assert type(observation['mean_time_to_peak']) is Quantity
            assert type(observation['time_to_peak_std']) is Quantity

        except Exception as e:
            raise ObservationError(("Observation must be of the form "
                                    "{'mean':float*mV,'std':float*mV}"))


    def generate_prediction(self, model, verbose=False):
        """Implementation of sciunit.Test.generate_prediction."""

        efel.reset()

        if self.base_directory:
            self.path_results = self.base_directory + 'results/' + 'oblique_integration/' + model.name + '/'
        else:
            self.path_results = model.base_directory + 'results/' + 'oblique_integration/'

        try:
            if not os.path.exists(self.path_results):
                os.makedirs(self.path_results)
        except OSError as e:
            if e.errno != 17:
                raise
            pass


        model.find_obliques_multiproc(self.trunk_origin)

        print('Dendrites and locations to be tested: ', model.dend_loc)

        traces = []

        global model_name_oblique
        model_name_oblique = model.name

        if not model.AMPA_name:
            print('')
            print('The built in Exp2Syn is used as the AMPA component. Tau1 =', model.AMPA_tau1, ',Tau2 =', model.AMPA_tau2 , '.')
            print('')
        if not model.NMDA_name:
            print('')
            print('The default NMDA model of HippoUnit is used with Jahr, Stevens voltage dependence.')
            print('')

        #pool0 = multiprocessing.pool.ThreadPool(self.npool)    # multiprocessing.pool.ThreadPool is used because a nested multiprocessing is used in the function called here (to keep every NEURON related task in independent processes)
        pool0 = NonDaemonPool(self.npool, maxtasksperchild = 1)

        print("Adjusting synaptic weights on all the locations ...")

        binsearch_ = functools.partial(self.binsearch, model)
        results0 = pool0.map(binsearch_, model.dend_loc, chunksize=1)  #model.dend_loc[0:2] - no need for info if it is distal or proximal
        #results0 = result0.get()

        pool0.terminate()
        pool0.join()
        del pool0


        num = numpy.arange(0,self.max_num_syn+1)
        dend_loc_num_weight=[]

        #dend0=[]
        #dend_loc00=[]
        indices_to_delete = []
        dend_loc000=list(model.dend_loc) #dend_loc000 will not contain the dendrites in which spike causes somatic AP or generates no spike or generates spike even for very small stimulus

        results00=list(results0)    #results00 will not contain the synaptic weights for dendrites in which spike causes somatic AP or generates no spike or generates spike even for very small stimulus

        for i in range(0, len(model.dend_loc)):

            if results0[i][0]==None:

                print('The dendritic spike at dendrite ' + str(model.dend_loc[i][0]) +'('+str(model.dend_loc[i][1])+')' + ' generated somatic AP - this location is not used in the test')

                indices_to_delete.append(i)

            elif results0[i][0]=='no spike':
                print('No dendritic spike could be generated at dendrite ' + str(model.dend_loc[i][0]) +'('+str(model.dend_loc[i][1])+')' +  ' - this location is not used in the test')

                indices_to_delete.append(i)

            elif results0[i][0]=='always spike':
                print('Dendrite ' + str(model.dend_loc[i][0]) +'('+str(model.dend_loc[i][1])+')' + ' generates dendritic spike even to smaller number of inputs - this location is not used in the test')

                indices_to_delete.append(i)

        for k in sorted(indices_to_delete, reverse=True):  #deleted in reverse order so that subsequent indices remains ok

            del dend_loc000[k]
            del results00[k]
        #print model.dend_loc
        #print dend_loc000

        if len(dend_loc000) > 0:        # if none of the dendrites was able to generate dendritic spike the list is empty
            for i in range(0, len(dend_loc000)):

                for j in num:

                    e=list(dend_loc000[i])
                    e.append(j)
                    e.append(results00[i][1])
                    dend_loc_num_weight.append(e)        #calculates, and adds the synaptic weights needed to a list
            #print dend_loc_num_weight
            interval_sync=0.1

            pool = multiprocessing.Pool(self.npool, maxtasksperchild=1)
            run_synapse_ = functools.partial(self.run_synapse, model, interval=interval_sync)
            results = pool.map(run_synapse_, dend_loc_num_weight, chunksize=1)
            #results = result.get()

            pool.terminate()
            pool.join()
            del pool


            pool1 = multiprocessing.Pool(self.npool, maxtasksperchild=1)

            interval_async=2

            run_synapse_ = functools.partial(self.run_synapse, model, interval=interval_async)
            results_async = pool1.map(run_synapse_,dend_loc_num_weight, chunksize=1)    # ordered results
            #results_async = result_async.get()

            pool1.terminate()
            pool1.join()
            del pool1

            plt.close('all') #needed to avoid overlapping of saved images when the test is run on multiple models in a for loop

            model_means, model_SDs, model_N, model_prox_N, model_dist_N, EPSPs_sync, sync_peak_derivatives = self.calcs_plots(model, results, dend_loc000, dend_loc_num_weight)

            mean_nonlin_at_th_asynch, SD_nonlin_at_th_asynch, EPSPs_async, async_peak_derivatives = self.calcs_plots_async(model, results_async, dend_loc000, dend_loc_num_weight)

            #errors = dict(sync_errors)
            #errors.update(async_errors)

            ''' printing to logFile'''
            filepath = self.path_results + self.test_log_filename
            self.logFile = open(filepath, 'w') # if it is opened before multiprocessing, the multiporeccing won't work under python3

            self.logFile.write('Dendrites and locations to be tested:\n'+ str(model.dend_loc)+'\n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

            if not model.AMPA_name:
                self.logFile.write('The built in Exp2Syn is used as the AMPA component. Tau1 = ' + str(model.AMPA_tau1) + ', Tau2 = ' + str(model.AMPA_tau2) + '.\n')
                self.logFile.write("---------------------------------------------------------------------------------------------------\n")
            if not model.NMDA_name:
                self.logFile.write('The default NMDA model of HippoUnit is used with Jahr, Stevens voltage dependence.\n')
                self.logFile.write("---------------------------------------------------------------------------------------------------\n")

            for i in range(0, len(model.dend_loc)):

                if results0[i][0]==None:
                    self.logFile.write('The dendritic spike at dendrite ' + str(model.dend_loc[i][0]) +'('+str(model.dend_loc[i][1])+')' + ' generated somatic AP - this location is not used in the test\n')
                    self.logFile.write("---------------------------------------------------------------------------------------------------\n")


                elif results0[i][0]=='no spike':
                    self.logFile.write('No dendritic spike could be generated at dendrite ' + str(model.dend_loc[i][0]) +'('+str(model.dend_loc[i][1])+')' +  ' - this location is not used in the test\n')
                    self.logFile.write("---------------------------------------------------------------------------------------------------\n")

                elif results0[i][0]=='always spike':
                    self.logFile.write('Dendrite ' + str(model.dend_loc[i][0]) +'('+str(model.dend_loc[i][1])+')' + ' generates dendritic spike even to smaller number of inputs - this location is not used in the test\n')
                    self.logFile.write("---------------------------------------------------------------------------------------------------\n")



            prediction = {'model_mean_threshold':model_means[0], 'model_threshold_std': model_SDs[0],
                            'model_mean_prox_threshold':model_means[1], 'model_prox_threshold_std': model_SDs[1],
                            'model_mean_dist_threshold':model_means[2], 'model_dist_threshold_std': model_SDs[2],
                            'model_mean_peak_deriv':model_means[3],'model_peak_deriv_std': model_SDs[3],
                            'model_mean_nonlin_at_th':model_means[4], 'model_nonlin_at_th_std': model_SDs[4],
                            'model_mean_nonlin_suprath':model_means[5], 'model_nonlin_suprath_std': model_SDs[5],
                            'model_mean_amp_at_th':model_means[6],'model_amp_at_th_std': model_SDs[6],
                            'model_mean_time_to_peak':model_means[7], 'model_time_to_peak_std': model_SDs[7],
                            'model_mean_async_nonlin':mean_nonlin_at_th_asynch, 'model_async_nonlin_std': SD_nonlin_at_th_asynch,
                            'model_n': model_N, 'model_prox_n' : model_prox_N, 'model_dist_n' : model_dist_N, }

        else:

            prediction = {'model_mean_threshold':float('NaN')*mV, 'model_threshold_std': float('NaN')*mV,
                            'model_mean_prox_threshold':float('NaN')*mV, 'model_prox_threshold_std': float('NaN')*mV,
                            'model_mean_dist_threshold':float('NaN')*mV, 'model_dist_threshold_std': float('NaN')*mV,
                            'model_mean_peak_deriv':float('NaN')*V / s,'model_peak_deriv_std': float('NaN')*V / s,
                            'model_mean_nonlin_at_th':float('NaN'), 'model_nonlin_at_th_std': float('NaN'),
                            'model_mean_nonlin_suprath':float('NaN'), 'model_nonlin_suprath_std': float('NaN'),
                            'model_mean_amp_at_th':float('NaN') *mV,'model_amp_at_th_std': float('NaN')*mV,
                            'model_mean_time_to_peak':float('NaN') *ms, 'model_time_to_peak_std': float('NaN')*ms,
                            'model_mean_async_nonlin':float('NaN'), 'model_async_nonlin_std': float('NaN'),
                            'model_n': 0.0 }

            EPSPs_sync = None
            sync_peak_derivatives = None
            EPSPs_async = None
            async_peak_derivatives = None

            filepath = self.path_results + self.test_log_filename
            self.logFile = open(filepath, 'w') # if it is opened before multiprocessing, the multiporeccing won't work under python3

            print("There isn\'t any appropriately behaving oblique dendrite to be tested")
            self.logFile.write("There isn\'t any appropriately behaving oblique dendrite to be tested\n")
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        prediction_json = dict(prediction)
        '''
        for key, value in prediction_json .iteritems():
            prediction_json[key] = str(value)
        '''

        for key, val in list(prediction.items()):
            val = str(val)
            quantity_parts = val.split("*")
            prediction_json[key] = " ".join(quantity_parts)

        features = {'threshold' :{'mean': prediction_json['model_mean_threshold'] ,
                                          'sd':prediction_json['model_threshold_std']},
                             'proximal_threshold' : {'mean':prediction_json['model_mean_prox_threshold'],
                                                     'sd':prediction_json['model_prox_threshold_std']},
                             'distal_threshold' : {'mean':prediction_json['model_mean_dist_threshold'],
                                                   'sd':prediction_json['model_dist_threshold_std']},
                             'peak_dV/dt_at_th' : {'mean':prediction_json['model_mean_peak_deriv'],
                                                   'sd':prediction_json['model_peak_deriv_std']},
                             'degree_of_nonlinearity_at_th' : {'mean':prediction_json['model_mean_nonlin_at_th'],
                                                               'sd':prediction_json['model_nonlin_at_th_std']},
                             'suprath_degree_of_nonlinearity' : {'mean':prediction_json['model_mean_nonlin_suprath'],
                                                                 'sd':prediction_json['model_nonlin_suprath_std']},
                             'amplitude_at_th' : {'mean':prediction_json['model_mean_amp_at_th'],
                                                  'sd':prediction_json['model_amp_at_th_std']},
                             'time_to_peak_at_th' : {'mean':prediction_json['model_mean_time_to_peak'],
                                                     'sd':prediction_json['model_time_to_peak_std']},
                             'asynch_degree_of_nonlin_at_th': {'mean':prediction_json['model_mean_async_nonlin'],
                                                               'sd':prediction_json['model_async_nonlin_std']},
                             'num_of_locations': prediction_json['model_n']
                             }

        '''
        file_name_json = self.path_results + 'oblique_model_features.json'

        json.dump(prediction_json, open(file_name_json, "wb"), indent=4)
        '''

        file_name_json = self.path_results + 'oblique_model_features.json'

        json.dump(features, open(file_name_json, "w"), indent=4)


        #file_name_json_errors = self.path_results + 'oblique_model_errors.json'

        #json.dump(errors, open(file_name_json_errors, "wb"), indent=4)

        if self.save_all:
            file_name_epsps_sync = self.path_results + 'oblique_model_epsp_amps_sync.p'

            pickle.dump(EPSPs_sync, gzip.GzipFile(file_name_epsps_sync, "wb"))


            file_name_epsps_async = self.path_results + 'oblique_model_epsp_amps_async.p'

            pickle.dump(EPSPs_async, gzip.GzipFile(file_name_epsps_async, "wb"))


            file_name_derivs_sync = self.path_results + 'oblique_model_mean_peak_derivs_sync.p'

            pickle.dump(sync_peak_derivatives, gzip.GzipFile(file_name_derivs_sync, "wb"))


            file_name_derivs_async = self.path_results + 'oblique_model_mean_peak_derivs_async.p'

            pickle.dump(async_peak_derivatives, gzip.GzipFile(file_name_derivs_async, "wb"))


        print("Results are saved in the directory: ", self.path_results)

        efel.reset()

        return prediction

    def compute_score(self, observation, prediction, verbose=False):
        """Implementation of sciunit.Test.score_prediction."""


        file_name = self.path_results + 'oblique_features.p'

        results=[]
        results.append(observation)
        results.append(prediction)
        if self.save_all:
            pickle.dump(results, gzip.GzipFile(file_name, "wb"))

        score0 = scores.P_Value_ObliqueIntegration.ttest_calc(observation,prediction)
        score1, errors_dict = scores.ZScore_ObliqueIntegration.compute(observation,prediction)

        score=scores.ZScore_ObliqueIntegration(score1)

        p_values = {'threshold' : score0[0],
                     'proximal_threshold' : score0[1],
                     'distal_threshold' : score0[2],
                     'peak_dV/dt_at_th' : score0[3],
                     'degree_of_nonlinearity_at_th' : score0[4],
                     'suprath_degree_of_nonlinearity' : score0[5],
                     'amplitude_at_th' : score0[6],
                     'time_to_peak_at_th' : score0[7],
                     'asynch_degree_of_nonlin_at_th': score0[8]}

        file_name_p = self.path_results + 'p_values.json'
        json.dump(p_values, open(file_name_p, "w"), indent=4)

        file_name_errors = self.path_results + 'oblique_model_errors.json'
        json.dump(errors_dict, open(file_name_errors, "w"), indent=4)

        plt.figure()
        x =numpy.arange(1,10)
        labels=['threshold', 'proximal threshold', 'distal threshold', 'peak dV/dt at th.','degree of nonlinearity at th.', 'suprath. degree of nonlinearity', 'amplitude at th.', 'time to peak at th.', 'asynch. degree of nonlin. at th.']
        plt.semilogy(x, score0, linestyle='None', marker='o')
        plt.xticks(x, labels, rotation=20)
        plt.tick_params(labelsize=11)
        plt.axhline(y=0.05, label='0.05', color='red')
        plt.legend()
        plt.margins(0.1)
        plt.ylabel("p values")
        fig = plt.gcf()
        fig.set_size_inches(12, 10)
        if self.save_all and self.path_figs is not None:
            plt.savefig(self.path_figs + 'p_values' + '.pdf', dpi=600, bbox_inches='tight')

        plt.figure(figsize = (210/25.4, 210/25.4))
        y=0
        labels=[]
        for key, value in errors_dict.items():
            labels.append(key)
            plt.plot(value, y, linestyle='None', marker='o', color = 'b' )
            y+=1
        plt.yticks(list(range(len(list(errors_dict.keys())))), labels)
        plt.title('Errors')
        plt.xlabel('error (# sd)')
        if self.save_all and self.path_figs is not None:
            plt.savefig(self.path_figs + 'mean_errors' + '.pdf', dpi=600, bbox_inches='tight')


        if self.show_plot:
            plt.show()

        final_score={'score' : str(score1)}
        file_name_score= self.path_results + 'final_score.json'
        json.dump(final_score, open(file_name_score, "w"), indent=4)

        self.logFile.write(str(score)+'\n')
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        self.logFile.close()

        self.logFile = self.path_results + self.test_log_filename

        return score

    def bind_score(self, score, model, observation, prediction):

        if self.path_figs is not None:
             score.related_data["figures"] = [self.path_figs + 'errors_sync.pdf', self.path_figs + 'input_output_curves_async.pdf',
                                        self.path_figs + 'input_output_curves_sync.pdf', self.path_figs + 'mean_errors.pdf',
                                        self.path_figs + 'mean_values_sync.pdf', self.path_figs + 'nonlin_errors_async.pdf',
                                        self.path_figs + 'nonlin_values_async.pdf', self.path_figs + 'peak_derivative_plots_async.pdf',
                                        self.path_figs + 'peak_derivative_plots_sync.pdf', self.path_figs + 'p_values.pdf',
                                        self.path_figs + 'somatic_traces_async.pdf', self.path_figs + 'somatic_traces_sync.pdf',
                                        self.path_figs + 'summary_input_output_curve_sync.pdf', self.path_figs + 'traces_async.pdf',
                                        self.path_figs + 'traces_sync.pdf', self.path_figs + 'values_sync.pdf', self.path_results + 'oblique_model_features.json',
                                        self.path_results + 'p_values.json', self.path_results + 'oblique_model_errors.json',
                                        self.path_results + 'final_score.json', self.path_results + self.test_log_filename]
        score.related_data["results"] = [self.path_results + 'oblique_model_features.json', self.path_results + 'oblique_model_epsp_amps_sync.p',
                                        self.path_results + 'oblique_model_mean_peak_derivs_sync.p', self.path_results + 'oblique_model_epsp_amps_async.p',
                                        self.path_results + 'oblique_model_mean_peak_derivs_async.p', self.path_results + 'oblique_features.p',
                                        self.path_results + 'p_values.json', self.path_results + 'oblique_model_errors.json',
                                        self.path_results + 'final_score.json']
        return score
