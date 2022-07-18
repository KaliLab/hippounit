from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
#from builtins import str
from quantities.quantity import Quantity
from sciunit import Test
try:
    from sciunit import ObservationError
except:
    from sciunit.errors import ObservationError
import hippounit.capabilities as cap

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

from quantities import mV, nA, ms, V, s, Hz

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


class APFrequencyTest(Test):
    """
    Evaluates the AP frequency over a range of injected current stimuli.

    Parameters
    ----------
    observation : dict
        JSON file containing the experimental mean and std values of AP frequency for each level of stimulus
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
                 observation = [],
                 name="AP Frequency test" ,
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

        description = "Evaluates the AP frequency over a range of injected current stimuli."

    score_type = scores.ZScore_APfreq

    def format_data(self, observation):
        
        for idx, entry in enumerate(observation):
          for key, val in entry.items():
                try:
                    assert type(observation[key]) is Quantity
                except Exception as e:
                    quantity_parts = val.split(" ")
                    number = float(quantity_parts[0])
                    units = " ".join(quantity_parts[1:])
                    observation[idx][key] = Quantity(number, units)
        return observation

    def validate_observation(self, observation):

        try:
            assert type(observation) is list
            for entry in observation:
                assert type(entry) is dict
                assert all(key in entry.keys() for key in ["i_inj", "mean", "std"])
                for key in entry.keys():
                    assert type(entry[key]) is Quantity
        except Exception as e:
            raise ObservationError(("Observation must be of the form "
                                    "[{'i_inj':float*mV, 'mean':float*mV, 'std':float*mV}, ...]"))

    def cclamp(self, model, amp, delay, dur):

        if self.base_directory:
            self.path_temp_data = self.base_directory + 'temp_data/' + 'ap_freq/' + model.name + '/'
        else:
            self.path_temp_data = model.base_directory + 'temp_data/' + 'ap_freq/'

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

            trace['T'] = t
            trace['V'] = v
            trace['stim_start'] = [delay]
            trace['stim_end'] = [delay + dur]
            trace['stim_amp'] = [amp]
            traces.append(trace)

            traces_results = efel.getFeatureValues(traces, ['Spikecount'])
            traces.append(traces_results)
            if self.save_all:
                pickle.dump(traces, gzip.GzipFile(file_name, "wb"))
        else:
            traces = pickle.load(gzip.GzipFile(file_name, "rb"))

        return traces

    def generate_prediction(self, model, verbose=False):
        """Implementation of sciunit.Test.generate_prediction."""

        efel.reset()

        pool = multiprocessing.Pool(self.npool, maxtasksperchild=1)

        # stimulus levels (current injection) extracted from observation
        amps = [x["i_inj"] for x in self.observation]

        cclamp_ = functools.partial(self.cclamp, model, delay = 500, dur = 1000)
        results = pool.map(cclamp_, amps, chunksize=1)

        # Generate prediction
        # frequency is in Hz since stimulus and evaluation is for 1000 ms (1s)
        prediction = []
        for entry in results:
            prediction.append({"i_inj": entry[0]['stim_amp'][0], "freq": entry[1][0]["Spikecount"][0] * Hz})
        sorted(prediction, key=lambda d: d['i_inj']) 

        pool.terminate()
        pool.join()
        del pool

        plt.close('all') #needed to avoid overlapping of saved images when the test is run on multiple models in a for loop

        efel.reset()

        # Create required output directories - 1
        if self.base_directory:
            self.path_figs = self.base_directory + 'figures/' + 'ap_freq/' + model.name + '/'
        else:
            self.path_figs = model.base_directory + 'figures/' + 'ap_freq/'

        try:
            if not os.path.exists(self.path_figs) and self.save_all:
                os.makedirs(self.path_figs)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        # Create required output directories - 2
        if self.base_directory:
            self.path_results = self.base_directory + 'results/' + 'ap_freq/' + model.name + '/'
        else:
            self.path_results = model.base_directory + 'results/' + 'ap_freq/'
        
        try:
            if not os.path.exists(self.path_results):
                os.makedirs(self.path_results)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        # Generate ap_freq_stim_X.pdf figures
        plt.figure()
        for entry in results:
            plt.plot(entry[0]['T'], entry[0]['V'])
            plt.title("Somatic response to stimulus = " + str(entry[0]['stim_amp'][0]))
            plt.xlabel("Time (ms)")
            plt.ylabel("Somatic voltage (mV)")
            #plt.tick_params(labelsize=18)
            if self.save_all:
                fig_name = "ap_freq_stim_" + str(entry[0]['stim_amp'][0]) + '.pdf'
                plt.savefig(self.path_figs + fig_name, dpi=600, bbox_inches='tight')
                self.figures.append(self.path_figs + fig_name)
            plt.close('all') 

        self.logFile = open(self.path_results + self.test_log_filename, 'w')

        return prediction

    def compute_score(self, observation, prediction, verbose=False):
        """Implementation of sciunit.Test.score_prediction."""

        # Following result related files will be generated by this test:
        #   JSON
        #       - current_amps_spikecounts.json (for each simulated i_inj freq, exp mean, std, score)
        #       - ap_freq_summary.json  (observation, prediction, final score)
        #   Logs
        #       - test_log.txt
        #   Figures
        #       - ap_freq_fI_plot.pdf  (f-I relationship plot; model vs exp Z-score for each level of i_inj)
        #       - ap_freq_stim_X.pdf (see generate_prediction(); multiple Vm vs t plots; one per stimulus level 'X') 

        # Evaluate the score
        score, stim_spikecounts = scores.ZScore_APfreq.compute(observation, prediction)

        # Generate current_amps_spikecounts.json
        file_name_sc = self.path_results + 'current_amps_spikecounts.json'
        json.dump(stim_spikecounts, open(file_name_sc, "w"), default=str, indent=4)

        # Generate ap_freq_summary.json
        summary = {
            "observation": observation,
            "prediction": prediction,
            "score": score
        }
        file_name_summary = self.path_results + 'ap_freq_summary.json'
        json.dump(summary, open(file_name_summary, "w"), default=str, indent=4)

        # Generate ap_freq_fI_plot.pdf
        amps = [ float(x["i_inj"]) for x in stim_spikecounts ]
        freqs = [ float(x["freq"]) for x in stim_spikecounts ]
        means = [ float(x["mean"]) for x in stim_spikecounts ]
        stds = [ float(x["std"]) for x in stim_spikecounts ]
        
        plt.figure()
        fig = plt.gcf()
        plt.plot(amps, freqs, 'o-r')
        plt.errorbar(amps, means, stds, linestyle='None', marker='o', capsize=4, color='black')
        plt.title("f-I relationship")
        plt.legend(['Model simulation', 'Experimental data'], loc='best')
        plt.xlabel("$I_{inj} (nA)$")
        plt.ylabel("AP freq (Hz)")
        plt.margins(0.01)
        if self.save_all:
            fig_name = self.path_figs + "ap_freq_fI_plot.pdf"
            plt.savefig(fig_name, dpi=600, bbox_inches='tight')
            self.figures.append(fig_name)

        if self.show_plot:
            plt.show()

        self.logFile.write("Overall score: " + str(score) + "\n")
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")
        self.logFile.close()

        return score

    def bind_score(self, score, model, observation, prediction):

        self.figures.append(self.path_results + 'current_amps_spikecounts.json')
        self.figures.append(self.path_results + 'ap_freq_summary.json')
        self.figures.append(self.path_results + self.test_log_filename)
        score.related_data["figures"] = self.figures
        return score
