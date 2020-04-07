from __future__ import print_function
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
import pkg_resources

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
    import copyreg
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


class SomaticFeaturesTest(Test):
    """
    Tests some somatic features under current injection of increasing amplitudes.

    Parameters
    ----------
    config : dict
        dictionary loaded from a JSON file, containing the parameters of the simulation
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
    specify_data_set : str
        When set to a string, output will be saved into subdirectory (within the model_name subderotory) named like this. This makes it possible to run the validation on a specific model, against different data sets, and save the results separately.
    """

    def __init__(self,
                 observation = {}  ,
                 config = {},
                 name="Somatic features test" ,
                 force_run=False,
                 base_directory=None,
                 show_plot=True,
                 save_all = True,
                 specify_data_set = ''):


        Test.__init__(self,observation,name)

        self.required_capabilities += (cap.ReceivesSquareCurrent_ProvidesResponse,)

        self.force_run = force_run
        self.show_plot = show_plot
        self.save_all = save_all
        self.config = config

        self.base_directory = base_directory

        self.path_temp_data = None #added later, because model name is needed
        self.path_figs = None
        self.path_results = None
        self.npool = multiprocessing.cpu_count() - 1

        self.logFile = None
        self.test_log_filename = 'test_log.txt'
        self.specify_data_set = specify_data_set  #this is added to the name of the directory (somaticfeat), so tests runs using different data sets can be saved into different directories


        plt.close('all') #needed to avoid overlapping of saved images when the test is run on multiple models in a for loop
        #with open('./stimfeat/PC_newfeat_No14112401_15012303-m990803_stimfeat.json') as f:
            #self.config = json.load(f, object_pairs_hook=collections.OrderedDict)

        description = "Tests some somatic features under current injection of increasing amplitudes."

    score_type = scores.ZScore_somaticSpiking

    def create_stimuli_list(self):

        #with open('./stimfeat/PC_newfeat_No14112401_15012303-m990803_stimfeat.json') as f:
            #config = json.load(f, object_pairs_hook=collections.OrderedDict)

        stimulus_list=[]
        stimuli_list=[]
        stimuli_names = []

        stimuli_names=list(self.config['stimuli'].keys())



        for i in range (0, len(stimuli_names)):
            stimulus_list.append(stimuli_names[i])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['Amplitude'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['Delay'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['Duration'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['StimSectionName'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['StimLocationX'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['Type'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['RecSectionName'])
            stimulus_list.append(self.config['stimuli'][stimuli_names[i]]['RecLocationX'])
            stimuli_list.append(stimulus_list)
            stimulus_list=[]

        return stimuli_list

    def create_features_list(self, observation):

        feature_list=[]
        features_list=[]
        features_names=(list(observation.keys()))


        for i in range (0, len(features_names)):
            feature_list.append(features_names[i])
            feature_list.append(observation[features_names[i]]['Std'])
            feature_list.append(observation[features_names[i]]['Mean'])
            feature_list.append(observation[features_names[i]]['Stimulus'])
            feature_list.append(observation[features_names[i]]['Type'])
            features_list.append(feature_list)
            feature_list=[]

        return features_names, features_list

    def run_stim(self, model, stimuli_list):


        stimulus_name, amplitude, delay, duration, stim_section_name, stim_location_x, stim_type, rec_section_name, rec_location_x = stimuli_list

        traces_result={}

        if self.specify_data_set != '':
            specify_data_set = '_' + self.specify_data_set
        else:
            specify_data_set = self.specify_data_set

        if self.base_directory:
            self.path_temp_data = self.base_directory + 'temp_data/' + 'somaticfeat' + specify_data_set + '/' + model.name + '/'
        else:
            self.path_temp_data = model.base_directory + 'temp_data/' + 'somaticfeat' + specify_data_set + '/'

        try:
            if not os.path.exists(self.path_temp_data) and self.save_all:
                os.makedirs(self.path_temp_data)
        except OSError as e:
            if e.errno != 17:
                raise
            pass


        if stim_type == "SquarePulse":
            file_name = self.path_temp_data + stimulus_name + '.p'

            if self.force_run or (os.path.isfile(file_name) is False):


                t, v = model.get_vm(float(amplitude), float(delay), float(duration), stim_section_name, stim_location_x, rec_section_name, rec_location_x)

                traces_result[stimulus_name]=[t,v]
                if self.save_all:
                    pickle.dump(traces_result, gzip.GzipFile(file_name, "wb"))

            else:
                traces_result = pickle.load(gzip.GzipFile(file_name, "rb"))

        else:
            traces_result=None

        return traces_result

    def analyse_traces(self, stimuli_list, traces_results, features_list):

        feature_name, target_sd, target_mean, stimulus, feature_type = features_list
        target_sd=float(target_sd)
        target_mean=float(target_mean)

        feature_result={}
        trace = {}
        for i in range (0, len(traces_results)):
            for key, value in traces_results[i].items():
                stim_name=key
            if stim_name == stimulus:

                trace['T'] = traces_results[i][stim_name][0]
                trace['V'] = traces_results[i][stim_name][1]


        for i in range (0, len(stimuli_list)):
            if stimuli_list[i][0]==stimulus:

                trace['stim_start'] = [float(stimuli_list[i][2])]
                trace['stim_end'] = [float(stimuli_list[i][2])+float(stimuli_list[i][3])]

        traces = [trace]
        #print traces

        efel_results = efel.getFeatureValues(traces,[feature_type])

        feature_values=efel_results[0][feature_type]

        if feature_values is not None and feature_values.size != 0:

            if (feature_type == 'AP_rise_time' or feature_type == 'AP_amplitude' or feature_type == 'AP_duration_half_width' or feature_type == 'AP_begin_voltage'
                or feature_type == 'AP_rise_rate' or feature_type == 'fast_AHP' or feature_type == 'AP_begin_time' or feature_type == 'AP_begin_width' or feature_type == 'AP_duration'
                or feature_type == 'AP_duration_change' or feature_type == 'AP_duration_half_width_change' or feature_type == 'fast_AHP_change' or feature_type == 'AP_rise_rate_change' or feature_type == 'AP_width'):
               """
               In case of features that are AP_begin_time/AP_begin_index, the 1st element of the resulting vector, which corresponds to AP1, is ignored
               This is because the AP_begin_time/AP_begin_index feature often detects the start of the stimuli instead of the actual beginning of AP1
               """
               feature_mean=numpy.mean(feature_values[1:])
               feature_sd=numpy.std(feature_values[1:])
            else:
               feature_mean=numpy.mean(feature_values)
               feature_sd=numpy.std(feature_values)


        else:
            feature_mean = float('nan')
            feature_sd = float('nan')

        #feature_mean=numpy.mean(feature_values)
        #feature_sd=numpy.std(feature_values)

        feature_result={feature_name:{'feature values': feature_values,
                                      'feature mean': feature_mean,
                                      'feature sd': feature_sd}}
        return feature_result

    def create_figs(self, model, traces_results, features_names, feature_results_dict, observation):


        if self.specify_data_set != '':
            specify_data_set = '_' + self.specify_data_set
        else:
            specify_data_set = self.specify_data_set

        if self.base_directory:
            self.path_figs = self.base_directory + 'figs/' + 'somaticfeat' + specify_data_set + '/' + model.name + '/'
        else:
            self.path_figs = model.base_directory + 'figs/' + 'somaticfeat' + specify_data_set + '/'

        try:
            if not os.path.exists(self.path_figs) and self.save_all:
                os.makedirs(self.path_figs)
        except OSError as e:
            if e.errno != 17:
                raise
            pass
        if self.save_all:
            print("The figures are saved in the directory: ", self.path_figs)

        plt.figure(1)
        #key=sorted()
        for i in range (0, len(traces_results)):
            for key, value in traces_results[i].items():
                plt.plot(traces_results[i][key][0], traces_results[i][key][1], label=key)
        plt.legend(loc=2)
        if self.save_all:
            plt.savefig(self.path_figs + 'traces' + '.pdf', dpi=600,)


        columns = 2
        width_ratios=[1]*columns
        frames = len(traces_results)
        rows = int(numpy.ceil(frames/float(columns)))
        height_ratios=[1]*rows
        #axs=[]

        fig = plt.figure(figsize = (210/25.4, 297/25.4))
        gs = matplotlib.gridspec.GridSpec(rows, columns, height_ratios=height_ratios, width_ratios=width_ratios)
        gs.update(top=0.97, bottom=0.04, left=0.07, right=0.97, hspace=0.75, wspace=0.3)
        #fig, axes = plt.subplots(nrows=int(round(len(traces_results)/2.0)), ncols=2)
        #fig.tight_layout()
        for i in range (0, len(traces_results)):

            for key, value in traces_results[i].items():

                #plt.subplot(round(len(traces_results)/2.0),2,i+1)
                plt.subplot(gs[i])
                #axs.append(fig.add_subplot(gs[i]))
                plt.plot(traces_results[i][key][0], traces_results[i][key][1])
                plt.title(key)
                plt.xlabel("ms")
                plt.ylabel("mV")
                minx = float(self.config['stimuli'][key]['Delay']) - 200
                maxx = float(self.config['stimuli'][key]['Delay']) + float(self.config['stimuli'][key]['Duration']) + 200
                plt.xlim(minx, maxx)
                #plt.tick_params(labelsize=15)
        #gs.tight_layout(fig)
        #fig = plt.gcf()
        #fig.set_size_inches(12, 10)
        if self.save_all:
            plt.savefig(self.path_figs + 'traces_subplots' + '.pdf', dpi=600, bbox_inches='tight')

        axs = plottools.tiled_figure("absolute features", figs={}, frames=1, columns=1, orientation='page',
                                height_ratios=None, top=0.97, bottom=0.05, left=0.25, right=0.97, hspace=0.1, wspace=0.2)

        plt.gcf().set_size_inches(210/25.4, 297/25.4*2 )

        label_added = False

        for i in range (len(features_names)):
            feature_name=features_names[i]
            y=i
            if not label_added:
                axs[0].errorbar(feature_results_dict[feature_name]['feature mean'], y, xerr=feature_results_dict[feature_name]['feature sd'], marker='o', color='blue', label = model.name)
                axs[0].errorbar(float(observation[feature_name]['Mean']), y, xerr=float(observation[feature_name]['Std']), marker='o', color='red', label = 'experiment')
                label_added = True
            else:
                axs[0].errorbar(feature_results_dict[feature_name]['feature mean'], y, xerr=feature_results_dict[feature_name]['feature sd'], marker='o', color='blue')
                axs[0].errorbar(float(observation[feature_name]['Mean']), y, xerr=float(observation[feature_name]['Std']), marker='o', color='red')
        axs[0].yaxis.set_ticks(list(range(len(features_names))))
        axs[0].set_yticklabels(features_names)
        axs[0].set_ylim(-1, len(features_names))
        axs[0].set_title('Absolute Features')
        lgd=plt.legend(bbox_to_anchor=(1.0, 1.0), loc = 'upper left')
        if self.save_all:
            plt.savefig(self.path_figs + 'absolute_features' + '.pdf', dpi=600, bbox_extra_artists=(lgd,), bbox_inches='tight')

    def generate_prediction(self, model, verbose=False):
        """Implementation of sciunit.Test.generate_prediction."""

        efel.reset()

        self.observation = collections.OrderedDict(sorted(self.observation.items()))

        global model_name_soma
        model_name_soma = model.name

        pool = multiprocessing.Pool(self.npool, maxtasksperchild=1)

        stimuli_list=self.create_stimuli_list()

        run_stim_ = functools.partial(self.run_stim, model)
        traces_results = pool.map(run_stim_, stimuli_list, chunksize=1)
        #traces_results = traces_result.get()


        pool.terminate()
        pool.join()
        del pool

        pool2 = multiprocessing.Pool(self.npool, maxtasksperchild=1)

        features_names, features_list = self.create_features_list(self.observation)

        analyse_traces_ = functools.partial(self.analyse_traces, stimuli_list, traces_results)
        feature_results = pool2.map(analyse_traces_, features_list, chunksize=1)
        #feature_results = feature_result.get()

        pool2.terminate()
        pool2.join()
        del pool2

        feature_results_dict={}
        for i in range (0,len(feature_results)):
            feature_results_dict.update(feature_results[i])  #concatenate dictionaries

        if self.specify_data_set != '':
            specify_data_set = '_' + self.specify_data_set
        else:
            specify_data_set = self.specify_data_set
        if self.base_directory:
            self.path_results = self.base_directory + 'results/' + 'somaticfeat' + specify_data_set + '/' + model.name + '/'
        else:
            self.path_results = model.base_directory + 'results/' + 'somaticfeat' + specify_data_set + '/'

        try:
            if not os.path.exists(self.path_results):
                os.makedirs(self.path_results)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        file_name=self.path_results+'soma_features.p'

        SomaFeaturesDict={}
        SomaFeaturesDict['traces_results']=traces_results
        SomaFeaturesDict['features_names']=features_names
        SomaFeaturesDict['feature_results_dict']=feature_results_dict
        SomaFeaturesDict['observation']=self.observation
        if self.save_all:
            pickle.dump(SomaFeaturesDict, gzip.GzipFile(file_name, "wb"))

        plt.close('all') #needed to avoid overlapping of saved images when the test is run on multiple models in a for loop

        self.create_figs(model, traces_results, features_names, feature_results_dict, self.observation)

        #prediction = feature_results_dict

        soma_features={}
        needed_keys = { 'feature mean', 'feature sd'}
        for i in range(len(SomaFeaturesDict['features_names'])):
            feature_name = SomaFeaturesDict['features_names'][i]
            soma_features[feature_name] = { key:value for key,value in list(feature_results_dict[feature_name].items()) if key in needed_keys }

        file_name_json = self.path_results + 'somatic_model_features.json'

        json.dump(soma_features, open(file_name_json, "w"), indent=4)

        prediction=soma_features

        efel.reset()

        return prediction

    def compute_score(self, observation, prediction, verbose=False):
        """Implementation of sciunit.Test.score_prediction."""

        try:
            if not os.path.exists(self.path_figs) and self.save_all:
                os.makedirs(self.path_figs)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        try:
            if not os.path.exists(self.path_results):
                os.makedirs(self.path_results)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        filepath = self.path_results + self.test_log_filename
        self.logFile = open(filepath, 'w')

        score_avg, feature_results_dict, features_names, bad_features  = scores.ZScore_somaticSpiking.compute(observation,prediction)


        if len(bad_features) > 0:
            self.logFile.write('Features excluded (due to invalid values):\n' + ', '.join(str(f) for f in bad_features) + '\n')
            self.logFile.write("---------------------------------------------------------------------------------------------------\n")

            print('Features excluded (due to invalid values):', ', '.join(str(f) for f in bad_features))

        self.logFile.write('Number of features succesfully evaluated: ' + str(len(list(feature_results_dict.keys())) - len(bad_features)) +'/' + str(len(list(feature_results_dict.keys())))+ '\n')
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        print('Number of features succesfully evaluated: ' + str(len(list(feature_results_dict.keys())) - len(bad_features)) +'/' + str(len(list(feature_results_dict.keys()))))



        file_name=self.path_results+'soma_errors.p'

        SomaErrorsDict={}
        SomaErrorsDict['features_names']=features_names
        SomaErrorsDict['feature_results_dict']=feature_results_dict
        if self.save_all:
            pickle.dump(SomaErrorsDict, gzip.GzipFile(file_name, "wb"))

        file_name_json = self.path_results + 'somatic_model_errors.json'
        json.dump(SomaErrorsDict['feature_results_dict'], open(file_name_json, "w"), indent=4)

        print("Results are saved in the directory: ", self.path_results)

        axs2 = plottools.tiled_figure("features", figs={}, frames=1, columns=1, orientation='page',
                                      height_ratios=None, top=0.97, bottom=0.05, left=0.25, right=0.97, hspace=0.1, wspace=0.2)
        plt.gcf().set_size_inches(210/25.4, 297/25.4*2 )

        for i in range (len(features_names)):
            feature_name=features_names[i]
            y=i
            axs2[0].plot(feature_results_dict[feature_name], y, marker='o', color='blue')
        axs2[0].yaxis.set_ticks(list(range(len(features_names))))
        axs2[0].set_yticklabels(features_names)
        axs2[0].set_ylim(-1, len(features_names))
        axs2[0].set_title('Feature errors')
        if self.save_all:
            plt.savefig(self.path_figs + 'Feature_errors' + '.pdf', dpi=600,)

        if self.show_plot:
            plt.show()

        final_score={'score' : str(score_avg)}
        file_name_score= self.path_results + 'final_score.json'
        json.dump(final_score, open(file_name_score, "w"), indent=4)

        score=scores.ZScore_somaticSpiking(score_avg)

        self.logFile.write(str(score)+'\n')
        self.logFile.write("---------------------------------------------------------------------------------------------------\n")

        self.logFile.close()

        self.logFile = self.path_results + self.test_log_filename

        return score

    def bind_score(self, score, model, observation, prediction):
        score.related_data["figures"] = [self.path_figs + 'traces.pdf', self.path_figs + 'absolute_features.pdf',
                                        self.path_figs + 'Feature_errors.pdf', self.path_figs + 'traces_subplots.pdf',
                                        self.path_results + 'somatic_model_features.json', self.path_results + 'somatic_model_errors.json',
                                        self.path_results + 'final_score.json', self.path_results + self.test_log_filename]
        score.related_data["results"] = [self.path_results + 'somatic_model_features.json', self.path_results + 'somatic_model_errors.json', self.path_results+'soma_errors.p', self.path_results+'soma_features.p', self.path_results + 'final_score.json']
        return score
