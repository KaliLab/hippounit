from quantities.quantity import Quantity
import sciunit
from sciunit import Test,Score,ObservationError
import hippounit.capabilities as cap
from sciunit.utils import assert_dimensionless# Converters.
from sciunit.scores import BooleanScore,ZScore # Scores.

try:
	import numpy
except:
	print("NumPy not loaded.")

import matplotlib.pyplot as plt
import matplotlib
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
    import cPickle as pickle
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
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
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
	copyreg.pickle(MethodType, _pickle_method, _unpickle_method)
except:
	copy_reg.pickle(MethodType, _pickle_method, _unpickle_method)


class DepolarizationBlockTest(Test):
	"""Tests if the model enters depolarization block under current injection of increasing amplitudes."""

	def __init__(self,
			     observation = {'mean_Ith':None, 'Ith_std':None, 'mean_Veq': None, 'Veq_std': None}  ,
			     name="Depolarization block test" ,
				 force_run=False,
				 base_directory= '/home/osboxes/BBP_project/150904_neuronunit/neuronunit/',
				show_plot=True):

		Test.__init__(self,observation,name)

		self.required_capabilities += (cap.ReceivesSquareCurrent_ProvidesResponse,)

		self.force_run = force_run
		self.directory = base_directory + 'temp_data/'
		self.directory_results = base_directory + 'results/'
		self.directory_figs = base_directory + 'figs/'

		self.show_plot = show_plot

		self.path_temp_data = None #added later, because model name is needed
		self.path_figs = None
		self.path_results = None

		self.npool = 4



		description = "Tests if the model enters depolarization block under current injection of increasing amplitudes."

	score_type = scores.ZScore_depolblock


	def cclamp(self, model, amp, delay, dur):

		#path = self.directory + model.name + '/depol_block/'
		self.path_temp_data = self.directory + model.name + '/depol_block/'

		try:
			if not os.path.exists(self.path_temp_data):
				os.makedirs(self.path_temp_data)
		except OSError, e:
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
			pickle.dump(traces, gzip.GzipFile(file_name, "wb"))

		else:
		    traces = pickle.load(gzip.GzipFile(file_name, "rb"))

		return traces

	def find_Ith_Veq(self, model, results, amps):


		#path_figs = self.directory_figs + 'depol_block/' + model.name + '/'
		self.path_figs = self.directory_figs + 'depol_block/' + model.name + '/'

		try:
			if not os.path.exists(self.path_figs):
				os.makedirs(self.path_figs)
		except OSError, e:
			if e.errno != 17:
				raise
			pass

		print "The figures are saved in the directory: ", self.path_figs

		spikecount_array=numpy.array([])

		for i, amp in enumerate(amps):

		    spikecount_array=numpy.append(spikecount_array, results[i][1][0]['Spikecount'])


		max=numpy.amax(spikecount_array)


		Ith_index = numpy.where(spikecount_array==max)[0]		# this is an array if there are a lot of same values in spike_count_array

		if Ith_index.size > 1 or Ith_index == (spikecount_array.size) -1:     # If the max num AP is the last element, it didn`t enter depol. block
		    Ith=float('NaN')                                 # If Ith == None, it means it didn`t enter depol block!!!
		    Veq=float('NaN')
		    Veq_index=Ith_index
		    Veq_index = int(Veq_index[0])
		    plt.figure(1)
		    plt.plot(results[spikecount_array.size-1][0]['T'],results[spikecount_array.size-1][0]['V'])
		    plt.title("somatic response to the highest current intensity\n (The model did not enter depol. block.)")
		    print " the model did not enter depolarization block"
		    plt.savefig(self.path_figs + 'somatic_resp_at_depol_block' + '.pdf', dpi=300)


		else:
		    Ith=amps[Ith_index]
		    Ith=Ith[0]
		    Veq_index=Ith_index+1
		    Veq_index = int(Veq_index[0])

		    plt.figure(1)
		    plt.plot(results[Veq_index][0]['T'],results[Veq_index][0]['V'])
		    plt.title("somatic response at Ith + 0.05 nA", fontsize=20)
		    plt.xlabel("time (ms)", fontsize=20)
		    plt.ylabel("Somatic voltage (mV)", fontsize=20)
		    plt.tick_params(labelsize=18)
		    plt.savefig(self.path_figs + 'somatic_resp_at_depol_block' + '.pdf', dpi=300)

		    Veq_trace=results[Veq_index][0]['V']
		    time=numpy.array(results[Veq_index][0]['T'])
		    indices1 = numpy.where(1400<=time)[0]           # we need the last 100ms of the current pulse
		    indices2 = numpy.where(1500>=time)[0]
		    trace_end_index_beginning = numpy.amin(indices1)
		    trace_end_index_end=numpy.amax(indices2)
		    trace_end=Veq_trace[trace_end_index_beginning:trace_end_index_end]      # THIS CONTAINS the last 100ms of the current pulse

		    Veq=numpy.average(trace_end)


		print "Ith (the current intensity for which the model exhibited the maximum number of APs):", Ith *nA
		print "Veq (the equilibrium value during the depolarization block):", Veq * mV


		plt.figure(2)
		fig = plt.gcf()
		#fig.set_size_inches(14, 12)
		plt.plot(amps,spikecount_array,'o-', markersize=10)
		plt.tick_params(labelsize=20)
		plt.xlabel("I (nA)",fontsize=20)
		#plt.ylabel("number of APs")
		plt.ylabel("num. of APs",fontsize=20)
		plt.margins(0.01)
		plt.savefig(self.path_figs + 'number_of_APs' + '.pdf', dpi=600)
		#plt.savefig(self.path_figs + 'num. of Aps' + '.pdf')

		if Ith_index.size > 1:
			Ith_index = int(Ith_index[-1])
			#Veq_index = int(Veq_index[-1])
		elif Ith_index.size == 1:
			Ith_index = int(Ith_index[0])
			#Veq_index = int(Veq_index[0])

		plt.figure(3)
		plt.plot(results[Ith_index][0]['T'],results[Ith_index][0]['V'])
		plt.title("somatic response at Ith", fontsize=20)
		plt.xlabel("time (ms)", fontsize=20)
		plt.ylabel("Somatic voltage (mV)", fontsize=20)
		plt.tick_params(labelsize=18)
		plt.savefig(self.path_figs + 'somatic_resp_at_Ith' + '.pdf', dpi=600)

		x =numpy.array([1, 2])
		Ith_array = numpy.array([self.observation['mean_Ith'], Ith])
		labels = ['Target Ith with SD', model.name]
		e = numpy.array([self.observation['Ith_std'], 0.0])

		x2 =numpy.array([1])
		y2 = numpy.array([self.observation['mean_Ith']])
		e = numpy.array([self.observation['Ith_std']])

		plt.figure(4)
		plt.plot(x, Ith_array, 'o')
		plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
		plt.xticks(x, labels, rotation=10)
		plt.margins(0.2)
		plt.ylabel("Ith (nA)")
		plt.savefig(self.path_figs + 'Ith' + '.pdf', dpi=600)

		x =numpy.array([1, 2])
		Veq_array = numpy.array([self.observation['mean_Veq'], Veq])
		labels = ['Target Veq with SD',model.name]
		e = numpy.array([self.observation['Veq_std'], 0.0])

		x2 =numpy.array([1])
		y2 = numpy.array([self.observation['mean_Veq']])
		e = numpy.array([self.observation['Veq_std']])

		plt.figure(5)
		plt.plot(x, Veq_array, 'o')
		plt.errorbar(x2, y2, e, linestyle='None', marker='o', color='blue')
		plt.xticks(x, labels, rotation=10)
		plt.margins(0.2)
		plt.ylabel("Veq (mV)")
		plt.savefig(self.path_figs + 'Veq' + '.pdf', dpi=600)

    	#errors
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


		num_AP_min=13
		num_AP_max=82
		labels2 = [model.name]

		plt.figure(7)
		plt.axhline(y=num_AP_min, label='Min. num.of AP\n observed\n experimentally', color='green')
		plt.axhline(y=num_AP_max, label='Max. num.of AP\n observed\n experimentally', color='red')
		plt.legend()

		plt.plot([1], spikecount_array[Ith_index], 'o')
		plt.title("For the models that doesn't enter depol block,\n the plot shows the num of AP-s for the highest current intensity")
		plt.xticks([1], labels2)
		plt.margins(0.2)
		plt.ylabel("Number of AP at Ith")
		plt.savefig(self.path_figs + 'num_of_APs_at_Ith' + '.pdf', dpi=600)

		#path_dir = self.directory_results + model.name + '/'
		self.path_results = self.directory_results + model.name + '/'

		try:
			if not os.path.exists(self.path_results):
				os.makedirs(self.path_results)
		except OSError, e:
			if e.errno != 17:
				raise
			pass

		file_name_f = self.path_results + 'depol_block_features_traces.p'
		file_name_json = self.path_results + 'depol_block_model_features.json'
		file_name_json_err = self.path_results + 'depol_block_model_errors.json'

		errors={}
		errors['Ith_error']=float(Ith_error)
		errors['Veq_error']=float(Veq_error)

		json.dump(errors, open(file_name_json_err, "wb"), indent=4)

		features_json={}
		features_json['Ith']=str(float(Ith) * nA)
		features_json['Veq']=str(float(Veq) * mV)
		json.dump(features_json, open(file_name_json, "wb"), indent=4)

		features={}
		features['Ith']=[Ith]
		features['Veq']=[Veq]
		features['Ith_error']=[float(Ith_error)]
		features['Veq_error']=[float(Veq_error)]
		features['Spikecount']=spikecount_array
		features['Ith_trace']=results[Ith_index][0]['V']  # trace where Ith is measured (max num of APs)
		features['Ith_time']=results[Ith_index][0]['T']
		features['Veq_trace']=results[Veq_index][0]['V']  # trace where Veq is measured
		features['Veq_time']=results[Veq_index][0]['T']

		pickle.dump(features, gzip.GzipFile(file_name_f, "wb"))

		print "Results are saved in the directory: ", self.path_results


		return Ith, Veq

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

		Ith, Veq = self.find_Ith_Veq(model, results, amps)

		prediction = {'model_Ith':float(Ith)*nA,'model_Veq': float(Veq)*mV}

		return prediction

	def compute_score(self, observation, prediction, verbose=False):
		"""Implementation of sciunit.Test.score_prediction."""
		score0 = scores.ZScore_depolblock.compute(observation,prediction)
		score=scores.ZScore_depolblock(score0)

		if self.show_plot:
			plt.show()

		return score

	def bind_score(self, score, model, observation, prediction):

		score.related_data["figures"] = [self.path_figs + 'Ith.pdf', self.path_figs + 'Veq.pdf', self.path_figs + 'number_of_APs.pdf', self.path_figs + 'num_of_APs_at_Ith.pdf', self.path_figs + 'somatic_resp_at_depol_block.pdf', self.path_figs + 'somatic_resp_at_Ith.pdf']
		score.related_data["results"] = [self.path_results + 'depol_block_model_errors.json', self.path_results + 'depol_block_model_features.json']
		return score
