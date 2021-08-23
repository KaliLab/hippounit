from __future__ import print_function
from __future__ import division
# from builtins import str
from builtins import range
import os
import numpy
import sciunit
import hippounit.capabilities as cap
from quantities import ms,mV,Hz
from neuron import h

import multiprocessing
import zipfile
import collections

import collections

import json

import pkg_resources
import sys



class ModelLoader(sciunit.Model,
                 cap.ProvidesGoodObliques,
                 cap.ReceivesSquareCurrent_ProvidesResponse,
                 cap.ReceivesSynapse,
                 cap.ReceivesMultipleSynapses,
                 cap.ReceivesSquareCurrent_ProvidesResponse_MultipleLocations,
                 cap.ProvidesRecordingLocationsOnTrunk,
                 cap.ProvidesRandomDendriticLocations,
                 cap.ReceivesEPSCstim):

    def __init__(self, name="model", mod_files_path=None):
        """ Constructor. """

        """ This class should be used with Jupyter notebooks"""

        self.modelpath = mod_files_path
        self.libpath = 'x86_64/.libs/libnrnmech.so'
        self.hocpath = None

        self.cvode_active = False

        self.template_name = None
        self.SomaSecList_name = None
        self.max_dist_from_soma = 150
        self.v_init = -70
        self.celsius = 34

        self.name = name
        self.threshold = -20
        self.stim = None
        self.soma = None
        sciunit.Model.__init__(self, name=name)

        self.c_step_start = 0.00004
        self.c_step_stop = 0.000004
        self.c_minmax = numpy.array([0.00004, 0.04])

        self.ObliqueSecList_name = None
        self.TrunkSecList_name = None
        self.dend_loc = []  #self.dend_loc = [['dendrite[80]',0.27],['dendrite[80]',0.83],['dendrite[54]',0.16],['dendrite[54]',0.95],['dendrite[52]',0.38],['dendrite[52]',0.83],['dendrite[53]',0.17],['dendrite[53]',0.7],['dendrite[28]',0.35],['dendrite[28]',0.78]]
        self.dend_locations = collections.OrderedDict()
        self.NMDA_name = None
        self.default_NMDA_name = 'NMDA_CA1_pyr_SC'
        self.default_NMDA_path = pkg_resources.resource_filename("hippounit", "tests/default_NMDAr/")

        self.AMPA_name = None
        self.AMPA_NMDA_ratio = 2.0

        self.AMPA_tau1 = 0.1
        self.AMPA_tau2 = 2.0
        self.start=150

        self.ns = None
        self.ampa = None
        self.nmda = None
        self.ampa_nc = None
        self.nmda_nc = None

        self.ampa_list = []
        self.nmda_list = []
        self.ns_list = []
        self.ampa_nc_list = []
        self.nmda_nc_list = []

        self.ndend = None
        self.xloc = None

        self.base_directory = './validation_results/'   # inside current directory

        self.find_section_lists = False

        self.compile_mod_files()
        self.compile_default_NMDA()

    def translate(self, sectiontype, distance=0):

        if "soma" in sectiontype:
            return self.soma
        else:
            return False

    def compile_mod_files(self):

        if self.modelpath is None:
            raise Exception("Please give the path to the mod files (eg. mod_files_path = \'/home/models/CA1_pyr/mechanisms/\') as an argument to the ModelLoader class")

        if os.path.isfile(self.modelpath + self.libpath) is False:
            os.system("cd " + "\'" + self.modelpath + "\'" + "; nrnivmodl")

    def compile_default_NMDA(self):
        if os.path.isfile(self.default_NMDA_path + self.libpath) is False:
            os.system("cd " + "\'" + self.default_NMDA_path  + "\'" + "; nrnivmodl")

    def load_mod_files(self):

        h.nrn_load_dll(str(self.modelpath + self.libpath))


    def initialise(self):

        save_stdout=sys.stdout                   #To supress hoc output from Jupyter notebook 
        # sys.stdout=open("trash","w")
        sys.stdout=open('/dev/stdout', 'w')      #rather print it to the console 

        self.load_mod_files()

        if self.hocpath is None:
            raise Exception("Please give the path to the hoc file (eg. model.modelpath = \"/home/models/CA1_pyr/CA1_pyr_model.hoc\")")


        h.load_file("stdrun.hoc")
        h.load_file(str(self.hocpath))

        if self.soma is None and self.SomaSecList_name is None:
            raise Exception("Please give the name of the soma (eg. model.soma=\"soma[0]\"), or the name of the somatic section list (eg. model.SomaSecList_name=\"somatic\")")

        try:
            if self.template_name is not None and self.SomaSecList_name is not None:

                h('objref testcell')
                h('testcell = new ' + self.template_name)

                exec('self.soma_ = h.testcell.'+ self.SomaSecList_name)

                for s in self.soma_ :
                    self.soma = h.secname()

            elif self.template_name is not None and self.SomaSecList_name is None:
                h('objref testcell')
                h('testcell = new ' + self.template_name)
                # in this case self.soma is set in the jupyter notebook
            elif self.template_name is None and self.SomaSecList_name is not None:
                exec('self.soma_ = h.' +  self.SomaSecList_name)
                for s in self.soma_ :
                    self.soma = h.secname()
            # if both is None, the model is loaded, self.soma will be used
        except AttributeError:
            print ("The provided model template is not accurate. Please verify!")
        except Exception:
            print ("If a model template is used, please give the name of the template to be instantiated (with parameters, if any). Eg. model.template_name=CCell(\"morph_path\")")
            raise


        sys.stdout=save_stdout    #setting output back to normal 

    def inject_current(self, amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec):

        self.initialise()

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        stim_section_name = self.translate(section_stim, distance=0)
        rec_section_name = self.translate(section_rec, distance=0)
        #exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        exec("self.sect_loc_stim=h." + str(stim_section_name)+"("+str(loc_stim)+")")

        print("- running amplitude: " + str(amp)  + " on model: " + self.name + " at: " + stim_section_name + "(" + str(loc_stim) + ")")

        self.stim = h.IClamp(self.sect_loc_stim)

        self.stim.amp = amp
        self.stim.delay = delay
        self.stim.dur = dur

        #print "- running model", self.name, "stimulus at: ", str(self.soma), "(", str(0.5), ")"

        exec("self.sect_loc_rec=h." + str(rec_section_name)+"("+str(loc_rec)+")")

        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc_rec._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init#-65

        h.celsius = self.celsius
        h.init()
        h.tstop = delay + dur + 200
        h.run()

        t = numpy.array(rec_t)
        v = numpy.array(rec_v)

        return t, v

    def inject_current_record_respons_multiple_loc(self, amp, delay, dur, section_stim, loc_stim, dend_locations):

        self.initialise()

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        stim_section_name = self.translate(section_stim, distance=0)
        #rec_section_name = self.translate(section_rec, distance=0)
        #exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        exec("self.sect_loc_stim=h." + str(stim_section_name)+"("+str(loc_stim)+")")
        exec("self.sect_loc_rec=h." + str(stim_section_name)+"("+str(loc_stim)+")")

        print("- running amplitude: " + str(amp)  + " on model: " + self.name + " at: " + stim_section_name + "(" + str(loc_stim) + ")")

        self.stim = h.IClamp(self.sect_loc_stim)

        self.stim.amp = amp
        self.stim.delay = delay
        self.stim.dur = dur

        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v_stim = h.Vector()
        rec_v_stim.record(self.sect_loc_rec._ref_v)

        rec_v = []
        v = collections.OrderedDict()
        self.dend_loc_rec =[]

        '''
        for i in range(0,len(dend_loc)):

            exec("self.dend_loc_rec.append(h." + str(dend_loc[i][0])+"("+str(dend_loc[i][1])+"))")
            rec_v.append(h.Vector())
            rec_v[i].record(self.dend_loc_rec[i]._ref_v)
            #print self.dend_loc[i]
        '''
        #print dend_locations
        for key, value in dend_locations.items():
            for i in range(len(dend_locations[key])):
                exec("self.dend_loc_rec.append(h." + str(dend_locations[key][i][0])+"("+str(dend_locations[key][i][1])+"))")
                rec_v.append(h.Vector())

        for i in range(len(self.dend_loc_rec)):
            rec_v[i].record(self.dend_loc_rec[i]._ref_v)
            #print self.dend_loc[i]

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init#-65

        h.celsius = self.celsius
        h.init()
        h.tstop = delay + dur + 200
        h.run()

        t = numpy.array(rec_t)
        v_stim = numpy.array(rec_v_stim)

        '''
        for i in range(0,len(dend_loc)):
            v.append(numpy.array(rec_v[i]))
        '''

        i = 0
        for key, value in dend_locations.items():
            v[key] = collections.OrderedDict()
            for j in range(len(dend_locations[key])):
                loc_key = (dend_locations[key][j][0],dend_locations[key][j][1]) # list can not be a key, but tuple can
                v[key][loc_key] = numpy.array(rec_v[i])     # the list that specifies dendritic location will be a key too.
                i+=1

        return t, v_stim, v

    def classify_apical_point_sections(self, icell):

        import os
        import neurom as nm
        from hippounit import classify_apical_sections as cas

        '''
        for file_name in os.listdir(self.morph_path[1:-1]):
            filename = self.morph_path[1:-1]+ '/' + file_name
            break
        '''

        morph = nm.load_neuron(self.morph_full_path)

        apical_point_sections = cas.multiple_apical_points(morph)

        sections = cas.get_list_of_diff_section_types(morph, apical_point_sections)

        apical_trunk_isections = cas.get_neuron_isections(icell, sections['trunk'])
        #print sorted(apical_trunk_isections)

        apical_tuft_isections = cas.get_neuron_isections(icell, sections['tuft'])
        #print sorted(apical_tuft_isections)

        oblique_isections = cas.get_neuron_isections(icell, sections['obliques'])
        #print sorted(oblique_isections)

        return apical_trunk_isections, apical_tuft_isections, oblique_isections

    def find_trunk_locations(self, distances, tolerance, trunk_origin):

        if self.TrunkSecList_name is None and not self.find_section_lists:
            raise NotImplementedError("Please give the name of the section list containing the trunk sections. (eg. model.TrunkSecList_name=\"trunk\" or set model.find_section_lists to True)")

        #locations={}
        locations=collections.OrderedDict()
        actual_distances ={}
        dend_loc=[]

        if self.TrunkSecList_name is not None:
            self.initialise()

            if self.template_name is not None:
                exec('self.trunk=h.testcell.' + self.TrunkSecList_name)
            else:
                exec('self.trunk=h.' + self.TrunkSecList_name)


        if self.find_section_lists:

            self.initialise()

            if self.template_name is not None:
                exec('self.icell=h.testcell')

            apical_trunk_isections, apical_tuft_isections, oblique_isections = self.classify_apical_point_sections(self.icell)

            self.trunk = []
            for i in range(len(apical_trunk_isections)):
                exec('self.sec = h.testcell.apic[' + str(apical_trunk_isections[i]) + ']')
                self.trunk.append(self.sec)

        for sec in self.trunk:
            #for seg in sec:
            if not trunk_origin:
                h(self.soma + ' ' +'distance(0,1)') # For apical dendrites the default reference point is the end of the soma (point 1)
            elif len(trunk_origin) == 1:
                h(self.soma + ' ' +'distance(0,'+str(trunk_origin[0]) + ')') # Trunk origin point (reference for distance measurement) can be
            elif len(trunk_origin) == 2:
                h(trunk_origin[0] + ' ' +'distance(0,'+str(trunk_origin[1]) + ')') # Trunk origin point (reference for distance measurement) can be added by the user as an argument to the test
            #print sec.name()
            if self.find_section_lists:
                h('access ' + sec.name())

            for seg in sec:
                #print('SEC: ', sec.name())
                #print('SEG.X', seg.x)
                #print('DIST', h.distance(seg.x, sec=sec))
                #print('DIST0', h.distance(0, sec=sec))
                #print('DIST1', h.distance(1, sec=sec))
                for i in range(0, len(distances)):
                    locations.setdefault(distances[i], []) # if this key doesn't exist it is added with the value: [], if exists, value not altered
                    if h.distance(seg.x, sec=sec) < (distances[i] + tolerance) and h.distance(seg.x, sec=sec) > (distances[i]- tolerance): # if the seq is between distance +- 20
                        #print 'SEC: ', sec.name()
                        #print 'seg.x: ', seg.x
                        #print 'DIST: ', h.distance(seg.x)
                        locations[distances[i]].append([sec.name(), seg.x])
                        actual_distances[sec.name(), seg.x] = h.distance(seg.x, sec=sec)

        #print actual_distances
        return locations, actual_distances

    def get_random_locations(self, num, seed, dist_range, trunk_origin):

        if self.TrunkSecList_name is None and not self.find_section_lists:
            raise NotImplementedError("Please give the name of the section list containing the trunk sections. (eg. model.TrunkSecList_name=\"trunk\" or set model.find_section_lists to True)")

        locations=[]
        locations_distances = {}

        if self.TrunkSecList_name is not None:
            self.initialise()

            if self.template_name is not None:
                exec('self.trunk=h.testcell.' + self.TrunkSecList_name)

            else:
                exec('self.trunk=h.' + self.TrunkSecList_name)

        if self.find_section_lists:

            self.initialise()

            if self.template_name is not None:
                exec('self.icell=h.testcell')

            apical_trunk_isections, apical_tuft_isections, oblique_isections = self.classify_apical_point_sections(self.icell)
            apical_trunk_isections = sorted(apical_trunk_isections) # important to keep reproducability

            self.trunk = []
            for i in range(len(apical_trunk_isections)):
                exec('self.sec = h.testcell.apic[' + str(apical_trunk_isections[i]) + ']')
                self.trunk.append(self.sec)
        else:
            self.trunk = list(self.trunk)

        kumm_length_list = []
        kumm_length = 0
        num_of_secs = 0


        for sec in self.trunk:
            #print sec.L
            num_of_secs += sec.nseg
            kumm_length += sec.L
            kumm_length_list.append(kumm_length)
        #print 'kumm' ,kumm_length_list
        #print num_of_secs

        if num > num_of_secs:
            for sec in self.trunk:
                if not trunk_origin:
                    h(self.soma + ' ' +'distance(0,1)') # For apical dendrites the default reference point is the end of the soma (point 1)
                elif len(trunk_origin) == 1:
                    h(self.soma + ' ' +'distance(0,'+str(trunk_origin[0]) + ')') # Trunk origin point (reference for distance measurement) can be
                elif len(trunk_origin) == 2:
                    h(trunk_origin[0] + ' ' +'distance(0,'+str(trunk_origin[1]) + ')') # Trunk origin point (reference for distance measurement) can be added by the user as an argument to the test
                h('access ' + sec.name())
                for seg in sec:
                    if h.distance(seg.x, sec=sec) > dist_range[0] and h.distance(seg.x, sec=sec) < dist_range[1]:     # if they are out of the distance range they wont be used
                        locations.append([sec.name(), seg.x])
                        locations_distances[sec.name(), seg.x] = h.distance(seg.x, sec=sec)
            #print 'Dendritic locations to be tested (with their actual distances):', locations_distances

        else:

            norm_kumm_length_list = [i/kumm_length_list[-1] for i in kumm_length_list]
            #print 'norm kumm',  norm_kumm_length_list

            import random

            _num_ = num  # _num_ will be changed
            num_iterations = 0

            while len(locations) < num and num_iterations < 50 :
                #print 'seed ', seed
                random.seed(seed)
                rand_list = [random.random() for j in range(_num_)]
                #print rand_list

                for rand in rand_list:
                    #print 'RAND', rand
                    for i in range(len(norm_kumm_length_list)):
                        if rand <= norm_kumm_length_list[i] and (rand > norm_kumm_length_list[i-1] or i==0):
                            #print norm_kumm_length_list[i-1]
                            #print norm_kumm_length_list[i]
                            seg_loc = (rand - norm_kumm_length_list[i-1]) / (norm_kumm_length_list[i] - norm_kumm_length_list[i-1])
                            #print 'seg_loc', seg_loc
                            segs = [seg.x for seg in self.trunk[i]]
                            d_seg = [abs(seg.x - seg_loc) for seg in self.trunk[i]]
                            min_d_seg = numpy.argmin(d_seg)
                            segment = segs[min_d_seg]
                            #print 'segment', segment
                            if not trunk_origin:
                                h(self.soma + ' ' +'distance(0,1)') # For apical dendrites the default reference point is the end of the soma (point 1)
                            elif len(trunk_origin) == 1:
                                h(self.soma + ' ' +'distance(0,'+str(trunk_origin[0]) + ')') # Trunk origin point (reference for distance measurement) can be
                            elif len(trunk_origin) == 2:
                                h(trunk_origin[0] + ' ' +'distance(0,'+str(trunk_origin[1]) + ')') # Trunk origin point (reference for distance measurement) can be added by the user as an argument to the test
                            h('access ' + self.trunk[i].name())
                            if [self.trunk[i].name(), segment] not in locations and h.distance(segment) >= dist_range[0] and h.distance(segment) < dist_range[1]:
                                locations.append([self.trunk[i].name(), segment])
                                locations_distances[self.trunk[i].name(), segment] = h.distance(segment)
                _num_ = num - len(locations)
                #print '_num_', _num_
                seed += 10
                num_iterations += 1
                #print len(locations)
        #print 'Dendritic locations to be tested (with their actual distances):', locations_distances

        return locations, locations_distances

    def find_good_obliques(self, trunk_origin):
        """Used in ObliqueIntegrationTest"""

        if (self.ObliqueSecList_name is None or self.TrunkSecList_name is None) and not self.find_section_lists:
            raise NotImplementedError("Please give the names of the section lists containing the oblique dendrites and the trunk sections. (eg. model.ObliqueSecList_name=\"obliques\", model.TrunkSecList_name=\"trunk\" or set model.find_section_lists to True)")


        #self.initialise()

        good_obliques = h.SectionList()
        dend_loc=[]

        if self.TrunkSecList_name is not None and self.ObliqueSecList_name is not None:
            self.initialise()

            if self.template_name is not None:

                exec('self.oblique_dendrites=h.testcell.' + self.ObliqueSecList_name)   # so we can have the name of the section list as a string given by the user
                #exec('oblique_dendrites = h.' + oblique_seclist_name)
                exec('self.trunk=h.testcell.' + self.TrunkSecList_name)
            else:
                exec('self.oblique_dendrites=h.' + self.ObliqueSecList_name)   # so we can have the name of the section list as a string given by the user
                #exec('oblique_dendrites = h.' + oblique_seclist_name)
                exec('self.trunk=h.' + self.TrunkSecList_name)

        if self.find_section_lists:

            self.initialise()

            if self.template_name is not None:
                exec('self.icell=h.testcell')

            apical_trunk_isections, apical_tuft_isections, oblique_isections = self.classify_apical_point_sections(self.icell)

            self.trunk = []
            for i in range(len(apical_trunk_isections)):
                exec('self.sec = h.testcell.apic[' + str(apical_trunk_isections[i]) + ']')
                self.trunk.append(self.sec)

            self.oblique_dendrites = []
            for i in range(len(oblique_isections)):
                exec('self.sec = h.testcell.apic[' + str(oblique_isections[i]) + ']')
                self.oblique_dendrites.append(self.sec)

        good_obliques_added = 0

        while good_obliques_added == 0 and self.max_dist_from_soma <= 190:
            for sec in self.oblique_dendrites:
                if not trunk_origin:
                    h(self.soma + ' ' +'distance(0,1)') # For apical dendrites the default reference point is the end of the soma (point 1)
                elif len(trunk_origin) == 1:
                    h(self.soma + ' ' +'distance(0,'+str(trunk_origin[0]) + ')') # Trunk origin point (reference for distance measurement) can be
                elif len(trunk_origin) == 2:
                    h(trunk_origin[0] + ' ' +'distance(0,'+str(trunk_origin[1]) + ')') # Trunk origin point (reference for distance measurement) can be added by the user as an argument to the test
                if self.find_section_lists:
                    h('access ' + sec.name())
                parent = h.SectionRef(sec).parent
                child_num = h.SectionRef(sec).nchild()
                dist = h.distance(0, sec=sec)
                #print 'SEC: ', sec.name()
                #print 'NCHILD: ', child_num
                #print 'PARENT: ', parent.name()
                #print 'DIST: ', h.distance(0)
                """
                for trunk_sec in trunk:
                    if self.find_section_lists:
                        h('access ' + trunk_sec.name())
                    if h.issection(parent.name()) and dist < self.max_dist_from_soma and child_num == 0:   # true if string (parent.name()) is contained in the name of the currently accessed section.trunk_sec is the accessed section,
                        #print sec.name(), parent.name()
                        h('access ' + sec.name())         # only currently accessed section can be added to hoc SectionList
                        good_obliques.append(sec.name())
                        good_obliques_added += 1
                """
                if dist < self.max_dist_from_soma and child_num == 0:   # now the oblique section can branch from another oblique section, but it has to be a tip (terminal) section
                    #print(sec.name(), parent.name())
                    #print(sec.name(), dist)
                    h('access ' + sec.name())         # only currently accessed section can be added to hoc SectionList
                    good_obliques.append()
                    good_obliques_added += 1
            if good_obliques_added == 0:
                self.max_dist_from_soma += 15
                print("Maximum distance from soma was increased by 15 um, new value: " + str(self.max_dist_from_soma))

        for sec in good_obliques:

            dend_loc_prox=[]
            dend_loc_dist=[]
            seg_list_prox=[]
            seg_list_dist=[]

            h(sec.name() + ' ' +'distance()')  #set the 0 point of the section as the origin
            # print(sec.name())


            for seg in sec:
                # print(seg.x, h.distance(seg.x))
                if h.distance(seg.x, sec=sec) > 5 and h.distance(seg.x, sec=sec) < 50:
                    seg_list_prox.append(seg.x)
                if h.distance(seg.x, sec=sec) > 60 and h.distance(seg.x, sec=sec) < 126:
                    seg_list_dist.append(seg.x)

            #print seg_list_prox
            #print seg_list_dist

            if len(seg_list_prox) > 1:
                s = int(numpy.ceil(len(seg_list_prox)/2.0))
                dend_loc_prox.append(sec.name())
                dend_loc_prox.append(seg_list_prox[s])
                dend_loc_prox.append('prox')
            elif len(seg_list_prox) == 1:
                dend_loc_prox.append(sec.name())
                dend_loc_prox.append(seg_list_prox[0])
                dend_loc_prox.append('prox')

            if len(seg_list_dist) > 1:
                s = int(numpy.ceil(len(seg_list_dist)/2.0)-1)
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(seg_list_dist[s])
                dend_loc_dist.append('dist')
            elif len(seg_list_dist) == 1:
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(seg_list_dist[0])
                dend_loc_dist.append('dist')
            elif len(seg_list_dist) == 0:                # if the dendrite is not long enough to meet the criteria, we stimulate its end
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(0.9)
                dend_loc_dist.append('dist')

            if dend_loc_prox:
                dend_loc.append(dend_loc_prox)
            if dend_loc_dist:
                dend_loc.append(dend_loc_dist)

        #print 'Dendrites and locations to be tested: ', dend_loc

        return dend_loc


    def set_ampa_nmda(self, dend_loc):
        """Currently not used - Used to be used in ObliqueIntegrationTest"""

        ndend, xloc, loc_type = dend_loc

        exec("self.dendrite=h." + ndend)

        self.ampa = h.Exp2Syn(xloc, sec=self.dendrite)
        self.ampa.tau1 = self.AMPA_tau1
        self.ampa.tau2 = self.AMPA_tau2

        exec("self.nmda = h."+self.NMDA_name+"(xloc, sec=self.dendrite)")

        self.ndend = ndend
        self.xloc = xloc


    def set_netstim_netcon(self, interval):
        """Currently not used - Used to be used in ObliqueIntegrationTest"""

        self.ns = h.NetStim()
        self.ns.interval = interval
        self.ns.number = 0
        self.ns.start = self.start

        self.ampa_nc = h.NetCon(self.ns, self.ampa, 0, 0, 0)
        self.nmda_nc = h.NetCon(self.ns, self.nmda, 0, 0, 0)


    def set_num_weight(self, number, AMPA_weight):
        """Currently not used - Used to be used in ObliqueIntegrationTest"""

        self.ns.number = number
        self.ampa_nc.weight[0] = AMPA_weight
        self.nmda_nc.weight[0] =AMPA_weight/self.AMPA_NMDA_ratio

    def run_syn(self, dend_loc, interval, number, AMPA_weight):
        """Currently not used - Used to be used in ObliqueIntegrationTest"""

        self.initialise()

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        self.set_ampa_nmda(dend_loc)
        self.set_netstim_netcon(interval)
        self.set_num_weight(number, AMPA_weight)

        exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/ dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop = 500
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend

    def set_multiple_ampa_nmda(self, dend_loc, number):
        """Used in ObliqueIntegrationTest"""

        ndend, xloc, loc_type = dend_loc

        exec("self.dendrite=h." + ndend)

        for i in range(number):

            if self.AMPA_name: # if this is given, the AMPA model defined in a mod file is used, else the built in Exp2Syn
                exec("self.ampa_list[i] = h."+self.AMPA_name+"(xloc, sec=self.dendrite)")
            else:
                self.ampa_list[i] = h.Exp2Syn(xloc, sec=self.dendrite)
                self.ampa_list[i].tau1 = self.AMPA_tau1
                self.ampa_list[i].tau2 = self.AMPA_tau2
                #print 'The built in Exp2Syn is used as the AMPA component. Tau1 = ', self.AMPA_tau1, ', Tau2 = ', self.AMPA_tau2 , '.'

            if self.NMDA_name: # if this is given, the NMDA model defined in a mod file is used, else the default NMDA model of HippoUnit
                exec("self.nmda_list[i] = h."+self.NMDA_name+"(xloc, sec=self.dendrite)")
            else:
                try:
                    exec("self.nmda_list[i] = h."+self.default_NMDA_name+"(xloc, sec=self.dendrite)")
                except:
                    h.nrn_load_dll(self.default_NMDA_path + self.libpath)
                    exec("self.nmda_list[i] = h."+self.default_NMDA_name+"(xloc, sec=self.dendrite)")

        self.ndend = ndend
        self.xloc = xloc


    def set_multiple_netstim_netcon(self, interval, number, AMPA_weight):
        """Used in ObliqueIntegrationTest"""

        for i in range(number):
            self.ns_list[i] = h.NetStim()
            self.ns_list[i].number = 1
            self.ns_list[i].start = self.start + (i*interval)

            self.ampa_nc_list[i] = h.NetCon(self.ns_list[i], self.ampa_list[i], 0, 0, 0)
            self.nmda_nc_list[i] = h.NetCon(self.ns_list[i], self.nmda_list[i], 0, 0, 0)

            self.ampa_nc_list[i].weight[0] = AMPA_weight
            self.nmda_nc_list[i].weight[0] =AMPA_weight/self.AMPA_NMDA_ratio


    def run_multiple_syn(self, dend_loc, interval, number, weight):
        """Used in ObliqueIntegrationTest"""

        self.ampa_list = [None] * number
        self.nmda_list = [None] * number
        self.ns_list = [None] * number
        self.ampa_nc_list = [None] * number
        self.nmda_nc_list = [None] * number


        self.initialise()

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        self.set_multiple_ampa_nmda(dend_loc, number)

        self.set_multiple_netstim_netcon(interval, number, weight)


        exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop =500
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend



    def set_Exp2Syn(self, dend_loc, tau1, tau2):
        """Used in PSPAttenuationTest"""

        ndend, xloc = dend_loc

        exec("self.dendrite=h." + ndend)

        self.ampa = h.Exp2Syn(xloc, sec=self.dendrite)
        self.ampa.tau1 = tau1
        self.ampa.tau2 = tau2

        self.ndend = ndend
        self.xloc = xloc


    def set_netstim_netcon_Exp2Syn(self):
        """Used in PSPAttenuationTest"""
        self.start = 300

        self.ns = h.NetStim()
        #self.ns.interval = interval
        #self.ns.number = 0
        self.ns.start = self.start

        self.ampa_nc = h.NetCon(self.ns, self.ampa, 0, 0, 0)

    def set_weight_Exp2Syn(self, weight):
        """Used in PSPAttenuationTest"""

        self.ns.number = 1
        self.ampa_nc.weight[0] = weight

    def run_EPSCstim(self, dend_loc, weight, tau1, tau2):
        """Used in PSPAttenuationTest"""

        self.initialise()

        if self.cvode_active:
            h.cvode_active(1)
        else:
            h.cvode_active(0)

        self.set_Exp2Syn(dend_loc, tau1, tau2)
        self.set_netstim_netcon_Exp2Syn()
        self.set_weight_Exp2Syn(weight)

        exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        # initiate recording
        rec_t = h.Vector()
        rec_t.record(h._ref_t)

        rec_v = h.Vector()
        rec_v.record(self.sect_loc._ref_v)

        rec_v_dend = h.Vector()
        rec_v_dend.record(self.dendrite(self.xloc)._ref_v)

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1/dt
        h.v_init = self.v_init #-80

        h.celsius = self.celsius
        h.init()
        h.tstop = 450
        h.run()

        # get recordings
        t = numpy.array(rec_t)
        v = numpy.array(rec_v)
        v_dend = numpy.array(rec_v_dend)

        return t, v, v_dend

class ModelLoader_BPO(ModelLoader):

    def __init__(self, name="model", model_dir=None, SomaSecList_name=None):
        """ Constructor. """
        """ This class should be used with Jupyter notebooks"""
        super(ModelLoader_BPO, self).__init__(name=name)
        self.SomaSecList_name = SomaSecList_name
        self.morph_full_path = None
        self.find_section_lists = True

        self.setup_dirs(model_dir)
        self.setup_values()
        self.compile_mod_files_BPO()

    def compile_mod_files(self):
        """This method is called by the parent class (ModelLoader), but as the path to the mod files is unknown at this point, this is not used. compile_mode_files_BPO is used instead to compile the mod files."""
        pass

    def compile_mod_files_BPO(self):

        if self.modelpath is None:
            raise Exception("Please give the path to the mod files (eg. model.modelpath = \"/home/models/CA1_pyr/mechanisms/\")")

        if os.path.isfile(self.modelpath + self.libpath) is False:
            os.system("cd " + self.modelpath + "; nrnivmodl")

    def load_mod_files(self):

        h.nrn_load_dll(str(self.modelpath + self.libpath))

    def setup_dirs(self, model_dir=""):

        '''
        split_dir = model_dir.split('/')
        del split_dir[-1]
        outer_dir = '/'.join(split_dir)

        if not os.path.exists(model_dir):
            try:

                #split_dir = model_dir.split('/')
                #del split_dir[-1]
                #outer_dir = '/'.join(split_dir)

                zip_ref = zipfile.ZipFile(model_dir + '.zip', 'r')
                zip_ref.extractall(outer_dir)
            except IOError:
                print "Error accessing directory/zipfile named: ", model_dir
        '''

        base_path = os.path.join(model_dir, self.name)
        if os.path.exists(base_path) or os.path.exists(base_path+".zip"):     # If the model_dir is the outer directory, that contains the zip
            self.base_path = base_path
            if not os.path.exists(self.base_path):
                file_ref = zipfile.ZipFile(self.base_path+".zip", 'r')
                file_ref.extractall(model_dir)
                file_ref.close()
            try:
                with open(self.base_path + '/' + self.name + '_meta.json') as f:
                    meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)
            except Exception as e1:
                try:
                    with open(model_dir + '/' + self.name + '_meta.json') as f:
                        meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)
                except Exception as e2:
                    print(e1,e2)
        else:                                                                   # If model_dir is the inner directory (already unzipped)
            self.base_path = model_dir
            split_dir = model_dir.split('/')
            del split_dir[-1]
            outer_dir = '/'.join(split_dir)

            try:
                with open(self.base_path + '/' + self.name + '_meta.json') as f:
                    meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)
            except Exception as e1:
                try:
                    with open(outer_dir + '/' + self.name + '_meta.json') as f:
                        meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)
                except Exception as e2:
                    print(e1,e2)
        '''
        try:
            with open(self.base_path + '/' + self.name + '_meta.json') as f:
                meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)
        except:
            try:
                with open(model_dir + '/' + self.name + '_meta.json') as f:
                    meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)
            except Exception as e:
                print e
        '''

        self.morph_path = "\"" + self.base_path + "/morphology\""

        for file_name in os.listdir(self.morph_path[1:-1]):
            self.morph_full_path = self.morph_path[1:-1]+ '/' + file_name
            break


        # path to mod files
        self.modelpath = self.base_path + "/mechanisms/"

        # if this doesn't exist mod files are automatically compiled
        self.libpath = "x86_64/.libs/libnrnmech.so"

        best_cell = meta_data["best_cell"]

        self.hocpath = self.base_path + "/checkpoints/" + str(best_cell)

        if not os.path.exists(self.hocpath):
            self.hocpath = None
            for file in os.listdir(self.base_path + "/checkpoints/"):
                if file.startswith("cell") and file.endswith(".hoc"):
                    self.hocpath = self.base_path + "/checkpoints/" + file
                    print("Model = " + self.name + ": cell.hoc not found in /checkpoints; using " + file)
                    break
            if not os.path.exists(self.hocpath):
                raise IOError("No appropriate .hoc file found in /checkpoints")

        self.base_directory = self.base_path +'/validation_results/'

    def setup_values(self):

        # get model template name
        # could also do this via other JSON, but morph.json seems dedicated for template info
        with open(os.path.join(self.base_path, "config", "morph.json")) as morph_file:
            template_name = list(json.load(morph_file, object_pairs_hook=collections.OrderedDict).keys())[0]

        self.template_name = template_name + "(" + self.morph_path+")"

        # access model config info
        with open(os.path.join(self.base_path, "config", "parameters.json")) as params_file:
            params_data = json.load(params_file, object_pairs_hook=collections.OrderedDict)

        # extract v_init and celsius (if available)
        v_init = None
        celsius = None
        try:
            for item in params_data[template_name]["fixed"]["global"]:
                # would have been better if info was stored inside a dict (rather than a list)
                if "v_init" in item:
                    item.remove("v_init")
                    v_init = float(item[0])
                if "celsius" in item:
                    item.remove("celsius")
                    celsius = float(item[0])
        except:
            pass
        if v_init == None:
            self.v_init = -70.0
            print("Could not find model specific info for `v_init`; using default value of {} mV".format(str(self.v_init)))
        else:
            self.v_init = v_init
        if celsius == None:
            self.celsius = 34.0
            print("Could not find model specific info for `celsius`; using default value of {} degrees Celsius".format(str(self.celsius)))
        else:
            self.celsius = celsius
        self.trunk_origin = [0.5]
