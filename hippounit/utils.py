import os
import numpy
import sciunit
from neuronunit.capabilities import ReceivesCurrent, ProducesMembranePotential
import hippounit.capabilities as cap
from quantities import ms,mV,Hz
from neuron import h

import multiprocessing
import zipfile
import collections

class ModelLoader(sciunit.Model,
                 cap.ProvidesGoodObliques,
                 cap.ReceivesSquareCurrent_ProvidesResponse,
                 cap.ReceivesSynapse,
                 cap.ReceivesSquareCurrent_ProvidesResponse_MultipleLocations,
                 cap.ProvidesRecordingLocationsOnTrunk):

    def __init__(self, name="model"):
        """ Constructor. """

        """ This class should be used with Jupyter notebooks"""

        self.modelpath = None
        self.libpath = None
        self.hocpath = None
        self.template_name = None
        self.SomaSecList_name = None
        self.max_dist_from_soma = None
        self.v_init = None

        self.name = name
        self.threshold = -20
        self.stim = None
        self.soma = None
        sciunit.Model.__init__(self, name=name)

        self.c_step_start = 0.00004
        self.c_step_stop = 0.000004
        self.c_minmax = numpy.array([0.00004, 0.004])

        self.ObliqueSecList_name = None
        self.TrunkSecList_name = None
        self.dend_loc = []  #self.dend_loc = [['dendrite[80]',0.27],['dendrite[80]',0.83],['dendrite[54]',0.16],['dendrite[54]',0.95],['dendrite[52]',0.38],['dendrite[52]',0.83],['dendrite[53]',0.17],['dendrite[53]',0.7],['dendrite[28]',0.35],['dendrite[28]',0.78]]
        self.dend_locations = collections.OrderedDict()
        self.NMDA_name = 'NMDA_JS'
        self.AMPA_NMDA_ratio = 0.4

        self.AMPA_tau1 = 0.1
        self.AMPA_tau2 = 2
        self.start=150

        self.ns = None
        self.ampa = None
        self.nmda = None
        self.ampa_nc = None
        self.nmda_nc = None

        self.ndend = None
        self.xloc = None

        self.base_directory = './validation_results/'   # inside current directory

    def translate(self, sectiontype, distance=0):

        if "soma" in sectiontype:
            return self.soma
        else:
            return False

    def load_mod_files(self):
        if os.path.isfile(self.modelpath + self.libpath) is False:
            os.system("cd " + self.modelpath + "; nrnivmodl")

        h.nrn_load_dll(self.modelpath + self.libpath)

    def initialise(self):

        self.load_mod_files()

        h.load_file("stdrun.hoc")
        h.load_file(self.hocpath)


        if self.template_name is not None and self.SomaSecList_name is not None:

            h('objref testcell')
            h('testcell = new ' + self.template_name)

            exec('soma = h.testcell.'+ self.SomaSecList_name)

            for s in soma :
                self.soma = h.secname()

        elif self.template_name is not None and self.SomaSecList_name is None:
            h('objref testcell')
            h('testcell = new ' + self.template_name)
            # in this case self.soma is set in the jupyter notebook
        elif self.template_name is None and self.SomaSecList_name is not None:
            exec('soma = h.' +  self.SomaSecList_name)
            for s in soma :
                self.soma = h.secname()
        # if both is None, the model is loaded, self.soma will be used

    def inject_current(self, amp, delay, dur, section_stim, loc_stim, section_rec, loc_rec):

        self.initialise()

        stim_section_name = self.translate(section_stim, distance=0)
        rec_section_name = self.translate(section_rec, distance=0)
        #exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        exec("self.sect_loc_stim=h." + str(stim_section_name)+"("+str(loc_stim)+")")

        print "- running amplitude: " + str(amp)  + " on model: " + self.name + " at: " + stim_section_name + "(" + str(loc_stim) + ")"

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
        h.steps_per_ms = 1 / dt
        h.v_init = self.v_init#-65

        h.celsius = 34
        h.init()
        h.tstop = delay + dur + 200
        h.run()

        t = numpy.array(rec_t)
        v = numpy.array(rec_v)

        return t, v

    def inject_current_record_respons_multiple_loc(self, amp, delay, dur, section_stim, loc_stim, dend_locations):

        self.initialise()


        stim_section_name = self.translate(section_stim, distance=0)
        #rec_section_name = self.translate(section_rec, distance=0)
        #exec("self.sect_loc=h." + str(self.soma)+"("+str(0.5)+")")

        exec("self.sect_loc_stim=h." + str(stim_section_name)+"("+str(loc_stim)+")")
        exec("self.sect_loc_rec=h." + str(stim_section_name)+"("+str(loc_stim)+")")

        print "- running amplitude: " + str(amp)  + " on model: " + self.name + " at: " + stim_section_name + "(" + str(loc_stim) + ")"

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
        for key, value in dend_locations.iteritems():
            for i in range(len(dend_locations[key])):
                exec("self.dend_loc_rec.append(h." + str(dend_locations[key][i][0])+"("+str(dend_locations[key][i][1])+"))")
                rec_v.append(h.Vector())

        for i in range(len(self.dend_loc_rec)):
            rec_v[i].record(self.dend_loc_rec[i]._ref_v)
            #print self.dend_loc[i]

        h.stdinit()

        dt = 0.025
        h.dt = dt
        h.steps_per_ms = 1 / dt
        h.v_init = self.v_init#-65

        h.celsius = 34
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
        for key, value in dend_locations.iteritems():
            v[key] = collections.OrderedDict()
            for j in range(len(dend_locations[key])):
                loc_key = (dend_locations[key][j][0],dend_locations[key][j][1]) # list can not be a key, but tuple can
                v[key][loc_key] = numpy.array(rec_v[i])     # the list that specifies dendritic location will be a key too.
                i+=1

        return t, v_stim, v


    def find_trunk_locations(self, distances):

        self.initialise()

        #locations={}
        locations=collections.OrderedDict()
        actual_distances ={}
        dend_loc=[]

        if self.template_name is not None:
            exec('trunk=h.testcell.' + self.TrunkSecList_name)
        else:
            exec('trunk=h.' + self.TrunkSecList_name)

        for sec in trunk:
            #for seg in sec:
            h(self.soma + ' ' +'distance()') #set soma as the origin
            for seg in sec:
                #print 'SEC: ', sec.name(),
                #print 'SEG.X', seg.x
                #print 'DIST', h.distance(seg.x)
                #print 'DIST0', h.distance(0)
                #print 'DIST1', h.distance(1)
                for i in range(0, len(distances)):
                    locations.setdefault(distances[i], []) # if this key doesn't exist it is added with the value: [], if exists, value not altered
                    if h.distance(seg.x) < distances[i] +20 and h.distance(seg.x) > distances[i]-20: # if the seq is between distance +- 20
                        #print 'SEC: ', sec.name()
                        #print 'seg.x: ', seg.x
                        #print 'DIST: ', h.distance(seg.x)
                        locations[distances[i]].append([sec.name(), seg.x])
                        actual_distances[sec.name(), seg.x] = h.distance(seg.x)

        return locations, actual_distances


    def find_good_obliques(self):
        """Used in ObliqueIntegrationTest"""


        self.initialise()

        good_obliques = h.SectionList()
        dend_loc=[]

        if self.template_name is not None:

            exec('oblique_dendrites=h.testcell.' + self.ObliqueSecList_name)   # so we can have the name of the section list as a string given by the user
            #exec('oblique_dendrites = h.' + oblique_seclist_name)
            exec('trunk=h.testcell.' + self.TrunkSecList_name)
        else:
            exec('oblique_dendrites=h.' + self.ObliqueSecList_name)   # so we can have the name of the section list as a string given by the user
            #exec('oblique_dendrites = h.' + oblique_seclist_name)
            exec('trunk=h.' + self.TrunkSecList_name)

        for sec in oblique_dendrites:
            h(self.soma + ' ' +'distance()') #set soma as the origin
            parent = h.SectionRef(sec).parent
            child_num = h.SectionRef(sec).nchild()
            dist = h.distance(0)
            #print 'SEC: ', sec.name()
            #print 'NCHILD: ', child_num
            #print 'PARENT: ', parent.name()
            #print 'DIST: ', h.distance(0)

            for trunk_sec in trunk:
                if h.issection(parent.name()) and dist < self.max_dist_from_soma and child_num == 0:   # true if string (parent.name()) is contained in the name of the currently accessed section.trunk_sec is the accessed section,
                    #print sec.name(), parent.name()
                    h('access ' + sec.name())         # only currently accessed section can be added to hoc SectionList
                    good_obliques.append(sec.name())


        for sec in good_obliques:

            dend_loc_prox=[]
            dend_loc_dist=[]
            seg_list_prox=[]
            seg_list_dist=[]

            h(sec.name() + ' ' +'distance()')  #set the 0 point of the section as the origin
            #print sec.name()


            for seg in sec:
                #print seg.x, h.distance(seg.x)
                if h.distance(seg.x) > 5 and h.distance(seg.x) < 50:
                    seg_list_prox.append(seg.x)
                if h.distance(seg.x) > 60 and h.distance(seg.x) < 126:
                    seg_list_dist.append(seg.x)

            #print seg_list_prox
            #print seg_list_dist

            if len(seg_list_prox) > 1:
                s = int(round(len(seg_list_prox)/2.0))
                dend_loc_prox.append(sec.name())
                dend_loc_prox.append(seg_list_prox[s])
            else:
                dend_loc_prox.append(sec.name())
                dend_loc_prox.append(seg_list_prox[0])

            if len(seg_list_dist) > 1:
                s = int(round(len(seg_list_dist)/2.0)-1)
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(seg_list_dist[s])
            elif len(seg_list_dist) == 1:
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(seg_list_dist[0])
            elif len(seg_list_dist) == 0:                # if the dendrite is not long enough to meet the criteria, we stimulate its end
                dend_loc_dist.append(sec.name())
                dend_loc_dist.append(0.9)


            dend_loc.append(dend_loc_prox)
            dend_loc.append(dend_loc_dist)

        print 'Dendrites and locations to be tested: ', dend_loc
        return dend_loc


    def set_ampa_nmda(self, dend_loc):
        """Used in ObliqueIntegrationTest"""

        ndend, xloc = dend_loc

        exec("self.dendrite=h." + ndend)

        self.ampa = h.Exp2Syn(xloc, sec=self.dendrite)
        self.ampa.tau1 = self.AMPA_tau1
        self.ampa.tau2 = self.AMPA_tau2

        exec("self.nmda = h."+self.NMDA_name+"(xloc, sec=self.dendrite)")

        self.ndend = ndend
        self.xloc = xloc


    def set_netstim_netcon(self, interval):
        """Used in ObliqueIntegrationTest"""

        self.ns = h.NetStim()
        self.ns.interval = interval
        self.ns.number = 0
        self.ns.start = self.start

        self.ampa_nc = h.NetCon(self.ns, self.ampa, 0, 0, 0)
        self.nmda_nc = h.NetCon(self.ns, self.nmda, 0, 0, 0)


    def set_num_weight(self, number, AMPA_weight):
        """Used in ObliqueIntegrationTest"""

        self.ns.number = number
        self.ampa_nc.weight[0] = AMPA_weight
        self.nmda_nc.weight[0] =AMPA_weight/self.AMPA_NMDA_ratio

    def run_syn(self, dend_loc, interval, number, AMPA_weight):
        """Used in ObliqueIntegrationTest"""

        self.initialise()
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
        h.steps_per_ms = 1 / dt
        h.v_init = self.v_init #-80

        h.celsius = 34
        h.init()
        h.tstop = 300
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
        self.setup_dirs(model_dir)

    def setup_dirs(self, model_dir=""):
        if not os.path.exists(model_dir):
            try:
                zip_ref = zipfile.ZipFile(model_dir + '.zip', 'r')
                zip_ref.extractall(dir)
            except IOError:
                print "Error accessing directory/zipfile named: ", model_dir

        self.morph_path = "\"" + model_dir + "/morphology\""

        # path to mod files
        self.modelpath = model_dir + "/mechanisms/"

        # if this doesn't exist mod files are automatically compiled
        self.libpath = "x86_64/.libs/libnrnmech.so.0"
        self.hocpath = model_dir + "/checkpoints/cell.hoc"
        if not os.path.exists(self.hocpath):
            self.hocpath = None
            for file in os.listdir(model_dir + "/checkpoints/"):
                if file.startswith("cell") and file.endswith(".hoc"):
                    self.hocpath = model_dir + "/checkpoints/" + file
                    print "Model = " + self.name + ": cell.hoc not found in /checkpoints; using " + file
                    break
            if not os.path.exists(self.hocpath):
                raise IOError("No appropriate .hoc file found in /checkpoints")

        self.base_directory = model_dir +'/validation_results/'
