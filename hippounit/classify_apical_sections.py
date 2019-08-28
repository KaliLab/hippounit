from builtins import range
import neurom as nm
from neurom.core.dataformat import COLS
from neurom.morphmath import point_dist2
from neurom.viewer import draw
from neuron import h
import neuron
import numpy
import matplotlib.pyplot as plt


""" Authors: Michael Gevaert, Armando Romani, Christian Rossert, Sara Saray"""

def get_apical_point(neurite, morph, tuft_percent=27):
    '''Attempt to find the apical point in 'tufted' neurons
    Consider a neuron:
        |   /    | Tuft = 20%
        |--/     |
        |   /
        |--/
        |
    ----.-----
    All endpoints in the top 'tuft_percent' are found, then their common
    branch segment, furthest from the soma, is identified.
    Args:
        morph: neurom.fst._core.Neuron
        tuft_percent: percentage of the 'height' of the apical dendrite that
        would enclose the tuft, only leaves in this volume are considered as
        endpoints.  Note that this a spherical shell centered at the soma
    Returns:
        Section whose *end* is where the apical branching begins, or None if there
        is a problem
    '''

    apical = neurite
    max_distance2 = -float('inf')

    for leaf in apical.root_node.ileaf():
        point = leaf.points[-1, COLS.XYZ]
        max_distance2 = max(max_distance2,
                            point_dist2(point, morph.soma.center))


    min_distance2 = max_distance2 * (1 - tuft_percent / 100.) ** 2

    common_parents = set(nm.iter_sections(apical))
    #Iterator to the sections in a neurite, neuron or neuron population.
    all_parents = set([])

    for leaf in apical.root_node.ileaf():
        point = leaf.points[-1, COLS.XYZ]
        if min_distance2 <= point_dist2(point, morph.soma.center):
            parents = leaf.iupstream()
            set_parents =  set(parents)
            common_parents &= set_parents
            all_parents |= set_parents

    apical_point_section = None
    for parent_section in nm.iter_sections(apical):
        if parent_section in common_parents:
            common_parents.remove(parent_section)
            if not common_parents:
                #print parent_section
                if parent_section in all_parents:
                    #return parent_section
                    apical_point_section = parent_section
                else:
                    apical_point_section = None

    #return None

    return apical_point_section

def multiple_apical_points(morphology):

    morph = morphology

    apical = [neurite for neurite in morph.neurites
       if nm.NeuriteType.apical_dendrite == neurite.type]


    if not apical:
        L.warning('No apical found')
        return None
    elif 1 < len(apical):
        L.warning('Too many apical dendrites')
        return None

    apical = apical[0]

    dist_apical_point = []
    apical_points_and_distances = []

    closer_apical_points =[]

    point = get_apical_point(apical, morph)

    if point is not None:
        dist = nm.morphmath.point_dist(morph.soma.center, point.points[-1, COLS.XYZ])

        apical_points_and_distances.append([point, dist])


    while any(point[1] < 200 for point in apical_points_and_distances):
        #point = get_apical_point(apical, morph)
        #apical_point.append(point)
        #dist_apical_point.append(nm.morphmath.point_dist(morph.soma.center, point.points[-1, COLS.XYZ]))
        #a, b = min(point[1] for idx, point in enumerate(apical_points_and_distances))
        mn,idx = min( (apical_points_and_distances[i][1],i) for i in range(len(apical_points_and_distances)) )
        #print mn, idx
        #print apical_points_and_distances
        new_apical_points = []
        new_apical_point_added = []

        for child in apical_points_and_distances[idx][0].children:

            root_node = child
            tree = nm.core._neuron.Neurite(root_node)
            point = get_apical_point(tree, morph)
            if point is not None:

                dist = nm.morphmath.point_dist(morph.soma.center, point.points[-1, COLS.XYZ])

                new_apical_point_added.append(True)

                apical_points_and_distances.append([point, dist])
            else:
                new_apical_point_added.append(False)

        """ avoid infinite loop"""
        #print new_apical_point_added
        if any(added == True for added in new_apical_point_added):
            del apical_points_and_distances[idx]
        else:
            closer_apical_points.append(apical_points_and_distances[idx])
            del apical_points_and_distances[idx]

    apical_points_and_distances = apical_points_and_distances + closer_apical_points
    #print apical_points_and_distances

    apical_points = []

    for point in apical_points_and_distances:
        apical_points.append(point[0])

    _, ax = draw(morph)
    for point in apical_points:
        pos = point.points[-1, COLS.XYZ]
        ax.scatter([pos[0], ], [pos[1], ])
    plt.show()

    return apical_points

def get_list_of_diff_section_types(morphology, apical_point_sections):

    morph = morphology

    _, ax = draw(morph)

    apical_trunk_sections = []

    for apical_point_section in apical_point_sections:
        some_trunk_sections = list(apical_point_section.iupstream())
        apical_trunk_sections = list(set(apical_trunk_sections) | set(some_trunk_sections))

    apical_trunk_ids = [section.id for section in apical_trunk_sections]
    apical_trunk_points = [section.points[-1, :2] for section in apical_trunk_sections]
    for point in apical_trunk_points:
            ax.plot(point[0], point[1], 'bo', ms=5)

    apical_tuft_sections = []

    for apical_point_section in apical_point_sections:
        some_tuft_sections = list(apical_point_section.ipreorder())
        #print some_tuft_sections
        apical_tuft_sections = apical_tuft_sections + some_tuft_sections
    """!!!!!!"""
    apical_tuft_sections = [sec for sec in apical_tuft_sections if sec not in apical_point_sections]
    """!!!!!!"""

    apical_tuft_ids = [section.id for section in apical_tuft_sections]
    apical_tuft_points = [section.points[-1, :2] for section in apical_tuft_sections]
    for point in apical_tuft_points:
            ax.plot(point[0], point[1], 'go', ms=5)

    filter1 = lambda n : n.type == nm.APICAL_DENDRITE
    all_sections = list(nm.iter_sections(morph, neurite_filter=filter1))
    all_ids = [section.id for section in all_sections]
    oblique_ids = set(all_ids) - (set(apical_trunk_ids) | set(apical_tuft_ids))
    oblique_sections = [section for section in all_sections if section.id in oblique_ids]
    oblique_points = [section.points[-1, :2] for section in oblique_sections]
    for point in oblique_points:
            ax.plot(point[0], point[1], 'yo', ms=5)

    plt.show()

    section_types = {'tuft' : apical_tuft_sections, 'trunk' : apical_trunk_sections, 'obliques' : oblique_sections}

    return section_types

def get_neuron_isections(icell, list_neurom_sections):
    #icell=h.CCell
    #seclist=icell.apical

    seclist=list(icell.apical)


    xyz = numpy.zeros((len(seclist),3))

    for i, sec in enumerate(seclist):
        #print sec.name()
        iend = int(neuron.h.n3d(sec=sec))-1 # apical point always at endpoint of section
        x = neuron.h.x3d(iend, sec=sec)
        y = neuron.h.y3d(iend, sec=sec)
        z = neuron.h.z3d(iend, sec=sec)
        xyz[i] = [x,y,z]

    neuron_section_list = []

    for section in list_neurom_sections:
        #print section.id

        coord_xyz = section.points[-1, COLS.XYZ]
        #print coord_xyz

        d = (xyz-coord_xyz)**2 # search for minimal distance here
        isec = numpy.argmin(numpy.sum(d, axis=1))
        #print "ISEC", isec

        neuron_section_list.append(isec)

    return neuron_section_list
