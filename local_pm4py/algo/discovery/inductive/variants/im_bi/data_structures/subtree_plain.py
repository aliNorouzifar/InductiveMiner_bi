'''
    This file is part of PM4Py (More Info: https://pm4py.fit.fraunhofer.de).

    PM4Py is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PM4Py is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PM4Py.  If not, see <https://www.gnu.org/licenses/>.
'''
from copy import copy
import pkgutil
import time

from pm4py.algo.discovery.dfg.utils.dfg_utils import get_activities_from_dfg, \
    infer_start_activities, infer_end_activities
from pm4py.algo.discovery.dfg.utils.dfg_utils import get_ingoing_edges, get_outgoing_edges
from pm4py.algo.discovery.dfg.utils.dfg_utils import negate, get_activities_self_loop, transform_dfg_to_directed_nx_graph
from pm4py.algo.discovery.dfg.variants import native as dfg_inst
from pm4py.algo.filtering.dfg.dfg_filtering import clean_dfg_based_on_noise_thresh
from local_pm4py.algo.discovery.inductive.variants.im.util import base_case, fall_through
from pm4py import util as pmutil
from local_pm4py.algo.discovery.inductive.variants.im_bi.util import splitting as split
from pm4py.algo.discovery.inductive.util import parallel_cut_utils, detection_utils, cut_detection
from pm4py.statistics.attributes.log import get as attributes_get
from pm4py.statistics.end_activities.log import get as end_activities_get
from pm4py.statistics.start_activities.log import get as start_activities_get
from pm4py.util import exec_utils
from pm4py.objects.log.util import filtering_utils
import logging
from pm4py.util import constants
from enum import Enum
from pm4py.objects.log import obj as log_instance
from pm4py.util import xes_constants
import datetime
from local_pm4py.algo.discovery.dfg import algorithm as dfg_discovery
import networkx as nx
import re
from pm4py.algo.filtering.log.start_activities import start_activities_filter
from pm4py.algo.filtering.log.end_activities import end_activities_filter
from pm4py.algo.conformance.alignments.dfg import algorithm as dfg_alignment
from pm4py.algo.evaluation.replay_fitness import algorithm as replay_fitness
from pm4py.visualization.dfg import visualizer as dfg_visualization
import pandas as pd
from pm4py.objects.dfg.filtering import dfg_filtering


class Parameters(Enum):
    ACTIVITY_KEY = constants.PARAMETER_CONSTANT_ACTIVITY_KEY
    START_TIMESTAMP_KEY = constants.PARAMETER_CONSTANT_START_TIMESTAMP_KEY
    TIMESTAMP_KEY = constants.PARAMETER_CONSTANT_TIMESTAMP_KEY
    CASE_ID_KEY = constants.PARAMETER_CONSTANT_CASEID_KEY
    NOISE_THRESHOLD = "noiseThreshold"
    EMPTY_TRACE_KEY = "empty_trace"
    ONCE_PER_TRACE_KEY = "once_per_trace"
    CONCURRENT_KEY = "concurrent"
    STRICT_TAU_LOOP_KEY = "strict_tau_loop"
    TAU_LOOP_KEY = "tau_loop"


from itertools import chain, combinations

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def artificial_start_end(log):
    st = 'start'
    en = 'end'
    activity_key = xes_constants.DEFAULT_NAME_KEY
    timestamp_key = xes_constants.DEFAULT_TIMESTAMP_KEY
    start_event = log_instance.Event()
    start_event[activity_key] = st

    end_event = log_instance.Event()
    end_event[activity_key] = en

    for trace in log:
        # start_event[timestamp_key] = trace[0]['time:timestamp'] - datetime.timedelta(microseconds=10)
        trace.insert(0, start_event)
        # end_event[timestamp_key] = trace[-1]['time:timestamp'] + datetime.timedelta(microseconds=10)
        trace.append(end_event)
    return log

def remove_empty(log):
    import copy
    new_log = copy.deepcopy(log)
    for tr in new_log:
        if len(tr)==0:
            new_log.remove(tr)
    return new_log

def generate_nx_graph_from_dfg(dfg):
    dfg_acts = set()
    for x in dfg:
        dfg_acts.add(x[0])
        dfg_acts.add(x[1])
    G = nx.DiGraph()
    for act in dfg_acts:
        G.add_node(act)
    for edge in dfg:
        G.add_edge(edge[0], edge[1])
    return G


def filter_dfg_on_threshold(dfg,inf_thr):
    filtered_dfg = {}
    for element in dfg:
        threshold = max([dfg[x] for x in dfg if x[0] == element[0]]) * inf_thr
        if dfg[element] == threshold or dfg[element] > threshold:
            filtered_dfg[element] = dfg[element]
    return filtered_dfg

class SubtreePlain(object):
    def __init__(self, logp,logm, dfg, master_dfg, initial_dfg, activities, counts, rec_depth, noise_threshold=0,
                 start_activities=None, end_activities=None, initial_start_activities=None,
                 initial_end_activities=None, parameters=None, real_init=True, sup= None, ratio = None, size_par = None):
        """
        Constructor

        Parameters
        -----------
        dfg
            Directly follows graph of this subtree
        master_dfg
            Original DFG
        initial_dfg
            Referral directly follows graph that should be taken in account adding hidden/loop transitions
        activities
            Activities of this subtree
        counts
            Shared variable
        rec_depth
            Current recursion depth
        """
        if real_init:
            self.master_dfg = copy(master_dfg)
            self.initial_dfg = copy(initial_dfg)
            self.counts = counts
            self.rec_depth = rec_depth
            self.noise_threshold = noise_threshold
            self.start_activities = start_activities_filter.get_start_activities(logp)
            self.start_activitiesM = start_activities_filter.get_start_activities(logm)
            self.end_activities = end_activities_filter.get_end_activities(logp)
            self.end_activitiesM = end_activities_filter.get_end_activities(logm)
            self.initial_start_activities = initial_start_activities
            if self.initial_start_activities is None:
                self.initial_start_activities = infer_start_activities(master_dfg)
            self.initial_end_activities = initial_end_activities
            if self.initial_end_activities is None:
                self.initial_end_activities = infer_end_activities(master_dfg)

            self.second_iteration = None
            self.activities = None
            self.activitiesM = None
            self.dfg = None
            self.outgoing = None
            self.ingoing = None
            self.self_loop_activities = None
            self.initial_ingoing = None
            self.initial_outgoing = None
            self.activities_direction = None
            self.activities_dir_list = None
            self.negated_dfg = None
            self.negated_activities = None
            self.negated_outgoing = None
            self.negated_ingoing = None
            self.detected_cut = None
            self.children = None
            self.must_insert_skip = False
            self.log = logp
            self.log_art = artificial_start_end(logp.__deepcopy__())
            self.logM = logm
            self.logM_art = artificial_start_end(logm.__deepcopy__())
            self.inverted_dfg = None
            self.original_log = logp

            self.initialize_tree(dfg, logp,logm, initial_dfg, activities, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)

    def __deepcopy__(self, memodict={}):
        """
            def __init__(self, log, dfg, master_dfg, initial_dfg, activities, counts, rec_depth, noise_threshold=0,
                 start_activities=None, end_activities=None, initial_start_activities=None,
                 initial_end_activities=None, parameters=None, real_init=False):
        :param memodict:
        :return:
        """
        S = SubtreePlain(None,None, None, None, None, None, None, None, real_init=False, sup= None, ratio = None, size_par = None)
        S.master_dfg = self.master_dfg
        S.initial_dfg = self.initial_dfg
        S.counts = self.counts
        S.rec_depth = self.rec_depth
        S.noise_threshold = self.noise_threshold
        S.start_activities = self.start_activities
        S.start_activitiesM = self.start_activitiesM
        S.end_activities = self.end_activities
        S.end_activitiesM = self.end_activitiesM
        S.initial_start_activities = self.initial_start_activities
        S.initial_end_activities = self.initial_end_activities
        S.second_iteration = self.second_iteration
        S.activities = self.activities
        S.activitiesM = self.activitiesM
        S.dfg = self.dfg
        S.outgoing = self.outgoing
        S.ingoing = self.ingoing
        S.self_loop_activities = self.self_loop_activities
        S.initial_ingoing = self.initial_ingoing
        S.initial_outgoing = self.initial_outgoing
        S.activities_direction = self.activities_direction
        S.activities_dir_list = self.activities_dir_list
        S.negated_dfg = self.negated_dfg
        S.negated_activities = self.negated_activities
        S.negated_outgoing = self.negated_outgoing
        S.negated_ingoing = self.negated_ingoing
        S.detected_cut = self.detected_cut
        S.children = self.children
        S.must_insert_skip = self.must_insert_skip
        S.log = self.log
        S.log_art = self.log_art
        S.logM = self.logM
        S.logM_art = self.logM_art
        S.inverted_dfg = self.inverted_dfg
        S.original_log = self.original_log
        try:
            S.parameters = self.parameters
        except:
            pass
        return S

    def initialize_tree(self, dfg, logp,logm, initial_dfg, activities, second_iteration=False, end_call=True,
                        parameters=None, sup= None, ratio = None, size_par = None):
        """
            Initialize the tree


            Parameters
            -----------
            dfg
                Directly follows graph of this subtree
            log
                the event log
            initial_dfg
                Referral directly follows graph that should be taken in account adding hidden/loop transitions
            activities
                Activities of this subtree
            second_iteration
                Boolean that indicates if we are executing this method for the second time
            """

        self.second_iteration = second_iteration

        if activities is None:
            self.activities = get_activities_from_dfg(dfg)
        else:
            self.activities = copy(activities)

        if second_iteration:
            self.dfg = clean_dfg_based_on_noise_thresh(self.dfg, self.activities, self.noise_threshold)
        else:
            self.dfg = copy(dfg)

        self.initial_dfg = initial_dfg

        self.outgoing = get_outgoing_edges(self.dfg)
        self.ingoing = get_ingoing_edges(self.dfg)
        self.self_loop_activities = get_activities_self_loop(self.dfg)
        self.initial_outgoing = get_outgoing_edges(self.initial_dfg)
        self.initial_ingoing = get_ingoing_edges(self.initial_dfg)
        self.negated_dfg = negate(self.dfg)
        self.negated_activities = get_activities_from_dfg(self.negated_dfg)
        self.negated_outgoing = get_outgoing_edges(self.negated_dfg)
        self.negated_ingoing = get_ingoing_edges(self.negated_dfg)
        self.detected_cut = None
        self.children = []
        self.log = logp
        self.log_art = artificial_start_end(logp.__deepcopy__())
        self.logM = logm
        self.logM_art = artificial_start_end(logm.__deepcopy__())
        self.original_log = logp
        self.parameters = parameters

        self.detect_cut(second_iteration=False, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)

    def create_dfg(self, parameters=None):
        if parameters is None:
            parameters = {}

        dfg = [(k, v) for k, v in dfg_inst.apply(self.log, parameters=parameters).items() if v > 0]

        return dfg


    def contains_empty_trace(self):
        contains = False
        for trace in self.log:
            if len(trace) == 0:
                contains = True
        return contains


    def find_possible_partitions(self, dfg2):
        def adj(node_set, net):
            adj_set = set()
            for node in node_set:
                adj_set = adj_set.union(set(net.neighbors(node)))
            return adj_set

        activity_list = set()
        for trace in self.log:
            for event in trace:
                activity_list.add(event['concept:name'])
        activity_list = list(activity_list)

        net = generate_nx_graph_from_dfg(dfg2)
        queue = []
        queue.append((set(), {'start'}))
        visited = []
        valid = []
        while len(queue) != 0:
            current = queue.pop()
            for x in powerset(current[1]):
                if len(x) != 0:
                    new_state = current[0].union(set(x))
                    if len(new_state & self.end_activities.keys()) != 0:
                        new_state.add('end')

                    if new_state not in visited:
                        new_adj = current[1].union(adj(set(x),net)) - new_state
                        queue.append((new_state, new_adj))
                        visited.append(new_state)
                        B = set(activity_list) - new_state
                        if (len(B) == 0) or (len(B) == len(activity_list)):
                            continue
                        if len(B & self.end_activities.keys()) != 0:
                            B.add('end')
                        if len(B & self.start_activities.keys()) != 0:
                            B.add('start')
                        BB = net.subgraph(B)
                        if 'end' in B:
                            disc_nodes_BB = set(BB.nodes) - set(nx.ancestors(BB, 'end')) - {'end'}
                            if len(disc_nodes_BB) == 0:
                                if ('end' in new_state) and ('start' in B) and (B not in visited):
                                    valid.append((new_state, B, {"seq", "exc", "par", "loop"}))
                                elif 'end' in new_state:
                                    valid.append((new_state, B, {"loop", "seq"}))
                                else:
                                    valid.append((new_state, B, {"seq"}))

                            else:
                                AA = net.subgraph(new_state.union(disc_nodes_BB))
                                if (set(nx.descendants(AA, 'start')) == set(AA.nodes) - {'start'}) and (
                                        new_state.union(disc_nodes_BB) not in visited):
                                    visited.append(new_state.union(disc_nodes_BB))
                                    if ('end' in new_state) and ('start' in B - disc_nodes_BB) and (
                                            B - disc_nodes_BB not in visited):
                                        valid.append((new_state.union(disc_nodes_BB), B - disc_nodes_BB,
                                                      {"seq", "exc", "par", "loop"}))
                                    elif 'end' in new_state:
                                        valid.append(
                                            (new_state.union(disc_nodes_BB), B - disc_nodes_BB, {"loop", "seq"}))
                                    else:
                                        valid.append((new_state.union(disc_nodes_BB), B - disc_nodes_BB, {"seq"}))
                                    queue.append((new_state.union(disc_nodes_BB),
                                                  new_adj.union(adj(disc_nodes_BB, net)) - new_state.union(disc_nodes_BB)))

                        if ('end' not in BB) and ('start' not in BB):
                            if nx.is_weakly_connected(BB):
                                valid.append((new_state, B, {"loop"}))
                                # queue.append((new_state, new_adj))
        return valid

    def project_a_b(self, dfg, A, B):
        new_dfg = {('start', 'a'): 0, ('start', 'b'): 0, ('a', 'a'): 0, ('b', 'b'): 0, ('a', 'b'): 0, ('b', 'a'): 0,
                       ('a', 'end'): 0, ('b', 'end'): 0, ('start', 'end'): 0}
        for edg in dfg:
            if edg[0] == 'start':
                if edg[1] in A - {'start', 'end'}:
                    new_dfg[('start', 'a')] += dfg[edg]
                elif edg[1] in B - {'start', 'end'}:
                    new_dfg[('start', 'b')] += dfg[edg]
                else:
                    new_dfg[('start', 'end')] += dfg[edg]
            elif edg[1] == 'end':
                if edg[0] in A - {'start', 'end'}:
                    new_dfg[('a', 'end')] += dfg[edg]
                elif edg[0] in B - {'start', 'end'}:
                    new_dfg[('b', 'end')] += dfg[edg]
                else:
                    new_dfg[('start', 'end')] += dfg[edg]
            # elif edg[0] in pp[0] and edg[1] in pp[0]:
            #     newlogP_dfg[('a', 'a')] += dfgP[edg]
            elif edg[0] in A - {'start', 'end'} and edg[1] in B - {'start', 'end'}:
                new_dfg[('a', 'b')] += dfg[edg]
            elif edg[0] in B - {'start', 'end'} and edg[1] in A - {'start', 'end'}:
                new_dfg[('b', 'a')] += dfg[edg]
            # elif edg[0] in pp[1] and edg[1] in pp[1]:
            #     newlogP_dfg[('b', 'b')] += dfgP[edg]
        return new_dfg

    def check_base_case(self,sup_thr, ratio, size_par):
        base_check = True
        # for tr in self.log:
        #     if len(set([ev['concept:name'] for ev in tr])) > 1:
        #         base_check = False
        #         cut = "not_base"
        #         break

        if len(self.activities.keys()) > 1:
            base_check = False
            cut = "not_base"


        if base_check==True:
            counter = 0
            counter_loop = 0
            visited_events = set()
            for tr in self.log:
                # if (len(set([ev['concept:name'] for ev in tr])) == 1) and len(tr) > 1:
                #     counter_loop += len(tr) - 1
                if len(tr) == 0:
                    counter += 1
                visited_events = visited_events.union(set([ev['concept:name'] for ev in tr]))
            if sum(self.activities.values())>0:
                p_prime_Lp = (len(self.log)-counter)/((len(self.log)-counter)+sum(self.activities.values()))
            else:
                p_prime_Lp = 'nd'

            counter_loopM = 0
            counterM = 0
            for tr in self.logM:
                # if (len(set([ev['concept:name'] for ev in tr])) == 1) and len(tr) > 1:
                #     counter_loopM += len(tr) - 1
                if len(tr) == 0:
                    counterM += 1
            if sum([len(x) for x in self.logM]) > 0:
                p_prime_Lm = (len(self.logM)-counterM) / ((len(self.logM)-counterM) + sum([len(x) for x in self.logM]))
            else:
                p_prime_Lm = 'nd'

            if (counter == len(self.log)) or (len(self.log) == 0):
                empty_log = True
                self.detected_cut = 'empty_log'
                cut = ('none', 'empty_log', 'none', 'none')
            else:
                empty_log = False
                if (counter > (sup_thr/2) * len(self.log)):
                    cost_single_exc = max(0, (sup_thr) * len(self.log) - counter) - ratio * size_par * max(0,(sup_thr) * len(self.logM) - counterM)
                    if cost_single_exc <= 0:
                        tau_xor, first_activity = True, visited_events.pop()
                    else:
                        tau_xor = False
                        new_log = log_instance.EventLog()
                        for n, tr in enumerate(self.log):
                            if len(tr) != 0:
                                new_log.append(tr)
                        self.log = new_log
                else:
                    if counter>0:
                        new_log = log_instance.EventLog()
                        for n, tr in enumerate(self.log):
                            if len(tr) != 0:
                                new_log.append(tr)
                        self.log = new_log
                    tau_xor = False
                if tau_xor==True:
                    cut = (({first_activity}, set()), 'exc', 'none', 'none')
                else:
                    # if (counter_loop > 0):
                    #     cost_single_loop = max(0, (sup_thr) * len(self.log) - counter_loop) - ratio * size_par * max(0,(sup_thr) * len(self.logM) - counter_loopM)
                    #     if cost_single_loop <= 0:
                    #         tau_loop, first_activity = True, visited_events.pop()
                    #     else:
                    #         tau_loop = False
                    # else:
                    #     tau_loop = False
                    if abs(p_prime_Lp - 0.5) > sup_thr/2:
                        if p_prime_Lm != 'nd':
                            cost_single_loop = max(0,sup_thr - abs(p_prime_Lp - 0.5)) - ratio * size_par * max(0, sup_thr - abs(p_prime_Lm - 0.5))
                        else:
                            cost_single_loop = max(0,sup_thr - ratio * size_par * abs(p_prime_Lp - 0.5))
                        if cost_single_loop <= 0:
                            tau_loop, first_activity = True, visited_events.pop()
                        else:
                            tau_loop = False
                    else:
                        tau_loop = False
                    if tau_loop==True:
                        cut = (({first_activity}, set()), 'loop1', 'none', 'none')
                    else:
                        single_activity = True
                        self.detected_cut = 'single_activity'
                        cut = ('none', 'single_activity', 'none', 'none')
        return base_check, cut

    def max_flow_graph(self, dfg):
        net = nx.DiGraph()
        net.add_weighted_edges_from([(x[0], x[1], dfg[x]) for x in dfg], weight='weight')
        flow_graph = {}
        for x in net.nodes:
            for y in net.nodes:
                if (x != y):
                    flow_graph[(x, y)] = nx.algorithms.flow.maximum_flow(net, x, y, capacity='weight')[0]
        return flow_graph

    def detect_cut(self, second_iteration=False, parameters=None, sup= None, ratio = None, size_par = None):
        ratio = ratio
        sup_thr = sup
        if pkgutil.find_loader("networkx"):
            start3 = time.time()
            import networkx as nx

            if parameters is None:
                parameters = {}
            activity_key = exec_utils.get_param_value(Parameters.ACTIVITY_KEY, parameters,
                                                      pmutil.xes_constants.DEFAULT_NAME_KEY)

            # check base cases:
            isbase, cut = SubtreePlain.check_base_case(self, sup, ratio, size_par)

            if isbase==False:
                dfg2 = dfg_discovery.apply(self.log_art, variant=dfg_discovery.Variants.FREQUENCY)
                del dfg2[('start', 'end')]
                possible_partitions = SubtreePlain.find_possible_partitions(self, dfg2)

                # seq, exc, par
                start_pre = time.time()
                cut = []

                dfgP = dfg_discovery.apply(self.log_art, variant=dfg_discovery.Variants.FREQUENCY)
                # gviz = dfg_visualization.apply(dfgP)
                # dfg_visualization.view(gviz)
                number_of_edges_P = sum([dfgP[x] for x in dfgP])
                dfgM = dfg_discovery.apply(self.logM_art, variant=dfg_discovery.Variants.FREQUENCY)
                activitiesM = attributes_get.get_attribute_values(self.logM, activity_key)
                # gviz = dfg_visualization.apply(dfgM)
                # dfg_visualization.view(gviz)

                end_pre1 = time.time()
                # print("pre_part1: " + str(end_pre1 - start_pre))
                #########################
                fP = SubtreePlain.max_flow_graph(self, dfgP)
                fM = SubtreePlain.max_flow_graph(self, dfgM)

                start_acts_P = set(self.start_activities.keys())
                start_acts_M = set(self.start_activitiesM.keys())
                end_acts_P = set(self.end_activities.keys())
                end_acts_M = set(self.end_activitiesM.keys())

                if True:
                    missing_loopP = 0
                    missing_loopM = 0
                    rej_tau_loop = False
                    c_rec = 0

                    if len(start_acts_P.intersection(end_acts_P)) !=0:
                        rej_tau_loop = True
                    total_loop_P = sum([dfgP[x] for x in dfgP if (x[0] in end_acts_P) and (x[1] in start_acts_P)])
                    for x in start_acts_P:
                        for y in end_acts_P:
                            # L1P = max(0, total_loop_P * sup_thr * (self.start_activities[x] / (sum(self.start_activities.values()))) * (self.end_activities[y] / (sum(self.end_activities.values()))) - dfgP[(y, x)])
                            L1P = max(0, len(self.log) * sup_thr * (self.start_activities[x] / (sum(self.start_activities.values()))) * (self.end_activities[y] / (sum(self.end_activities.values()))) - dfgP[(y, x)])
                            missing_loopP += L1P
                            c_rec += dfgP[(y, x)]

                    total_loop_M = sum([dfgM[x] for x in dfgM if (x[0] in end_acts_P) and (x[1] in start_acts_P)])
                    for x in start_acts_P.intersection(self.start_activitiesM.keys()):
                        for y in end_acts_P.intersection(self.end_activitiesM.keys()):
                            # L1M = max(0, total_loop_M * sup_thr * (self.start_activitiesM[x] / (sum(self.start_activitiesM.values()))) * (self.end_activitiesM[y] / (sum(self.end_activitiesM.values()))) - dfgM[(y, x)])
                            L1M = max(0, len(self.logM) * sup_thr * (self.start_activitiesM[x] / (sum(self.start_activitiesM.values()))) * (self.end_activitiesM[y] / (sum(self.end_activitiesM.values()))) - dfgM[(y, x)])
                            missing_loopM += L1M

                    cost_loop_P = missing_loopP
                    cost_loop_M = missing_loopM
                    # if c_rec >= 0.3*len(self.log):
                    # if c_rec >0:
                    #     print(c_rec==total_loop_P)
                    if rej_tau_loop == False and c_rec >0:
                        cut.append(((start_acts_P, end_acts_P), 'loop_tau', cost_loop_P,  cost_loop_P - ratio * size_par * cost_loop_M))

                end_pre = time.time()
                # print("pre: " + str(end_pre - start_pre))
                start_cut = time.time()
                ratio_backup = ratio
                for pp in possible_partitions:
                    A = pp[0]
                    B = pp[1]

                    type = pp[2]
                    if len(set(activitiesM.keys()).intersection(A))==0 or len(set(activitiesM.keys()).intersection(B))==0:
                        ratio = 0
                    else:
                        ratio = ratio_backup

                    newlogP_dfg = SubtreePlain.project_a_b(self, dfgP, A, B)
                    newlogM_dfg = SubtreePlain.project_a_b(self, dfgM, A, B)


                    #####################################################################
                    # seq check
                    if "seq" in type:
                        ba_count_P = newlogP_dfg[('b', 'a')]
                        ba_count_M = newlogM_dfg[('b', 'a')]
                        OE_A_P = set([x[0] for x in dfgP if (x[0] in A and (x[1] in B or x[1]=='end'))])
                        IS_B_P = set([x[1] for x in dfgP if (x[1] in B and (x[0] in A or x[0]=='start'))])
                        OE_A_M = set([x[0] for x in dfgM if (x[0] in A and (x[1] in B or x[1]=='end'))])
                        IS_B_M = set([x[1] for x in dfgM if (x[1] in B and (x[0] in A or x[0]=='start'))])
                        missing_seq_P = 0

                        if A=={'end', 'A_ACCEPTED', 'start'}:
                            print('wait')

                        for x in A-{'start','end'}:
                            for y in B-{'start','end'}:
                                if x!=y:
                                    missing_seq_P += max(0, self.activities[x] * sup_thr * (self.activities[y]/(sum([self.activities[p] for p in B-{'start','end'}])+sum([self.activities[p] for p in A-{'start','end'}]))) - fP[(x,y)])

                        for x in OE_A_P-{'start','end'}:
                            for y in IS_B_P-{'start','end'}:
                                missing_seq_P += max(0,sum([dfgP[z] for z in dfgP if (z[0]==x and (z[1] in B or z[1]=='end'))]) * sup * (sum([dfgP[z] for z in dfgP if (z[1]==y and (z[0] in A or z[0]=='start'))])/(sum([dfgP[z] for z in dfgP if (z[1] in B and (z[0] in A or z[0]=='start'))])+sum([dfgP[z] for z in dfgP if (z[0] in A and (z[1] in B or z[1]=='end'))]))) - dfgP[(x,y)])

                        missing_seq_M = 0
                        for x in A.intersection(activitiesM.keys())-{'start','end'}:
                            for y in B.intersection(activitiesM.keys())-{'start','end'}:
                                if ((x,y) in fM) and (x != y):
                                    missing_seq_M += max(0, activitiesM[x] * sup_thr * (activitiesM[y]/(sum([activitiesM[p] for p in B.intersection(activitiesM.keys())-{'start','end'}])+sum([activitiesM[p] for p in A.intersection(activitiesM.keys())-{'start','end'}]))) - fM[(x, y)])
                                    # print(str(x)+".."+str(y)+"::"+str(missing_seq_M))

                        for x in OE_A_M.intersection(activitiesM.keys())-{'start','end'}:
                            for y in IS_B_M.intersection(activitiesM.keys())-{'start','end'}:
                                missing_seq_M += max(0,sum([dfgM[z] for z in dfgM if (z[0]==x and (z[1] in B or z[1]=='end'))]) * sup * (sum([dfgM[z] for z in dfgM if (z[1]==y and (z[0] in A or z[0]=='start'))])/(sum([dfgM[z] for z in dfgM if (z[1] in B and (z[0] in A or z[0]=='start'))])++sum([dfgM[z] for z in dfgM if (z[0] in A and (z[1] in B or z[1]=='end'))]))) - dfgM[(x,y)])

                        cost_seq_P = ba_count_P + missing_seq_P
                        cost_seq_M = ba_count_M + missing_seq_M

                        cut.append(((A, B), 'seq', cost_seq_P, cost_seq_P - ratio* size_par * cost_seq_M))
                    #####################################################################

                    #####################################################################
                    # xor check
                    # if A_reachable_from_start and B_reachable_from_start and A_reachable_to_end and B_reachable_to_end:
                    if "exc" in type:
                        ab_count_P = newlogP_dfg[('a', 'b')]
                        ab_count_M = newlogM_dfg[('a', 'b')]
                        ba_count_P = newlogP_dfg[('b', 'a')]
                        ba_count_M = newlogM_dfg[('b', 'a')]

                        cost_exc_P = ba_count_P + ab_count_P
                        cost_exc_M = ba_count_M + ab_count_M
                        cut.append(((A, B), 'exc', cost_exc_P, cost_exc_P - ratio* size_par * cost_exc_M))
                    #####################################################################

                    #####################################################################
                    # xor-tau check
                    # if A_reachable_from_start and B_reachable_from_start and A_reachable_to_end and B_reachable_to_end and (
                    #         newlogP_dfg[('start', 'end')] > (sup_thr/2 * len(self.log))):
                    # if "exc" in type and (newlogP_dfg[('start', 'end')] > (sup_thr * len(self.log))):
                    if newlogP_dfg[('start', 'end')]>0:
                        missing_exc_tau_P = 0
                        missing_exc_tau_P += max(0, sup_thr * len(self.log) - newlogP_dfg[('start', 'end')])

                        missing_exc_tau_M = 0
                        missing_exc_tau_M += max(0, sup_thr * len(self.logM) - newlogM_dfg[('start', 'end')])

                        cost_exc_tau_P = missing_exc_tau_P
                        cost_exc_tau_M = missing_exc_tau_M
                        cut.append(((A.union(B), set()), 'exc2', cost_exc_tau_P,cost_exc_tau_P - ratio * size_par * cost_exc_tau_M))
                    #####################################################################

                    #####################################################################
                    # parallel check
                    # if A_reachable_from_start and B_reachable_from_start and A_reachable_to_end and B_reachable_to_end:
                    if "par" in type:
                        missing_par_P = 0
                        for a in A-{'start','end'}:
                            for b in B-{'start','end'}:
                                missing_par_P += max(0, (self.activities[a] * sup_thr * self.activities[b]) / ((sum([self.activities[p] for p in B - {'start', 'end'}])) + (sum([self.activities[p] for p in A - {'start', 'end'}])))- dfgP[(a,b)])
                                missing_par_P += max(0, (self.activities[b] * sup_thr * self.activities[a]) / ((sum([self.activities[p] for p in B - {'start', 'end'}])) + (sum([self.activities[p] for p in A - {'start', 'end'}])))- dfgP[(b, a)])

                        missing_par_M = 0
                        for a in A.intersection(activitiesM.keys()) - {'start', 'end'}:
                            for b in B.intersection(activitiesM.keys()) - {'start', 'end'}:
                                if len([dfgM[x] for x in dfgM if x[0] == a])!=0:
                                    # missing_par_M += max(0, activitiesM[a] * sup_thr * (activitiesM[b] / (activitiesM[a]-sum([dfgM[e] for e in dfgM if (e[0]==a and (e[1] in B))]) +sum([activitiesM[p] for p in B.intersection(activitiesM.keys()) - {'start', 'end'}]))) - dfgM[(a,b)])
                                    missing_par_M += max(0, (activitiesM[a] * sup_thr * activitiesM[b]) / (sum([activitiesM[p] for p in B.intersection(activitiesM.keys()) - {'start', 'end'}])+sum([activitiesM[p] for p in A.intersection(activitiesM.keys()) - {'start', 'end'}])) - dfgM[(a, b)])
                                if len([dfgM[x] for x in dfgM if x[0] == b])!=0:
                                    # missing_par_M += max(0, activitiesM[b] * sup_thr * (activitiesM[a] / (activitiesM[b]-sum([dfgM[e] for e in dfgM if (e[0]==b and (e[1] in A))]) +sum([activitiesM[p] for p in A.intersection(activitiesM.keys()) - {'start', 'end'}]))) - dfgM[(b,a)])
                                    missing_par_M += max(0, (activitiesM[a] * sup_thr * activitiesM[b]) / (sum([activitiesM[p] for p in B.intersection(activitiesM.keys()) - {'start', 'end'}])+sum([activitiesM[p] for p in A.intersection(activitiesM.keys()) - {'start', 'end'}])) - dfgM[(b, a)])

                        cost_par_P = missing_par_P
                        cost_par_M = missing_par_M
                        cut.append(((A, B), 'par', cost_par_P, cost_par_P - ratio * size_par * cost_par_M))
                    #####################################################################

                    #####################################################################
                    # loop check
                    # if A_reachable_from_start and A_reachable_to_end and (B_reachable_from_start==False) and (B_reachable_to_end==False):
                    # if A_reachable_from_start and A_reachable_to_end and B_weakly_connected:
                    if "loop" in type:
                        flag_loop_valid1 = False
                        flag_loop_valid2 = False
                        # flag_loop_valid1 = True
                        # flag_loop_valid2 = True


                        for a in start_acts_P:
                            for b in B - {'start', 'end'}:
                                if dfgP[(b, a)]>0:
                                    flag_loop_valid1 = True
                        for a in end_acts_P:
                            for b in B - {'start', 'end'}:
                                if dfgP[(a, b)]>0:
                                    flag_loop_valid2 = True

                        bad_edges_loop_P = 0
                        bad_edges_loop_P += newlogP_dfg[('start', 'b')]
                        bad_edges_loop_P += newlogP_dfg[('b', 'end')]
                        for a in (A-({'start', 'end'}.union(start_acts_P))):
                            for b in B-{'start', 'end'}:
                                bad_edges_loop_P += dfgP[(b, a)]
                        for a in (A-({'start', 'end'}.union(end_acts_P))):
                            for b in B-{'start', 'end'}:
                                bad_edges_loop_P += dfgP[(a, b)]



                        bad_edges_loop_M = 0
                        bad_edges_loop_M += newlogM_dfg[('start', 'b')]
                        bad_edges_loop_M += newlogM_dfg[('b', 'end')]
                        for a in (A - ({'start', 'end'}.union(start_acts_M))):
                            for b in B - {'start', 'end'}:
                                bad_edges_loop_M += dfgM[(b, a)]
                        for a in (A - ({'start', 'end'}.union(end_acts_M))):
                            for b in B - {'start', 'end'}:
                                bad_edges_loop_M += dfgM[(a, b)]


                        b_out_P = set()
                        b_in_P = set()
                        BotoAs_P = 0
                        dict_BotoAs_P = {}
                        AetoBi_P = 0
                        dict_AetoBi_P = {}
                        for a in A.intersection(start_acts_P):
                            for b in B-{'start','end'}:
                                if dfgP[(b,a)]>0:
                                    BotoAs_P += dfgP[(b,a)]
                                    b_out_P.add(b)
                                    if b not in dict_BotoAs_P:
                                        dict_BotoAs_P[b] = dfgP[(b,a)]
                                    else:
                                        dict_BotoAs_P[b] += dfgP[(b, a)]
                        for a in A.intersection(end_acts_P):
                            for b in B-{'start','end'}:
                                if dfgP[(a,b)]>0:
                                    AetoBi_P += dfgP[(a,b)]
                                    b_in_P.add(b)
                                    # if a not in dict_AetoBi_P:
                                    #     dict_AetoBi_P[a] = dfgP[(a,b)]
                                    # else:
                                    #     dict_AetoBi_P[a] += dfgP[(a,b)]
                                    if b not in dict_AetoBi_P:
                                        dict_AetoBi_P[b] = dfgP[(a,b)]
                                    else:
                                        dict_AetoBi_P[b] += dfgP[(a,b)]

                        b_out_M = set()
                        b_in_M = set()
                        BotoAs_M = 0
                        AetoBi_M = 0
                        dict_BotoAs_M = {}
                        dict_AetoBi_M = {}
                        for a in A.intersection(start_acts_M):
                            for b in B-{'start','end'}:
                                if dfgM[(b, a)] > 0:
                                    BotoAs_M += dfgM[(b, a)]
                                    b_out_M.add(b)
                                    # if a not in dict_BotoAs_M:
                                    #     dict_BotoAs_M[a] = dfgM[(b,a)]
                                    # else:
                                    #     dict_BotoAs_M[a] += dfgM[(b, a)]
                                    if b not in dict_BotoAs_M:
                                        dict_BotoAs_M[b] = dfgM[(b,a)]
                                    else:
                                        dict_BotoAs_M[b] += dfgM[(b, a)]

                        for a in A.intersection(end_acts_M):
                            for b in B-{'start','end'}:
                                if dfgM[(a, b)] > 0:
                                    AetoBi_M += dfgM[(a, b)]
                                    b_in_M.add(b)
                                    # if a not in dict_AetoBi_M:
                                    #     dict_AetoBi_M[a] = dfgM[(a,b)]
                                    # else:
                                    #     dict_AetoBi_M[a] += dfgM[(a,b)]
                                    if b not in dict_AetoBi_M:
                                        dict_AetoBi_M[b] = dfgM[(a,b)]
                                    else:
                                        dict_AetoBi_M[b] += dfgM[(a,b)]

                        missing_loop_P = 0
                        M_P = max(BotoAs_P, AetoBi_P)
                        if M_P == 0:
                            t_P = 0
                        else:
                            t_P = abs(BotoAs_P - AetoBi_P)/M_P
                        # print(str(B)+': '+ str(BotoAs_P) + "--" + str(AetoBi_P))
                        if len(b_out_P)!=0:
                            for a in A.intersection(start_acts_P):
                                for b in b_out_P:
                                    # missing_loop_P += max(0, sum([dfgP[x] for x in dfgP if (x[0]==b and x[1] in A)]) * sup_thr * (self.start_activities[a] / sum([self.start_activities[p] for p in start_acts_P])) - dfgP[(b,a)])
                                    missing_loop_P += max(0, M_P * sup * (self.start_activities[a]/sum(self.start_activities.values())) * (dict_BotoAs_P[b]/sum(dict_BotoAs_P.values()))- dfgP[(b,a)])

                        if len(b_in_P) !=0:
                            for a in A.intersection(end_acts_P):
                                for b in b_in_P:
                                    missing_loop_P += max(0, M_P * sup * (self.end_activities[a]/sum(self.end_activities.values())) * (dict_AetoBi_P[b]/sum(dict_AetoBi_P.values()))- dfgP[(a,b)])

                        missing_loop_M = 0
                        M_M = max(BotoAs_M, AetoBi_M)
                        if len(b_out_M) != 0:
                            if len((A.intersection(start_acts_M)).intersection(activitiesM.keys())) != 0:
                                for a in (A.intersection(start_acts_M)).intersection(self.start_activitiesM.keys()):
                                    for b in b_out_M.intersection(activitiesM.keys()):
                                        if len([dfgM[x] for x in dfgM if x[0] == b]) != 0:
                                            if a in start_acts_M:
                                                missing_loop_M += max(0, M_M * sup * (self.start_activitiesM[a]/sum(self.start_activitiesM.values())) * (dict_BotoAs_M[b]/sum(dict_BotoAs_M.values()))- dfgM[(b,a)])

                        if len(b_in_M) != 0:
                            for a in (A.intersection(end_acts_M)).intersection(self.end_activitiesM.keys()):
                                for b in b_in_M.intersection(activitiesM.keys()):
                                    if len([dfgM[x] for x in dfgM if x[0] == a]) != 0:
                                        if a in end_acts_M:
                                            missing_loop_M += max(0, M_M * sup * (self.end_activitiesM[a]/sum(self.end_activitiesM.values())) * (dict_AetoBi_M[b]/sum(dict_AetoBi_M.values()))- dfgM[(a,b)])

                        cost_loop_P = bad_edges_loop_P + missing_loop_P
                        cost_loop_M = bad_edges_loop_M + missing_loop_M
                        if (flag_loop_valid1 is True) and (flag_loop_valid2 is True) and t_P<0.3:
                            cut.append(((A, B), 'loop', cost_loop_P, cost_loop_P - ratio * size_par *  cost_loop_M))

                        ratio = ratio_backup

                # print(len(cut))
                cut = [x for x in cut if x[2]<0.3*number_of_edges_P]
                # cut = [x for x in cut if x[2] < sup_thr * len(self.log)]

                sorted_cuts = sorted(cut, key=lambda x: (x[3], x[2],['exc','exc2','seq','par','loop','loop_tau'].index(x[1]), (len(x[0][0]) * len(x[0][1]) / (len(x[0][0]) + len(x[0][1])))))
                # sorted_cuts = sorted(cut, key=lambda x: (x[3], x[2],['exc','exc2','seq','par','loop','loop_tau'].index(x[1]), -(len(x[0][0]) * len(x[0][1]) / (len(x[0][0]) + len(x[0][1])))))
                if len(sorted_cuts) != 0:
                    cut = sorted_cuts[0]
                else:
                    cut= ('none', 'none','none','none')

                # the following part searches for a cut in the current log
                # if a cut is found, the log is split according to the cut, the resulting logs are saved in new_logs
                # recursion is used on all the logs in new_logs

                end_cut = time.time()
                # print("checking cuts: " + str(end_cut - start_cut))

            # print(cut)

            end3 = time.time()
            # print("all: " + str(end3-start3))


            if cut[1] == 'par':
                self.detected_cut = 'parallel'
                LAP,LBP = split.split('par', [cut[0][0], cut[0][1]], self.log, activity_key)
                LAM, LBM = split.split('par', [cut[0][0], cut[0][1]], self.logM, activity_key)
                new_logs = [[LAP,LAM],[LBP,LBM]]
                for l in new_logs:
                    new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                    activities = attributes_get.get_attribute_values(l[0], activity_key)
                    start_activities = list(
                        start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                    end_activities = list(
                        end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                    self.children.append(
                        SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
                                     self.rec_depth + 1,
                                     noise_threshold=self.noise_threshold, start_activities=start_activities,
                                     end_activities=end_activities,
                                     initial_start_activities=self.initial_start_activities,
                                     initial_end_activities=self.initial_end_activities,
                                     parameters=parameters, sup= sup, ratio = ratio, size_par = size_par))
            elif cut[1] == 'seq':
                self.detected_cut = 'sequential'
                LAP,LBP = split.split('seq', [cut[0][0], cut[0][1]], self.log, activity_key)
                LAM, LBM = split.split('seq', [cut[0][0], cut[0][1]], self.logM, activity_key)
                new_logs = [[LAP,LAM],[LBP,LBM]]
                for l in new_logs:
                    new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                    activities = attributes_get.get_attribute_values(l[0], activity_key)
                    start_activities = list(
                        start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                    end_activities = list(
                        end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                    self.children.append(
                        SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
                                     self.rec_depth + 1,
                                     noise_threshold=self.noise_threshold, start_activities=start_activities,
                                     end_activities=end_activities,
                                     initial_start_activities=self.initial_start_activities,
                                     initial_end_activities=self.initial_end_activities,
                                     parameters=parameters, sup= sup, ratio = ratio, size_par = size_par))
            elif (cut[1] == 'exc') or (cut[1] == 'exc2'):
                self.detected_cut = 'concurrent'
                LAP,LBP = split.split('exc', [cut[0][0], cut[0][1]], self.log, activity_key)
                LAM, LBM = split.split('exc', [cut[0][0], cut[0][1]], self.logM, activity_key)
                new_logs = [[LAP,LAM],[LBP,LBM]]
                for l in new_logs:
                    new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                    activities = attributes_get.get_attribute_values(l[0], activity_key)
                    start_activities = list(
                        start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                    end_activities = list(
                        end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                    self.children.append(
                        SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
                                     self.rec_depth + 1,
                                     noise_threshold=self.noise_threshold, start_activities=start_activities,
                                     end_activities=end_activities,
                                     initial_start_activities=self.initial_start_activities,
                                     initial_end_activities=self.initial_end_activities,
                                     parameters=parameters, sup= sup, ratio = ratio, size_par = size_par))

            elif cut[1] == 'loop':
                self.detected_cut = 'loopCut'
                LAP,LBP = split.split('loop', [cut[0][0], cut[0][1]], self.log, activity_key)
                LAM, LBM = split.split('loop', [cut[0][0], cut[0][1]], self.logM, activity_key)
                new_logs = [[LAP,LAM],[LBP,LBM]]
                for l in new_logs:
                    new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                    activities = attributes_get.get_attribute_values(l[0], activity_key)
                    start_activities = list(
                        start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                    end_activities = list(
                        end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                    self.children.append(
                        SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
                                     self.rec_depth + 1,
                                     noise_threshold=self.noise_threshold, start_activities=start_activities,
                                     end_activities=end_activities,
                                     initial_start_activities=self.initial_start_activities,
                                     initial_end_activities=self.initial_end_activities,
                                     parameters=parameters, sup= sup, ratio = ratio, size_par = size_par))

            elif cut[1] == 'loop1':
                self.detected_cut = 'loopCut'
                LAP,LBP = split.split('loop1', [cut[0][0], cut[0][1]], self.log, activity_key)
                LAM, LBM = split.split('loop1', [cut[0][0], cut[0][1]], self.logM, activity_key)
                new_logs = [[LAP,LAM],[LBP,LBM]]
                for l in new_logs:
                    new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                    activities = attributes_get.get_attribute_values(l[0], activity_key)
                    start_activities = list(
                        start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                    end_activities = list(
                        end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                    self.children.append(
                        SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
                                     self.rec_depth + 1,
                                     noise_threshold=self.noise_threshold, start_activities=start_activities,
                                     end_activities=end_activities,
                                     initial_start_activities=self.initial_start_activities,
                                     initial_end_activities=self.initial_end_activities,
                                     parameters=parameters, sup= sup, ratio = ratio, size_par = size_par))

            elif cut[1] == 'loop_tau':
                self.detected_cut = 'loopCut'
                LAP,LBP = split.split('loop_tau', [cut[0][0], cut[0][1]], self.log, activity_key)
                LAM, LBM = split.split('loop_tau', [cut[0][0], cut[0][1]], self.logM, activity_key)
                new_logs = [[LAP,LAM],[LBP,LBM]]
                for l in new_logs:
                    new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                    activities = attributes_get.get_attribute_values(l[0], activity_key)
                    start_activities = list(
                        start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                    end_activities = list(
                        end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                    self.children.append(
                        SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
                                     self.rec_depth + 1,
                                     noise_threshold=self.noise_threshold, start_activities=start_activities,
                                     end_activities=end_activities,
                                     initial_start_activities=self.initial_start_activities,
                                     initial_end_activities=self.initial_end_activities,
                                     parameters=parameters, sup= sup, ratio = ratio, size_par = size_par))

            elif cut[1] == 'none':
                self.detected_cut = 'flower'

        else:
            msg = "networkx is not available. inductive miner cannot be used!"
            logging.error(msg)
            raise Exception(msg)



def make_tree(logp, logm, dfg, master_dfg, initial_dfg, activities, c, recursion_depth, noise_threshold, start_activities,
              end_activities, initial_start_activities, initial_end_activities, parameters=None, sup= None, ratio = None, size_par = None):
    tree = SubtreePlain(logp,logm, dfg, master_dfg, initial_dfg, activities, c, recursion_depth, noise_threshold,
                        start_activities,
                        end_activities, initial_start_activities, initial_end_activities, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)

    return tree
