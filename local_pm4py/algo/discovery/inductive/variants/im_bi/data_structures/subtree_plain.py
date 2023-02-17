import collections
from copy import copy
import time
from pm4py.algo.discovery.dfg.utils.dfg_utils import infer_start_activities, infer_end_activities

from pm4py.algo.discovery.dfg.variants import native as dfg_inst
from pm4py import util as pmutil
from local_pm4py.algo.discovery.inductive.variants.im_bi.util import splitting as split
from pm4py.statistics.attributes.log import get as attributes_get
from pm4py.statistics.end_activities.log import get as end_activities_get
from pm4py.statistics.start_activities.log import get as start_activities_get
from pm4py.util import exec_utils
from pm4py.util import constants
from enum import Enum
from pm4py.objects.log import obj as log_instance
from pm4py.util import xes_constants
from local_pm4py.algo.discovery.dfg import algorithm as dfg_discovery
import networkx as nx
from pm4py.algo.filtering.log.start_activities import start_activities_filter
from pm4py.algo.filtering.log.end_activities import end_activities_filter
from pm4py.algo.discovery.dfg.utils.dfg_utils import get_activities_from_dfg
from local_pm4py.algo.analysis import dfg_functions
import copy
from collections import Counter

def artificial_start_end(log):
    st = 'start'
    en = 'end'
    activity_key = xes_constants.DEFAULT_NAME_KEY
    start_event = log_instance.Event()
    start_event[activity_key] = st

    end_event = log_instance.Event()
    end_event[activity_key] = en

    for trace in log:
        trace.insert(0, start_event)
        trace.append(end_event)
    return log

def generate_nx_graph_from_dfg(dfg):
    dfg_acts = set()
    for x in dfg:
        dfg_acts.add(x[0])
        dfg_acts.add(x[1])
    G = nx.DiGraph()
    for act in dfg_acts:
        G.add_node(act)
    for edge in dfg:
        G.add_edge(edge[0], edge[1], weight=dfg[edge])
    return G


class SubtreePlain(object):
    def __init__(self, logp,logm, dfg, master_dfg, initial_dfg, activities, counts, rec_depth, noise_threshold=0,
                 start_activities=None, end_activities=None, initial_start_activities=None,
                 initial_end_activities=None, parameters=None, real_init=True, sup= None, ratio = None, size_par = None):

        if real_init:
            self.master_dfg = copy.copy(master_dfg)
            self.initial_dfg = copy.copy(initial_dfg)
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

            self.log = logp
            self.log_art = artificial_start_end(copy.deepcopy(logp))
            self.logM = logm
            self.logM_art = artificial_start_end(copy.deepcopy(logm))
            self.inverted_dfg = None
            self.original_log = logp
            self.activities = None

            self.initialize_tree(dfg, logp,logm, initial_dfg, activities, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)


    def initialize_tree(self, dfg, logp,logm, initial_dfg, activities, second_iteration=False, end_call=True,
                        parameters=None, sup= None, ratio = None, size_par = None):


        if activities is None:
            self.activities = get_activities_from_dfg(dfg)
        else:
            self.activities = copy.copy(activities)
        self.detected_cut = None
        self.children = []
        self.log = logp
        self.log_art = artificial_start_end(logp.__deepcopy__())
        self.logM = logm
        self.logM_art = artificial_start_end(logm.__deepcopy__())
        self.original_log = logp
        self.parameters = parameters

        self.detect_cut(second_iteration=False, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)


    def detect_cut(self, second_iteration=False, parameters=None, sup= None, ratio = None, size_par = None):
        ratio = ratio
        sup_thr = sup

        logP_var = Counter([tuple([x['concept:name'] for x in t]) for t in self.log])
        logM_var = Counter([tuple([x['concept:name'] for x in t]) for t in self.logM])


        if parameters is None:
            parameters = {}
        activity_key = exec_utils.get_param_value(constants.PARAMETER_CONSTANT_ACTIVITY_KEY, parameters,
                                                  pmutil.xes_constants.DEFAULT_NAME_KEY)

        # check base cases:
        isbase, cut = dfg_functions.check_base_case(self, logP_var,logM_var, sup, ratio, size_par)

        if isbase==False:
            dfg2 = dfg_discovery.apply(self.log_art, variant=dfg_discovery.Variants.FREQUENCY)
            netP = generate_nx_graph_from_dfg(dfg2)
            del dfg2[('start', 'end')]

            dfg2M = dfg_discovery.apply(self.logM_art, variant=dfg_discovery.Variants.FREQUENCY)
            netM = generate_nx_graph_from_dfg(dfg2M)
            del dfg2M[('start', 'end')]

            if parameters == {}:
                feat_scores_togg = collections.defaultdict(lambda: 1, {})
                feat_scores = collections.defaultdict(lambda: 1, {})
                for x in dfg2.keys():
                    feat_scores[x] = 1
                    feat_scores_togg[x] = 1
                for y in dfg2M.keys():
                    feat_scores[y] = 1
                    feat_scores_togg[y] = 1


            possible_partitions = dfg_functions.find_possible_partitions(netP)

            cut = []

            dfgP = dfg_discovery.apply(self.log_art, variant=dfg_discovery.Variants.FREQUENCY)
            dfgM = dfg_discovery.apply(self.logM_art, variant=dfg_discovery.Variants.FREQUENCY)
            activitiesM = set(a for x in logM_var.keys() for a in x)


            #########################
            fP = dfg_functions.max_flow_graph(netP)
            fM = dfg_functions.max_flow_graph(netM)


            start_acts_P = set([x[1] for x in dfgP if (x[0] == 'start')])-{'end'}
            end_acts_P = set([x[0] for x in dfgP if (x[1] == 'end')])-{'start'}

            if True:
                missing_loopP = 0
                missing_loopM = 0
                rej_tau_loop = False
                c_rec = 0

                if len(start_acts_P.intersection(end_acts_P)) !=0:
                    rej_tau_loop = True
                for x in start_acts_P:
                    for y in end_acts_P:
                        L1P = max(0, len(self.log) * sup_thr * (self.start_activities[x] / (sum(self.start_activities.values()))) * (self.end_activities[y] / (sum(self.end_activities.values()))) - dfgP[(y, x)])
                        missing_loopP += L1P
                        c_rec += dfgP[(y, x)]

                for x in start_acts_P.intersection(self.start_activitiesM.keys()):
                    for y in end_acts_P.intersection(self.end_activitiesM.keys()):
                        L1M = max(0, len(self.logM) * sup_thr * (self.start_activitiesM[x] / (sum(self.start_activitiesM.values()))) * (self.end_activitiesM[y] / (sum(self.end_activitiesM.values()))) - dfgM[(y, x)])
                        missing_loopM += L1M

                cost_loop_P = missing_loopP
                cost_loop_M = missing_loopM

                if rej_tau_loop == False and c_rec >0:
                    cut.append(((start_acts_P, end_acts_P), 'loop_tau', cost_loop_P, cost_loop_M,  cost_loop_P - ratio * size_par * cost_loop_M,1))
            ratio_backup = ratio

            for pp in possible_partitions:
                A = pp[0] - {'start', 'end'}
                B = pp[1] - {'start', 'end'}

                start_A_P = set([x[1] for x in dfgP if ((x[0] == 'start') and (x[1] in A))])
                end_A_P = set([x[0] for x in dfgP if (x[0] in A and (x[1] == 'end'))])
                start_B_P = set([x[1] for x in dfgP if (x[1] in B and (x[0] == 'start'))])
                input_B_P = set([x[1] for x in dfgP if ((x[0] not in B) and (x[1] in B))])
                output_B_P = set([x[0] for x in dfgP if ((x[0] in B) and (x[1] not in B))])

                start_A_M = set([x[1] for x in dfgM if ((x[0] == 'start') and (x[1] in A))])
                end_A_M = set([x[0] for x in dfgM if (x[0] in A and (x[1] == 'end'))])
                start_B_M = set([x[1] for x in dfgM if (x[1] in B and (x[0] == 'start'))])
                input_B_M = set([x[1] for x in dfgM if ((x[0] not in B) and (x[1] in B))])
                output_B_M = set([x[0] for x in dfgM if ((x[0] in B) and (x[1] not in B))])

                type = pp[2]
                if len(set(activitiesM).intersection(A))==0 or len(set(activitiesM).intersection(B))==0:
                    ratio = 0
                else:
                    ratio = ratio_backup

                #####################################################################
                # seq check

                    fit_seq = dfg_functions.fit_seq(logP_var, A, B)
                    if fit_seq > 0.0:
                        cost_seq_P = dfg_functions.cost_seq(netP, A, B, start_B_P, end_A_P, sup, fP, feat_scores)
                        cost_seq_M = dfg_functions.cost_seq(netM, A.intersection(activitiesM), B.intersection(activitiesM), start_B_M.intersection(activitiesM), end_A_M.intersection(activitiesM), sup, fM, feat_scores_togg)
                        cut.append(((A, B), 'seq', cost_seq_P, cost_seq_M, cost_seq_P - ratio* size_par * cost_seq_M, fit_seq))
                #####################################################################

                #####################################################################
                # xor check
                if "exc" in type:
                    fit_exc = dfg_functions.fit_exc(logP_var, A, B)
                    if fit_exc > 0.0:
                        cost_exc_P = dfg_functions.cost_exc(netP, A, B, feat_scores)
                        cost_exc_M = dfg_functions.cost_exc(netM, A.intersection(activitiesM), B.intersection(activitiesM), feat_scores)
                        cut.append(((A, B), 'exc', cost_exc_P, cost_exc_M, cost_exc_P - ratio* size_par * cost_exc_M, fit_exc))
                #####################################################################

                #####################################################################
                # xor-tau check
                if dfg_functions.n_edges(netP,{'start'},{'end'})>0:
                    missing_exc_tau_P = 0
                    missing_exc_tau_P += max(0, sup_thr * len(self.log) - dfg_functions.n_edges(netP,{'start'},{'end'}))


                    missing_exc_tau_M = 0
                    missing_exc_tau_M += max(0, sup_thr * len(self.logM) - dfg_functions.n_edges(netM, {'start'}, {'end'}))


                    cost_exc_tau_P = missing_exc_tau_P
                    cost_exc_tau_M = missing_exc_tau_M
                    cut.append(((A.union(B), set()), 'exc2', cost_exc_tau_P, cost_exc_tau_M,cost_exc_tau_P - ratio * size_par * cost_exc_tau_M,1))
                #####################################################################

                #####################################################################
                # parallel check
                if "par" in type:
                    cost_par_P = dfg_functions.cost_par(netP, A.intersection(activitiesM), B.intersection(activitiesM), sup, feat_scores)
                    cost_par_M = dfg_functions.cost_par(netM, A.intersection(activitiesM), B.intersection(activitiesM), sup, feat_scores)
                    cut.append(((A, B), 'par', cost_par_P, cost_par_M, cost_par_P - ratio * size_par * cost_par_M,1))
                #####################################################################

                #####################################################################
                # loop check
                if "loop" in type:
                    fit_loop = dfg_functions.fit_loop(logP_var, A, B, end_A_P, start_A_P)
                    if (fit_loop > 0.0):
                        cost_loop_P = dfg_functions.cost_loop(netP, A, B, sup, start_A_P, end_A_P, input_B_P, output_B_P, feat_scores)
                        cost_loop_M = dfg_functions.cost_loop(netM, A, B, sup, start_A_M, end_A_M, input_B_M, output_B_M, feat_scores)

                        if cost_loop_P is not False:
                            cut.append(((A, B), 'loop', cost_loop_P, cost_loop_M, cost_loop_P - ratio * size_par * cost_loop_M, fit_loop))



            sorted_cuts = sorted(cut, key=lambda x: (x[4], x[2],['exc','exc2','seq','par','loop','loop_tau'].index(x[1]), -(len(x[0][0]) * len(x[0][1]) / (len(x[0][0]) + len(x[0][1])))))
            if len(sorted_cuts) != 0:
                cut = sorted_cuts[0]
            else:
                cut = ('none', 'none', 'none','none','none', 'none')

        # print(cut)

        if cut[1] == 'par':
            self.detected_cut = 'parallel'
            LAP,LBP = split.split('par', [cut[0][0], cut[0][1]], self.log, activity_key)
            LAM, LBM = split.split('par', [cut[0][0], cut[0][1]], self.logM, activity_key)
            new_logs = [[LAP,LAM],[LBP,LBM]]
            for l in new_logs:
                new_dfg = [(k, v) for k, v in dfg_inst.apply(l[0], parameters=parameters).items() if v > 0]
                activities = attributes_get.get_attribute_values(l[0], activity_key)
                start_activities = list(start_activities_get.get_start_activities(l[0], parameters=parameters).keys())
                end_activities = list(end_activities_get.get_end_activities(l[0], parameters=parameters).keys())
                self.children.append(SubtreePlain(l[0],l[1], new_dfg, self.master_dfg, self.initial_dfg, activities, self.counts,
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


def make_tree(logp, logm, dfg, master_dfg, initial_dfg, activities, c, recursion_depth, noise_threshold, start_activities,
              end_activities, initial_start_activities, initial_end_activities, parameters=None, sup= None, ratio = None, size_par = None):

    tree = SubtreePlain(logp,logm, dfg, master_dfg, initial_dfg, activities, c, recursion_depth, noise_threshold,
                        start_activities,
                        end_activities, initial_start_activities, initial_end_activities, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)

    return tree
