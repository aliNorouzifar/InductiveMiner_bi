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
import pkgutil

from pm4py import util as pmutil
from pm4py.algo.discovery.dfg.variants import native as dfg_inst
from local_pm4py.algo.discovery.inductive.util import tree_consistency
from local_pm4py.algo.discovery.inductive.util.petri_el_count import Counts
from local_pm4py.algo.discovery.inductive.variants.im_bi.data_structures import subtree_plain as subtree
from local_pm4py.algo.discovery.inductive.variants.im_bi.util import get_tree_repr_implain
from pm4py.objects.conversion.log import converter
from pm4py.objects.conversion.process_tree import converter as tree_to_petri
from pm4py.objects.log.obj import EventLog, Trace, Event
from pm4py.objects.log.util import filtering_utils
from pm4py.objects.process_tree.utils import generic
from pm4py.objects.process_tree.utils.generic import tree_sort
from pm4py.statistics.attributes.log import get as attributes_get
from pm4py.statistics.end_activities.log import get as end_activities_get
from pm4py.statistics.start_activities.log import get as start_activities_get
from pm4py.util import exec_utils
from pm4py.util import variants_util
from pm4py.util import xes_constants
from pm4py.util import constants
from enum import Enum
import deprecation

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


def apply(logp,logm, parameters=None,sup= None, ratio = None, size_par = None):
    """
    Apply the IM algorithm to a log obtaining a Petri net along with an initial and final marking

    Parameters
    -----------
    log
        Log
    parameters
        Parameters of the algorithm, including:
            Parameters.ACTIVITY_KEY -> attribute of the log to use as activity name
            (default concept:name)

    Returns
    -----------
    net
        Petri net
    initial_marking
        Initial marking
    final_marking
        Final marking
    """
    if pkgutil.find_loader("pandas"):
        import pandas as pd
        from pm4py.statistics.variants.pandas import get as variants_get


    net, initial_marking, final_marking = tree_to_petri.apply(apply_tree(logp,logm, parameters,sup= sup, ratio = ratio, size_par = size_par))
    return net, initial_marking, final_marking


def apply_variants(variants, parameters=None):
    """
    Apply the IM algorithm to a dictionary of variants, obtaining a Petri net along with an initial and final marking

    Parameters
    -----------
    variants
        Variants
    parameters
        Parameters of the algorithm, including:
            Parameters.ACTIVITY_KEY -> attribute of the log to use as activity name
            (default concept:name)

    Returns
    -----------
    net
        Petri net
    initial_marking
        Initial marking
    final_marking
        Final marking
    """
    net, im, fm = tree_to_petri.apply(apply_tree_variants(variants, parameters=parameters))
    return net, im, fm


@deprecation.deprecated('2.2.10', '3.0.0', details='use newer IM implementation (IM_CLEAN)')
def apply_tree(logp,logm, parameters=None, sup= None, ratio = None, size_par = None):
    """
    Apply the IM algorithm to a log obtaining a process tree

    Parameters
    ----------
    log
        Log
    parameters
        Parameters of the algorithm, including:
            Parameters.ACTIVITY_KEY -> attribute of the log to use as activity name
            (default concept:name)

    Returns
    ----------
    process_tree
        Process tree
    """
    if parameters is None:
        parameters = {}

    if pkgutil.find_loader("pandas"):
        import pandas as pd
        from pm4py.statistics.variants.pandas import get as variants_get


    activity_key = exec_utils.get_param_value(Parameters.ACTIVITY_KEY, parameters,
                                              pmutil.xes_constants.DEFAULT_NAME_KEY)

    # # since basic IM is influenced once per variant, it makes sense to keep one trace per variant
    # log = filtering_utils.keep_one_trace_per_variant(log, parameters=parameters)

    # keep only the activity attribute (since the others are not used)
    # logp = filtering_utils.keep_only_one_attribute_per_event(logp, activity_key)
    # logm = filtering_utils.keep_only_one_attribute_per_event(logm, activity_key)

    dfgp = [(k, v) for k, v in dfg_inst.apply(logp, parameters=parameters).items() if v > 0]
    dfgm = [(k, v) for k, v in dfg_inst.apply(logm, parameters=parameters).items() if v > 0]

    c = Counts()
    activitiesp = attributes_get.get_attribute_values(logp, activity_key)
    start_activitiesp = list(start_activities_get.get_start_activities(logp, parameters=parameters).keys())
    end_activitiesp = list(end_activities_get.get_end_activities(logp, parameters=parameters).keys())
    contains_empty_traces = False
    traces_length = [len(trace) for trace in logp]
    if traces_length:
        contains_empty_traces = min([len(trace) for trace in logp]) == 0

    recursion_depth = 0
    sub = subtree.make_tree(logp,logm, dfgp, dfgp, dfgp, activitiesp, c, recursion_depth, 0.0, start_activitiesp,
                            end_activitiesp,
                            start_activitiesp, end_activitiesp, parameters, sup= sup, ratio = ratio, size_par = size_par)

    process_tree = get_tree_repr_implain.get_repr(sub, 0, contains_empty_traces=contains_empty_traces)
    # Ensures consistency to the parent pointers in the process tree
    tree_consistency.fix_parent_pointers(process_tree)
    # Fixes a 1 child XOR that is added when single-activities flowers are found
    tree_consistency.fix_one_child_xor_flower(process_tree)
    # folds the process tree (to simplify it in case fallthroughs/filtering is applied)
    process_tree = generic.fold(process_tree)
    # sorts the process tree to ensure consistency in different executions of the algorithm
    tree_sort(process_tree)

    return process_tree


def apply_tree_variants(variants, parameters=None):
    """
    Apply the IM algorithm to a dictionary of variants obtaining a process tree

    Parameters
    ----------
    variants
        Variants
    parameters
        Parameters of the algorithm, including:
            Parameters.ACTIVITY_KEY -> attribute of the log to use as activity name
            (default concept:name)

    Returns
    ----------
    process_tree
        Process tree
    """
    log = EventLog()
    activity_key = exec_utils.get_param_value(Parameters.ACTIVITY_KEY, parameters, xes_constants.DEFAULT_NAME_KEY)

    var_keys = list(variants.keys())
    for var in var_keys:
        trace = Trace()
        activities = variants_util.get_activities_from_variant(var)
        for act in activities:
            trace.append(Event({activity_key: act}))
        log.append(trace)

    return apply_tree(log, parameters=parameters)
