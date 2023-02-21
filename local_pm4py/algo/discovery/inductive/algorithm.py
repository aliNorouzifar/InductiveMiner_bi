from enum import Enum

from local_pm4py.algo.discovery.inductive.variants.im_bi import algorithm as im_bi_algo
from pm4py.util import exec_utils
from typing import Optional, Dict, Any, Union, Tuple, List
from pm4py.objects.log.obj import EventLog, EventStream
import pandas as pd
from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.process_tree.obj import ProcessTree



class Variants(Enum):
    IMbi = im_bi_algo


IMbi = Variants.IMbi

DEFAULT_VARIANT_LOG = IMbi


def apply(log: Union[EventLog, EventStream, pd.DataFrame], parameters: Optional[Dict[Any, Any]] = None, variant=DEFAULT_VARIANT_LOG) -> Tuple[PetriNet, Marking, Marking]:
    """
    Apply the chosen IM algorithm to a log obtaining a Petri net along with an initial and final marking

    Parameters
    -------------
    log
        Log
    variant
        Variant of the algorithm to apply, possible values:
        - Variants.IMd
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
    print('hi')
    return exec_utils.get_variant(variant).apply(log, parameters=parameters)


def apply_bi(logp: Union[EventLog, EventStream, pd.DataFrame], logm: Union[EventLog, EventStream, pd.DataFrame], parameters: Optional[Dict[Any, Any]] = None, variant=DEFAULT_VARIANT_LOG, sup= None, ratio = None, size_par = None) -> Tuple[PetriNet, Marking, Marking]:
    """
    Apply the chosen IM algorithm to a log obtaining a Petri net along with an initial and final marking

    Parameters
    -------------
    log
        Log
    variant
        Variant of the algorithm to apply, possible values:
        - Variants.IMd
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
    return exec_utils.get_variant(variant).apply(logp,logm, parameters=parameters, sup= sup, ratio = ratio, size_par = size_par)



def apply_tree(log: Union[EventLog, EventStream, pd.DataFrame], parameters: Optional[Dict[Any, Any]] = None, variant=DEFAULT_VARIANT_LOG) -> ProcessTree:
    """
    Apply the chosen IM algorithm to a log obtaining a process tree

    Parameters
    ----------
    log
        Log
    variant
        Variant of the algorithm to apply, possible values:
        - Variants.IMd
    parameters
        Parameters of the algorithm, including:
            Parameters.ACTIVITY_KEY -> attribute of the log to use as activity name
            (default concept:name)

    Returns
    ----------
    tree
        Process tree
    """
    return exec_utils.get_variant(variant).apply_tree(log, parameters=parameters)


