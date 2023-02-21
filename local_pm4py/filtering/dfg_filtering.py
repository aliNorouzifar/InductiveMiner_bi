
from copy import deepcopy
import networkx as nx

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


def noise_filtering(dfg0, nt):
    dfg = deepcopy(dfg0)
    log_size = sum([dfg[x] for x in dfg if x[0] == 'start'])
    noisy_edges = sorted([(x,dfg[x]) for x in dfg if (dfg[x]/log_size) < nt], key=lambda z:z[1])
    net = generate_nx_graph_from_dfg(dfg0)
    for ne in noisy_edges:
        net_copy = deepcopy(net)
        nodes_set = set(net_copy.nodes)
        net_copy.remove_edge(ne[0][0],ne[0][1])
        if (set(nx.ancestors(net_copy, 'end')) == nodes_set-{'end'}):
            if(set(nx.descendants(net_copy, 'start')) == nodes_set-{'start'}):
                # print(ne)
                del dfg[ne[0]]
                net = net_copy
    return dfg