from pm4py.objects.log.importer.xes import importer as xes_importer
from local_pm4py.algo.discovery.inductive import algorithm as inductive_miner
from pm4py.visualization.petri_net import visualizer as pn_visualizer
from local_pm4py.algo.analysis import Optimzation_Goals
from local_pm4py.algo.analysis import gui
import pm4py
from pm4py.objects.petri_net.exporter import exporter as pnml_exporter
import time

support, ratio, LPlus_LogFile, LMinus_LogFile= gui.input()
logP = xes_importer.apply(LPlus_LogFile)
logM = xes_importer.apply(LMinus_LogFile)

start = time.time()
net, initial_marking, final_marking = inductive_miner.apply_bi(logP,logM, variant= inductive_miner.Variants.IMbi, sup=support, ratio=ratio, size_par=len(logP)/len(logM))
end = time.time()

parameters = {pn_visualizer.Variants.WO_DECORATION.value.Parameters.FORMAT:"pdf"}
gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters)
pn_visualizer.view(gviz)

# file_name = "petri_r"+str(ratio)+"_s"+str(support)
# pm4py.write_pnml(net, initial_marking, final_marking, file_name)

mes = Optimzation_Goals.apply_petri(logP,logM,net,initial_marking,final_marking)
print(mes)

gui.output(str(mes['acc']), str(mes['F1']), str(mes['acc_perf']), str(mes['F1_perf']), str(mes['fitP']), str(mes['precision']), str(round(end - start,2)))

