import yt
#from yt.mods import *           # deprecated
#yt.enable_parallelism()
import ytree

import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg

from yt import YTArray
#from yt import YTRegionBase
from yt.visualization.volume_rendering.api import Scene, Camera
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource
import datetime
import unyt
from unyt import Mpc, kpc
from matplotlib import rc_context

import warnings
import yt
# from yt.utilities.parallel_tools.parallel_analysis_interface \
#     import communication_system

quantity = 'mass'#str(input('Maximum quantity to search for in halo catalogue: '))


a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")
# a.add_analysis_field("significance", "Msun*Myr")

# ap = ytree.AnalysisPipeline()
# ap.add_operation(calc_significance)


i = np.argmax(a[quantity])
    
# do it for one specific tree
my_tree = a[i]

my_ancestors = list(my_tree.ancestors)
print(my_tree.descendent)
#print (my_ancestors)

def calc_significance(node):
    if node.descendent is None:
        dt = 0. * node["time"]
    else:
        dt = node.descendent["time"] - node["time"]
        
    sig = node["mass"] * dt
    ancestor_list = list(node.ancestors)
    sig_ancestors = []

    
    if node.ancestors is not None:
        sig_list = []
        anc_list = []
        next_anc = 0
        # get ancestors with greatest significance
        
    # for halo in a.select_halos("(tree['forest', 'mass'].to('Msun') > 6e14) & (tree['forest', 'mass'].to('Msun') > 690072000000000.0)"):        # exclude the most massive halo
        # pair_pos = halo["position"].to("unitary")
        # distance = np.linalg.norm(np.array(pair_pos) - np.array(initial_position))
        # positions.append(np.array(pair_pos))
        # distances.append(distance)
    
    
        for anc in ancestor_list:
            sig += calc_significance(anc)[0]
            sig_list.append(sig.value)
            anc_list.append(anc)
        
            if len(sig_list) != 0:
                i = np.argmax(sig_list)
                next_anc = anc_list[i]
            else:
                next_anc = anc
            print(next_anc['position'])
            sig_ancestors.append(next_anc['position'])
        return sig
    return sig_ancestors


print(calc_significance(my_tree))


#for tree in ytree.parallel_nodes(trees):
# for tree in trees:
# ap.process_target(trees)
#print(trees.descendent)

# fn = a.save_arbor(filename="significance", trees=trees)
# a2 = ytree.load(fn)
# a2.set_selector("max_field_value", "significance")
# prog = list(a2[0]["prog"])
# print (prog)

# p = ytree.TreePlot(a[0], dot_kwargs={'rankdir': 'LR', 'size': '"12,4"'})
# p.save('tree.png')

