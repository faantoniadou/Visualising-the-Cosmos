import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
import yt
from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource
import datetime
import unyt
from unyt import Mpc, kpc
from matplotlib import rc_context
#from tqdm import tqdm
import warnings
import yt
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system

#yt.enable_parallelism()

comm = communication_system.communicators[-1]


# this makes the file name for the covering grids directory
def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name

# this makes file name for halo trees 
def get_halos():
    halo_dir = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5"
    return halo_dir



class Camera():
    
    def get_halo_centre():
        # Finds halo with the maximum of the quantity (field) the user provides.
        self.halo_dir = get_halos()
        self.a = ytree.load("/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/mergertree_h5/rockstar/rockstar.h5")

        # quantity maximum halo to find:
        self.quantity = str(input('Maximum quantity to search for in halo catalogue: '))
        
        # find most massive halo
        self.i = np.argmax(a[quantity])             #['mass'])

        my_tree = self.a[self.i]

        # some other useful stuff:
        # print(f'Position of most massive progenitor: {my_tree["prog", "position"]}') # give position of most massive progenitor
        # print(f'Redshifts of progenitors: {my_tree["prog", "redshift"]}')

        halo_position = my_tree["position"]
        radius = my_tree["virial_radius"].to("unitary")

        return halo_position, radius
    
    
    