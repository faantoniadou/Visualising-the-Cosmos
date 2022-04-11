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

yt.enable_parallelism()

def make_grid_name(fnumber):
    grid_name = f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5"
    return grid_name

def main():
    fnumber = int(input('Grid number: '))
    frame = make_grid_name(fnumber)
    ds = yt.load(frame)
    initialWidth = float(ds.domain_width.in_units(kpc)[0]) 
    print(initialWidth)
    
    rho, c = ds.find_max("nbody_mass") # find the highest density peak (value and location)
    print(rho, c)
    
    sc = yt.create_scene(ds)
    sc.camera.resolution = (1920, 1080)
    
    sc.camera.set_focus(c) # focus on the highest density peak
    source, bounds = sc[0], (2e-28, 1e-2) # very large range of densities 
    source.set_field(("grid", "nbody_mass")) # field to render
    tf = yt.ColorTransferFunction(x_bounds = np.log10(bounds))
    tf.add_layers(N=20, w=0.03, colormap='plasma') # add 20 Gaussians filters 
    source.tfh.tf, source.tfh.bounds = tf, bounds # tfh stands for TransferFunctionHelper 
    source.tfh.plot('transferFunction.png', profile_field='nbody_mass')
    
    # in 1795 log steps change the window from 97.8 kpc down to 9.78e-11 kpc = 0.0202 AU = 3.02e6 km 
    for i, coef in enumerate(np.linspace(start=0,stop=12,num=1795)): # i=0..1794
        width = initialWidth/10.**coef # width of the visualization window 
        sc.camera.set_width(((width*192/108,'kpc'),(width,'kpc'),(width,'kpc'))) 
        sc.save('frame%04d.png' % (i+1), sigma_clip=4)

main()