import yt
#yt.enable_parallelism()
import ytree
from yt import YTArray
from yt.visualization.volume_rendering.api import Scene, Camera, create_volume_source
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import VolumeSource

import unyt
from unyt import Mpc, kpc
from yt.utilities.parallel_tools.parallel_analysis_interface \
    import communication_system
from yt.utilities.cosmology import Cosmology

import time
import pandas as pd
import numpy as np
from numpy.polynomial import polynomial
import datetime
import math
import warnings

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.animation as animation
import matplotlib.image as mpimg
from matplotlib.pyplot import figure
from matplotlib import rc_context, cm
from matplotlib.colors import Normalize 
    
from tqdm import tqdm
# to create smoother camera path:
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline, interpn
from scipy.stats import gaussian_kde

from IPython.display import Image

from main import ExploreHalo


co = Cosmology()


today = datetime.datetime.now()


def make_grid_name(fnumber=int(101)):
    # this makes the file name for the covering grids directory
    return yt.load(f"/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_{fnumber:03}_covering_grid.h5")


class DatasetBasics:
    def __init__(self, snap_number=101):
        self.ds = make_grid_name(snap_number)
        self.ad = self.ds.all_data()                                # this is a region describing the entire box
        self.redshift = self.ds.current_redshift
        self.time = self.ds.current_time.in_units("Gyr")
        
    def plot_quant(self):
        snaps = np.arange(102)
        redshifts = np.empty(len(snaps))
        times = np.empty(len(snaps))
        
        for i, snap in enumerate(snaps):
            now_ds = make_grid_name(snap)
            times[i] = now_ds.current_time.in_units("Gyr")
            redshifts[i] =  now_ds.current_redshift
        
        fig = plt.figure(figsize=(9, 5))
        plt.title('Redshift corresponding to each Snapshot Number in the dataset')
        plt.xlabel('Snapshot Number')
        plt.ylabel('Redshift')
        # plt.xlim(np.min(snaps))
        # plt.ylim(np.min(redshifts))
        plt.scatter(snaps, redshifts, marker='.', c='black')#, s = 0.6)
        plt.plot(snaps, redshifts, linewidth=0.5, linestyle='--', c='k', alpha=0.6)
        plt.grid(which='minor', alpha=0.4)
        plt.minorticks_on()
        plt.xlim(0)
        plt.ylim(0)
        plt.savefig('SnapRed')

        plt.clf()

        fig = plt.figure(figsize=(9, 5))
        plt.title('Cosmic Time corresponding to each Snapshot Number in the dataset')
        plt.xlabel('Snapshot Number')
        plt.ylabel('Cosmic Time (Gyr)')
        # plt.xlim(np.min(snaps))
        # plt.ylim(np.min(times))
        plt.scatter(snaps, times, marker='.', c='black')#, s = 0.6)
        plt.plot(snaps, times, linewidth=0.5, linestyle='--', c='k', alpha=0.6)
        plt.grid(which='minor', alpha=0.4)
        plt.minorticks_on()
        plt.xlim(0)
        plt.ylim(0)
        plt.savefig('SnapTime')
        plt.clf()