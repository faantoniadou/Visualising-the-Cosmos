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


class SimplePlots(ExploreHalo):
    '''
    Here we include projections and slices 
    that can be either a movie or a multiplot or a single plot
    we can also play around with focus but by default that is the origin if not we need input
    '''
    def __init__(self, fnumber=101):           ### check if units are right
        self.fnumber = fnumber
        ExploreHalo.__init__(self, self.fnumber)
    

    def time_slices(self): 
        fns = [
        "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_000_covering_grid.h5",
        "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_033_covering_grid.h5",
        "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_067_covering_grid.h5",
        "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/covering_grids/snapshot_101_covering_grid.h5",
        ]

        fig = plt.figure()

        # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
        # These choices of keyword arguments produce a four panel plot with a single
        # shared narrow colorbar on the right hand side of the multipanel plot. Axes
        # labels are drawn for all plots since we're slicing along different directions
        # for each plot.
        grid = AxesGrid(
            fig,
            (0.075, 0.075, 0.85, 0.85),
            nrows_ncols=(2, 2),
            axes_pad=0.05,
            label_mode="L",
            share_all=True,
            cbar_location="right",
            cbar_mode="single",
            cbar_size="3%",
            cbar_pad="0%",
        )


        for i, fn in enumerate(fns):
            # Load the data and create a single plot
            self.ds = make_grid_name(fn)
            p = yt.SlicePlot(ds, "z", self.field, width=(80, "Mpc"))
            p.set_cmap(field=self.field, cmap="octar")
            # Ensure the colorbar limits match for all plots
            p.set_zlim(self.field, 1e36, 1e46)
            # nicen up the plot by setting the background color to the minimum of the colorbar
            p.set_background_color(self.field)
            # hide the colorbar:
            p.hide_colorbar()

            # hide the axes, while still keeping the background color correct:
            p.hide_axes(draw_frame=True)
            p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
            p.annotate_scale(corner="upper_right")

            # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            plot = p.plots[self.field]
            plot.figure = fig
            plot.axes = grid[i].axes
            plot.cax = grid.cbar_axes[i]

            # Finally, this actually redraws the plot.
            p._setup_plots()
            plt.savefig("{:%Y_%m_%d_%H_%M_%S}".format(today)+"_multi.png")

