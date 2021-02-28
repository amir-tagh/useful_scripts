#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches  as patches
from   matplotlib import colors, ticker, cm, rc
from   matplotlib.colors import ListedColormap
from   scipy.interpolate import spline
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   mpl_toolkits.axes_grid1 import ImageGrid
import sys
import glob
import time

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

font1 = 12
font2 = 16

import seaborn as sns
sns.set(style='ticks')
sns.set(rc={"figure.figsize": (6, 8)})
plt.clf()
f=plt.figure(1, figsize=(6, 8))

#f,  ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
#axarr = [ax1, ax2, ax3, ax4]


axarr = ImageGrid(f, 111,
                  nrows_ncols=(3,2),
                  axes_pad=0.275, aspect=False,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="4%",
                  cbar_pad=0.12)
                   


threeSpace = cm.get_cmap('rainbow', 3)  

x     = np.linspace(1, 24, 24)


i_ax   = 0
csList = []
for stem in [ 'arg', 'plain', 'et']:
    for seq in [ 'gggggg', 'ggcgac' ]:
        if seq == "gggggg":
            fullSeq   = 'GGGGGGGGGGGGCCCCCCCCCCCC'
            fullSeq_2 = 'CCCCCCCCCCCCGGGGGGGGGGGG'
        elif seq == 'ggcgac':
            fullSeq    = 'GGCGGCGGCGGCGACGACGACGAC'
            fullSeq_2  = 'CCGCCGCCGCCGCTGCTGCTGCTG'
        
        #work_distance_file="%s_%s_swails.dat" % (seq, stem) 
        #work = np.loadtxt(work_distance_file)
        #for i in range(len(work)):
        #    if force[i] 

        #dist_vs_t_stretch7_joined.dat
        ##extension per frame? 200 frames per sub-run.
        extnFile  = "%s_%s/dist_vs_t_stretch15_joined.dat" % (seq,stem)
        print ('file name', extnFile)
 	
        #extnFile_1  = "%s_%s/dist_vs_t_stretch15_14.dat"
        #print (extnFile_1,extnFile)
        #time.sleep(3)       
        xtnseries = np.loadtxt(extnFile)
        delta_xtn = -1*(xtnseries[-1,0] - xtnseries[0,0])/200.
        print('delta xtn:', delta_xtn)
        base_xtn  = xtnseries[0,0]
        print('base extension:',base_xtn)
        time.sleep(5)
        print("calculated delta_xtn: %.9e Angstrom/frame base_ext: %.9e A" %\
              ( delta_xtn, base_xtn ) )
        n_frames  = 15*200
        max_extn  = n_frames * delta_xtn 
        print('max extn:',max_extn)
        time.sleep(5)
        mapFile   = "curves/%s_%s/meanRise.dat" % (seq,stem)
        popMap    = np.loadtxt(mapFile)
        print("expect %i frames and 23 bp steps" % n_frames)
        print(np.shape(popMap))

        ##for quicker plot
        popMap = popMap[::10,:]

        n_tpoints = np.shape(popMap)[0]
        n_steps   = np.shape(popMap)[1]
        
#        y     = np.linspace(0., max_extn/23., n_tpoints)
        y     = np.linspace(1., (max_extn+base_xtn)/base_xtn, n_tpoints)
        X, Y  = np.meshgrid(x, y)

        axarr[i_ax].set_xlim([1,24])
#        axarr[i_ax].set_ylim([0,max_extn/23.])
        axarr[i_ax].set_ylim([1.,(max_extn+base_xtn)/base_xtn])
        ##axarr[0].set_xlabel(r"$L_z$")

        if i_ax % 2 == 0:
            axarr[i_ax].set_ylabel(r"Extension [\AA~per bp step]")
        if i_ax >= 2:
            ticList = list(fullSeq)
            for i in range(len(ticList)):
                ticList[i] = "%s\n%s" % (fullSeq[i], fullSeq_2[i])

            axarr[i_ax].set_xticks(x)
            axarr[i_ax].set_xticklabels(ticList)

        if stem == 'arg':
            axarr[i_ax].set_title(r"Arginine")
        elif stem == 'plain':
            axarr[i_ax].set_title(r"No Intercalator")
        elif stem == 'et':
            axarr[i_ax].set_title(r"EtBr")
            

        cs = axarr[i_ax].pcolor( X[:], Y[:], popMap[:,:],rasterized=True,
                                     linewidth=0, cmap="gnuplot2")
        cs.set_edgecolor('face')
        csList.append(cs)

        i_ax += 1

#divider = make_axes_locatable(axarr[1])
#cax     = divider.append_axes("right", size="4%", pad=0.1)
#cbar    = f.colorbar(csList[-1], cax=cax,  label=r'Rise [\AA]')
#cbar    = f.colorbar(cs, ax=axarr, label=r'Rise [\AA]')
#cax, kw  = mpl.colorbar.make_axes([ax for ax in axarr])
#plt.colorbar(cs, cax=cax, **kw)

bar = axarr[-1].cax.colorbar(cs)  
axarr[-1].cax.toggle_label(True)
axarr[-1].cax.set_label(r'Rise [\AA]')
axarr[-1].cax.set_title(r'Rise [\AA]')
bar.solids.set_edgecolor("face")
bar.solids.set_rasterized(True)

chars='abcdef'
for i in range(len(axarr)):
    if i % 2 == 1:
        x = -1.5
    else:
        x = -3.25
    axarr[i].text(x, -0.2, r'({\bf %s})'%chars[i])

#plt.tight_layout()
plt.savefig("kymo_extension_TEST.pdf")
