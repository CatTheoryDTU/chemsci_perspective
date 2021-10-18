import sys, os
#from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import numpy as np
from scipy.optimize import curve_fit, leastsq
from copy import deepcopy
import pickle as pkl
import matplotlib.patheffects as path_effects
from ase.units import kB

from catmaphelpers.catmap_plotter import _make_plots, _plot_model, _plot_CV, post_process_CV, _compute_tafel_prate, plot_tafel_screening, _get_unit_cell_area, _plot_pH, _plot_map, _plot_pH_bar_coverage
# from catmaphelpers.catmap_output_handler import _get_data, _sum_catmap_data, _get_rate_control, read_catmap_wdir
# from catmaphelpers.extra_tools import make_FEDs



#plt.style.use(['science','ieee'])
from scipy import stats
import numpy as np
from matplotlib import rc
plt.rc('axes', labelsize=32)    # fontsize of the x and y labels
rc('font',**{'family':'sans-serif','sans-serif':['Palatino']})
# plt.rcParams['figure.figsize'] = (8,6)
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

# plt.rcParams["figure.figsize"] = (7,7)
# plt.rcParams["font.family"] = "Times New Roman"
# plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
# plt.rcParams['xtick.labelsize'] = 18
# plt.rcParams['ytick.labelsize'] = 18
# plt.rc('text', usetex = True)
# markersize=10

#sys.path.append("../.")
#from post_catmap import read_catmap_wdir
surf='211'
AperSite = _get_unit_cell_area(a_input=3.74)
ec = 1.602176e-19; A2cm = 1e-8
T=300
#tof2cur = (1e3*ec/(A2cm*A2cm))/AperSite[str(surf)]
#tof2cur = {'C2+':8*tof2cur, 'Ac':4*tof2cur, 'tot':1.0,'EtOH':8*tof2cur,'C1':6*tof2cur,'HER':2*tof2cur,'C2H4':8*tof2cur}

if __name__ == "__main__":
    data_in = pkl.load(open('model_py3.pkl','rb'),encoding='latin1')
    print(data_in.keys())
    data={'rate':{},'cov':{}}
    pots,phs=[],[]
    nsteps=1
    irate=[1,2]
    iads=[0,1]

    for dat in data_in['rate_map']:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data['rate']:
            data['rate'][pot] = {}
        data['rate'][pot][ph] = dat[1]
        if pot not in pots:
            #pots.append(np.around(pot,2))
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)

    X=np.array(sorted(np.unique(pots)))
    Y=np.array(sorted(phs))
    if nsteps == 1:cols=1
    else: cols=2
    cols=1

    fig,ax=plt.subplots(1,cols)

    #Analytical solution
    selanal=[]
    xanal=np.linspace(-0.1,0.1,100)
    for x in xanal:
        selanal.append(1/(np.exp(-x/(kB*T))+1))
    ax.plot(xanal,selanal,'-k',lw=3,label='Analytical')
    ax.fill_between(xanal,selanal,1,color='tab:red',alpha=0.0)
    ax.fill_between(xanal,selanal,0,color='tab:blue',alpha=0.0)

    sel=[]
    #Numerical solution
    for ix,x in enumerate(X):
        for iy,y in enumerate(Y):
            sel.append(data['rate'][x][y][irate[0]]/np.sum([data['rate'][x][y][i] for i in irate])*1)

    ax.plot(X/10.,sel,'o',color='cyan',markersize=8,markeredgecolor='Black',alpha=1.0,label='MKM')


    # Text
#     text=ax.annotate(r'P$_2$',(-0.05,0.7),color='r',fontweight='bold',
#            fontsize=30,ha='center',bbox=dict(facecolor='w', edgecolor=None))
#     text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),#
#                   path_effects.Normal()])
#     text=ax.annotate(r'P$_1$',(0.05,0.3),color='b',fontweight='bold',
#            fontsize=30,ha='center',bbox=dict(facecolor='w',edgecolor=None))
#     text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
#                   path_effects.Normal()])
    # Reaction scheme
    text=ax.annotate(r'R',(-0.08,0.8),color='k',fontweight='bold',
           fontsize=10,ha='left',va='center')
    text=ax.annotate(r'RDS',(-0.07,0.75),color='k',fontweight='bold',
           fontsize=7,ha='left',va='center')
    text=ax.annotate(r'SDS',(-0.025,0.8),color='k',fontweight='bold',
           fontsize=7,ha='left',va='center')
    text=ax.annotate(r'A',(-0.05,0.8),color='k',fontweight='bold',
           fontsize=10,ha='left',va='center')
    text=ax.annotate(r'P$_1$',(-0.01,0.85),color='k',fontweight='bold',
           fontsize=10,ha='left',va='center')
    text=ax.annotate(r'P$_2$',(-0.01,0.7),color='k',fontweight='bold',
           fontsize=10,ha='left',va='center')
    ax.arrow(-0.072,0.8,0.02,0.0,color='k',head_length=0.005,length_includes_head=True,
           head_width=0.01)
    ax.arrow(-0.04,0.8,0.02,0.05,color='k',head_length=0.02,length_includes_head=True,
           head_width=0.005)
    ax.arrow(-0.04,0.8,0.02,-0.05,color='k',head_length=0.02,length_includes_head=True,
           head_width=0.005)

    #Labels
    ax.set_xlabel('$\Delta G^\ddagger_{P_2}-\Delta G^\ddagger_{P_1}$(eV)',fontsize=16)
    ax.set_ylabel('Selectivity towards P$_1$',fontsize=16)
    ax.set_xlim([-0.10,0.10])
    ax.legend(loc='lower right')
    ax.set_ylim([0,1])
    fig.savefig('selectivity.png',dpi=1000)
