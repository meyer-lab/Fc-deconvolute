from common import subplotLabel, getSetup
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

def makeFigure():
    #set up legend elements
       legend_elements1 = [mlines.Line2D([], [], marker='^', markeredgecolor='k', label='Complement C1q',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='d', markeredgecolor='k', label='Complement C4',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='s', markeredgecolor='k', label='FcγRIa',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='v', markeredgecolor='k', label='FcγRII',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='o', markeredgecolor='k', label='FcγRIII',markersize=18, ls = 'none', markerfacecolor=  'none')]
                   

       legend_elements2 = [mpatches.Patch(facecolor = 'lightcoral', label = 'ADCC'),
                   mpatches.Patch(facecolor = 'gold', label = 'Complement Activation'),
                   mpatches.Patch(facecolor = 'mediumturquoise', label = 'Binding')]

       #set up plot
       ax, f = getSetup((7, 7), (1, 1))
       legend1 = ax[0].legend(handles= legend_elements1, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=1, borderaxespad=0., prop= {'size': 19})
       legend2 = ax[0].legend(handles= legend_elements2, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',
                      ncol=1, borderaxespad=0., prop= {'size': 19}, handleheight = 2)

       ax[0].add_artist(legend1)
       ax[0].add_artist(legend2)

       return(f)


def makeFigure2():
    #set up legend elements
       legend_elements1 = [mlines.Line2D([], [], marker='^', markeredgecolor='k', label='F',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='o', markeredgecolor='k', label='FN',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='d', markeredgecolor='k', label='FS',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='D', markeredgecolor='k', label='FS2',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='s', markeredgecolor='k', label='FNS',markersize=18, ls = 'none', markerfacecolor=  'none'),
                   mlines.Line2D([], [], marker='X', markeredgecolor='k', label='FNS2',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='x', markeredgecolor='k', label='none',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='+', markeredgecolor='k', label='N',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='1', markeredgecolor='k', label='S',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='2', markeredgecolor='k', label='S2',markersize=18, ls = 'none', markerfacecolor=  'none'),
                   mlines.Line2D([], [], marker='3', markeredgecolor='k', label='NS',markersize=18, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='4', markeredgecolor='k', label='NS2',markersize=18, ls = 'none', markerfacecolor=  'none'),]
                   

       legend_elements2 = [mpatches.Patch(facecolor = 'orchid', label = 'G0'),
                   mpatches.Patch(facecolor = 'darkblue', label = 'G1'),
                   mpatches.Patch(facecolor = 'palegreen', label = 'G2')]

       #set up plot
       ax, f = getSetup((7, 7), (1, 1))
       legend1 = ax[0].legend(handles= legend_elements1, bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                      ncol=2, borderaxespad=0., prop= {'size': 19}, handlelength = 2)
       legend2 = ax[0].legend(handles= legend_elements2, bbox_to_anchor=(0., 1.02, .7, .102), loc='lower right',
                      ncol=1, borderaxespad=0., prop= {'size': 19}, handleheight = 2)

       ax[0].add_artist(legend1)
       ax[0].add_artist(legend2)

       return(f)