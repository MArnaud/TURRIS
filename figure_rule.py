import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib

#customize font
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

#customize axis
def ax_settings_plot(ax):
    #ax.set_yticks([])
    
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.spines['bottom'].set_edgecolor('#444444')
    ax.spines['left'].set_edgecolor('#444444')
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    
    return None


def ax_settings_plot2(ax):
    ax.set_yticks([])
    
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    ax.spines['top'].set_edgecolor('#444444')
    ax.spines['top'].set_linewidth(1.5)
    ax.xaxis.tick_top()


def box_style(bp):
	for box in bp['boxes']:
		box.set( facecolor = 'silver' )	
	for median in bp['medians']:
		median.set(color='k', linewidth=2)
	for flier in bp['fliers']:
		flier.set(marker='o', color='darkgray')
	for means in bp['means']:
		means.set(color='b', linewidth=2)


