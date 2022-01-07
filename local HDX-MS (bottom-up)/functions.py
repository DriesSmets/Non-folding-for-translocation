import proplot as pplt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pyhdx.support import rgb_to_hex

settings = {
    'font.family': 'Arial',
    'font.size': 7}

pplt.rc.update(settings)


#node_pos = [10, 25, 40]  # in kJ/mol
node_pos = [19, 24, 33]
linear_colors = ['#ff0000', '#00ff00', '#0000ff']  # red, green, blue
#linear_colors = ['#ee3377', '#009988', '#33bbee']  # magenta, teal, cyan
linear_colors = ['#8c8c8c', '#fdae61', '#d73027'][::-1] # srinath's strange color preference

no_coverage = '#8c8c8c'

rgb_norm = plt.Normalize(node_pos[0], node_pos[-1], clip=True)
rgb_cmap = mpl.colors.LinearSegmentedColormap.from_list("rgb_cmap", list(zip(rgb_norm(node_pos), linear_colors)))
rgb_cmap.set_bad(color=no_coverage)


diff_colors = ['#54278e', '#ffffff', '#006d2c'][::-1]
diff_node_pos = [-10, 0, 10]
diff_norm = plt.Normalize(diff_node_pos[0], diff_node_pos[-1], clip=True)
diff_cmap = mpl.colors.LinearSegmentedColormap.from_list("diff_cmap", list(zip(diff_norm(diff_node_pos), diff_colors)))
diff_cmap.set_bad(color=no_coverage)

cbar_width = 0.075

dG_ylabel = 'ΔG (kJ/mol)'
ddG_ylabel = 'ΔΔG (kJ/mol)'

r_xlabel = 'Residue Number'


errorbar_kwargs = {
    'fmt': 'o',
    'ecolor': 'k',
    'elinewidth': 0.2,
    'markersize': 0,
    'alpha': 0.5
}

scatter_kwargs = {
    's': 7
}


def get_colors(array, colors='abs', no_coverage=no_coverage):
    """
    Not sure if coverage kwarg does anything since cmap has set_bad set to no coverage
    Takes df object and adds a colors column
    (candidate for adding to PyHDX)

    """

    if colors == 'abs':
        cmap = rgb_cmap
        norm = rgb_norm
    elif colors == 'diff':
        cmap = diff_cmap
        norm = diff_norm

    yvals = array
    rgba_colors = cmap(norm(yvals), bytes=True)
    hex_colors = rgb_to_hex(rgba_colors)
    hex_colors[np.isnan(yvals)] = no_coverage

    return hex_colors

def single_deltaG_scatter(ax, protein):
    yvals = protein['deltaG'] * 1e-3
    ax.errorbar(protein.index, yvals, yerr=protein['covariance']*1e-3, **errorbar_kwargs, zorder=-1)
    ax.scatter(protein.index, yvals, c=get_colors(yvals), **scatter_kwargs)

    ax.format(ylabel='', xlabel=r_xlabel)


def diff_deltaG_scatter(ax, p1, p2):
    diff = p1 - p2
    yvals = diff['deltaG'] * 1e-3
    cov = np.sqrt(p1['covariance'] ** 2 + p2['covariance'])

    ax.errorbar(diff.index, yvals, yerr=cov*1e-3, **errorbar_kwargs, zorder=-1)
    ax.scatter(diff.index, yvals, c=get_colors(yvals, colors='diff'), **scatter_kwargs)

    ax.format(ylabel='', xlabel=r_xlabel)


def single_deltaG_colorbar(fig, ax, ticks, colors='abs', **kwargs):
    ymin, ymax = ax.get_ylim()

    if colors == 'abs':
        normal_colors = linear_colors
        positions = node_pos
        label = dG_ylabel
    elif colors == 'diff':
        normal_colors = diff_colors
        positions = diff_node_pos
        label = ddG_ylabel

    assert ymin <= min(positions)
    assert ymax >= max(positions)
    extended_colors = [normal_colors[0]] + normal_colors + [normal_colors[-1]]
    extended_pos = [ymin] + positions + [ymax]


    extended_norm = plt.Normalize(ymin, ymax)
    extended_cmap = mpl.colors.LinearSegmentedColormap.from_list("extended_cmap", list(zip(extended_norm(extended_pos), extended_colors)))

    tick_pos = extended_norm(ticks)

    cbar = fig.colorbar(extended_cmap, width=cbar_width, ticks=tick_pos, label=label, space=0, tickminor=True, **kwargs)
    cbar.ax.set_yticklabels(ticks)

    ax.format(yticklabelloc='None', ytickloc='None')


def single_deltaG_inverted_colorbar(fig, ax, ticks, colors='abs', **kwargs):
    ymax, ymin = ax.get_ylim()

    if colors == 'abs':
        normal_colors = linear_colors
        positions = node_pos
        label = dG_ylabel
    elif colors == 'diff':
        normal_colors = diff_colors
        positions = diff_node_pos
        label = ddG_ylabel

    assert ymin <= min(positions)
    assert ymax >= max(positions)

    extended_colors = [normal_colors[0]] + normal_colors + [normal_colors[-1]]
    extended_pos = [ymin] + positions + [ymax]
    extended_norm = plt.Normalize(ymin, ymax)
    normed_positions = np.array(extended_norm(extended_pos))
    normed_positions = 1 - normed_positions[::-1]
    extended_cmap = mpl.colors.LinearSegmentedColormap.from_list("extended_cmap", list(zip(normed_positions, extended_colors[::-1])))

    tick_pos = 1 - np.array(extended_norm(ticks))[::-1]

    cbar = fig.colorbar(extended_cmap, width=cbar_width, ticks=tick_pos, label=label, space=0, tickminor=True, **kwargs)
    cbar.ax.set_yticklabels(ticks[::-1])

    ax.format(yticklabelloc='None', ytickloc='None')