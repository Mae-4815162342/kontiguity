"""Plots for prediciton results"""
import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon
from matplotlib import gridspec
import cmasher as cmr
from scipy import interpolate
from scipy.signal import find_peaks
import os


def compute_contact(c, chrom, plasmid):
    d=c.info
    total_reads = d['sum']

    mat = c.matrix().fetch(chrom, plasmid)
    mat[np.isnan(mat)] = 0
    coverage_plasmid = mat.sum(axis=1)
    coverage_plasmid[coverage_plasmid==0] = np.nan

    matscn = c.matrix().fetch(chrom, chrom)
    matscn[np.isnan(matscn)] = 0
    coverage = matscn.sum(axis=0)
    coverage[coverage==0] = np.nan

    m1 = c.matrix(balance=False).fetch(chrom, chrom)
    reads_chrom= np.nansum(m1)
    m2 = c.matrix(balance=True).fetch(plasmid, plasmid)
    reads_plasmid= np.nansum(m2)
    m12 = c.matrix(balance=False).fetch(chrom, plasmid)
    reads_chr12= np.nansum(m12)

    coverage_plasmid = (coverage_plasmid/reads_chrom)/(reads_plasmid/reads_chr12)

    return coverage_plasmid

def compute_4C(c, chromosomes, plasmid):
    computed_4C_signal = []
    positions = {}
    i = 0
    for chromosome in chromosomes:
        if chromosome != plasmid:
            signal = compute_contact(c, chromosome, plasmid)
            computed_4C_signal += list(signal)
            positions[chromosome] = (i, i + len(signal))
            i += len(signal)

    return np.array(computed_4C_signal), positions#np.array(pd.DataFrame(np.array(computed_4C_signal).T, columns=["Signal"]).rolling(10).median()["Signal"]), positions

def plot_4C(signal, positions, contig_name, save_to=""):
    plt.figure(figsize=(30,3))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax_signal = plt.subplot(gs[0])
    ax_signal.set_title(f"{contig_name} contact with the chromosomes")
    ax_signal.plot(signal, color='black')
    ax_signal.set_ylabel("Contact signal")
    index = []
    label = []
    colors = cmr.get_sub_cmap('viridis', 0.25, 1, N=len(positions)).colors
    i = 0
    for chromosome in positions.keys():
        start, stop = positions[chromosome]
        rec = Rectangle((start, -1), stop - start, 5, color = colors[i])
        ax_signal.add_patch(rec)
        index += [start + ((stop - start) // 2)]
        label += [chromosome]
        i += 1
    ax_signal.set_xticks(index, label, rotation=45)
    if not np.isnan(np.nanmax(signal)) and not np.isinf(np.nanmax(signal)):
        ax_signal.set_ylim(- np.nanmax(signal)*0.1, np.nanmax(signal)*1.25)
    if len(save_to) > 0:
        plt.savefig(save_to + '.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(save_to + '.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()
    

def plot_contact(ax, contact, centromere=None, size = 1):
    to_plot = [contact for i in range(int(7 * size))]
    ax.axis("off")
    ax.set_ylim((7*size, -size))

    if centromere == None:
        patch = FancyBboxPatch((0, 0), len(contact) - 1, len(to_plot) - 1, boxstyle=f'round,rounding_size={3*size}', fill=None)
        ax.add_patch(patch)
        image = ax.imshow(to_plot, cmap='YlOrBr', clip_path=patch, clip_on=True)#, interpolation="gaussian", vmin=0, vmax=0.001)
    else:
    
        patch = FancyBboxPatch((0, 0), centromere, len(to_plot) - 1, boxstyle=f'round,rounding_size={3*size}', fill=None)
        ax.add_patch(patch)
        ax.imshow(to_plot, cmap='YlOrBr', clip_path=patch, clip_on=True, interpolation="gaussian")
        patch2 = FancyBboxPatch((centromere, 0), len(contact) - centromere - 1, len(to_plot) - 1, boxstyle=f'round,rounding_size={3*size}', fill=None)
        ax.add_patch(patch2)
        image = ax.imshow(to_plot, cmap='YlOrBr', clip_path=patch2, clip_on=True, interpolation="gaussian")

        # central bin to cover boxes joining
        cmap = image.get_cmap()
        r,g,b,_ = cmap(image.norm(contact[centromere]))
        patch3 = Rectangle((centromere - 0.45, 1.25 * size), 0.9, len(to_plot) - (3 * size), color=(r,g,b))
        ax.add_patch(patch3)

        # triangles
        decalage = 0.25 * size
        patch4 = Polygon([[centromere - decalage, 1.25 * size], [centromere + decalage, 1.25 * size], [centromere, len(to_plot)/2]], color="w")
        ax.add_patch(patch4)
        patch5 = Polygon([[centromere - decalage, len(to_plot) - (1.5 * size)], [centromere + decalage, len(to_plot) - (1.5 * size)], [centromere, len(to_plot)/2]], color="w")
        ax.add_patch(patch5)

        # triangle borders
        ax.plot([centromere - decalage, centromere + decalage], [1.25 * size, len(to_plot) - (1.5 * size)], color="black")
        ax.plot([centromere + decalage, centromere - decalage], [1.25 * size, len(to_plot) - (1.5* size)], color="black")
    return image

def plot_all_contacts(contacts, save_to=""):
    gs = gridspec.GridSpec(len(contacts), 3, height_ratios=[1]*len(contacts), width_ratios=[1, 50, 0.5])
    chromosomes = list(contacts.keys())
    max_size = np.max([len(contacts[chrom]) for chrom in chromosomes])
    ax0 = None
    for i in range(len(chromosomes)):
        if ax0 == None:
            axplot = plt.subplot(gs[i, 1])
            ax0 = axplot
        else:
            axplot = plt.subplot(gs[i, 1], sharex=ax0)
        axplot.set_xlim(-1, max_size)
        axtext = plt.subplot(gs[i, 0])
        axtext.text(0.5,0.5,chromosomes[i])
        axtext.axis("off")
        image = plot_contact(axplot, contacts[chromosomes[i]], size = 4)
    ax_colorbar = plt.subplot(gs[:, -1])
    ax_colorbar.axis("off")
    plt.colorbar(image, ax=ax_colorbar)
    if len(save_to) > 0:
        plt.savefig(save_to + '.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(save_to + '.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()

def plot_4C_signals(signal, pos, select_name, save_to=""):
    """Plots the signals """
    pass

def plot_plasmids_4C(cool_path, binning, chromosomes, contigs_selection, outpath = ""):
    cool_file = f"{cool_path}::resolutions/{binning}"
    file_name = cool_path.split('/')[-1].split('.')[0]

    data = cooler.Cooler(cool_file)
    for select_name in contigs_selection:
        signal, pos = compute_4C(data, chromosomes, select_name)
        if not np.isinf(np.nansum(signal)) and not np.isnan(np.nansum(signal)) and np.nansum(signal) > 0:
            save_to = f'{outpath}/{file_name}_{select_name}' if len(outpath) > 0 else ""
            plot_4C(signal, pos, select_name, save_to=save_to)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

def get_ax_size(ax):
    fig = ax.get_figure()
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

def draw_centromere(ax_chromosome, centromere, signal, index_value):
    """Draws a triangle to annotate centromeric region."""
    # computing centromere index
    k = 0
    while k < len(index_value) - 1:
        if centromere < index_value[k + 1]:
            break
        k += 1
            
    # computing y coordinates
    yvalue = signal[k]
    ymin, ymax = ax_chromosome.get_ylim()
    ydif = ymax - ymin
    yvalue_norm = (yvalue - ymin)/ydif
    y1 = (yvalue_norm + 0.2) * ydif + ymin
    y2 = (yvalue_norm + 0.4) * ydif + ymin
    
    while y2 > ymax:
        ax_chromosome.set_ylim(ymin, y2)
        ydif = y2 - ymin
        ymax = y2
        yvalue_norm = (yvalue - ymin)/ydif
        y1 = (yvalue_norm + 0.2) * ydif + ymin
        y2 = (yvalue_norm + 0.4) * ydif + ymin
    
    # computing x coordinates
    xsize, ysize = get_ax_size(ax_chromosome)
    xydif = xsize / ysize
    xmin, xmax = ax_chromosome.get_xlim()
    xdif = xmax - xmin
    to_add = (0.1 / xydif) * xdif
    x1 = k
    x2, x3 = x1 - to_add, x1 + to_add
    
    arrow = Polygon([[x1, y1], [x2, y2], [x3, y2]], color = "black")
    ax_chromosome.add_patch(arrow)

def plot_signals(signals, signal_names, signal_limits, centromeres = [], title = "", save_to = "", xlabel = "Genomic coordinates", ylabel = "Contact value"):

    nb_signals = len(signals)
    max_size = np.max([len(signal) for signal in signals])
    fig = plt.figure(figsize=(max_size // 5, 1.5 * nb_signals))
    outer_gs = GridSpec(nb_signals, 1, hspace=1)
    colors = cmr.get_sub_cmap('viridis', 0.25, 1, N=nb_signals).colors
    global_min = np.nanmin([np.nanmin(signal) for signal in signals])
    global_max = np.nanmax([np.nanmax(signal) for signal in signals])
    
    for i in range(nb_signals):
        sig_size = len(signals[i])
        sig_coordinates = signal_limits[i]
        min_val = global_min # np.nanmin(signals[i])
        
        figure_proportion = int((sig_size/max_size) * 100)
        inner_gs = GridSpecFromSubplotSpec(
            1, 2,
            subplot_spec=outer_gs[i],
            width_ratios=[figure_proportion, 100 - figure_proportion],
            wspace=0.25
        )

        ax = fig.add_subplot(inner_gs[0, 0])
        indexes = np.arange(sig_size)
        ax.fill_between(indexes, signals[i], min_val, color = colors[i])
        
        # indexing xaxis
        indexes_coordinates = np.round(np.linspace(sig_coordinates[0], sig_coordinates[1], sig_size), 2)
        step = sig_size//5 if sig_size//5 > 0 else 1
        labels_indexes = np.arange(0, sig_size, step)
        ax.set_xticks(indexes[labels_indexes], indexes_coordinates[labels_indexes])
        
        # adding centromere display
        if len(centromeres) > 0:
            draw_centromere(ax, centromeres[i], signals[i], indexes_coordinates)
        
        # setting axes labels
        ax.set_xlabel(signal_names[i] + " " + xlabel)
        ax.set_ylabel(ylabel)
        ax.set_ylim(global_min, global_max)
        
        ax.spines[['right', 'top']].set_visible(False)
        
        if len(title) > 0 and len(signals[i]) == max_size:
            ax.get_figure().suptitle(title)

    if len(save_to) > 0:
        plt.savefig(save_to + '.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(save_to + '.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()

# # test   
# signals = [
#     np.random.rand((i+1) * 10)
#     for i in range(5)
# ]
# rand_starts = [np.random.randint(10) for i in range(5)]
# coordinates = [
#     (rand_starts[i], rand_starts[i] + np.random.randint(1, 10))
#     for i in range(5)
# ]
# chromnames = [
#     "chrom" + str(np.random.randint(10))
#     for i in range(5)
# ]
# centromeres = [
#     np.random.randint(coordinates[i][0], coordinates[i][1])
#     for i in range(5)
# ]

# plot_signals(signals, chromnames, coordinates, centromeres = centromeres, title = "Contig contact with chromosomes")

def compute_contact(c, chrom, plasmid, epsilon = 1e-5):
    d=c.info
    total_reads = d['sum']

    mat = c.matrix().fetch(chrom, plasmid)
    mat[np.isnan(mat)] = 0
    coverage_plasmid = mat.sum(axis=1)
    coverage_plasmid[coverage_plasmid==0] = np.nan

    matscn = c.matrix().fetch(chrom, chrom)
    matscn[np.isnan(matscn)] = 0
    coverage = matscn.sum(axis=0)
    coverage[coverage==0] = np.nan

    m1 = c.matrix(balance=False).fetch(chrom, chrom)
    reads_chrom= np.nansum(m1) + epsilon
    m2 = c.matrix(balance=True).fetch(plasmid, plasmid)
    reads_plasmid= np.nansum(m2) + epsilon
    m12 = c.matrix(balance=False).fetch(chrom, plasmid)
    reads_chr12= np.nansum(m12) + epsilon
    
    coverage_plasmid = (coverage_plasmid/reads_chrom)/(reads_plasmid/reads_chr12)

    return coverage_plasmid

def plot_signals(signals, signal_names, signal_limits, reference_signal = "", centromeres = [], title = "", save_to = "", xlabel = "Genomic coordinates", ylabel = "Contact value", mitos = []):

    nb_signals = len(signals)
    max_size = np.max([len(signal) for signal in signals])
    fig = plt.figure(figsize=(max_size // 5, 1.5 * nb_signals))
    outer_gs = GridSpec(nb_signals, 1, hspace=1)
    colors = cmr.get_sub_cmap('viridis', 0.25, 1, N=nb_signals).colors
    global_min = np.nanmin([np.nanmin(signals[i]) for i in range(nb_signals) if signal_names[i] != reference_signal])
    global_max = np.nanmax([np.nanmax(signals[i]) for i in range(nb_signals) if signal_names[i] != reference_signal])
    
    for i in range(nb_signals):
        sig_size = len(signals[i])
        sig_coordinates = signal_limits[i]
        min_val = global_min # np.nanmin(signals[i])
        
        figure_proportion = int((sig_size/max_size) * 100)
        inner_gs = GridSpecFromSubplotSpec(
            1, 2,
            subplot_spec=outer_gs[i],
            width_ratios=[figure_proportion, 100 - figure_proportion],
            wspace=0.25
        )

        ax = fig.add_subplot(inner_gs[0, 0])
        indexes = np.arange(sig_size)
        ax.fill_between(indexes, signals[i], min_val, color = colors[i])
        
        # indexing xaxis
        indexes_coordinates = np.round(np.linspace(sig_coordinates[0], sig_coordinates[1], sig_size), 2)
        step = sig_size//5 if sig_size//5 > 0 else 1
        labels_indexes = np.arange(0, sig_size, step)
        ax.set_xticks(indexes[labels_indexes], indexes_coordinates[labels_indexes])
        ax.set_ylim(global_min, global_max)
        
        # adding centromere display
        if len(centromeres) > 0:
            draw_centromere(ax, centromeres[i], signals[i], indexes_coordinates)
        
        # setting axes labels
        name = signal_names[i] if signal_names[i] not in mitos else "Mitochondria"
        ax.set_xlabel(name + " " + xlabel)
        ax.set_ylabel(ylabel)
        
        ax.spines[['right', 'top']].set_visible(False)
        
        if len(title) > 0 and len(signals[i]) == max_size:
            ax.get_figure().suptitle(title)

    if len(save_to) > 0:
        plt.savefig(save_to + '.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(save_to + '.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()

def plot_contig_4C(cool_path, binning, chromosomes_list, contigs_selection, centromeres = [], outpath = "", mitochondrias = []):
    cool_file = f"{cool_path}::resolutions/{binning}"
    file_name = cool_path.split('/')[-1].split('.')[0]

    c = cooler.Cooler(cool_file)
    for select_name in contigs_selection:
        chromosomes = np.append(select_name, chromosomes_list)
        signals = [
            compute_contact(c, chromosome, select_name)
            for chromosome in chromosomes
        ]
        indexes_to_keep = np.array([not np.isnan(np.nanmax(signal)) and not np.isinf(np.nanmax(signal)) and np.nansum(signal) > 0 and len(signal) > 0 for signal in signals])
        signals_to_keep = [signals[i] for i in range(len(signals)) if indexes_to_keep[i]]
        chroms_to_keep = chromosomes[indexes_to_keep]
        chroms_coordinates = [
            [0, c.chromsizes[chrom]]
            for chrom in chroms_to_keep
        ]
        title = f"{select_name} contacts with chromosomes ({binning}kb)"
        
        if len(signals_to_keep) > 0:
            save_to = f'{outpath}/{file_name}_{select_name}' if len(outpath) > 0 else ""
            plot_signals(signals_to_keep, chroms_to_keep, chroms_coordinates, reference_signal = select_name, centromeres = centromeres, save_to = save_to, title = title, mitos = mitochondrias)

def make_mini_matrix(cool_file, balance = True):
    chroms = cool_file.chromnames
    matrix = cool_file.matrix(balance= balance)
    mini_mat = np.zeros((len(chroms), len(chroms)))
    
    for i in range(len(chroms)):
        for j in range(i, len(chroms)):
            if i == j:
                mini_mat[i, i] = 0
                continue
            
            value = np.nansum(matrix.fetch(chroms[i], chroms[j]))   
            value = value if not np.isinf(value) else np.nan
            mini_mat[i,j] = value
            mini_mat[j,i] = value

    return mini_mat

def plot_mini_matrix(matrix, cool_file, outpath, name, nb_chroms, min_coverage = 1):
    plt.figure(figsize=(len(matrix), nb_chroms))
    
    selected_indexes = [i for i in range(len(matrix)) if np.nansum(matrix[i]) > min_coverage]
    
    mat = plt.matshow(np.log10(matrix[selected_indexes, :nb_chroms]))
    labels = np.array(cool_file.chromnames)
    plt.title(f"{name} mini matrix")
    indexes = [i for i in range(len(selected_indexes))]
    plt.yticks(indexes, labels[selected_indexes],  fontsize=(11 if len(selected_indexes) < 100 else 5))
    chroms_indexes= [i for i in range(nb_chroms) if i in selected_indexes]
    indexes = [i for i in range(len(chroms_indexes))]
    plt.xticks(indexes, labels[chroms_indexes], rotation = 90, fontsize=(11 if len(selected_indexes) < 100 else 4))
    cbar = plt.colorbar(mat, fraction=0.05, label="Log10 of contact signal")

    plt.savefig(f'{outpath}/{name}_mito_small_matrix.pdf', bbox_inches='tight')
    plt.savefig(f'{outpath}/{name}_mito_small_matrix.png', bbox_inches='tight')
    #plt.show()

# dataset
folder_path = "/media/sardine/data_2/Eukaroytic_plasmids_manual_research/mcools/"
mcools = [file.replace(".mcool", "") for file in os.listdir(folder_path) if ".mcool" in file]
mitochondrias = [
    "ENA|OX637640|OX637640.1",
    "ENA|OY756929|OY756929.1",
    "ENA|OX637793|OX637793.1",
    "ENA|OX344747|OX344747.1", 
    "ENA|OX344748|OX344748.1",
    "ENA|OY730176|OY730176.1",
    "ENA|OY745785|OY745785.1",
    "ENA|OY752154|OY752154.1",
    "ENA|OX596347|OX596347.1",
    "ENA|OZ124174|OZ124174.1",
    "ENA|OZ066578|OZ066578.1",
    "ENA|OZ076499|OZ076499.1",
    "ENA|OZ203567|OZ203567.1",
    "ENA|OW971919|OW971919.1",
    "ENA|OW971919|OW971919.2",
    "ENA|OZ026805|OZ026805.1",
    "ENA|OX596135|OX596135.1",
    "ENA|OZ035952|OZ035952.1",
    "ENA|OY756152|OY756152.1",
    "ENA|OZ204873|OZ204873.1",
    "ENA|OY725478|OY725478.1"
]
binning = 50000

for library in mcools:
    cool_file = f"/media/sardine/data_2/Eukaroytic_plasmids_manual_research/mcools/{library}.mcool"
    plot_outpath = f"/media/sardine/data_2/Eukaroytic_plasmids_manual_research/contact_plots/{library}"
    
    if os.path.exists(f'{plot_outpath}/{library}_mito_small_matrix.pdf'):
        continue
    cool_file = f"/media/sardine/data_2/Eukaroytic_plasmids_manual_research/mcools/{library}.mcool"
    plot_outpath = f"/media/sardine/data_2/Eukaroytic_plasmids_manual_research/contact_plots/{library}"
    if not os.path.exists(plot_outpath):
        os.mkdir(plot_outpath)
    try:
        cool = cooler.Cooler(f"{cool_file}::resolutions/{binning}")
    except:
        continue
    chromosomes = np.array([ chrom for chrom in list(cool.chromsizes.keys()) if chrom[:3] == "ENA" or (chrom[:3] != "ERR" and chrom[0] != "k" and chrom[:3] != "tig")])
    print(f"Plotting {library} ({len(cool.chromnames) - len(chromosomes)} contigs)")
    mini_matrix = make_mini_matrix(cool, balance = False)
    plot_mini_matrix(mini_matrix, cool, plot_outpath, library, len(chromosomes))

for library in mcools:
    cool_file = f"/media/sardine/data_2/Eukaroytic_plasmids_manual_research/mcools/{library}.mcool"
    plot_outpath = f"/media/sardine/data_2/Eukaroytic_plasmids_manual_research/contact_plots/{library}"
    if os.path.exists(f"{plot_outpath}/matrix_{binning}.png"):
        continue
    if not os.path.exists(plot_outpath):
        os.mkdir(plot_outpath)
    try:
        cool = cooler.Cooler(f"{cool_file}::resolutions/{binning}")
    except:
        continue
    plt.figure(figsize=(20,20))
    plt.imshow(np.log10(cool.matrix()[:]))
    plt.savefig(f"{plot_outpath}/matrix_{binning}.png")
    plt.savefig(f"{plot_outpath}/matrix_{binning}.pdf")
    chromosomes = np.array([ chrom for chrom in list(cool.chromsizes.keys()) if chrom[:3] == "ENA" or (chrom[:3] != "ERR" and chrom[0] != "k" and chrom[:3] != "tig")])
    contigs_selection = [chrom for chrom in cool.chromsizes.keys() if chrom not in chromosomes]
    plot_contig_4C(cool_file, binning, chromosomes, contigs_selection, outpath = plot_outpath, mitochondrias = mitochondrias)