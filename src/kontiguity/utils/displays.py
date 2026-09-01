# display method for the describe function
from .imports import *

# local parameters
SIGNAL_TITLES = {
        "Raw":"Unbalanced",
        "Hi-C":"Average contact\nsignal (log10)",
        "Computed":"Computed\nsignal (log10)",
        "Normalized":"Normalized"
    }

def compute_log(array):
    """Computes log10 of array without errors. log10(0) and log10(NaN) are replaced by NaN."""
    log_array = np.copy(array)
    log_array[np.isnan(log_array)] = 0
    log_array[log_array > 0] = np.log10(log_array[log_array > 0])
    log_array[log_array == 0] = np.nan
    return log_array

def display_mini_matrix(matrix, ax, labels, nb_chroms,  title, left_labels = True, fontsize = 5):
    """Plots mini matrix on ax."""
    mat = ax.imshow(compute_log(matrix[:, :nb_chroms]))
    ax.set_title(f"{title}\nmini matrix")
    if left_labels:
        indexes = [i for i in range(len(labels))]
        ax.set_yticks(indexes, labels,  fontsize = fontsize)
    else:
        ax.set_yticks([], [])
    indexes = [i for i in range(nb_chroms)]
    ax.set_xticks(indexes, labels[:nb_chroms], rotation = 90, fontsize = fontsize)
    return mat
    
def display_mini_matrices(outer_gs, raw, norm, labels, nb_chrom, fontsize = 7):
    """Plots both mini-matrices with legends and colorbars."""
    gs = grid.GridSpecFromSubplotSpec(3, 4, subplot_spec=outer_gs, wspace=0.4, hspace=0.5, width_ratios = [1, 0.05, 1, 0.05], height_ratios = [1,1,1])

    ax_raw = plt.subplot(gs[0:3, 0])
    ax_cbar_raw = plt.subplot(gs[1, 1])
    ax_norm = plt.subplot(gs[0:3, 2])
    ax_cbar_norm = plt.subplot(gs[1, 3])
    mat = display_mini_matrix(raw, ax_raw, labels, nb_chrom, "Unbalanced")
    cbar = plt.colorbar(mat, cax = ax_cbar_raw, shrink = 0.1)
    ax_cbar_raw.set_title("Contact\nintensity\n(log10)", fontsize = fontsize)
    mat = display_mini_matrix(norm, ax_norm, labels, nb_chrom, "Normalized", left_labels = False)
    cbar = plt.colorbar(mat, cax = ax_cbar_norm, shrink = 0.1)
    ax_cbar_norm.set_title("Contact\nintensity\n(log10)", fontsize = fontsize)

def plot_mini_matrices(mini_raw, mini_norm, chromosomes, contigs, mapping, mitochondria = "", outpath = "", max_size = 20, formats = ["pdf"]):
    """Plots unbalanced and normalized mini matrices for all contigs. If the number of contigs to plot is above the max size, will divide in several figures, keeping for each the chromosomes and the mitochondria (if provided)."""
    for i in range(0, len(contigs) + max_size - 1, max_size):
        current_contigs_number = max_size if i + max_size < len(contigs) else len(contigs) - i
        selected_contigs = contigs[i: i + current_contigs_number]
        selected_sequences = chromosomes + [mitochondria + "(organelle)"] + selected_contigs if len(mitochondria) > 0 else chromosomes + selected_contigs
        nb_chrom = len(chromosomes) + (1 if len(mitochondria) > 0 else 0)
        
        selected_mini_raw = np.concatenate([mini_raw[:nb_chrom, :nb_chrom], mini_raw[nb_chrom + i: nb_chrom + i + current_contigs_number, :nb_chrom]], axis = 0)
        selected_mini_norm = np.concatenate([mini_norm[:nb_chrom, :nb_chrom], mini_norm[nb_chrom + i: nb_chrom + i + current_contigs_number, :nb_chrom]], axis = 0)

        plt.figure(figsize = (10, 10))
        outer_gs = grid.GridSpec(1, 1)
        display_mini_matrices(outer_gs[0], selected_mini_raw, selected_mini_norm, selected_sequences, nb_chrom, fontsize = 10)    

        figure_number = "" if len(contigs) <= max_size else f"{(i + current_contigs_number) % max_size}_out_of_{len(contigs % max_size) + 1}."
        for form in formats:
            plt.savefig(f"{outpath}/Mini_matrices{'.' + mapping if len(mapping) > 0 else ''}.{figure_number}{form}", dpi = 300, bbox_inches = "tight")
        if i + current_contigs_number >= len(contigs):
            break
    
def draw_info_box(ax, contig, species, binning, data, fontsize=11, title_fontsize=12):
    """Draw information rectangle on ax."""
    ax.axis("off")
    
    lines = [f"{key} : {value}" for key, value in data.items()]
    body_text = "\n".join(lines)
    box_props = dict(
        boxstyle="round,pad=0.6",
        facecolor="#f5f5f5",
        edgecolor="#333333",
        linewidth=1.2,
    )
    ax.set_title(f"Contig {contig} in {species.replace('_', ' ')} (binning: {int(binning // 1000)}kb)", fontsize=title_fontsize, fontweight="bold")
    ax.text(
        0.5, 0.75,
        body_text,
        transform=ax.transAxes,
        ha="center", va="top",
        fontsize=fontsize,
        linespacing=1.8,
        bbox=box_props,
    )

def display_contacts(outer_gs, contacts, size = 1):
    gs = grid.GridSpecFromSubplotSpec(len(contacts), 2, subplot_spec=outer_gs, wspace=0, hspace=0, width_ratios = [1, 0.025])
    contact_lists = np.concatenate([contact for contact in contacts.values()])
    vmin = np.nanmin(contact_lists)
    vmax = np.nanmax(contact_lists)
    
    shift = 0
    max_size = -np.inf
    axes = []
    for chromosome in contacts.keys():
        ax = plt.subplot(gs[shift, 0])
        ax.axis("off")
        
        contact = contacts[chromosome]
        if len(contact) > max_size:
            max_size = len(contact)
        to_plot = [contact for i in range(int(7 * size))]
        
        patch = patches.FancyBboxPatch((0, 0), len(contact) - 1, len(to_plot) - 1, boxstyle=f'round,rounding_size={3*size}', fill=None)
        ax.add_patch(patch)
        
        image = ax.imshow(to_plot, cmap='YlOrBr', clip_path=patch, clip_on=True, vmin = vmin, vmax = vmax)
        ax.text(len(contact) + 1, len(to_plot)/2, chromosome, ha = "left", va = "center", fontsize = size * 6)
        
        if shift == 0:
            ax.set_title("Contact profile")
        shift += 1
        
        axes.append(ax)
    for ax in axes:
        ax.set_xlim(0 - 0.05 * max_size, max_size * 1.05)
        ax.set_ylim((7*size, -size))
    
    ax_colorbar = plt.subplot(gs[0:2, 1])
    plt.colorbar(image, cax = ax_colorbar)
    ax_colorbar.set_title("Contact\nlevel", fontsize = 7)

def display_signal(ax, signals, binning):
    """Displays signals as dispensed in the signals dictionnary."""
    SIGNAL_COLORS = {
        "Contig":"red",
        "Chromosome":"blue",
        "Mitochondria":"green"
    }
    for entity, signal in signals.items():
        ax.plot(compute_log(signal), label = entity, color=SIGNAL_COLORS[entity])
        index = range( 0, len(signal), int(len(signal) // 5))
        ax.set_xticks(index, [i * binning for i in index])
    
def display_hic(outer_gs, contig, data):
    """Displays the data provided in data from Hi-C experiment."""
    gs = grid.GridSpecFromSubplotSpec(8, 4, subplot_spec=outer_gs,
                                             wspace=0.3, hspace=0.7)
    
    # Hi-C matrix
    matrix = data["Chrom_matrix"]
    mat_size = len(matrix) * data["Binning"]
    ax_matrix = plt.subplot(gs[:3, :2])
    ax_matrix.set_title("Hi-C contact matrix\n(chromosomes and contigs)")
    im = ax_matrix.imshow(compute_log(matrix), cmap = "afmhot_r", extent = [0, mat_size, mat_size, 0])
    ax_matrix.set_ylabel("Genomic coordinates (in bp)", fontsize = 7)
    ax_matrix.set_xlabel("Genomic coordinates (in bp)", fontsize = 7)
    plt.colorbar(im, ax = ax_matrix, shrink = 0.5)
    
    # Mini-matrices
    mini_raw, mini_norm = data["Mini_matrices"]
    selection = data["Chromosomes"] + [data["Contig"]]
    nb_chrom = len(data["Chromosomes"])
    display_mini_matrices(gs[3:6,0:2], mini_raw, mini_norm, selection, nb_chrom)
    
    # Summary statistics
    size = data["Length"]
    coverage = data["Coverage"]
    copies = data["Estimated_copies"]
    GC = data["GC"]
    circularity = data["Circularity"]
    infos = {"Length":size, "Coverage":coverage, "Copy number": round(copies, 2), "Circularity": circularity} | ({"GC content": GC} if len(GC) > 0 else {})
    ax_summary = plt.subplot(gs[0, 2:4])
    draw_info_box(ax_summary, contig, data["Species"], data["Binning"], infos, fontsize = 11)
    
    # Contact plots
    contacts = data["Computed_contacts"]
    display_contacts(gs[2:6, 2:4], contacts)
    
    # Signals
    start_signal = 6
    signals = data["Signals"]
    binning = data["Binning"]
    raw_hic = ["Raw", "Hi-C"], [start_signal, 0]
    raw_computed = ["Raw","Computed"], [start_signal + 1, 0]
    norm_hic = ["Normalized","Hi-C"], [start_signal, 2]
    norm_computed = ["Normalized","Computed"], [start_signal + 1, 2]
    
    for keys, ax_coord in [raw_hic, raw_computed, norm_hic, norm_computed]:
        ax_signal = plt.subplot(gs[ax_coord[0], ax_coord[1]:ax_coord[1]+2])
        current_signal = signals[keys[0]][keys[1]]
        display_signal(ax_signal, current_signal, binning)
        if ax_coord[0] == start_signal:
            ax_signal.set_title(SIGNAL_TITLES[keys[0]])
        else:
            ax_signal.set_xlabel("Genomic coordinates (in bp, all chromosomes)", fontsize = 7)
        if ax_coord[1] == 0:
            ax_signal.set_ylabel(SIGNAL_TITLES[keys[1]])
        elif ax_coord[0] == start_signal:
            ax_signal.legend(loc = "lower right", bbox_to_anchor = (1.0, 1.1), fontsize = 7)
    
def build_display(contig, hic_data, tracks_data, sequence_data, outpath = "", formats = ["png"]):
    """Builds display from data."""
    has_hic = len(hic_data) > 0
    has_tracks = len(tracks_data) > 0
    has_sequence = len(sequence_data) > 0
    
    nb_subfig = int(has_hic) + int(has_tracks) + int(has_sequence)
    
    plt.figure(figsize = (nb_subfig * 10, 10))
    plt.tight_layout()
    gs = grid.GridSpec(nb_subfig, 1)
    subfig_count = 0
    
    plt.rcParams['xtick.labelsize'] = 7
    plt.rcParams['ytick.labelsize'] = 7
    plt.rcParams['axes.titlesize'] = "medium"
    plt.rcParams['axes.labelsize'] = "small"
    
    if has_hic:
        outer_gs = gs[subfig_count, 0]
        subfig_count += 1
        display_hic(outer_gs, contig, hic_data)
    
    # if has_tracks:
    #     outer_gs = gs[subfig_count, 0]
    #     subfig_count += 1
        
    # if has_sequence:
    #     outer_gs = gs[subfig_count, 0]
    #     subfig_count += 1

    for fmt in formats:
        to_append = f"{'.' + hic_data['Mapping'] if len(hic_data['Mapping']) > 0 else ''}.{hic_data['Binning']}." if has_hic else "."
        plt.savefig(f"{outpath}/{contig}{to_append}summary.{fmt}", dpi = 300, bbox_inches = "tight")