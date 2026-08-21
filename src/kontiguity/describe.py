from kontiguity.utils.circularity import *
from kontiguity.utils.functions import *
from kontiguity.utils.displays import *
from kontiguity.utils.imports import *

def get_trans_cov(cool, contig):
    """Returns the trans coverage of a contig in the cool matrix."""
    coverage = 0
    for chrom in cool.chromsizes.keys():
        if chrom != contig:
            coverage += np.nansum(cool.matrix(balance = False).fetch(chrom, contig))
    return coverage

# Hi-C data methods
def get_chrom_matrix(cool, balance = True):
    """Returns the numpy aray representing the contact matrix contained in the cool file."""
    return np.array(cool.matrix(balance = balance)[:])

def get_mini_matrices(cool_file, chrom_list = []):
    """Returns a tuple (raw, normalized) mini-matrices (chromosome-size binned) for the provided cool."""
    chroms = cool_file.chromnames if len(chrom_list) == 0 else chrom_list
    matrix_raw = cool_file.matrix(balance = False)
    matrix_normalized = cool_file.matrix(balance = True)
    mini_mat_raw = np.zeros((len(chroms), len(chroms)))
    mini_mat_normalized = np.zeros((len(chroms), len(chroms)))
    
    for i in range(len(chroms)):
        for j in range(i, len(chroms)):
            if i == j:
                mini_mat_raw[i, i] = 0
                mini_mat_normalized[i, i] = 0
                continue
            
            # raw
            value = np.nansum(matrix_raw.fetch(chroms[i], chroms[j]))   
            value = value if not np.isinf(value) else np.nan
            mini_mat_raw[i,j] = value
            mini_mat_raw[j,i] = value
            
            # normalized
            value = np.nansum(matrix_normalized.fetch(chroms[i], chroms[j]))   
            value = value if not np.isinf(value) else np.nan
            mini_mat_normalized[i,j] = value
            mini_mat_normalized[j,i] = value

    return (mini_mat_raw, mini_mat_normalized)

def compute_contact(c, chrom, contig, balance = True, coords = [], epsilon = 1e-5):
    """Returns computed contact between a contig and a selected chromosome in the cool file.
    If coords is provided as the start and end coordinates delimiting a region of the contig, contact is only computed for said region."""
    d=c.info
    total_reads = d['sum']
    
    contig = f"{contig}:{coords[0]}-{coords[1]}" if len(coords) > 0 else contig

    mat = c.matrix(balance = balance).fetch(chrom, contig)
    mat[np.isnan(mat)] = 0
    coverage_contig = mat.sum(axis=1).astype("float64")
    coverage_contig[coverage_contig==0] = np.nan

    matscn = c.matrix(balance = balance).fetch(chrom, chrom)
    matscn[np.isnan(matscn)] = 0
    coverage = matscn.sum(axis=0).astype("float64")
    coverage[coverage==0] = np.nan

    m1 = c.matrix(balance=False).fetch(chrom, chrom)
    reads_chrom= np.nansum(m1) + epsilon
    m2 = c.matrix(balance=False).fetch(contig, contig)
    reads_contig= np.nansum(m2) + epsilon
    m12 = c.matrix(balance=False).fetch(chrom, contig)
    reads_chr12= np.nansum(m12) + epsilon
    
    coverage_contig = (coverage_contig/reads_chrom)/(reads_contig/reads_chr12)

    return coverage_contig

def get_computed_contacts(cool, chromosomes, contig, coords = [], balance = True):
    """Returns the computed contacts between a contig and the chromosome selection."""
    return { chromosome: compute_contact(cool, chromosome, contig, coords = coords, balance = balance) for chromosome in chromosomes }

def get_contig_signal(cool, chromosomes, contig, coords = [], balance = True):
    """Returns continuous signal of the selected contig with the entire genome."""
    return np.concatenate([compute_contact(cool, chromosome, contig, coords = coords, balance = balance) for chromosome in chromosomes])

def get_contig_contact(cool, chromosomes, contig, coords = [], balance = True):
    """Returns continuous average contact of the selected contig with the entire genome."""
    contig = f"{contig}:{coords[0]}-{coords[1]}" if len(coords) > 0 else contig
    matrix = cool.matrix(balance = balance)
    return np.concatenate([np.nanmedian(matrix.fetch(chromosome, contig), axis = 1) for chromosome in chromosomes])

# def get_ps(cool, chrom, coords = [], ignore_diags=1):
#     """Computes P(s) for selected region""" 
#     if len(coords) == 2:
#         view_df = bioframe.from_any([[chrom, coords[0], coords[1], chrom]])
#     else:
#         view_df = bioframe.from_any([[chrom, 0, cool.chromsizes[chrom], chrom]])
#     view_df.columns = ["chrom", "start", "end", "name"]
 
#     expected = cooltools.expected_cis(
#         cool, view_df=view_df, ignore_diags=ignore_diags, smooth=True
#     )
#     expected["dist_bp"] = expected["dist"] * cool.binsize
#     return expected

def get_random_coords(chromsize, size, nb_rand = 2):
    """Returns a selection of [start, stop] random regions of size."""
    selection = []
    for _ in range(nb_rand):
        random_chromosome_start = np.random.randint(1, chromsize - size)
        selection.append([random_chromosome_start, random_chromosome_start + size])
    return selection

def get_contact_signals(cool, chromosomes, selected_chromosome, contig, mitochondria = None, size = 100000, nb_rand = 5, average_method = np.nanmean):
    """Computes the contact signals (Hi-C, computed, raw, normalized) of a random sequence in each provided contig.
    If the mitochondria is provided, it will be included. 
    The used size will either be the minimum between the provided size and the smallest of the two (or three) provided sequences."""
    
    used_size = min([size, cool.chromsizes[selected_chromosome], cool.chromsizes[contig]] + ([cool.chromsizes[mitochondria]] if mitochondria in cool.chromsizes else []))
    chromosome_coords = get_random_coords(cool.chromsizes[selected_chromosome], used_size, nb_rand = nb_rand)
    mitochondria_coords = get_random_coords(cool.chromsizes[mitochondria], used_size, nb_rand = nb_rand) if mitochondria in cool.chromsizes else []

    signals = {}
    
    # raw signals
    signals["Raw"] = {
        "Hi-C":{
            "Chromosome": average_method([get_contig_contact(cool, chromosomes, selected_chromosome, coords = coords, balance = False) for coords in chromosome_coords], axis = 0),
            "Contig": get_contig_contact(cool, chromosomes, contig, balance = False)
        },
        "Computed":{
            "Chromosome": average_method([get_contig_signal(cool, chromosomes, selected_chromosome, coords = coords, balance = False) for coords in chromosome_coords], axis = 0),
            "Contig": get_contig_signal(cool, chromosomes, contig, balance = False)
        },
    }
    if mitochondria in cool.chromsizes:
        signals["Raw"]["Hi-C"]["Mitochondria"] = average_method([get_contig_contact(cool, chromosomes, mitochondria, coords = coords, balance = False) for coords in mitochondria_coords], axis = 0)
        signals["Raw"]["Computed"]["Mitochondria"] = average_method([get_contig_signal(cool, chromosomes, mitochondria, coords = coords, balance = False) for coords in mitochondria_coords], axis = 0)
        
    # normalized signals
    signals["Normalized"] = {
        "Hi-C":{
            "Chromosome": average_method([get_contig_contact(cool, chromosomes, selected_chromosome, coords = coords, balance = True) for coords in chromosome_coords], axis = 0),
            "Contig": get_contig_contact(cool, chromosomes, contig, balance = True)
        },
        "Computed":{
            "Chromosome": average_method([get_contig_signal(cool, chromosomes, selected_chromosome, coords = coords, balance = True) for coords in chromosome_coords], axis = 0),
            "Contig": get_contig_signal(cool, chromosomes, contig, balance = True)
        },
    }
    if mitochondria in cool.chromsizes:
        signals["Normalized"]["Hi-C"]["Mitochondria"] = average_method([get_contig_contact(cool, chromosomes, mitochondria, coords = coords, balance = True) for coords in mitochondria_coords], axis = 0)
        signals["Normalized"]["Computed"]["Mitochondria"] = average_method([get_contig_signal(cool, chromosomes, mitochondria, coords = coords, balance = True) for coords in mitochondria_coords], axis = 0)
    
    # # P(s)
    # signals["P(s)"] = {
    #     "Chromosome": average_method([get_ps(cool, selected_chromosome, coords = coords),
    #     "Contig": get_ps(cool, contig)
    # }
    # if mitochondria in cool.chromsizes:
    #     signals["P(s)"]["Mitochondria"] = average_method([get_ps(cool, mitochondria, coords = coords)
    
    
    return signals
    
def get_coverage_hic(cool, contig):
    """Returns total coverage of given contig in the cool file."""
    coverage = cool.matrix(balance = False).fetch(contig).sum()
    for chrom in cool.chromsizes.keys():
        coverage += cool.matrix(balance = False).fetch(chrom, contig).sum()
    return coverage

def get_copies(cool, reference_chromosome, contig, ploidy = 1):
    """Return the estimated number of contig's copies in the cool file, normalized by the reference_chromosome."""
    chrom_cov = get_coverage_hic(cool, reference_chromosome)
    contig_cov = get_coverage_hic(cool, contig)
    
    return (contig_cov/cool.chromsizes[contig]) / (ploidy * (chrom_cov/cool.chromsizes[reference_chromosome]))

def get_GC(sequence):
    """Computes GC content of a sequence"""
    return round(np.sum([1 if s in ["G", "C"] else 0 for s in sequence])/len(sequence), 2)

def get_GCs(fasta, chromosomes, contigs):
    """Returns the individual GC content of each contig present in the fasta file, as well as the average genome GC content."""
    sequences = parse_fasta(fasta)
    average_genome = []
    GC_contigs = {}
    for seq_id in sequences:
        if seq_id not in chromosomes and seq_id not in contigs:
            continue
        sequence = sequences[seq_id]
        GC = get_GC(sequence)
        if seq_id in contigs:
            GC_contigs[seq_id] = GC
        else:
            average_genome.append(GC)
    return GC_contigs, np.nanmedian(average_genome)

def describe(
    name = "",
    outpath = "",
    chroms = "",
    mitochondria = "",
    chromstart = "NC_",
    min_chrom_size = 100000,
    contigs = "",
    fasta = None,
    gff = None,
    tracks = None,
    mcool = None,
    binning = "10000",
    cool = None,
    table = None,
    formats = "pdf",
    sbatch = False,
    sbtach_partition = 'dedicated',
    sbtach_qos = 'fast',
    sbtach_mem = '40G',
    sbatch_ncpus = 30
):
    sbatch_params = {
        '--partition': sbtach_partition,
        '--qos': sbtach_qos,
        '--mem': sbtach_mem,
        '-c': sbatch_ncpus
    }

    outfolder = f"{outpath}/{name.replace(' ', '_')}/describe"
    build_arborescence(outfolder)

    cool = cooler.Cooler(cool if cool is not None else f"{mcool}::resolutions/{binning}")
    contigs = contigs.split(",")
    if len(chroms) == 0:
        chromosomes = [ chrom for chrom in list(cool.chromsizes.keys()) if chrom[:3] == chromstart and chrom not in contigs and cool.chromsizes[chrom] >= min_chrom_size] 
    else: 
        chromosomes = [ chrom for chrom in chroms.split(",") if cool.chromsizes[chrom] >= min_chrom_size ]

    ## selecting contigs interacting with the genome (more likely to be intra-nuclear)
    contig_selection = []
    for contig in cool.chromsizes.keys():
        if contig in chromosomes or contig not in contigs:
            continue
        trans_coverage = get_trans_cov(cool, contig)
        if trans_coverage > 0:
            contig_selection.append(contig)
    # sequence_selection = chromosomes + contig_selection

    ## computing circularity of contigs if the fasta file is provided
    circulars = None
    GC_contigs = {}
    GC_global = np.nan
    if os.path.exists(fasta):
        circulars = get_circulars(fasta, selection = contig_selection, min_overlap = 20, max_overlap = 100, max_mismatch_rate =  0.2)
        GC_contigs, GC_global = get_GCs(fasta, chromosomes, contig_selection)

    for contig in contig_selection:
        global_data = {
            "Contig":contig,
            "Species": name,
            "Chromosome":chromosomes[0],
            "Contigs": contig_selection,
            "Chromosomes": chromosomes, 
            "Mitochondria": mitochondria,
            "Length": cool.chromsizes[contig],
            "GC": f"{GC_contigs[contig]} (average: {GC_global})" if contig in GC_contigs else "",
            "global_GC": GC_global, 
            "Binning": cool.binsize,
            "Circularity": "Not computed" if circulars is None else circulars[contig] if contig in circulars else "No overlap detected"
        }
        
        hic_data = global_data | {
            "Chrom_matrix": get_chrom_matrix(cool, binning),
            "Mini_matrices": get_mini_matrices(cool, chrom_list = chromosomes + [contig]),
            "Computed_contacts": get_computed_contacts(cool, chromosomes, contig),
            "Signals": get_contact_signals(cool, chromosomes, chromosomes[0], contig, mitochondria=mitochondria),
            "Coverage": get_coverage_hic(cool, contig),
            "Estimated_copies": get_copies(cool, chromosomes[0], contig)
        }
        tracks_data = {}
        sequence_data = {}

        ## global display
        build_display(contig, hic_data, tracks_data, sequence_data, outpath = outfolder, formats = formats.split(","))