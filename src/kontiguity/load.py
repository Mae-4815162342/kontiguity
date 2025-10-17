from kontiguity.utils.functions import *
from kontiguity.workers import ScriptExecutor as scexec

def create_table(name, ref, WGSs, HICs):
    """Builds a kontiguity dataset table. At least one wgs and hic must be provided."""
    table_dict = []
    for wgs in WGSs:
        for hic in HICs:
            table_dict.append({
                "name":name,
                "ref":ref,
                "wgs":wgs,
                "hic":hic
            })
    if len(table_dict) == 0:
        return None
    return pd.DataFrame.from_dict(table_dict)

def load_dtol(outpath):
    """Retrieves datas from the Darwin Tree of Life and builds a kontiguity dataset table from it."""
    pass

def load_ref(ref_dict, outpath, chroms = None, threads = 8, sequence_types = 'chromosome,organelle', sbatch = False, **sbatch_params):
    """Builds the arborescence and eventualy retrieve fastas from the GCA database before building bowtie indexes."""
    # building script
    header = get_header(sbtach = sbatch, outpath = outpath, **sbatch_params)
    path_to_scripts = "/".join(__file__.split('/')[:-1])
    load_ref_script = path_to_scripts + "/scripts/loaders/load_ref.sh"
    format_ref_script = path_to_scripts + "/scripts/loaders/format_ref.sh"
    scripts_outpath = f"{outpath}/scripts"
    build_arborescence(scripts_outpath)
    script = write_script(header, {"load":(load_ref_script, 3), "format":(format_ref_script, 4)}, scripts_outpath, name = "load_ref")

    # initialisations
    genome_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(genome_queue, script) for _ in range(threads)
    ]

    # feading fastas to the input queue
    for ref in ref_dict:
        species = ref.replace(" ","_")

        # creating outfolder
        outfolder = f'{outpath}/{species}/dataset/genomes'
        build_arborescence(outfolder)

        # adding chromosome information to the outfolder if provided
        if chroms is not None and os.path.isfile(chroms):
            shutil.copyfile(chroms, f"{outfolder}/chromosomes.csv")

        # parameters
        for k in range(len(ref_dict[ref])):
            # scripts to call
            to_load = "true" if not os.path.isfile(ref_dict[ref][k]) else 'none'
            to_format = "true"

            # parameters
            fasta_path = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{ref_dict[ref][k]}?download=true&gzip=true" if not os.path.exists(ref_dict[ref][k]) else ref_dict[ref][k]
            loaded_fasta_path = f"{outfolder}/genome.fa" if not os.path.exists(ref_dict[ref][k]) else ref_dict[ref][k]
            genome_name = species + f"_{k + 1}"

            # queuing
            genome_queue.put(([
                to_load,
                to_format,
                fasta_path,
                outfolder,
                species,
                loaded_fasta_path,
                outfolder,
                genome_name,
                sequence_types
            ], sbatch))

    # closing queue and workers
    for _ in range(len(workers)):
        genome_queue.put("DONE")

    return workers

def load_fastqs(fastq_dict, outpath, experiment, threads = 8, sbatch = False, **sbatch_params):
    """Builds the arborescence and eventualy retrieve fastq files with fasterq-dump before splitting in subfastqs."""
    # building script
    header = get_header(sbtach = sbatch, outpath = outpath, **sbatch_params)
    path_to_scripts = "/".join(__file__.split('/')[:-1])
    load_fastq_script = path_to_scripts + "/scripts/loaders/load_fastq.sh"
    scripts_outpath = f"{outpath}/scripts"
    build_arborescence(scripts_outpath)
    script = write_script(header, {"load":(load_fastq_script, 3)}, scripts_outpath, name = "load_fastq")

    # initialisations
    fastq_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(fastq_queue, script) for _ in range(threads)
    ]

    # feading fastas to the input queue
    for ref in fastq_dict:
        species = ref.replace(" ","_")

        # creating outfolder
        outfolder = f'{outpath}/{species}/dataset/{experiment}'
        build_arborescence(outfolder)

        # parameters
        for k in range(len(fastq_dict[ref])):
            fastq1, fastq2, paired = fastq_dict[ref][k]

            # scripts to call
            to_load1 = "true" if not os.path.isfile(fastq1) else 'none'
            to_load2 = "true" if paired == "paired" and not os.path.isfile(fastq2) else 'none'

            # queuing first fastq
            fastq_queue.put(([
                to_load1,
                fastq1,
                outfolder,
                str(threads)
            ], sbatch))

            # if paired, queuing second fastq
            if paired == "paired":
                fastq_queue.put(([
                    to_load2,
                    fastq2,
                    outfolder,
                    str(threads)
                ], sbatch))

    return workers

def load(
    name = "",
    outpath = None,
    ref = None,
    chroms = None,
    wgs = [],
    hic = [],
    table = None,
    dtol = False,
    threads =  8,
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

    outfolder = f"{outpath}/{name.replace(' ', '_')}"
    build_arborescence(outfolder)

    ## building data tables
    if dtol:
        data = load_dtol(outfolder)
    elif table is not None:
        data = pd.read_csv(table)
    else:
        data = create_table(name, ref, wgs, hic)
        if data is None:
            print("Error: missing WGS or HIC input.")
            return
    data.to_csv(f"{outfolder}/samples_data.csv", index=False)

    ## retrieving unique name subset
    subset_ref = {}
    subset_wgs = {}
    subset_hic = {}
    for subname in np.unique(data['name']):
        subset = data[data['name'] == subname]
        subset_ref[name] = np.unique(subset['ref'])
        subset_wgs[name] = np.unique(subset['wgs'])
        subset_hic[name] = np.unique(subset['hic'])

    ## launching loaders
    ref_loaders = []#load_ref(subset_ref, outpath = outpath, chroms = chroms, threads = threads, sbatch = sbatch, **sbatch_params)
    wgs_loaders = load_fastqs(subset_wgs, outpath = outpath, experiment='wgs', threads = threads, sbatch = sbatch, **sbatch_params)
    hic_loaders = load_fastqs(subset_hic, outpath = outpath, experiment='hic', threads = threads, sbatch = sbatch, **sbatch_params)

    ### joining loaders
    join_workers(ref_loaders)
    join_workers(wgs_loaders)
    join_workers(hic_loaders)