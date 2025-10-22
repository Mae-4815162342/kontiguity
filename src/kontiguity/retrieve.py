from kontiguity.utils.functions import *
from kontiguity.workers import ScriptExecutor as scexec

def create_table(name, index, wgs, min_size = 1000):
    """Creates the table with the elements."""
    pass

def retrieve_contigs(retrieval_dict, outpath, experiment, logan = False, threads = 8, sbatch = False, **sbatch_params):
    """Launches contigs retrieval."""

    # building script
    header = get_header(sbtach = sbatch, outpath = outpath, **sbatch_params)
    path_to_scripts = "/".join(__file__.split('/')[:-1])
    logan_script = path_to_scripts + "/scripts/retrieve/logan.sh"
    build_contigs_script = path_to_scripts + "/scripts/retrieve/build_contigs.sh"
    filter_script = path_to_scripts + "/scripts/retrieve/filter.sh"
    scripts_outpath = f"{outpath}/scripts"
    build_arborescence(scripts_outpath)
    script = write_script(header, {
        "logan":(logan_script, 2),
        "build":(build_contigs_script, 3),
        "filter":(filter_script, 3),
        }, scripts_outpath, name = "retrieve_contigs")

    # initialisations
    contigs_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(contigs_queue, script) for _ in range(threads)
    ]
    
    # feading fastas to the input queue
    for ref in retrieval_dict[]:
        species = ref.replace(" ","_")

        # creating outfolder
        outfolder = f'{outpath}/{species}/contigs'
        build_arborescence(outfolder)

        # parameters
        for k in retrieval_dict:
            # scripts to call
            to_build = "true"
            to_filter = "true"

            # parameters
            fasta_path = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{ref_dict[ref][k]}?download=true&gzip=true" if not os.path.exists(ref_dict[ref][k]) else ref_dict[ref][k]
            genome_name = species + f"_{k + 1}"
            loaded_fasta_path = f"{outfolder}/{genome_name}.all_seqs.fa" if not os.path.exists(ref_dict[ref][k]) else ref_dict[ref][k]

            # queuing
            genome_queue.put(([
                to_load,
                to_format,
                fasta_path,
                outfolder,
                genome_name,
                loaded_fasta_path,
                outfolder,
                genome_name,
                sequence_types
            ], sbatch))

    # closing queue and workers
    for _ in range(len(workers)):
        genome_queue.put("DONE")

    return workers

def retrieve(
    name = "",
    outpath = "",
    index = "",
    min_size = 1000,
    wgs = "",
    table = None,
    logan = False,
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
    if table is not None:
        data = pd.read_csv(table)
    else:
        data = create_table(name, index, wgs, min_size = min_size)
        if data is None:
            print("Error: missing WGS or HIC input.")
            return
    data.to_csv(f"{outfolder}/contigs_data.csv", index=False)
    
    ## retrieving unique name subset
    subset_retrieval = {}
    is_single_name = True
    for subname in np.unique(data['name']):
        if is_single_name and name != subname:
            is_single_name = False
        subset_retrieval[subname] = data[data['name'] == subname]

    ## launching retrievers
    out_tmp = outpath if is_single_name else outfolder
    retrievers = retrieve_contigs(subset_retrieval, outpath = out_tmp, logan = logan, threads = threads, sbatch = sbatch, **sbatch_params)

    ### joining loaders
    join_workers(retrievers)