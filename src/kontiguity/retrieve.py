from kontiguity.utils.functions import *
from kontiguity.workers import ScriptExecutor as scexec

def create_table(name, index, WGSs, min_size = 1000):
    """Creates the table with the elements."""
    table_els = [
        {
            "name":name,
            "index":index,
            "fastq1_wgs":wgs[0],
            "fastq2_wgs":wgs[1] if paired else "",
            "is_paired_wgs":paired,
            "min_size":min_size
        }
        for wgs, paired in WGSs
    ]

    return pd.DataFrame.from_dict(table_els)

def retrieve_contigs(retrieval_dict, outpath, logan = False, no_tmp = False, threads = 8, sbatch = False, **sbatch_params):
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
        "logan":(logan_script, 3),
        "build":(build_contigs_script, 10),
        "filter":(filter_script, 6),
        }, scripts_outpath, name = "retrieve_contigs")

    # initialisations
    contigs_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(contigs_queue, script) for _ in range(threads)
    ]
    
    # feading fastas to the input queue
    for ref in retrieval_dict:
        species = ref.replace(" ","_")
        table_data = retrieval_dict[ref]

        # creating outfolder
        outfolder = f'{outpath}/{species}/contigs'
        build_arborescence(outfolder)

        # parameters
        for k in range(len(retrieval_dict[ref])):
            # parameters
            row_data = table_data.iloc[k]
            local_outpath = f"{outfolder}/{row_data['contigs']}"
            build_arborescence(outpath)
            process_stmp = f"{row_data['contigs'].split('_')[-1]}"
            min_len = f"{row_data['min_size']}"

            # index is looked for in an eventual dataset previously loaded with the load command if a path is not provided
            index = row_data["index"]
            if not os.path.isfile(index + '.fna') and not os.path.isfile(index + '.fa'):
                index = f"{outpath}/{species}/dataset/genomes/{index}"

            # logan params
            to_logan = "true" if logan else "false"
            accession = row_data["fastq1_wgs"].split("/")[-1].split('.')[0].split('_')[0] # retrieving the SRA id in the fastq name if provided

            # building params
            to_build = "true"

            # fastq check
            ## if it is not a path to a file, will be looked for in the kontiguity arborescence. If still not found, will return the ref which will be considered as a SRA id for Logan.
            paired_fastq = f"{outpath}/{species}/dataset/wgs/{row_data['fastq1_wgs']}_1.fastq"
            single_fastq = f"{outpath}/{species}/dataset/wgs/{row_data['fastq1_wgs']}.fastq"

            if os.path.exists(paired_fastq):
                is_paired = "true"
                fastq_R1  = paired_fastq
                fastq_R2 = f"{outpath}/{species}/dataset/wgs/{row_data['fastq1_wgs']}_2.fastq"
                fastq = "."
            elif os.path.exists(single_fastq):
                is_paired = "false"
                fastq_R1  = "."
                fastq_R2 = "."
                fastq = single_fastq
            else:
                fastq_R1= row_data["fastq1_wgs"] if row_data["is_paired_wgs"] else "."
                fastq_R2= row_data["fastq2_wgs"] if row_data["is_paired_wgs"] else "."
                fastq = row_data["fastq1_wgs"] if not row_data["is_paired_wgs"] else "."
                is_paired = "true" if row_data["is_paired_wgs"] else "false"

            # filtering params
            to_filter = "true"

            # queuing
            contigs_queue.put(([
                to_logan,
                to_build,
                to_filter,
                local_outpath, # logan parameters
                accession,
                process_stmp,
                index, # build params
                local_outpath,
                fastq_R1,
                fastq_R2,
                fastq,
                is_paired,
                min_len,
                f"{threads}",
                process_stmp,
                f'{"true" if no_tmp else "false"}',
                index, # filter params
                local_outpath,
                min_len,
                f"{threads}",
                process_stmp,
                f'{"true" if no_tmp else "false"}',
            ], sbatch))

    # closing queue and workers
    for _ in range(len(workers)):
        contigs_queue.put("DONE")

    return workers

def retrieve(
    name = "",
    outpath = "",
    index = "",
    min_size = 1000,
    wgs = "",
    table = None,
    logan = False,
    no_tmp = False,
    threads =  8,
    sbatch = False,
    sbtach_partition = 'dedicated',
    sbtach_qos = 'fast',
    sbtach_mem = '40G',
    sbatch_ncpus = 30,
    **kwargs
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
    data["contigs"] = [f"contigs_{k + 1}" for k in range(len(data))] if "contigs" not in data.columns else data["contigs"]
    data.to_csv(f"{outfolder}/contigs_data.csv", index=False)
    
    ## retrieving unique name subset
    subset_retrieval = {}
    is_single_name = len(np.unique(data['name'])) == 1
    for subname in np.unique(data['name']):
        subset_retrieval[subname] = data[data['name'] == subname]

    ## launching retrievers
    out_tmp = outpath if is_single_name else outfolder
    retrievers = retrieve_contigs(subset_retrieval, outpath = out_tmp, logan = logan, no_tmp = no_tmp, threads = threads, sbatch = sbatch, **sbatch_params)

    ### joining loaders
    join_workers(retrievers)