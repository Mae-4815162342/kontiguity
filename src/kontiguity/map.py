from kontiguity.utils.functions import *
from kontiguity.workers import ScriptExecutor as scexec

def create_table(name, index, HICs, enzymes = "HinfI,DpnII", binnings = "5000", format = "cool"):
    """Creates the table with the elements."""
    table_els = [
        {
            "name":name,
            "index":index,
            "fastq1_hic":hic[0],
            "fastq2_hic":hic[1] if paired else "",
            "enzymes":enzymes,
            "binnings":binnings,
            "format":format
        }
        for hic, paired in HICs
    ]

    return pd.DataFrame.from_dict(table_els)

def map_hics(hic_dict, outpath, no_tmp = False, threads = 8, sbatch = False, **sbatch_params):
    """Launches Hi-C mapping."""

    # building script
    header = get_header(sbtach = sbatch, outpath = outpath, **sbatch_params)

    path_to_scripts = "/".join(__file__.split('/')[:-1])
    generate_cool_script = path_to_scripts + "/scripts/map/generate_cool.sh"
    rebin_script = path_to_scripts + "/scripts/map/rebin.sh"
    zoomify_script = path_to_scripts + "/scripts/map/zoomify.sh"

    scripts_outpath = f"{outpath}/scripts"
    build_arborescence(scripts_outpath)

    script = write_script(header, {
        "generate":(generate_cool_script, 8),
        "rebin":(rebin_script, 6),
        "zoomify":(zoomify_script, 6),
        }, scripts_outpath, name = "map_hic")
    
    # initialisations
    hic_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(hic_queue, script) for _ in range(threads)
    ]
    
    # feading fastas to the input queue
    for ref in hic_dict:
        species = ref.replace(" ","_")
        table_data = hic_dict[ref]

        # creating outfolder
        outfolder = f'{outpath}/{species}/maps'
        build_arborescence(outfolder)

        # parameters
        for k in range(len(table_data)):
            # parameters
            row_data = table_data.iloc[k]
            local_outpath = f"{outfolder}/{row_data['mapping']}"
            build_arborescence(outpath)
            hic_name = row_data['mapping']
            binnings = sorted(list(str(row_data['binnings']).split(';')))
            min_binning = str(reduce(math.gcd, [int(value) for value in binnings]))
            other_binnings = binnings[1:] if min_binning in binnings else binnings

            # index is looked for in an eventual dataset previously loaded with the load command if a path is not provided
            index = row_data["index"] if "contigs" not in row_data or len(row_data["contigs"]) == 0 else f'{outpath}/{species}/contigs/{row_data["contigs"]}/genome'
            if not os.path.isfile(index + '.fna') and not os.path.isfile(index + '.fa'):
                index = f"{outpath}/{species}/dataset/genomes/{index}"

            # generate map params
            to_map = "true"
            enzymes=row_data['enzymes']

            # fastq check
            ## if it is not a path to a file, will be looked for in the kontiguity arborescence.
            paired_fastq = f"{outpath}/{species}/dataset/hic/{row_data['fastq1']}_1.fastq.gz"
            if os.path.exists(paired_fastq):
                hic_R1  = paired_fastq
                hic_R2 = f"{outpath}/{species}/dataset/hic/{row_data['fastq1']}_2.fastq.gz"
            else:
                hic_R1= "." if len(row_data["fastq1_hic"]) == 0 else row_data["fastq1_hic"]
                hic_R2= "." if len(row_data["fastq2_hic"]) == 0 else row_data["fastq2_hic"]

            # rebin params
            to_rebin = "true" if row_data['format'] == "cool" and len(other_binnings) > 0 else "false"
            binnings_formated = " ".join(other_binnings)

            # zoomify params
            to_zoomify = "true"
            zoomings = ",".join(binnings)

            # queuing
            hic_queue.put(([
                to_map,
                to_rebin,
                to_zoomify,
                index, # generate map
                local_outpath,
                hic_name,
                hic_R1,
                hic_R2,
                enzymes,
                min_binning,
                str(threads),
                local_outpath, # rebin
                hic_name,
                min_binning,
                binnings_formated,
                str(threads),
                "true" if no_tmp else "false",
                local_outpath, # zoomify
                hic_name,
                min_binning,
                zoomings,
                str(threads),
                "true" if no_tmp else "false"
            ], sbatch))

    # closing queue and workers
    for _ in range(len(workers)):
        hic_queue.put("DONE")

    return workers

def map(
    name = "",
    outpath = "",
    index = "",
    hic = "",
    enzymes = "HinfI,DpnII",
    binnings = [5000],
    table = None,
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
        data = create_table(name, index, hic, enzymes = enzymes, binnings = binnings, format = format)
        if data is None:
            print("Error: missing HIC input.")
            return

    data["mapping"] = [f"mapping_{k + 1}" for k in range(len(data))] if "mapping" not in data.columns else data["mapping"]
    data.to_csv(f"{outfolder}/mapping_data.csv", index=False)
    
    ## retrieving unique name subset
    subset_hic = {}
    is_single_name = len(np.unique(data['name'])) == 1
    for subname in np.unique(data['name']):
        subset_hic[subname] = data[data['name'] == subname]

    ## launching retrievers
    out_tmp = outpath if is_single_name else outfolder
    retrievers = map_hics(subset_hic, outpath = out_tmp, no_tmp = no_tmp, threads = threads, sbatch = sbatch, **sbatch_params)

    ### joining loaders
    join_workers(retrievers)