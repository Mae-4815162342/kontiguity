from kontiguity.utils.functions import *
from kontiguity.utils.logging_setup import get_logger
from kontiguity.workers import ScriptExecutor as scexec
from kontiguity.workers.DToLScrapper import DToLScrapperScheduler
from kontiguity.workers.DToLFormater import DToLFormaterScheduler

def create_table(name, ref, WGSs, HICs, no_hic = False, no_wgs = False):
    """Builds a kontiguity dataset table. At least one wgs and hic must be provided."""
    table_dict = []
    if no_hic and not no_wgs:
        for wgs in WGSs:
            table_dict.append({
                "name":name,
                "ref":ref,
                "wgs":wgs,
                "hic":"."
            })
    elif no_wgs and not no_hic:
        for hic in HICs:
            table_dict.append({
                "name":name,
                "ref":ref,
                "wgs":".",
                "hic":hic
            })
    elif no_wgs and no_hic:
        table_dict.append({
            "name":name,
            "ref":ref,
            "wgs":".",
            "hic":"."
        })
    else:
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

def build_dataset(names, indexes, fastqs_wgs, fastqs_hic, **args):
    """Builds global parameter table for each """
    binnings = np.array(args["binnings"]).astype(str) if "binnings" in args else []
    rows = []
    for name in names:
        name = name.replace(" ","_")
        current_indexes = indexes[name] if name in indexes else [args["index"]] if "index" in args else None
        if current_indexes is None:
            print(f"ERROR: no reference or index provided for {name}")
        WGSs = fastqs_wgs[name]
        HICs = fastqs_hic[name]

        contigs = 1
        mapping = 1
        for index in current_indexes:
            for WGS in WGSs:
                wgs, wgs_pairing = WGS
                is_paired_load_wgs = wgs_pairing == "WAIT_LOAD" and os.path.exists(wgs[0] + "_2.fastq")
                wgs_1 = wgs[0] if wgs_pairing == "PAIRED" or wgs_pairing == "SINGLE" else wgs[0] + "_1.fastq" if is_paired_load_wgs else wgs[0] + ".fastq"
                wgs_2 = wgs[1] if wgs_pairing == "PAIRED" or wgs_pairing == "SINGLE" else wgs[0] + "_2.fastq" if is_paired_load_wgs else ""
                for HIC in HICs:
                    hic, hic_pairing = HIC
                    is_paired_load_hic = hic_pairing == "WAIT_LOAD" and os.path.exists(hic[0] + "_2.fastq")
                    hic_1 = hic[0] if hic_pairing == "PAIRED" or hic_pairing == "SINGLE" else hic[0] + "_1.fastq" if is_paired_load_hic else hic[0] + ".fastq"
                    hic_2 = hic[1] if hic_pairing == "PAIRED" or hic_pairing == "SINGLE" else hic[0] + "_2.fastq" if is_paired_load_hic else ""
                    rows.append({
                        "name":name,
                        "index":index,
                        "fastq1_wgs": wgs_1, # retrieve params
                        "fastq2_wgs": wgs_2,
                        "is_paired_wgs":wgs_pairing == "PAIRED" or is_paired_load_wgs,
                        "contigs":f"contigs_{contigs}",
                        "min_size":args["min_size"] if "min_size" in args else np.nan,
                        "fastq1_hic":hic_1, # map params
                        "fastq2_hic":hic_2 if hic_pairing else "",
                        "mapping":f"mapping_{mapping}",
                        "enzymes":args["enzymes"] if "enzymes" in args else np.nan,
                        "binnings":";".join(binnings) if "binnings" in args else np.nan,
                        "format":args["format"] if "format" in args else np.nan,
                        "chroms":args["chroms"] if "chroms" in args else np.nan, # describe params
                        "cool":args["cool"] if "cool" in args else "",
                        "mcool":f"mapping_{mapping}/mapping_{mapping}.mcool",
                        "formats":";".join(args["formats"]) if "formats" in args else np.nan
                    })
                    mapping += 1
                contigs += 1

    return pd.DataFrame.from_dict(rows)

def load_dtol(nb_per_page=100, threads=8):
    """Retrieves datas from the Darwin Tree of Life and builds a kontiguity dataset table from it."""

    # retrieving the count value
    URL_tmp = f"https://portal.darwintreeoflife.org/api/data_portal?limit=1&offset=0&sort=currentStatus:asc&current_class=kingdom"

    result = None
    try:
        result = requests.get(URL_tmp)
    except:
        print("Request failed")
    data = json.loads(result.text)
    count = data["count"]

    # initialisation
    offset = 0
    request_params_queue = Queue()
    result_queue = Queue()
    table_queue = Queue()

    # initialazing workers
    scrappers = [
        DToLScrapperScheduler(request_params_queue, result_queue)
        for _ in range(threads)
    ]
    formaters = [
        DToLFormaterScheduler(result_queue, table_queue)
        for _ in range(threads)
    ]

    # feeding input queue
    while offset <= count:
        request_params_queue.put((nb_per_page, offset))
        offset += nb_per_page
        
    # closing & joining
    for _ in range(threads):
        request_params_queue.put("DONE")
    for scrapper in scrappers:
        scrapper.join()
    for _ in range(threads):
        result_queue.put("DONE")
    for formater in formaters:
        formater.join()
    table_queue.put("DONE")

    # writing formated table
    table = []
    while True:
        try:
            res = table_queue.get(timeout=10)
        except Empty:
            break
        if res == "DONE":
            break
        table.append(res)

    species_df = pd.DataFrame.from_dict(table)    
    return species_df

def load_ref(ref_dict, outpath, chroms = None, threads = 8, sequence_types = 'chromosome,organelle', sbatch = False, verbose = False, **sbatch_params):
    """Builds the arborescence and eventualy retrieve fastas from the GCA database before building bowtie indexes."""
    logger = get_logger(outpath, verbose = verbose)

    # building script
    header = get_header(sbtach = sbatch, outpath = outpath, verbose = verbose, **sbatch_params)
    path_to_scripts = "/".join(__file__.split('/')[:-1])
    load_ref_script = path_to_scripts + "/scripts/loaders/load_ref.sh"
    format_ref_script = path_to_scripts + "/scripts/loaders/format_ref.sh"
    scripts_outpath = f"{outpath}/scripts"
    build_arborescence(scripts_outpath)
    script = write_script(header, {"load":(load_ref_script, 3), "format":(format_ref_script, 4)}, scripts_outpath, name = "load_ref")

    # initialisations
    genome_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(genome_queue, script, logger = logger) for _ in range(threads)
    ]

    indexes = {}

    # feading fastas to the input queue
    for ref in ref_dict:
        species = ref.replace(" ","_")
        indexes[species] = []

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
            genome_name = species + f"_{k + 1}"
            loaded_fasta_path = f"{outfolder}/{genome_name}.all_seqs.fa" if not os.path.exists(ref_dict[ref][k]) else ref_dict[ref][k]

            indexes[species].append(f"{outfolder}/{genome_name}")

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

    return workers, indexes

def load_fastqs(fastq_dict, outpath, experiment, threads = 8, sbatch = False, verbose = False, **sbatch_params):
    """Builds the arborescence and eventualy retrieve fastq files with fasterq-dump before splitting in subfastqs."""
    logger = get_logger(outpath, verbose = verbose)

    # building script
    header = get_header(sbtach = sbatch, outpath = outpath, verbose = verbose, **sbatch_params)
    path_to_scripts = "/".join(__file__.split('/')[:-1])
    load_fastq_script = path_to_scripts + "/scripts/loaders/load_fastq.sh"
    scripts_outpath = f"{outpath}/scripts"
    build_arborescence(scripts_outpath)
    script = write_script(header, {"load":(load_fastq_script, 3)}, scripts_outpath, name = "load_fastq")

    # initialisations
    fastq_queue = Queue()
    workers = [
        scexec.ScriptExecutorScheduler(fastq_queue, script, logger = logger) for _ in range(threads)
    ]
    fastqs_refs = {}

    # feading fastas to the input queue
    for ref in fastq_dict:
        species = ref.replace(" ","_")
        fastqs_refs[species] = []

        # creating outfolder
        outfolder = f'{outpath}/{species}/dataset/{experiment}'
        build_arborescence(outfolder)

        # parameters
        for k in range(len(fastq_dict[ref])):
            if not isinstance(fastq_dict[ref][k], str):
                fastqs, pairing = fastq_dict[ref][k]
                fastq1 = fastqs[0]
                fastq2 = fastqs[1] if len(fastqs) > 1 else ""
            else:
                fastq1 = fastq_dict[ref][k]
                fastq2 = ""

            # scripts to call
            to_load1 = "true" if not os.path.isfile(fastq1) else 'none'
            to_load2 = "true" if pairing == "PAIRED" and not os.path.isfile(fastq2) else 'none'

            fastq1_path = f"{outfolder}/{fastq1}"
            fastq2_path = ""
            match pairing:
                case "PAIRED":
                    fastq1_path = fastq1 if to_load1 == 'none' else f"{outfolder}/{fastq1}_1.fastq"
                    fastq2_path = fastq2 if to_load2 == 'none' else f"{outfolder}/{fastq1}_2.fastq"
                case "SINGLE":
                    fastq1_path = fastq1 if to_load1 == 'none' else f"{outfolder}/{fastq1}.fastq"
                    fastq2_path = ""
            
            fastqs_refs[species].append(([fastq1_path, fastq2_path], pairing))

            # queuing first fastq
            fastq_queue.put(([
                to_load1,
                fastq1,
                outfolder,
                str(threads)
            ], sbatch))

            # if paired, queuing second fastq
            if len(fastq2) > 0:
                fastq_queue.put(([
                    to_load2,
                    fastq2,
                    outfolder,
                    str(threads)
                ], sbatch))

    return workers, fastqs_refs

def load(
    name = "",
    outpath = None,
    ref = None,
    chroms = None,
    wgs = [],
    hic = [],
    table = None,
    dtol = False,
    no_wgs = False,
    no_hic = False,
    threads =  8,
    sbatch = False,
    sbtach_partition = 'dedicated',
    sbtach_qos = 'fast',
    sbtach_mem = '40G',
    sbatch_ncpus = 30,
    verbose = False,
    **kwargs
):
    sbatch_params = {
        '--partition': sbtach_partition,
        '--qos': sbtach_qos,
        '--mem': sbtach_mem,
        '-c': sbatch_ncpus
    }

    logger = get_logger(outpath, verbose = verbose)
    logger.info(f"load: starting (name={name!r}, sbatch={sbatch})")

    outfolder = f"{outpath}/{name.replace(' ', '_')}"
    build_arborescence(outfolder)

    ## building data tables
    if dtol:
        data = load_dtol(threads = threads)
    elif table is not None:
        data = pd.read_csv(table)
    else:
        data = create_table(name, ref, wgs, hic, no_hic = no_hic, no_wgs = no_wgs)
        if data is None:
            print("Error: missing WGS or HIC input.")
            return
    data.to_csv(f"{outfolder}/samples_data.csv", index=False)

    ## retrieving unique name subset
    subset_ref = {}
    subset_wgs = {}
    subset_hic = {}
    is_single_name = True
    names = np.unique(data['name'])
    for subname in names:
        if is_single_name and name != subname:
            is_single_name = False
        subset = data[data['name'] == subname]
        subset_ref[subname] = np.unique(subset['ref'])
        subset_wgs[subname] = np.unique(subset['wgs'])
        subset_hic[subname] = np.unique(subset['hic'])

    ## launching loaders
    outtmp = outpath if is_single_name else outfolder
    ref_loaders, indexes = load_ref(subset_ref, outpath = outtmp, chroms = chroms, threads = threads, sbatch = sbatch, verbose = verbose, **sbatch_params) if ref is not None or dtol or table is not None else ([], {})
    wgs_loaders, wgs_fastqs = load_fastqs(subset_wgs, outpath = outtmp, experiment='wgs', threads = threads, sbatch = sbatch, verbose = verbose, **sbatch_params) if not no_wgs else ([], {})
    hic_loaders, hic_fastqs = load_fastqs(subset_hic, outpath = outtmp, experiment='hic', threads = threads, sbatch = sbatch, verbose = verbose, **sbatch_params) if not no_hic else ([], {})

    ### joining loaders
    join_workers(ref_loaders)
    join_workers(wgs_loaders)
    join_workers(hic_loaders)

    ## creating dataset for other methods usage
    dataset = build_dataset(names, indexes, wgs_fastqs, hic_fastqs, **kwargs)
    dataset.to_csv(f"{outfolder}/dataset.csv")

    logger.info(f"load: done (name={name!r})")

    return f"{outfolder}/dataset.csv"