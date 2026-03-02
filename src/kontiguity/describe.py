from kontiguity.utils.functions import *
from kontiguity.workers import ScriptExecutor as scexec

def describe(
    name = "",
    outpath = "",
    index = "",
    hic = "",
    enzymes = "HinfI,DpnII",
    binnings = "5000",
    format = "cool",
    table = None,
    no_tmp = False,
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
    # if table is not None:
    #     data = pd.read_csv(table)
    # else:
    #     data = create_table(name, index, hic, enzymes = enzymes, binnings = binnings, format = format)
    #     if data is None:
    #         print("Error: missing HIC input.")
    #         return
    # data["mapping"] = [f"mapping_{k + 1}" for k in range(len(data))]
    # data.to_csv(f"{outfolder}/mapping_data.csv", index=False)
    
    # ## retrieving unique name subset
    # subset_hic = {}
    # is_single_name = True
    # for subname in np.unique(data['name']):
    #     if is_single_name and name != subname:
    #         is_single_name = False
    #     subset_hic[subname] = data[data['name'] == subname]

    # ## launching retrievers
    # out_tmp = outpath if is_single_name else outfolder
    # retrievers = map_hics(subset_hic, outpath = out_tmp, no_tmp = no_tmp, threads = threads, sbatch = sbatch, **sbatch_params)

    # ### joining loaders
    # join_workers(retrievers)