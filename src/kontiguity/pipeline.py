from .load import load
from .retrieve import retrieve
from .map import map
from .describe import describe
from kontiguity.utils.functions import *

def next(step):
    """Returns the next step of the pipeline."""
    STEPS = {
        "load":"retrieve",
        "retrieve":"map",
        "map":"describe"
    }
    return STEPS[step] if step in STEPS else None

def build_dataset(name, **args):
        """Builds global parameter table for each """
        binnings = np.array(args["binnings"]).astype("str")
        rows = []
        mapping_index = 1
        contig_index = 1
        for wgs, wgs_pairing in args['wgs']:
            for hic, hic_pairing in args["hic"]:
                rows.append({
                    "name":name,
                    "index":args['index'],
                    "fastq1_wgs": wgs[0], # retrieve params
                    "fastq2_wgs": wgs[1] if wgs_pairing == "PAIRED" else "",
                    "is_paired_wgs":wgs_pairing,
                    "contigs":f"contigs_{contig_index}",
                    "min_size":args["min_size"],
                    "fastq1_hic":hic[0], # map params
                    "fastq2_hic":hic[1] if hic_pairing =="PAIRED" else "",
                    "mapping":f"mapping_{mapping_index}",
                    "enzymes":args["enzymes"],
                    "binnings":";".join(binnings),
                    "format":args["format"],
                    "chroms":args["chroms"], # describe params
                    "cool":args["cool"],
                    "mcool":f"mapping_{mapping_index}/mapping_{mapping_index}.mcool",
                    "formats":";".join(args["formats"])
                })
                mapping_index += 1
            contig_index += 1
        return pd.DataFrame.from_dict(rows)

def pipeline(name, outpath, first_step = "load", last_step = "describe", **args):
    """
    Calls each method between the first and the last required step in the following order: load -> retrieve -> map -> describe.
    Provided **args must be the arguments of the first required step as described in the single commands.
    Will execute the commands sequentially, connecting each step with the next by the "table" argument. 
    """
    current_step = first_step
    dataset = ""
    build_arborescence(outpath)

    if first_step != "load":
        dataset = f"{outpath}/dataset.csv"
        dataset_df = build_dataset(name, **args)
        dataset_df.to_csv(dataset)

    while True:
        match current_step:
            case "load": # 1. Load dataset
                dataset = load(name, outpath, **args)
            case "retrieve": # 2. Retrieve contigs
                _ = retrieve(name, outpath, **args) if dataset is None else retrieve(name, outpath, table = dataset)
            case "map": # 3. Map contigs in Hi-C
                _ = map(name, outpath, **args) if dataset is None else map(name, outpath, table = dataset)
            case "describe": # 4. Descriptive statistics on selected contigs
                _ = describe(name, outpath, **args) if dataset is None else describe(name, outpath, table = dataset)
            case None:
                break
        if current_step == last_step or current_step is None:
            break
        current_step = next(current_step)