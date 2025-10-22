from kontiguity.utils.imports import *

class DToLFormaterScheduler(threading.Thread):
    def __init__(self, input_queue, output_queue, **args):
        super(DToLFormaterScheduler, self).__init__(**args)
        self._input_queue = input_queue
        self._output_queue = output_queue
        self.start()

    def run(self):
        formater = DToLFormater()
        while True:
            try:
                value = self._input_queue.get()
            except Empty:
                break
            if value == "DONE":
                break
            output = formater.extract_infos(value)
            if output is not None:
                species, genomes, HiCs, WGSs = output
                for row in formater.format(species, genomes, HiCs, WGSs):
                    self._output_queue.put(row)

class DToLFormater():
    def __init__(self):
        pass

    @staticmethod
    def format(species, genomes, HiCs, WGSs):
        """creates a table row for each combination of assembly, Hi-C and WGS for the species."""
        for i in range(len(genomes)):
            for j in range(len(WGSs)):
                for k in range(len(HiCs)):
                    yield {
                        "species": species,
                        "name": species,
                        **genomes[i],
                        **WGSs[j],
                        **HiCs[k]
                    }
    
    @staticmethod
    def extract_infos(species_data):
        species = species_data["_id"]

        # retrieving Hi-C and WGS. If no Hi-C or WGS is available, returns None
        HiCs = []
        WGSs = []
        column = {
            "Hi-C":"hic",
            "WGS":"wgs"
        }
        for experiment in species_data["_source"]["experiment"]:
            strategy = experiment["library_strategy"]
            if strategy in ["Hi-C", "WGS"]:
                exp_datas = {
                    f"{column[strategy]}": experiment["study_accession"],
                    f"protocol_{strategy}": experiment["library_construction_protocol"] if len(experiment["library_construction_protocol"]) > 0 else strategy,
                    f"run_acession_{strategy}": experiment["run_accession"],
                    f"layout_{strategy}": experiment["library_layout"]
                }
            match strategy:
                case "Hi-C":
                    HiCs.append(exp_datas)
                case "WGS":
                    WGSs.append(exp_datas)
        
        if len(HiCs) == 0 or len(WGSs) == 0:
            return None

        genomes = []
        if "annotation" in species_data["_source"]: #  if the species is annotated, the reference genome used will be the one corresponding to the annotation. 
            for annotation in species_data["_source"]["annotation"]:
                genomes.append({
                    "ref": annotation["accession"],
                    # "study_accession":"",
                    # "fasta":f"https://www.ebi.ac.uk/ena/browser/api/fasta/{annotation["accession"]}?download=true&gzip=true",
                })


        # Useless: only retrieving 20 more species for a several times longer execution time
        # elif "assemblies" in species_data["_source"]: # Otherwise, we retrieve all the available genomes from the assemblies.
        #     for assembly in species_data["_source"]["assemblies"]:
        #         genome_location = f"https://ftp.ebi.ac.uk/pub/ensemblorganisms/{species.replace(' ', '_')}/{assembly['accession']}.{assembly['version']}/genome"
        #         # checking the address exists
        #         try:
        #             time.sleep(0.5) # sleeping to avoid max retries connection refused from EBI
        #             response = requests.get(genome_location)
        #             if response.status_code != 200:
        #                 continue
        #         except Exception as e:
        #             continue
        #         genomes.append({
        #             "accession": assembly["accession"],
        #             "study_accession": assembly["study_accession"],
        #             "fasta": f"{genome_location}/softmasked.fa.gz",
        #             "chromosomes": f"{genome_location}/chromosomes.tsv.gz",
        #         })

        if len(genomes) == 0:
            return None
        
        return species, genomes, HiCs, WGSs
