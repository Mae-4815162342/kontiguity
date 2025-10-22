from kontiguity.utils.functions import *    

class FastaFormater():
    def __init__(self, sequence_types=['chromosome', 'organelle'], batch_size = 20):
        self._sequence_types = sequence_types
        self._batch_size = batch_size

    @staticmethod
    def parse_ENA_info(info_line):
        """Parses the ENA sequence info line from fasta file"""
        infos = info_line.replace('\n','').split(" ")
        id = infos[0][1:]
        leftover = " ".join(infos[1:]).split(',')
        description = leftover[0]
        sequence_type = leftover[1].replace(" ", "").split(":")
        return {
            'id':id,
            'description':description,
            'sequence_type':sequence_type[0],
            'sequence_name':sequence_type[1]
        }

    def format_fasta(self, fasta_path, outfolder):
        """Writes a new fasta file keeping only the sequences with the set sequence types. Writes the chromosomes infos in a tsv file."""
        chromosome_infos = []
        chrom_file = None
        chromosomes = []
        batched_lines = []
        fasta_name = fasta_path.split('/')[-1][:-len(".all_seqs.fa")]
        output_fasta = open(outfolder + f"/{fasta_name}.filtered.fa", "w")
        has_chrom_info = os.path.isfile(f"{outfolder}/{fasta_name}.chromosomes.csv") or os.path.isfile(f"{outfolder}/chromosomes.csv")

        if has_chrom_info:
            if os.path.isfile(f"{outfolder}/chromosomes.csv"):
                shutil.copyfile(f"{outfolder}/chromosomes.csv", f"{outfolder}/{fasta_name}.chromosomes.csv")
            chrom_file = pd.read_csv(f"{outfolder}/{fasta_name}.chromosomes.csv")
            chromosomes = np.unique(chrom_file['id'])

        with open(fasta_path, "r") as fasta:
            keep_current = False
            line = fasta.readline()
            while len(line) > 1:
                if line[0] == ">":
                    if not has_chrom_info:
                        infos = self.parse_ENA_info(line)
                        keep_current = infos['sequence_type'] in self._sequence_types
                        if keep_current:
                            chromosome_infos.append(infos)
                    else:
                        chromosome = line[1:].split(' ')[0]
                        keep_current = chromosome in chromosomes and list(chrom_file[chrom_file['id'] == chromosome]['sequence_type'])[0] in self._sequence_types

                if keep_current:
                    batched_lines.append(line)
                    if len(batched_lines) >= self._batch_size:
                        output_fasta.writelines(batched_lines)
                        batched_lines = []
                line = fasta.readline()

            if len(batched_lines) > 0:
                output_fasta.writelines(batched_lines)

        output_fasta.close()

        if not has_chrom_info:
            chromosome_infos_df = pd.DataFrame.from_dict(chromosome_infos)
            chromosome_infos_df.to_csv(outfolder + f"/{fasta_name}.chromosomes.csv", sep=',')

        return outfolder + f"/{fasta_name}.filtered.fa"

fasta = sys.argv[1]
sequence_types = sys.argv[2].split(',')
outfolder = sys.argv[3]

formater = FastaFormater(sequence_types = sequence_types)
formater.format_fasta(fasta, outfolder)