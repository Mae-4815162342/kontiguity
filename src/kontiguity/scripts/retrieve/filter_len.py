# filters fasta sequences per length and keeps only the ones above a given length in output_fasta
from kontiguity.utils.functions import *

def filter_fasta_per_length(fasta, length, output_fasta):

    output_folder = "/".join(output_fasta.split('/')[:-1])
    build_arborescence(output_folder)

    filtered = open(output_fasta, 'w')

    current_id = ""
    current_lines = []
    current_length = 0

    with open(fasta, 'r') as unfiltered:
        line = unfiltered.readline()
        while line is not None and len(line) > 0:
            if line[0] == '>':
                # writing current sequence
                if current_length >= length:
                    filtered.writelines([current_id] + current_lines)
                current_id = line
                current_lines = []
                current_length = 0
            else:
                current_lines.append(line)
                current_length += len(line.replace("\n", ""))

            line = unfiltered.readline()

    if current_length >= length:
        filtered.writelines([current_id] + current_lines)
    
    filtered.close()

if len(sys.argv) == 1: # local test
    local_path = "/".join(__file__.split('/')[:-5])
    fasta = local_path + "/test_data/retrieve/contigs.fa"
    min_length = 1000
    output_fasta = local_path + f"/test_results/filtering/contigs_filtered_{min_length}.fa"

else: # global usage
    fasta = sys.argv[1]
    min_length = int(sys.argv[2])
    output_fasta = sys.argv[3]

filter_fasta_per_length(fasta, min_length, output_fasta)

if len(sys.argv) == 1:
    assert(os.path.isfile(output_fasta))
