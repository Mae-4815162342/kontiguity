from kontiguity.utils.functions import *

def split_file(files, section_size, nb_split, batch_size=20):
    """Splits a file or a pair of files into a subset of smaller files"""
    out_files = [
        [open(".".join(file.split(".")[:-1]) + f".split_{i}." + file.split(".")[-1]) for file in files] for i in range(nb_split)
    ]
    for i in range(len(files)):
        batch = []
        current_file = 0
        with open(files[i]) as f:
            line = f.readline()
            k = 0
            while len(line) >0:
                batch.append(line)
                line = f.readline()
                k += 1
                if k >= batch_size and k % section_size == 0:
                    out_files[current_file % nb_split][i].writelines(batch)
                    current_file += 1
                    batch = []
            if len(batch) > 0:
                out_files[current_file % nb_split][i].writelines(batch)
            

# files = sys.argv[1]
# section_size = sys.argv[2] # for instance 4 for fastq files
# nb_split = sys.argv[3]
# paired = sys.argv[4]

files = ["/mnt/c/Users/maely/Desktop/Thèse/kontiguity/test_data/S_cerevisiae/dataset/WGS/SRR35504920_1.fast", 
         "/mnt/c/Users/maely/Desktop/Thèse/kontiguity/test_data/S_cerevisiae/dataset/WGS/SRR35504920_2.fast"]
section_size = 4
nb_split = 10
split_file(files, section_size, nb_split)

