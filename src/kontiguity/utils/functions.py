from .imports import *

def build_arborescence(path):
    """Creates the necessary folders for the provided path to be a valid folder path."""
    folders = path.split('/')
    
    for i in range(1, len(folders) + 1):
        current_path = "/".join(folders[:i])
        if not os.path.exists(current_path) and current_path != "":
            os.mkdir(current_path)

def join_workers(workers):
    """Joins each worker in the workers list."""
    for worker in workers:
        worker.join()

def get_header(sbtach = False, outpath = None, **sbatch_params):
    """Gets the header of the bash script to execute. Applies sbatch parameters if provided."""
    header = "#!/bin/bash\n"
    if sbtach:
        # adding sbatch defined parameters
        for param in sbatch_params:
            header += f"#SBATCH {param}={sbatch_params[param]}\n"

        # adding output path
        header += f"#SBATCH -o {outpath}/logs/R-%x-%j-reference_load.out.txt -e {outpath}/logs/R-%x-%j-reference_load.err.txt\n"
    return header

def write_script(header, subscripts, outpath, name = "script"):
    """Writes a script that executes the subscripts with condition (true, false). Parameters will be passed in the order of the scripts."""
    ### when calling the script, ALL the parameters must be provided, even for scripts set to false (none will do perfectly).
    ### script first arguments are the execution condition of each script, in order (set to "true" to execute, anything else not to).
    script_path = outpath + f"/{name}.sh"
    with open(script_path, 'w') as script:
        script.write(header + "\n")
        nb_args = 1 + len(subscripts)
        nb_scripts = 1
        for script_name in subscripts:
            subscript, args = subscripts[script_name]
            script.write(f'if [ ${nb_scripts} = "true" ];then\n')
            nb_scripts += 1
            arguments = ""
            for arg in range(1,args + 1):
                script.write(f"\tparam_{arg}=${nb_args}\n")
                arguments += f"\t$param_{arg} "
                nb_args += 1
            script.write(f"\tbash {subscript} {arguments}\n")
            script.write('fi\n\n')

    return script_path