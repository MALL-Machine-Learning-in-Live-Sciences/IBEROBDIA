setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/IBEROBDIA/')
source('03_training/config_file_mlr3.r')
wd = getwd()
# Execute in parallel from Cesga
input.paths = dir(path = input.dir.path)
input.algs = dir(path = path.algs, pattern = pattern.algs)
#input.algs = input.algs[grep('RF', input.algs)] #Comment/uncomment in order to select only one algorithm "RF", "glmnet"...

for (i in 1:length(input.paths)) {
  for (j in 1:length(input.algs)) {
    
    #Declare paths to sh scripts and outs from CESGA
    exec.dir = exec.dir.path
    outslurm.dir = outslurm.dir.path

    if (dir.exists(exec.dir) == FALSE) {
      dir.create(exec.dir)
      message('Creating directory --Exec-- ...')
    }

    if (dir.exists(outslurm.dir) == FALSE) {
      dir.create(outslurm.dir)
      message('Creating directory --Exec-- ...')
    }
    # METERLE un -4 para eliminar el ".rds" en los inputh pathts evidatndo el "all_28.rds.sh"
    sink(paste(exec.dir, substr(input.algs[j], 1, nchar(input.algs[j]) - 7),
               substr(input.paths[i], 1, nchar(input.paths[i])-4), '.sh', sep = ''))
    cat("#!/bin/bash \n")
    
    cat(paste("#SBATCH", "-p", part, '\n'))
    cat(paste("#SBATCH", "-t", time, '\n'))
    cat(paste("#SBATCH", paste0("--mem=",mem),'\n'))
    cat(paste("#SBATCH", "-N", nodes, '\n'))
    cat(paste("#SBATCH", "-n", ntasks, '\n'))
    cat(paste("SBATCH", paste0("--error=", wd, outslurm.dir.path, "error-", substr(input.algs[j], 1, nchar(input.algs[j]) - 7),
               substr(input.paths[i], 1, nchar(input.paths[i])-4), ".txt")))
    cat(paste("SBATCH", paste0("--output=", wd, outslurm.dir.path, "output-", substr(input.algs[j], 1, nchar(input.algs[j]) - 7),
               substr(input.paths[i], 1, nchar(input.paths[i])-4),".txt")))

    cat("module load cesga/2018 gcc/6.4.0 R/4.0.2", '\n')
    
    cat(paste("name=", input.paths[i], '\n', sep = ''))
    cat(paste("data=", input.dir.path, input.paths[i], '\n', sep = ''))
    cat(paste('Rscript /mnt/netapp2/Store_uni/home/ulc/co/jlb/git/DREAM-Microbiome/02_training/models/', input.algs[j],  " $data", " $name", sep=''))
    
    sink(file = NULL)
    
    system(paste('sbatch ', '/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/DREAM-Microbiome/02_training/Exec/', substr(input.algs[j], 1, nchar(input.algs[j])-7) ,substr(input.paths[i], 1, nchar(input.paths[i])), '.sh', sep = ''))
    
  } 
}
                    