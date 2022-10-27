setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/IBEROBDIA/')
source('03_training/config_file_mlr3.r')
source('03_training/rfmlr3_slurm.r')

args = commandArgs(trailingOnly = T)
data = readRDS(args[1])
names(data) = make.names(names(data))
name = substr(args[2], 1, nchar(args[2]))
set.seed(seed)

rf.bmr.slurm(data = data,
             name = name,
             path = out.path,
             filename = out.filename.rf,
             cv.inner = cv.inner, 
             cv.outer = cv.outer,
             rf.mtry = rf.mtry,
             rf.ntree = rf.ntree,
             rf.nodesize = rf.nodesize,
             workers = workers)
