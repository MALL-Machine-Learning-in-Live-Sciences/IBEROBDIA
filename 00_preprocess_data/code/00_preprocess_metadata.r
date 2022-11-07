# Clinical Data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
## Clinical Data Extraction
path = "../../extdata/FASTQ/"
names = list.files(path = path, pattern = ".fastq.gz")
library(stringr)
library(data.table)

# Clinic Data
vars = list()
for (i in seq_along(names)) {
  vars[i] = str_remove(names[i], ".fastq.gz")
  vars[i] = str_split(vars[i], "_", n = Inf, simplify = FALSE)
}
metalist = list()
for (i in seq_along(vars)) {
  if (length(vars[[i]]) == 3){
    metalist[[i]] = c(vars[[i]][1],vars[[i]][2],"Healthy",vars[[i]][3])
  }else{
    metalist[[i]] = vars[[i]]
  }
}

clinical = as.data.frame(do.call(rbind, metalist))
colnames(clinical)= c("SampleID","CI", "Status","Gender")
clinical <- data.frame(clinical[,-1], row.names = clinical[,1])
table(clinical$CI)
# NP OB 
# 23 56 
table(clinical$Status)
# DT2    GAA Health    TGA 
# 5     12     59      3 
saveRDS(object = clinical, file = "../data/clinical_data.rds")
