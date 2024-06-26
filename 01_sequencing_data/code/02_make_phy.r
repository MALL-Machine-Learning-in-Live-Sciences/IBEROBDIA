# Load packages
library(phyloseq)
library(dplyr)
# Load Clinical Data and prepare paths
metadata = readRDS("00_preprocess_metadata/data/clinical_data.rds") 
path ="01_sequencing_data/data/"
truncL <- sort(list.files(path, pattern = "TL_251"))
nams <- paste(truncL, "phy", sep="_")
phy_list = list()
# Phylo.object construction
for (i in seq_along(truncL)) {
  # Load and merge all data in phy.object
  tax = readRDS(file = paste(path,truncL[i],"tax_table.rds",sep = "/"))
  otu = readRDS(file = paste(path,truncL[i],"otu_table.rds",sep = "/"))
  otu = phyloseq::t(otu)
  ps = phyloseq(otu_table(otu, taxa_are_rows = TRUE), 
                 sample_data(metadata), 
                 tax_table(tax))
  
  # Shorten the name of our ASVs
  dna <- Biostrings::DNAStringSet(taxa_names(ps))
  names(dna) <- taxa_names(ps)
  ps <- merge_phyloseq(ps, dna)
  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
  
  # Correctly rename the species for a correct merge later on.
  # This is done in order to have a complete species name and
  #  not to make mistakes when agglomerating. 
  ## Replace NAs for 0 to filter spss in order to add genus name before
  tax <- data.frame(tax_table(ps))
  tax$Species[is.na(tax$Species)] <- 0
  ## Add genus name to Species names
  tax <- tax %>% 
    mutate(Species = ifelse(Species != 0, paste0(Genus,"_",Species), Species))
  ## Replace 0 for NA again
  tax <- tax %>% 
    mutate(Species = ifelse(Species == 0, NA, Species))
  tax = as.matrix(tax)
  tax_table(ps) <- tax
  phy_list[[i]] = ps
}

# Rename and save phy.objects
names(phy_list) = nams
for (i in seq_along(phy_list)) {
  saveRDS(object =phy_list[[i]],file = paste(path, truncL[i], paste0(names(phy_list[i]),".rds"),sep = "/" ))
}
