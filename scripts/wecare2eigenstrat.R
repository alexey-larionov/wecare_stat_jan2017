# Code fragment to illustrate export of wecare data to EIGENSTRAT format
# Alexey Larionov, 07Feb2017

# Notes:
# This is not executable sctipt! This is just a code fragment to illustrate the procedure. 
# The EIGENSTRAT format specifications are taken from EIGENSTRAT tool documentation.  

# -------------------- make snp file -------------------- #

# Extract data from variants table and emulate SNPs positions in "morgans". 
# Not sure whether it is or should be in morgans or santi-morgans? 
# It does not matter in our case: the "morgan" column is not used 
# for eigenvectors calculations anyway. 
library(dplyr)
wecare_snp.df <- wecare_variants.df %>% 
  mutate(morgan = POS / 1000000) %>%  
  select(SplitVarID, CHROM, morgan, POS, REF, ALT)

# Recode CRHOM data
wecare_snp.df$CHROM <- as.vector(wecare_snp.df$CHROM)

"23" -> wecare_snp.df[wecare_snp.df$CHROM == "X", "CHROM"]
"24" -> wecare_snp.df[wecare_snp.df$CHROM == "Y", "CHROM"]
"90" -> wecare_snp.df[wecare_snp.df$CHROM == "MT", "CHROM"]

wecare_snp.df$CHROM <- as.factor(wecare_snp.df$CHROM)

# Write file (tab-separated)
write.table(wecare_snp.df, 
            paste(interim_data_folder, "wecare.snp", sep="/"), 
            quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=FALSE) 

# -------------------- make geno file -------------------- #

# Get data from genotypes matrix
wecare_geno.mx <- wecare_genotypes.mx

# Recode NA
9 -> wecare_geno.mx[is.na(wecare_geno.mx)] # NA is coded as 9

# Write file (no delimiters)
write.table(wecare_geno.mx, 
            paste(interim_data_folder, "wecare.geno", sep="/"),             
            quote=FALSE, sep="", 
            row.names=FALSE, col.names=FALSE)

# -------------------- make ind file -------------------- #

# Samples column
samples <- as.vector(wecare_phenotypes.df$wes_id)

# Gender column
gender <- rep("F", length(samples)) # all females

# Population/Group column
group <- as.vector(wecare_phenotypes.df$cc)
"UBC" -> group[group == 0]
"CBC" -> group[group == 1]

# Get together
wecare_ind.mx <- cbind(samples, gender, group)

# Write file (tab-separated) 
write.table(wecare_ind.mx, 
            paste(interim_data_folder, "wecare.ind", sep="/"),
            quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=FALSE)

# -------------------- The End -------------------- #
