# Export the VCF from Plink or radiator with all the necessary filters 
# Edit the file with Textwrangler or Excel so there are only the rows with SNP information (nothing else above or below)

###################################################################################################################  
###########################  # Be careful and erase the # before CHROM in the vcf file. ###########################  
###########################   If not you will not be able to open it properly           ###########################  
###################################################################################################################


library(dplyr)
library(tidyr)


## list all edited vcf files into a feature
rm(list = ls())
vcf_files <- list.files(path = ".", pattern = "\\.vcf$", full.names = F)
vcf_files
## Create an empty dataframe 

missing_data_datasets <- data.frame(dataset = NA, num_loci = NA, SNP.numb = NA,  percent = NA)


## Loop all the datasets to find their missingness
for (j in 1:length(vcf_files)) {
  
  print(paste("starting dataset", vcf_files[j], sep = " ")) 
  
vcf <- read.table(vcf_files[j], sep = '\t', header = TRUE)

dt <- vcf %>% select(10:length(vcf))  # Be careful to select only the columns with individuals. That could vary between programs head(vcf)
#head(vcf)
#head(dt)

missing_data_table <- data.frame(sample = NA, num_loci = NA, SNP.numb = NA, missing_numb = NA, percent = NA)
name <- as.factor(names(dt))
number_of_sites<- NROW(dt)
for (i in 1:length(dt)) {
  snp <- as.character(name[i])
  Ind <- dt[,i]
  missing_data <- which(Ind == './.' | Ind == './1' | Ind == './0' | Ind == '0/.' | Ind == '1/.')
  loci_ind <- vcf[-missing_data,]
  missing_sites = length(missing_data)
  percent = missing_sites/number_of_sites
  missing_data_table[i,1] <- snp
  missing_data_table[i,2] <- length(unique(loci_ind[,1]))
  missing_data_table[i,3] <- number_of_sites - missing_sites
  missing_data_table[i,4] <- missing_sites
  missing_data_table[i,5] <- percent
  print(paste("specimen", i , "processed", sep = " "))
  
}

Total.num.Loci <- length(unique(vcf[,1]))
Total.SNP.numb <- sum(missing_data_table$SNP.numb)
Total_snps <- length(unique(vcf$ID))
Total.missing_numb <- sum(missing_data_table$missing_numb)
Total_missing_data_table <- data.frame(sample = 'TOTAL', num_loci = Total.num.Loci, SNP.numb = Total_snps, missing_numb = NA, percent = Total.missing_numb/(Total.SNP.numb+Total.missing_numb)*100)

missing_data_table <- rbind(missing_data_table, Total_missing_data_table)
write.table(missing_data_table, paste("RESULTS/",vcf_files[j], ".csv", sep = ""), sep = ';', quote = FALSE, row.names = FALSE)
missing_data_dataset <- data.frame(dataset = c(vcf_files[j]), num_loci = Total.num.Loci, SNP.numb = Total_snps, percent = Total.missing_numb/(Total.SNP.numb+Total.missing_numb)*100)
missing_data_datasets[j,] <- missing_data_dataset
}

write.table(missing_data_datasets, "RESULTS/missing_data_datastes.csv", sep = ';', quote = FALSE, row.names = FALSE) 




################################################
#### FOR NON POSTPROCESSED FILTERED DATASETS ###
################################################

rm(list = ls())
## list all edited vcf files into a feature
vcf_files <- list.files(path = ".", pattern = "\\.vcf$", full.names = F)
vcf_files
## Create an empty dataframe 

missing_data_non_processed <- data.frame(dataset = NA, num_loci = NA, SNP.numb = NA,  percent = NA)


## Loop all the datasets to find their missingness



for (j in 1:length(vcf_files)) {
  
  print(paste("starting dataset", vcf_files[j], sep = " "))
  
  vcf <- read.table(vcf_files[j], sep = '\t', header = TRUE)
  
  dt <- vcf %>% select(10:length(vcf))  # Be careful to select only the columns with individuals. That could vary between programs head(vcf)
  
  missing_data_table <- data.frame(sample = NA, num_loci = NA, SNP.numb = NA, missing_numb = NA, percent = NA)
  name <- as.factor(names(dt))
  number_of_sites<- NROW(dt)


for (i in 1:length(dt)) {
  snp <- as.character(name[i])
  dtt <- separate(data = dt, col = i, into = c('GENOTYPE','DEPTH', 'BASES'), sep = ":")
  Ind <- dtt$GENOTYPE
  missing_data <- which(Ind == './.' | Ind == './1' | Ind == './0' | Ind == '0/.' | Ind == '1/.')
  loci_ind <- vcf[-missing_data,]
  missing_sites = length(missing_data)
  percent = missing_sites/number_of_sites
  missing_data_table[i,1] <- snp
  missing_data_table[i,2] <- length(unique(loci_ind[,1]))
  missing_data_table[i,3] <- number_of_sites - missing_sites
  missing_data_table[i,4] <- missing_sites
  missing_data_table[i,5] <- percent
  print(paste("specimen", i , "processed", sep = " "))
}

  Total.num.Loci <- length(unique(vcf[,1]))
  Total_snps <- length(unique(vcf$ID))
  Total.SNP.numb <- sum(missing_data_table$SNP.numb)
  Total.missing_numb <- sum(missing_data_table$missing_numb)
  Total_missing_data_table <- data.frame(sample = 'TOTAL' , num_loci = Total.num.Loci, SNP.numb = Total_snps, missing_numb = NA, percent = Total.missing_numb/(Total.SNP.numb+Total.missing_numb)*100)
  
  missing_data_table <- rbind(missing_data_table, Total_missing_data_table)
  write.table(missing_data_table, paste("./",vcf_files[j], ".csv", sep = ""), sep = ';', quote = FALSE, row.names = FALSE)
   missing_data_process <- data.frame(dataset = c(vcf_files[j]), num_loci = Total.num.Loci, SNP.numb = Total_snps, percent = Total.missing_numb/(Total.SNP.numb+Total.missing_numb)*100)
  missing_data_non_processed[j,] <- missing_data_process
}

write.table(missing_data_non_processed, "./missing_data_non_processed.csv", sep = ';', quote = FALSE, row.names = FALSE) 
