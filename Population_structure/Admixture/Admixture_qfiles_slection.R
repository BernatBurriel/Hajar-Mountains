## Script to assign K values into a matrix for plotting in A MapThis script works together with the output generated with admixture wrapper
# https://github.com/dportik/admixture-wrapper

{rm(list = ls()) 
  library(dplyr)
  library(tidyr)
  library(pophelper)}

#setwd("PATH_TO_FILES")
setwd("PATH_TO_DIRECTORY")

# import the nosex file generated with plink
ped_files <- list.files(path = "PATH/TO/DIRECTORY", pattern = "\\.nosex$", full.names = TRUE)
ped_file <- read.table(ped_files, row.names = 2, header = FALSE, sep = "\t")
ped_file <- ped_file %>% dplyr::select("Samples"=V1)
ped_file$Samples <- row.names(ped_file)


all_cv_file <- read.table(list.files(path = "PATH/TO/ADMIXTURE_FOLDER", pattern = "\\All.txt$", full.names = TRUE), header = T)
clusters <- unique(all_cv_file$K)
list_of_qfiles <- data.frame(K = NA, Rep = NA, CV = NA)
for (i in 1:length(unique(all_cv_file$K))) {
  grp <- all_cv_file[all_cv_file$K == i,] %>%
    arrange(CV)
  list_of_qfiles[i,] = grp[1,]
}

list_of_qfiles
list_of_qfiles$pattern = paste(list_of_qfiles$K, ".", list_of_qfiles$Rep, ".Q$", sep = "")

folder <- list.files("PATH/TO/ADMIXTURE_FOLDER/", pattern = "Out", full.names = T)
qfiles <- c()
for (i in 1:NROW(list_of_qfiles)){
  qfiles <- c(qfiles, list.files(path = folder, pattern = list_of_qfiles$pattern[i], full.names = TRUE))  
}

slist <- readQ(qfiles, filetype = "auto", indlabfromfile = T)
tbl <- alignK(slist)

for (i in 1:length(tbl)) {
  name <- c(1:i)
  names_vec <- paste("K", i, ".", name, sep = "")
  tbl[[i]] <- tbl[[i]]/rowSums(tbl[[i]]) * 100
  tbl[[i]] <- round(tbl[[i]], digits = 0)
  colnames(tbl[[i]]) <- names_vec 
  ped_file <- cbind(ped_file, tbl[[i]])
}

# We will merge our admixture data with a dataset containing the coordinates of each individual 
sample <- list.files(path = "PATH/to/CSV", pattern = "\\.csv$", full.names = TRUE)
samples <- read.csv2(sample, header = TRUE)
head(samples)

# Check that all our samples are in the dataset
ped_file[!ped_file$Samples %in% samples$Specimen.Code,] ### MODIFY the column names to the ones that you have in your dataset

samples <- left_join(samples, ped_file, by = c("Specimen.Code" = "Samples"))
head(samples)
samples_admixture <- samples[!is.na(samples$K1.1),]
write.csv2(samples_admixture, "only_admixture_sample_localities.csv", quote = FALSE)