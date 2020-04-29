### user provided inputs ###
# supply one of phenofile + idcolumnname OR sample.id.file
phenofile <- NULL # 'path_to_phenotype_data.csv'
idcolumnname <- 'sample.id' # column name of sample IDs in pheno.file

#sample.id.file <- NULL # 'path_to_vector_of_sample_ids.RData'
sample.id.file <- "/proj/yunligrp/users/minzhi/asso/cohort/SOL/ready_data/common_sample_SOL.RData"

# path to directory where the matrix is stored
path <- './'

# name of output matrix
output.file <- 'common_kinship_SOL.RData'
#############################################################


library('data.table')
startTime <- Sys.time()

# load the dense kinship matrix
kmatr <- fread(cmd = paste0('zcat ', path, 'pcrelate_kinshipMatrix_dense_v2.txt.gz'), data.table = FALSE)
### NOTE: fread is going to write a temp file to unzip the .gz file. 
### If you get the error "gzip: stdout: No space left on device" you have a couple options:
### 1) set the environment variable TMPDIR to somewhere with enough disk space
### 2) uncompress the file on your system, then read in the uncompressed file with fread()

# set rownames
row.names(kmatr) <- colnames(kmatr)
dim(kmatr)

# get sample IDs that you want to keep
# if you have a CSV with phenotype data
if(!is.null(phenofile)){
    pheno <- read.csv(phenofile, as.is = TRUE)
    if(!idcolumnname %in% colnames(pheno)) stop('idcolumnname must be a column in phenofile')
    sample.keep <- pheno[[idcolumnname]]

# if you have an RData file with sample IDs
}else if(!is.null(sample.id.file)){
    sample.keep <- get(load(sample.id.file))
    
}else{
    stop('one of phenofile or sample.id.file must be specified at the top of the script')
}
message(length(sample.keep), ' sample IDs selected by user')
    
# subset matrix to just your samples
kmatrS = kmatr[row.names(kmatr) %in% sample.keep, colnames(kmatr) %in% sample.keep]
dim(kmatrS)

# memory intensive so free up space
rm(kmatr)
gc()

# convert to a matrix
kmatrS = as.matrix(kmatrS)

# copy upper triangle down to create symmetrical matrix
makeSymm <- function(m) {
   m[lower.tri(m)] <- t(m)[lower.tri(m)]
   return(m)
}
kmatrS = makeSymm(kmatrS)

# confirm there are no NA rows: result should be FALSE
if(any(is.na(rowMeans(kmatrS)))) warning('There are some NAs in the matrix!')

### if you want the matrix to be on the same scale as the standard GRM
### you need to double the values (entries are then 2*kinship)
# kmatrS <- 2*kmatrS

# save your subset kinship matrix
save(kmatrS, file = output.file)

endTime <- Sys.time()
message('Kinship matrix created and saved in ', round(difftime(endTime, startTime, units = 'mins'), 1), ' minutes')
