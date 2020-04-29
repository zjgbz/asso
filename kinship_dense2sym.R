library('data.table')
library(feather)

startTime <- Sys.time()

path <- '../raw_data/kinship/freeze8/'
# load the dense kinship matrix
kmatr <- fread(cmd = paste0('zcat ', path, 'pcrelate_dense_km.txt.gz'), data.table = FALSE)
### NOTE: fread is going to write a temp file to unzip the .gz file. 
### If you get the error "gzip: stdout: No space left on device" you have a couple options:
### 1) set the environment variable TMPDIR to somewhere with enough disk space
### 2) uncompress the file on your system, then read in the uncompressed file with fread()

# set rownames
row.names(kmatr) <- colnames(kmatr)
dim(kmatr)

kmatr = as.matrix(kmatr)

# copy upper triangle down to create symmetrical matrix
makeSymm <- function(m) {
   m[lower.tri(m)] <- t(m)[lower.tri(m)]
   return(m)
}
kmatr = makeSymm(kmatr)

# confirm there are no NA rows: result should be FALSE
if(any(is.na(rowMeans(kmatr)))) warning('There are some NAs in the matrix!')

### if you want the matrix to be on the same scale as the standard GRM
### you need to double the values (entries are then 2*kinship)
# kmatrS <- 2*kmatrS

# save your subset kinship matrix

kmatr = as.data.frame(kmatr)
write_feather(kmatr, "../prepro_data/kinship/freeze8_kinship.feather")

endTime <- Sys.time()
message('Kinship matrix created and saved in ', round(difftime(endTime, startTime, units = 'mins'), 1), ' minutes')
