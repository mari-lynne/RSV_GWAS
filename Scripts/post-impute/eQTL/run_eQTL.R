#Run matrix eQTL ####
library(MatrixEQTL)

#Set up file names and outputs ####

# Genotype file name
setwd("~/RNA/Daniel/mEQTL")
base.dir = getwd()
SNP_file_name = paste0(base.dir, "/Fc_genoD0.txt");
snps_location_file_name = paste0(base.dir, "/snp_loc.txt");

# Gene expression file name
expression_file_name = paste0(base.dir, "/Fc_exprsD0.txt");
gene_location_file_name = paste0(base.dir, "/gene_loc.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir,"/Fc_covarD0.txt");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

#Thresholds ####
# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-1;
pvOutputThreshold_tra = 1e-1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# Distance for local gene-SNP pairs
cisDist = 1e6;

#Load data in ####
#Slice reading is overkill for this small experiment but will be useful for GWAS

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      # Space delim
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " ";      # Space delim
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # Space delim
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Load location files
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);


## Run the analysis ###

#Define model
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name      = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis  = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results: ####

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

D0eqtl <- (me$trans$eqtls)

#Day0 summary:
D0eqtl <- filter(D0eqtl, pvalue<=0.05)

## Make the histogram of local and distant p-values
plot(me)


