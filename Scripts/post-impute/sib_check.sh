
name=RSV_imp_QC #name of study PLINK files #can we inherit args from r script?

wdir =/home/mari/RSV/data/post-imp/Assoc
mkdir -p $wdir/plink_log
log=$wdir/plink_log

mkdir -p $wdir/qc
qcdir=$wdir/qc


plink \
--bfile $name \
--genome \
--missing \
--out 1KG_merged

# Remove one sample from each pair with pi-hat (% IBD) above threshold (0.1875 below):
awk '$10 >= 0.1875 {print $1, $2}' $qcdir/1KG_merged.genome | uniq > $qcdir/1KG_merged.outliers.txt 
wc -l 1KG_merged.outliers.txt

echo outlier list - Done

#Use outlier list to filter samples in PLINK data 
plink2 \
--bfile $qcdir/1KG_merged.no_sib \
--remove $qcdir/1KG_merged.outliers.txt \
--make-bed \
--out $qcdir.1KG_merged.IBD