plink --bfile /home/mari/RSV/plink/topmed/rsv_remap --exclude /home/mari/RSV/plink/topmed/Exclude-rsv_remap-HRC.txt --make-bed --out /home/mari/RSV/plink/topmed/TEMP1
plink --bfile /home/mari/RSV/plink/topmed/TEMP1 --update-map /home/mari/RSV/plink/topmed/Chromosome-rsv_remap-HRC.txt --update-chr --make-bed --out /home/mari/RSV/plink/topmed/TEMP2
plink --bfile /home/mari/RSV/plink/topmed/TEMP2 --update-map /home/mari/RSV/plink/topmed/Position-rsv_remap-HRC.txt --make-bed --out /home/mari/RSV/plink/topmed/TEMP3
plink --bfile /home/mari/RSV/plink/topmed/TEMP3 --flip /home/mari/RSV/plink/topmed/Strand-Flip-rsv_remap-HRC.txt --make-bed --out /home/mari/RSV/plink/topmed/TEMP4
plink --bfile /home/mari/RSV/plink/topmed/TEMP4 --a2-allele /home/mari/RSV/plink/topmed/Force-Allele1-rsv_remap-HRC.txt --make-bed --out /home/mari/RSV/plink/topmed/rsv_remap-updated
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 1 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr1
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 1 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr1
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 2 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr2
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 2 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr2
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 3 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr3
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 3 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr3
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 4 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr4
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 4 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr4
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 5 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr5
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 5 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr5
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 6 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr6
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 6 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr6
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 7 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr7
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 7 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr7
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 8 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr8
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 8 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr8
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 9 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr9
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 9 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr9
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 10 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr10
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 10 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr10
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 11 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr11
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 11 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr11
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 12 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr12
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 12 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr12
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 13 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr13
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 13 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr13
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 14 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr14
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 14 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr14
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 15 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr15
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 15 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr15
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 16 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr16
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 16 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr16
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 17 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr17
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 17 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr17
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 18 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr18
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 18 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr18
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 19 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr19
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 19 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr19
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 20 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr20
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 20 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr20
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 21 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr21
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 21 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr21
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --make-bed --chr 22 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr22
plink --bfile /home/mari/RSV/plink/topmed/rsv_remap-updated --real-ref-alleles --recode vcf --chr 22 --out /home/mari/RSV/plink/topmed/rsv_remap-updated-chr22
rm TEMP*
