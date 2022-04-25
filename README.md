# RSV GWAS Analysis
Repositry containing scripts to perform GWAS and eQTL analysis using a Respiratory syncytial Virus, (RSV), case-control study

# Study Overview
This study aims to investigate potential genetic associations with severity of RSV disease, and characterise these associations using publically available data, in addition to running expression-quantative trait loci (eQTL) analysis.\
Participants are from a case-control study of infants presenting to hospital with RSV infection\
Part of the Rescue consortium https://resc-eu.org/ 

# Methods:
- Genotyping performed on Axiom chip, Microarray used for gene-expression analysis\
- Pre-processing of genotyping data performed using bash, perl and R v4.1
- Imputation performed using TOPMED imputation server
- https://imputation.biodatacatalyst.nhlbi.nih.gov/#!pages/home
- Quality control (QC) of genotyping data and association analyses were run in PLINK 1.9/2.0 and R v4.13
- https://www.cog-genomics.org/plink/2.0/
- SNP nexus was used to annotate SNPs of interest
- https://www.snp-nexus.org/index.html

# Git instructions

- To work on scripts please start a new Rstudio Project > version control > exisiting Repo https://github.com/mari-lynne/RSV_GWAS
- Check push/pull to origin is set up and working correctly, then in the R terminal use:

git branch user_dev #user_dev will be your new branch for developing code\
git checkout user_dev

- Within your branch you can pull from main and edit/add scripts (commit regularly).
- When happy with changes merge your branch to main 

# Directories

- **Scripts/pre-process** Contains scripts for formating and QC of genotyping data for upload to TOPMED
- **Scripts/post-impute folder** Contains scripts for data reformatting and Post-imputation QC steps
- **Scripts/post-impute/Assoc** Association analysis and SNP annotation scripts
- **Scripts/post-impute/eQTLeQTL** Contains scripts for performing eQTL analysis (currently copied from typhoid Fc analysis)
- **Scripts/meta** Contains details of metadata and scripts used to format data/IDs
