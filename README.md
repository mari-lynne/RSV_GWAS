# RSV GWAS Analysis
Scripts for processing RSV data for GWAS and eQTL analysis 

- Pre-processing performed using bash, perl and R
- Imputation performed using TOPMED imputation server
- Scripts/pre-impute

- Post-imputation QC steps in Scripts/post-impute folder

- Association analysis and SNP annotation can be found in Scripts/post-impute/Assoc
- eQTL analysis can be found in Scripts/post-impute/Assoc (currently copied from typhoid Fc analysis)
- Details of metadata and data formating are in Scripts/meta

# Git instructions

- To work on scripts please start a new Rstudio Project > version control > exisiting Repo https://github.com/mari-lynne/RSV_GWAS
- Check push/pull to origin is set up and working correctly, then in the R terminal use:

git branch user_dev #user_dev will be your new branch for developing code\
git checkout user_dev

- Within your branch you can pull from main and edit/add scripts (commit regularaly).
- When happy with changes merge your branch to main 
