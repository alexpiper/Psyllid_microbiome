# Psyllid_microbiome

This repository hosts the R based reproducible workflow that performed the microbiome analyses presented for the manuscript: _Martoni, F., Bulman, S. R., Piper, A. M., Pitman, A., Taylor, G. S., & Armstrong, K. F. (2023). Insect phylogeny structures the bacterial communities in the microbiome of psyllids (Hemiptera: Psylloidea) in Aotearoa New Zealand. PLoS One, 18(5), e0285587._

The reproducible workflow to conduct the analyses is seperated into 2 components, firstly the [bioinformatic](https://alexpiper.github.io/Psyllid_microbiome/Bioinformatics.html) processing of sequencing reads, and then the [statistical](https://alexpiper.github.io/Psyllid_microbiome/Statistics.html) analyses.

The input sequencing data are not included in the repository for size reasons, and are instead available from the SRA under accession: PENDING. However RDS files holding intermediate data objects such as the OTU and taxonomy tables suitable for performing the analyses are contained inside the output directory.

You can run these analyses on your own machine by (1) cloning the repository, (2) obtaining the raw sequencing data in fastq format from SRA, (3) Seperating the fastq files by run in the relevant subdirectory in the data directory, (4) installing required R libraries, and (5) pressing Run in the Rmarkdown file. Even without the sequencing data, the statistical analyses can be replciated using the stored rds files in the data directory.
