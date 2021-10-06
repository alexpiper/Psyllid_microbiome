#Set required packages
.cran_packages <- c("tidyverse", 
                    "patchwork", 
                    "vegan", 
                    "seqinr",
                    "ape", 
                    "sp",
                    "data.table", 
                    "RColorBrewer",
                    "ggtree", 
                    "castor", 
                    "picante", 
                    "phylosignal", 
                    "adephylo",
                    "paco",
                    "phytools")
.bioc_packages <- c("dada2",
                    "phyloseq", 
                    "DECIPHER",
                    "Biostrings",
                    "ShortRead", 
                    "philr",
                    "ALDEx2")

sapply(c(.cran_packages,.bioc_packages), require, character.only = TRUE)

# Github packages
library(taxreturn)
library(seqateurs)
library(speedyseq)

# read in ps --------------------------------------------------------------
print("conducting paco on rawdists")

ps2 <- readRDS("ps2.rds")
tree <- read.tree("phytree.nwk")

phy_tree(ps2) <- tree

#Rename taxa - only keep first 30 characters
taxa_names(ps2) <- substr(paste0("SV", seq(ntaxa(ps2)),"-",tax_table(ps2)[,7]), 1,30)

## Create species merged table

ps.sppmerged <- ps2 %>%
    merge_samples(group = "psyllid_spp", fun=mean)

#This loses the sample metadata - Need to add it agian
sample_data(ps.sppmerged) <- sample_data(ps2) %>%
  as("matrix") %>%
  as.data.frame() %>%
  filter(!duplicated(psyllid_spp)) %>%
  magrittr::set_rownames(.$psyllid_spp)

seqs <- refseq(ps2)
tree <- phy_tree(ps2)

#make new phyloseq object
ps3 <- phyloseq(tax_table(ps.sppmerged),
               sample_data(ps.sppmerged),
               otu_table(otu_table(ps.sppmerged), taxa_are_rows = FALSE),
               refseq(seqs),
               phy_tree(tree))
			   
# PACO Raw distances psyllid microbe ---------------------------------------------------

#Prepare co-occurance matrix
coocur <- ps3 %>%
  otu_table %>%
  as.matrix() %>% 
  apply(2, function(x) ifelse(x > 0, 1, 0))

colnames(coocur) <- colnames(coocur) %>%
  str_replace_all(" |-", "_")
rownames(coocur) <- rownames(coocur) %>%
  str_replace_all(" |-", "_")

# H distance
set.seed(666)
coi.dist <- read_csv("COI.csv", col_names = TRUE) %>%
  column_to_rownames("SampleID") 
coi.dist[upper.tri(coi.dist)] <- 0
coi.dist<- coi.dist + t(coi.dist)
coi.dist[is.na(coi.dist)] <- 0
h_dist <- coi.dist %>%
  magrittr::set_colnames(rownames(.)) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    sample_data(ps3) %>%
      as("matrix") %>%
      as.data.frame(stringsAsFactors =FALSE) %>%
      dplyr::select(SampleID, psyllid_spp)
  ) %>%
  mutate(psyllid_spp = na_if(psyllid_spp, "NA")) %>%
  filter(!is.na(psyllid_spp)) %>% 
  group_by(psyllid_spp) %>%
  dplyr::sample_n(1) %>%
  ungroup() %>%
  dplyr::select(.$SampleID, psyllid_spp) %>%
  magrittr::set_colnames(c(.$psyllid_spp, "psyllid_spp"))  %>%
  column_to_rownames("psyllid_spp") %>%
  as.matrix()

setdiff(rownames(h_dist), rownames(coocur))
# S distance
seqs <- readDNAStringSet("aligned_asv.fasta.gz")
names(seqs) <- names(seqs) %>%
  str_remove("\\ \\[.*$")

seqs <- seqs[names(seqs) %in% as.character(refseq(ps2))]

#Rename taxa - only keep first 30 characters
names(seqs) <- substr(paste0("SV", seq(ntaxa(ps2)),"-",tax_table(ps2)[,7]), 1,30) %>%
  str_replace_all(" |-", "_")

s_dist <- DECIPHER::DistanceMatrix(seqs)

coocur <- coocur[colnames(h_dist), colnames(s_dist)]

# prepare paco data
D <- prepare_paco_data(H=h_dist, P=s_dist, HP=coocur)
# Add pcord
D <- add_pcoord(D, correction='lingoes')

p_host <- ggplot(D$H_PCo[,c('Axis.1', 'Axis.2')] %>% as.data.frame, aes(Axis.1, Axis.2)) +
  geom_point() +
  theme_bw() +
  ggtitle("H PCA")

p_para <- ggplot(D$P_PCo[,c('Axis.1', 'Axis.2')] %>% as.data.frame, aes(Axis.1, Axis.2))  +
  geom_point() +
  theme_bw()+
  ggtitle("P PCA")

plot(p_host + p_para)

# run paco
print("running paco")
paco_run <- PACo(D, nperm=1000, seed=909, method='quasiswap', symmetric=TRUE) #Symetric - is one meant to track the evolution of another?
# Get interaction-specific cophylogenetic contributions using jacknifing
paco_run <- paco_links(paco_run,  .parallel = TRUE)

write_rds(paco_run, "paco_psyllid_microbe_rawdists.rds")

# Run Parafit
PF_run <- parafit(h_dist, s_dist, coocur, nperm=1000, test.links=TRUE, silent=TRUE)
write_rds(PF_run, "parafit_psyllid_microbe_rawdists.rds")