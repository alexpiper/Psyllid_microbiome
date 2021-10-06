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

ps2 <- readRDS("ps2.rds")
tree <- read.tree("phytree.nwk")

phy_tree(ps2) <- tree

#Rename taxa - only keep first 30 characters
taxa_names(ps2) <- substr(paste0("SV", seq(ntaxa(ps2)),"-",tax_table(ps2)[,7]), 1,30)

## Create species merged table

# Merge species for beta diversity
ps.sppmerged <- ps2 %>%
    merge_samples(group = "psyllid_spp", fun=mean)

#This loses the sample metadata - Need to add it agian
sample_data(ps.sppmerged) <- sample_data(ps2) %>%
  as("matrix") %>%
  as.data.frame() %>%
  filter(!duplicated(psyllid_spp)) %>%
  magrittr::set_rownames(.$psyllid_spp)

ps3 <- ps.sppmerged

# Read in tree ------------------------------------------------------------

psyllid_tree <- read.tree(text=readLines("psyllid_beast_tree.nwk"))

# Match names with sample sheet
psyllid_tree$tip.label <- psyllid_tree$tip.label %>%
  str_squish() %>%
  str_replace_all(pattern="\\.", replacement=" ") %>%
  str_replace_all(pattern="Acizzia hakae", replacement="Acizzia hakeae") %>%
  str_replace_all(pattern="POLLENISLAND", replacement="POLLEN ISLAND") %>%
  str_replace_all(pattern="Ctenarytaina fuchsiae$", replacement="Ctenarytaina fuchsia A") %>%
  str_replace_all(pattern="Ctenarytaina fuchsiaeB", replacement="Ctenarytaina fuchsia B") %>%
  str_replace_all(pattern="Ctenarytaina fuchsiaeC", replacement="Ctenarytaina fuchsia C") %>%
  str_replace_all(pattern="Ctenarytaina clavata", replacement="Ctenarytaina clavata sp ") %>%
  str_replace_all(pattern="Ctenarytaina clavata sp $", replacement="Ctenarytaina clavata sp A") %>%
  str_replace_all(pattern="Ctenarytaina sp$", replacement="Ctenarytaina sp ") %>%
  str_replace_all(pattern="Ctenarytaina spA", replacement="Ctenarytaina sp A") %>%
  str_replace_all(pattern="Ctenarytaina spB", replacement="Ctenarytaina sp B") %>%
  str_replace_all(pattern="Ctenarytaina unknown", replacement="Ctenarytaina insularis") %>%  
  str_replace_all(pattern="Psylla apicalisA", replacement="Psylla frodobagginsi") %>%
  str_replace_all(pattern="Psylla apicalisB", replacement="Psylla apicalis") %>%
  str_replace_all(pattern="carmichaeliae", replacement="carmichaeliae ") %>%
  str_replace_all(pattern="Trioza sp", replacement="Trioza sp ") %>%
  str_replace_all(pattern="Trioza acutaB", replacement="Trioza acuta B") %>%
  str_replace_all(pattern="Trioza gourlay", replacement="Trioza gourlayi") %>%
  str_replace_all(pattern="BRENDAMAY", replacement="BRENDA MAY") %>%
  str_replace_all(pattern="PRICES", replacement="PRICES VALLEY") %>%  
  str_replace_all(pattern="Acizzia sp", replacement="Acizzia errabunda") %>% 
  str_replace_all(pattern="Trioza ", replacement="Powellia ") %>%
  str_replace_all(pattern="Triozid sp", replacement="Casuarinicola sp") %>%
  str_replace_all(pattern="Powellia adventicia", replacement="Trioza adventicia") %>%
  str_replace_all(pattern="Powellia curta", replacement="Trioza curta") %>%
  str_replace_all(pattern=" ", replacement="_") %>%
  trimws(which="right")

# Subset to only those in sample data
psyllid_tree$tip.label[!psyllid_tree$tip.label %in% sample_data(ps2)$psyllid_spp]
pruned.tree <- drop.tip(psyllid_tree, psyllid_tree$tip.label[!psyllid_tree$tip.label %in% sample_data(ps2)$psyllid_spp] )



# paco Psylla microbe -----------------------------------------------------

#Prepare co-occurance matrix
coocur <- ps3 %>%
  subset_samples(psyllid_genus == "Psylla") %>%
  filter_taxa(function(x) mean(x) > 0, TRUE) %>%
  otu_table %>%
  as.matrix() %>% 
  apply(2, function(x) ifelse(x > 0, 1, 0))

colnames(coocur) <- colnames(coocur) %>%
  str_replace_all(" |-", "_")
rownames(coocur) <- rownames(coocur) %>%
  str_replace_all(" |-", "_")

# H cophenetic distance
h_tree <- pruned.tree
h_tree$tip.label <- h_tree$tip.label %>%
  str_replace_all(" |-", "_")
h_tree <- drop.tip(h_tree, setdiff(h_tree$tip.label, rownames(coocur)))
h_dist <- sqrt(cophenetic(h_tree))

# P cophenetic distance
s_tree <- phy_tree(ps3)
s_tree$tip.label <- s_tree$tip.label %>%
  str_replace_all(" |-", "_")
s_tree <- drop.tip(s_tree, setdiff(s_tree$tip.label, colnames(coocur)))
s_dist <- sqrt(cophenetic(s_tree) /1e+6) #convert to Mya so integers are small enough for PACO

coocur <- coocur[h_tree$tip.label, s_tree$tip.label]
# prepare paco data
D <- prepare_paco_data(H=h_dist, P=s_dist, HP=coocur)
# Add pcord
D <- add_pcoord(D, correction='none')

# run paco
print("running paco")
paco_run <- PACo(D, nperm=1000, seed=909, method='quasiswap', symmetric=TRUE)
paco_run <- paco_links(paco_run,  .parallel = TRUE)

write_rds(paco_run, "paco_psylla_microbe.rds")

# Run Parafit
PF_run <- parafit(h_dist, s_dist, coocur, nperm=1000, test.links=TRUE, silent=TRUE)
write_rds(PF_run, "parafit_psylla_microbe.rds")

