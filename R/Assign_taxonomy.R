
seqtab.nochim <- readRDS("output/rds/seqtab_final.rds")

# Assign Kingdom:Genus taxonomy using RDP classifier
tax <- assignTaxonomy(seqtab.nochim, "reference/merged_rdp_genus.fa.gz", multithread=TRUE, minBoot=60, outputBootstraps=FALSE)
colnames(tax) <- c("Root", "Phylum", "Class", "Order", "Family", "Genus")

##add species to taxtable using exact matching
tax_plus <- addSpecies(tax, "reference/merged_rdp_species.fa.gz", allowMultiple=TRUE)

test <- tax_plus
tax_plus <- propagate_tax(tax_plus,from="Family")

#add Genus_SPP
#for(col in seq(7,ncol(tax_plus))) { 
# propagate <- is.na(tax_plus[,col]) & !is.na(tax_plus[,col-1])
#  tax_plus[propagate,col:ncol(tax_plus)] <-  "spp."
#}

##join genus and species name in species rank column
#sptrue <- !is.na(tax_plus[,7])
#tax_plus[sptrue,7] <- paste(tax_plus[sptrue,6],tax_plus[sptrue,7], sep=" ")

#Check Output
taxa.print <- tax_plus # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Write taxonomy table to disk
saveRDS(tax_plus, "output/rds/tax_RDP_final.rds") 

##Plot assignment at each ranK!
