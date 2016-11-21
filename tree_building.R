library(seqinr)
library(Biostrings)
library(msa)
library(DECIPHER)
library(phangorn)
library(phytools)
library(ape)
library(phangorn)
library(dplyr)
library(tidyr)
library(beepr)


# load in my cleaned, unaligned sequences
spp_fasta <- read.fasta(file="phylo_data/all_species.fasta", seqtype = "DNA")
# GC content of my sequences
#allGC <- sapply(spp_fasta, GC)
#plot(allGC, xlab="Sequence", ylab="GC content", pch=19, col="tomato2")
#hist(allGC, xlab="Sequence GC Content", col="tomato", main=NA)
# normality check
#shapiro.test(allGC)

##--------------------------------------------------------------------##
#MULTIPLE SEQUENCE ALIGNMENT
spp_seq <- sapply(spp_fasta, getSequence)

spp_seq <- spp_seq[sapply(spp_seq, length) > 500]

spp_string <- sapply(spp_seq, c2s) #convert characters into string

spp_stringset <- DNAStringSet(spp_string) #convert into DNA string set

#ALIGN USING MUSCLE
spp_mus_al <- msaMuscle(spp_stringset)

#convert alignment back to string set
spp_mus_al_string <- as(spp_mus_al, "DNAStringSet")

#Browse alignment in the browser
#BrowseSeqs(spp_mus_al_string, highlight=1, htmlFile="spp_mus_al_string.html")

# trim 
spp_mus_al_string_cleaned <- sapply(spp_mus_al_string, c2s)

cleaned_strings <- lapply(spp_mus_al_string_cleaned, function(sequence) { 
                cleaned <- substr(sequence, 49, 731-(6 + 9))
                start_count <- attributes(gregexpr("^(-{1,})", cleaned)[[1]])$match.length
                if (start_count != -1) {
                  cleaned  <- gsub('^(-{1,})', replacement = paste(rep('N', start_count), collapse=''), cleaned)
                }
                end_count <- attributes(gregexpr("(-{1,})$", cleaned)[[1]])$match.length
                if (end_count != -1) {
                  cleaned  <- gsub('(-{1,})$', replacement = paste(rep('N', end_count), collapse=''), cleaned)
                }
                return(cleaned)
              })

for (i in seq_along(spp_mus_al_string_cleaned)) {
  spp_mus_al_string_cleaned[[i]] <- cleaned_strings[[i]]
}

#BrowseSeqs(as(spp_mus_al_string, 'DNAStringSet'), highlight=1, htmlFile='phylo_data/cleaned_alignment.html')

#write as FASTA
writeXStringSet(as(spp_mus_al_string_cleaned, 'DNAStringSet'), file="phylo_data/spp_mus_al_string.fasta") 


# Tree building!
aln <- read.phyDat(file="phylo_data/spp_mus_al_string.fasta", format="fasta", type="DNA")

seq_annots <- sapply(spp_fasta[sapply(spp_fasta, length) > 500], function(seq) {attributes(seq)$Annot})

seqid_lookup <- data_frame(id = names(seq_annots), annot = as.vector(seq_annots)) %>%
  extract(., annot, into = c('gi', 'species', 'gi2', 'species2'), regex = ">gi\\|(\\d*)\\|gb\\|.*\\| (\\w* \\w*)|>gi\\|(\\d*):\\d*-\\d* (\\w* \\w*)")
seqid_lookup[which(is.na(seqid_lookup$species)), c('gi', 'species')] <- 
seqid_lookup[which(!is.na(seqid_lookup$species2)), c('gi2', 'species2')]
seqid_lookup <- select(seqid_lookup, id, gi, species) %>% as_data_frame()

dm <- dist.ml(aln)

treeNJ <- NJ(dm)

fit <- pml(treeNJ, data = aln)

# With the default values pml will estimate a Jukes-Cantor model. The function
# update.pml allows to change parameters. We will change the model to
# the GTR + Γ(4) + I model and then optimize all the parameters.

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR

tips <- data_frame(id = fitGTR$tree$tip.label) %>% 
        left_join(., seqid_lookup) %>% 
        unite(col = gi_species, gi, species, sep = '-') %>% 
        select(gi_species)

fitGTR$tree$tip.label <- unlist(tips)

dev.new(height = 13, width = 8)
par(mar = c(1, 1, 1, 1))
plot(fitGTR$tree, cex = 0.3)

root(phy, outgroup, node = NULL, resolve.root = FALSE)

write.tree(fitGTR$tree, file="test_tree4.tre")


# randomly prune tree to single individual per species
random_tips <- 
  seqid_lookup %>% 
  filter(gi != '227937050') %>% 
  filter(!(species %in% c('Chaetodon ornatissimus', 'Chaetodon striatus', 
    'Chromis cyanea', 'Epinephelus hexagonatus', 'Grammistes sexlineatus', 
    'Mycteroperca microlepis', 'Parupeneus trifasciatus', 'Pseudupeneus maculatus'))) %>% 
  mutate(gi_species = paste(gi, species, sep = '-')) %>%
  group_by(species) %>% 
  sample_n(1)


tree1 <- drop.tip(fitGTR$tree, tip = unlist(anti_join(tips, random_tips, by = "gi_species")))

tree1$tip.label <- gsub('\\d*-', '', tree1$tip.label)
tree1 <- root(tree1, outgroup = 'Myxine glutinosa', node = NULL, resolve.root = TRUE)

chronos_tree1 <- chronos(tree1, model = "relaxed", lambda = 0)

plot(chronos_tree1)

chronos_tree1 <- drop.tip(chronos_tree1, "Myxine glutinosa")


tips

dev.copy2pdf(device = quartz, file = "jons_tree.pdf")

tips <- data_frame(id = treeNJ$tip.label) %>% 
        left_join(., seqid_lookup) %>% 
        unite(col = gi_species, gi, species, sep = '-') %>% 
        select(gi_species)

treeNJ$tip.label <- unlist(tips)


tree <- read.tree(text = '((((Acanthurus_nigricans:0.1091359933,Acanthurus_olivaceus:0.097197521):0.254431009,((Aphareus_furca:0.29883908,Lutjanus_bohar:0.47189922):0.03408358,(Caesio_teres:0.2462225709,Pterocaesio_tile:0.07438387):0.14789676):0.055139647):0.01510922,((((((Pseudanthias_dispar:0.63624986,Pseudanthias_olivaceus:0.042261708):0.862032437,(Cephalopholis_argus:0.26557206,Cephalopholis_urodeta:0.08370104):0.17164307):0.08437997,Variola_louti:0.297869):0.068573844,(Chromis_vanderbilti:0.49918127,Monotaxis_grandoculis:1.956863649):0):0,((Chlorurus_sordidus:0.16010243,(Scarus_frenatus:0.053014222,Scarus_rubroviolaceus:0.04210742):0.0462503305):0.360911637,(Paracirrhites_arcatus:0.6091793163,(Caranx_melampygus:0.36613484,Parupeneus_insularis:0.392413591):0.01419371):0):0):0.03287936,(Centropyge_flavissima:0.629484751,Chaetodon_ornatissimus:0.598817138):0):0):0.20922633,OUTGROUP_Myxine_glutinosa:3.275993);
')

dev.new()
plot(tree)

tips <- data_frame(id = treeUPGMA$tip.label) %>% 
        left_join(., seqid_lookup) %>% 
        unite(col = gi_species, gi, species, sep = '-') %>% 
        select(gi_species)

treeUPGMA$tip.label <- unlist(tips)



dev.new(height = 13, width = 8)
par(mar = c(1, 1, 1, 1))
plot(treeUPGMA, cex = 0.3)

treeNJ <- NJ(dm)

We can plot the trees treeUPGMA and treeNJ (figure 1) with the commands:
layout(matrix(c(1,2), 2, 1), height=c(1,2))
#par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="UPGMA", cex = 0.2)
plot(treeNJ, "unrooted", main="NJ")

# Aligning sequences
seqs <- readDNAStringSet("phylo_data/all_species.fasta")

# look at lengths of sequences
seq_df <- 
  data_frame(names = seqs@ranges@NAMES, length = seqs@ranges@width) %>% 
  extract(., names, into = c('gi', 'gb', 'species', 'desc', 'gi2', 'gb2', 'species2', 'desc2'), regex = "gi\\|(\\d*)\\|gb\\|(.*)\\| (\\w* \\w*) (.*)|gi\\|(\\d*):(\\d*-\\d*) (\\w* \\w*) (.*)")
seq_df[which(is.na(seq_df$gi)), c('gi', 'gb', 'species', 'desc')] <- 
seq_df[which(!is.na(seq_df$gi2)), c('gi2', 'gb2', 'species2', 'desc2')]

# eliminate any funky characters in names, but keep unique ids
unique_ids <- unite(seq_df, id, gi, species, sep = '|')
seqs@ranges@NAMES <- unique_ids[[1]]

#seq_df %>% 
#  filter(length > 600) %>%
#  group_by(species) %>% 
#  summarise(count = n(), min(length), max(length))

seqs <- seqs[seqs@ranges@width > 600]

alignment <- msaMuscle(seqs)

align_string <- as(alignment, "DNAStringSet")
writeXStringSet(align_string, filepath = 'test.fasta')

write.phylip(alignment, filepath = 'test.phy')

msaPrettyPrint(alignment, output='tex', file = 'phylo_data/alignment.tex', alFile='phylo_data/alignment.fasta', askForOverwrite=FALSE)


msaConvert(alignment, type=c("seqinr::alignment"))

alignments <- read.dna("phylo_data/alignment.fasta", format = "fasta")


alignments[1, 28:678]




test <- read.alignment("phylo_data/alignment.fasta", format = "fasta", forceToLower = TRUE)

x <- lapply(test$seq, function(sequence) {
            trimmed <- substr(sequence, 28, 678)
            replaced <- gsub('-', replacement = 'N', x = trimmed)})


test$seq[[1]]

spp_codes <- c('AC.NIGR', 'AC.OLIV', 'AP.FURC', 'LU.BOHA', 'CA.TERE', 
               'PT.TILE', 'PS.DISP', 'PS.OLIV', 'CE.ARGU', 'CE.UROD', 
               'VA.LOUT', 'CH.VAND', 'MO.GRAN', 'CH.SORD', 'SC.FREN', 
               'SC.RUBR', 'PA.ARCA', 'CA.MELA', 'PA.INSU', 'CE.FLAV')#, 
               #'CH.ORNA')

spp_list <- c("Acanthurus nigricans", "Acanthurus olivaceus", "Aphareus furca", "Lutjanus bohar", "Caesio teres", "Pterocaesio tile", "Pseudanthias dispar", "Pseudanthias olivaceus", "Cephalopholis argus", "Cephalopholis urodeta", "Variola louti", "Chromis vanderbilti", "Monotaxis grandoculis", "Chlorurus sordidus", "Scarus frenatus", "Scarus rubroviolaceus", "Paracirrhites arcatus", "Caranx melampygus", "Parupeneus insularis", "Centropyge flavissima"
)
