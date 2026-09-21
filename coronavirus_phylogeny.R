#################################################################
# Coronavirus Phylogenetic Analysis
# BIOL 3315: Evolution and Phylogenetics
# Langara College
#
# This script performs sequence alignment, phylogenetic analysis,
# model selection, bootstrap analysis, and visualization of
# coronavirus genome and spike protein sequences 
#
# Author: Akshita Sharma
# Term: FALL 2025
#
# Methods:
# - Multiple sequence alignment (ClustalW)
# - Substitution model selection (AIC)
# - NJ and ML phylogenetic inference
# - Bootstrapping and tree comparison (tanglegrams)
#
# Tools:
# R, ape, phangorn, msa, Biostrings, phytools, taxize
#
#################################################################

# ---------------------------
# Create results directory if it does not exist
# ---------------------------

if (!dir.exists("results")){
  dir.create("results")
}

# ---------------------------
# Load required packages
# ---------------------------

install.packages(c("seqinr", "ape", "phangorn", "phytools", "taxize"))
options(BioC_mirror = "https://bioconductor.statistik.tu-dortmund.de")
BiocManager::install("msa", ask = FALSE, update = FALSE)

# Sequence handling and alignment
# Installing the package msa (for performing multiple sequence alignment)
library(msa)        #Multiple sequence alignment
library(seqinr)     #Sequence format conversion
library(Biostrings) #FASTA sequence handling

# Phylogenetic analysis 
library(ape)        #Phylogenetic trees and utilities
library(phangorn)   #ML inference, model testing
library(phytools)   #Tree manipulation and comparison

# Taxonomic annotation
library(taxize)     #Retrieve organism metadata from NCBI

# ---------------------------
# Data acquisition
# ---------------------------

# Coronavirus genome nucleotide sequences (e.g., SARS, MERS)
genome_nt <- readDNAStringSet(
  filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/SARS_MERS_coronavirus.raw_sequence.fasta")
# Inspect
genome_nt

# Coronavirus spike protein amino acid sequences
spike_aa <- readAAStringSet("https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/spike.fa")
# Inspect
spike_aa

# Inspect sequence identifiers
names(genome_nt)
names(spike_aa)

# ---------------------------
# Multiple sequence alignment
# ---------------------------

# Align coronavirus genome nucleotide sequences
genome_nt_msa <- msa(genome_nt, method = "ClustalW")
# Inspect
genome_nt_msa

# Align coronavirus spike protein amino acid sequences
spike_aa_msa <- msa(spike_aa, method = "ClustalW")
# Inspect
spike_aa_msa

# ---------------------------
# Convert alignments for phylogenetic inference
# ---------------------------

# Convert nucleotide MSA objects to phyDat format
genome_nt_phy <- msaConvert(
  genome_nt_msa, 
  type = "phangorn::phyDat")

# Convert amino acid MSA to phyDat format
genome_aa_phy <- as.phyDat(
  msaConvert(
    spike_aa_msa, 
    "seqinr::alignment"),
  type = "AA")

# phangorn::phyDat stores aligned sequences in a format that is optimal 
# for distance calculations and likelihood-based phylogenetic inference

# ---------------------------
# Distance-based phylogenetic inference (Neighbor-Joining)
# ---------------------------

# Compute evolutionary distance matrices
# JC69 substitution model for nucleotide sequences
genome_nt_dist <- dist.ml(genome_nt_phy, model="JC69")
genome_nt_dist

# JTT model for amino acid sequences
spike_aa_dist <- dist.ml(genome_aa_phy, model="JTT")
spike_aa_dist

# Construct Neighbor-Joining trees
genome_nt_nj <- NJ(genome_nt_dist)
# Inspect NJ tree
genome_nt_nj

spike_aa_nj <- NJ(spike_aa_dist)
# Inspect NJ tree
spike_aa_nj

# ---------------------------
# Bootstrapping Neighbor-Joining trees
# ---------------------------

# Set seed for reproducibility
set.seed(123)

# Bootstrap NJ tree for genome nucleotide sequences
genome_nt_nj_bs <- bootstrap.phyDat(
  genome_nt_phy,
  FUN=function(x) NJ(dist.ml(x, model="JC69")),
  bs = 1000)

# Bootstrap NJ tree for spike protein amino acid sequences 
spike_aa_nj_bs <- bootstrap.phyDat(
  genome_aa_phy,
  FUN=function(x) NJ(dist.ml(x, model="JTT")),
  bs = 1000)

# ---------------------------
# Visualize NJ trees with bootstrap support
# ---------------------------

# Plot bootstrapped NJ tree for genome nucleotide sequences
pdf("results/genome_nt_nj_bs.pdf", width = 10, height =10)
plotBS(
  genome_nt_nj, 
  genome_nt_nj_bs, 
  p=50, 
  main = "Bootstrapped NJ Tree for Coronaviruses genome", 
  type = "phylogram"
)
dev.off()

# Plot bootstrapped NJ tree for spike protein sequences
pdf("results/spike_aa_nj_bs.pdf", width = 10, height =10)
plotBS(
  spike_aa_nj, 
  spike_aa_nj_bs, 
  p=50, 
  main = "Bootstrapped NJ Tree for spike protein", 
  type = "phylogram"
)
dev.off()

# Only bootstrap values greater than equal to 50% are displayed to 
# highlight moderately to strongly supported clades

#----------------------------
# Midpoint rooting of bootstrapped NJ trees
#----------------------------

# Rooting bootstrapped nucleotide NJ tree
genome_nt_nj_rooted <- midpoint(genome_nt_nj)

# Rooting bootstrapped amino acid NJ tree
spike_aa_nj_rooted <- midpoint(spike_aa_nj)

# Save rooted bootstrapped tree as a pdf
pdf("results/genome_nt_nj_bs_rooted.pdf", width = 10, height =10)
plotBS(genome_nt_nj_rooted, 
       genome_nt_nj_bs, 
       p=50, 
       main = "Bootstrapped And Rooted NJ Tree for Coronaviruses genome", 
       type = "phylogram"
)
dev.off()

pdf("results/spike_aa_nj_bs_rooted.pdf", width = 10, height =10)
plotBS(spike_aa_nj_rooted, 
       spike_aa_nj_bs, 
       p=50, 
       main = "Bootstrapped and Rooted NJ Tree for spike protein", 
       type = "phylogram"
)
dev.off()

# ---------------------------
# Substitution model selection
# ---------------------------

# Model testing for genome nucleotide sequences
genome_nt_models <- modelTest(genome_nt_phy, model="all")
# Inspect
genome_nt_models

# Model testing for spike protein sequences
spike_aa_models <- modelTest(genome_aa_phy, model="all")
# Inspect
spike_aa_models

# Identifying best-fitting models based on AIC
genome_nt_models$Model[genome_nt_models$AIC==min(genome_nt_models$AIC)]
spike_aa_models$Model[spike_aa_models$AIC==min(spike_aa_models$AIC)]

# ---------------------------
# Constructing ML trees and Bootstrapping
# ---------------------------

# Maximum likelihood tree for coronavirus nucleotide sequences
genome_nt_ml <- pml(
  genome_nt_nj, 
  data = genome_nt_phy
)

# Optimize model parameters
genome_nt_ml <- optim.pml(
  genome_nt_ml, 
  model = "GTR", 
  optGamma = TRUE, 
  optInv = TRUE, 
  k = 4)

# Root the ML genome tree
genome_nt_ml_rooted <- midpoint(genome_nt_ml$tree)

# Bootstrap ML genome tree
genome_nt_ml_bs <- bootstrap.pml(
  genome_nt_ml, 
  bs = 100, 
  optNni = TRUE
)


# Maximum likelihood tree for spike protein sequences
spike_aa_ml <- pml(
  spike_aa_nj, 
  data = genome_aa_phy,
  k=4,
  inv=0
)

# Optimize the model parameters
spike_aa_ml <- optim.pml(
  spike_aa_ml, 
  model="WAG", 
  optGamma = TRUE, 
  optInv = FALSE, 
  k = 4
)

# Root the ML spike protein tree
spike_aa_ml_rooted <- midpoint(spike_aa_ml$tree)

# Bootstrap ML spike protein tree
spike_aa_ml_bs <- bootstrap.pml(
  spike_aa_ml, 
  bs = 100, 
  optNni = TRUE
)


# Save bootstrapped maximum likelihood coronavirus genome tree as a pdf
pdf("results/genome_nt_ml_bs.pdf", width =10, height =10)
plotBS(genome_nt_ml_rooted, 
       genome_nt_ml_bs, 
       main = "Maximum Likelihood Bootstrapped tree of Coronaviruses genome", 
       type = "phylogram"
)
dev.off()


# Save bootstrapped maximum likelihood coronavirus genome tree as a pdf
pdf("results/spike_aa_ml_bs.pdf", width =10, height =10)
plotBS(spike_aa_ml_rooted, 
       spike_aa_ml_bs, 
       main = "Maximum Likelihood Bootstrapped tree of spike protein", 
       type = "phylogram"
)
dev.off()


# ---------------------------
# Helper function: relabel tree tips using NCBI common names
# ---------------------------

relabel_tree_with_common_names <- function(tree) {
  
  tip_labels <- tree$tip.label
  
  common_names <- tryCatch({
    df <- tax_name(
      query = tip_labels,
      get = "common",
      db = "ncbi"
    )
    df$common
  }, error = function(e) {
    message("NCBI lookup failed — using original tip labels")
    rep(NA, length(tip_labels))
  })
  
  common_names[is.na(common_names)] <- tip_labels[is.na(common_names)]
  tree$tip.label <- common_names
  
  return(tree)
}



# ---------------------------
# Relabel tree tips
# ---------------------------

# Genome trees
genome_nt_nj_labeled <- relabel_tree_with_common_names(genome_nt_nj)
genome_nt_ml_labeled <- relabel_tree_with_common_names(genome_nt_ml$tree)

# Spike protein trees
spike_aa_nj_labeled <- relabel_tree_with_common_names(spike_aa_nj)
spike_aa_ml_labeled <- relabel_tree_with_common_names(spike_aa_ml$tree)


# ---------------------------
# Generate tanglegrams
# ---------------------------

# NJ vs ML (Genome)
tangle_genome <- cophylo(
  genome_nt_nj_labeled,
  genome_nt_ml_labeled
)

# NJ vs ML (Spike protein)
tangle_spike <- cophylo(
  spike_aa_nj_labeled,
  spike_aa_ml_labeled
)

# ML Genome vs ML Spike protein
tangle_ml_compare <- cophylo(
  genome_nt_ml_labeled,
  spike_aa_ml_labeled
)

# Save each object as the pdf
pdf("results/tanglegram_genome_NJ_vs_ML.pdf", width = 10, height = 8)
plot(tangle_genome, link.type="curved", link.lwd=2, fsize=0.8, main="Genome Phylogeny: NJ vs ML")
dev.off()


# Save each object as the pdf
pdf("results/tanglegram_spike_NJ_vs_ML.pdf", width = 10, height = 8)
# mar sets the margins of the plot
plot(tangle_spike, link.type="curved", link.lwd=2, fsize=0.8, main="Spike Protein Phylogeny: NJ vs ML")
dev.off()


# Save each object as the pdf
pdf("results/tanglegram_genome_vs_spike_ML.pdf", width = 10, height = 8)
plot(tangle_ml_compare, link.type="curved", link.lwd=2, fsize=0.8, main="Genome vs Spike Protein: Maximum Likelihood Phylogenies")
dev.off()


