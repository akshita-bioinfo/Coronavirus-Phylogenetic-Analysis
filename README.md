# Coronavirus Phylogenetic Analysis

## Overview

This project investigates evolutionary relationships among selected coronaviruses using genome nucleotide sequences and spike protein amino acid sequences.

Two complementing tree-building techniques were used in R to conduct phylogenetic analyzes:

Neighbor-Joining (NJ), a distance-based method
Maximum Likelihood (ML), a model-based method

The generated phylogenies were compared to examine similarities and differences between inference methods and between genome-level and spike-protein analyses.

## Objectives

The main objectives of this project were to:

Perform multiple sequence alignment of coronavirus genome and spike protein sequences.
Construct Neighbor-Joining phylogenetic trees.
Select nucleotide and amino acid substitution models using Akaike Information Criterion (AIC).
Construct Maximum Likelihood phylogenetic trees using the selected models.
Assess phylogenetic support using bootstrap resampling.
Root trees using midpoint rooting.
Compare phylogenies using tanglegrams.
Compare genome-based and spike-protein-based evolutionary relationships.

## Data

The analysis used nine coronavirus sequences representing SARS-CoV-2, SARS-CoV, MERS-CoV, and related animal coronaviruses

Two sequence datasets were analyzed:
1. Genome nucleotide sequences
2. Spike protein amino acid sequences

During the study, the sequence files were obtained from the course GitHub repository and supplied as course materials. This repository does not distribute the original course data files.

## Methods

### Multiple Sequence Alignment

Genome nucleotide sequences and spike protein amino acid sequences were aligned separately using ClustalW through the msa R package.

The resulting alignments were converted into phyDat objects for downstream phylogenetic analysis.

### Neighbor-Joining Phylogenies

Neighbor-Joining trees were constructed from evolutionary distance matrices.

Genome nucleotide sequences: JC69 distance model
Spike protein amino acid sequences: JTT distance model

Bootstrap analysis was performed with 1,000 replicates for each NJ tree.

### Substitution Model Selection

Candidate evolutionary models were evaluated using modelTest() from the phangorn package.

The Akaike Information Criterion (AIC) was used to identify the selected substitution model for each dataset.

The selected models were:
1. GTR + Gamma (4 categories) + Invariant sites for Genome nucleotide sequences
2. WAG + Gamma (4 categories) for Spike protein amino acid sequences

### Maximum Likelihood Phylogenies

Maximum Likelihood trees were constructed using the selected substitution models.

Genome: GTR + Γ(4) + I
Spike protein: WAG + Γ(4)

Bootstrap analysis was performed with 100 replicates for each ML tree.

### Midpoint Rooting

The NJ and ML trees were midpoint-rooted to provide rooted representations of the inferred phylogenies.

### Tree Comparisons

Tanglegrams were generated using phytools to compare:

1. Genome NJ vs. Genome ML
2. Spike protein NJ vs. Spike protein ML
3. Genome ML vs. Spike protein ML

These comparisons were used to examine whether different datasets and phylogenetic inference methods produced similar or different evolutionary relationships.

## Results

The analysis generated nine PDF figures containing phylogenetic trees and tree comparisons.

### Phylogenetic Trees

Phylogenetic Trees
Bootstrapped genome NJ tree
Bootstrapped spike protein NJ tree
Rooted genome NJ tree
Rooted spike protein NJ tree
Bootstrapped genome ML tree
Bootstrapped spike protein ML tree

### Tanglegrams

Genome NJ vs. ML
Spike protein NJ vs. ML
Genome ML vs. spike protein ML

The tanglegrams provide visual comparisons of relationships among the same nine taxa across the different phylogenetic analyses

## Tools and R Packages

The analysis was performed using R and the following packages:

msa — multiple sequence alignment
Biostrings — biological sequence manipulation
seqinr — sequence analysis
ape — phylogenetic analysis
phangorn — phylogenetic inference and model selection
phytools — phylogenetic visualization and tree comparison
taxize — taxonomic name lookup

## Repository Structure

```
Coronavirus-Phylogenetic-Analysis/ 
│ 
├── README.md 
├── coronavirus_phylogeny.R 
├── results/ 
│ ├── genome_nt_nj_bs.pdf 
│ ├── spike_aa_nj_bs.pdf 
│ ├── genome_nt_nj_bs_rooted.pdf 
│ ├── spike_aa_nj_bs_rooted.pdf 
│ ├── genome_nt_ml_bs.pdf 
│ ├── spike_aa_ml_bs.pdf 
│ ├── tanglegram_genome_NJ_vs_ML.pdf 
│ ├── tanglegram_spike_NJ_vs_ML.pdf 
│ └── tanglegram_genome_vs_spike_ML.pdf 
│ 
└── .gitignore
```

## Reproducibility

The analysis is implemented in the coronavirus_phylogeny.R script.

The script performs the workflow from sequence retrieval through multiple sequence alignment, phylogenetic inference, bootstrap analysis, and figure generation.

The results/ directory is created automatically by the script if it does not already exist.

To reproduce the analysis, an R environment with the required packages installed is needed.


## Project Context

This repository contains an individual lab work completed for BIOL 3315: Evolution and Phylogenetics as part of the undergraduate Bioinformatics program at Langara College.

The sequence datasets and assignment framework were provided as course resources. I implemented the analysis workflow in R, including sequence alignment, evolutionary model selection using AIC, Neighbor-Joining and Maximum Likelihood phylogenetic inference, bootstrap analysis, midpoint rooting, and tanglegram-based tree comparisons.

This project provided hands-on experience with reproducible sequence analysis and phylogenetic workflows in R.
