# SomaticAnalysisPublic


Scripts for various kinds of somatic analyses

- AmpDelsCancerGenes.py - From FACETS copy number data, identify and plot genes that are amplified or deleted in multiple samples from the same patient
- CategorizeDeletions.py - Find whole intron deletions in the IMPACT cohort (deletions that coincide with introns in any gene)
- CheckMicrohomology.py - check sequence microhomology around a deletion, a tell-tale sign of HRD
- ClonalDecomposition.py - Perform clonal decomposition and identify various mutational clusters from Cancer Cell Fraction (CCF) values for each mutation in each sample for multiple samples from the same patient
- ComputeCTDNAFraction.py - compute ctDNA fraction in cfDNA (proportion of cfDNA originating from tumor cells) using known clonal mutations and their copy number status in corresponding tissue sample
- MutationPlotsSeaborn.py - Make longitudinal cfDNA allele fraction plots with mutation annotations
- MutationsPlotsAltair.py - Make interactive longitudinal cfDNA allele fraction plots, showing a rich array of information like treatment regimens, clinical response, hotspot/OncoKB/AlphaMissense status of mutations, read depths and supporting reads in cfDNA and in tissue, in normal etc. to help with interpretation of the cfDNA data
- PlotTumorVolume.py - make plots showing correspondence between ctDNA fraction and tumor volume as measured by CT scans
- SingleCellAnalyzeAndPlot.py - Analyze and visualize single cell DNA seq data from MissionBio, showing biologically relevant results
- WholeGenomeDoubling.py - Detect whole genome doubling from copy number data

 
