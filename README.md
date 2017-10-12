![supermerge_logo](https://cloud.githubusercontent.com/assets/9893806/24439486/7b1b26fc-141c-11e7-9a23-20b5b644f64f.png)

# SUPERmerge

## About

`SUPERmerge` is a ChIP-seq read pileup analysis and annotation algorithm for investigating alignment (`BAM`) files of diffuse histone modification ChIP-seq datasets with broad chromatin domains at a single base pair resolution level.  `SUPERmerge` allows flexible regulation of a variety of read pileup parameters, thereby revealing how read islands aggregate into areas of coverage across the genome and what annotation features they map to within individual biological replicates. `SUPERmerge` is especially useful for investigating low sample size ChIP-seq experiments in which epigenetic histone modifications (e.g., H3K9me1, H3K27me3) result in inherently broad peaks with a diffuse range of signal enrichment spanning multiple consecutive genomic loci and annotated features.

## Logo

`SUPERmerge` evaluates read coverage islands from histone modification ChIP-seq data with broad histone marks, which can often be wide in diameter and low in signal.  Just like having a cell phone on an island in the middle of the ocean.  

## Installation

`gcc -o supermerge supermerge.c`

## Usage

`./supermerge –d 20 –i 500 –g gencode.v19.annotation.gtf myFile.bam`

See manuscript for details: https://doi.org/10.1101/121897.  

## Funding

`SUPERmerge` is supported by the United States Department of Defense (DoD) through the National Defense Science and Engineering
Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under
and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering
Graduate (NDSEG) Fellowship, 32 CFR 168a.
