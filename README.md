# Antigenic variation by switching inter-chromosomal interactions with an RNA splicing locus in trypanosomes
Joana Faria, Vanessa Luzak, Laura S. M. Mueller, Benedikt G. Brink, Sebastian Hutchinson, Lucy Glover, David Horn, T. Nicolai Siegel
bioRxiv 2020.01.27.921452; doi: https://doi.org/10.1101/2020.01.27.921452
This repo contains the HiC analysis pipeline for this preprint. The raw data can be found at ENA under accession number PRJEB35632. Processed data and additional information can be found at Zenodo https://doi.org/10.5281/zenodo.3628213.

## Dependencies

The main dependency is the [mHiC](https://github.com/yezhengSTAT/mHiC) pipeline, which has been included as a submodule. Other dependencies are (modify the Snakefile if these are not in your PATH):

- [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
- g++ compiler
- [GEM-library](https://sourceforge.net/projects/gemlibrary/files/gem-library/) (to produce the mappability track required by mHiC)
- R with the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package
- [BEDOPS](https://bedops.readthedocs.io/en/latest/index.html) `wig2bed`
- [HiC-Pro Utilities](https://github.com/nservant/HiC-Pro/blob/master/doc/UTILS.md) (the `digest_genome.py` script for genome digestion is the only dependency)
- [HiC sunt draconis (HiCSD)](https://github.com/foerstner-lab/HiCsuntdracones) (>=0.2.0)


## Usage

Please edit the config.yaml to match your data and organism. For further information about snakemake, please refer to the official documentation. This pipeline has been adapted from the mHiC pipeline by the Keles Research Group, University of Wisconsin, so please check out their [Readme](https://github.com/yezhengSTAT/mHiC) for more details about mHiC and the parameters in this pipeline.
