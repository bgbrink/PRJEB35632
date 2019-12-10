# PRJEB35632
Further details to follow.


## Dependencies

The main dependency is the [mHiC](https://github.com/yezhengSTAT/mHiC) pipeline, which has been included as a submodule. Other dependencies are (modify the Snakefile if these are not in your PATH):

- [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
- g++ compiler
- [GEM-library](https://sourceforge.net/projects/gemlibrary/files/gem-library/) (we use it to produce the mappability track required by mHiC ICEing)
- R with the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package
- [BEDOPS](https://bedops.readthedocs.io/en/latest/index.html) [wig2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/wig2bed.html)
- [HiC-Pro Utilities](https://github.com/nservant/HiC-Pro/blob/master/doc/UTILS.md) (the `digest_genome.py` script for genome digestion is the only dependency)
- [HiC sunt draconis (HiCSD)](https://github.com/foerstner-lab/HiCsuntdracones)


## Usage

Please edit the config.yaml to match your data and organism. For further information about snakemake, please refer to the official documentation. This pipeline has been adapted from the mHiC pipeline by the Keles Research Group, University of Wisconsin, so please check out their [Readme](https://github.com/yezhengSTAT/mHiC) for more details about mHiC and the parameters in this pipeline.
