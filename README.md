bioConvert
==========

Format converters for Bioinformatics:

* _vcf2tplink.py_:  Transform a VCF with only SNPs to a TPED and TFAM PLINK format. More details running `./vcf2tplink.py --help`
* _shapeit2bgl.R_:  Transform the [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html "SHAPEIT") output (extensions .haps and .sample assumed) to BGL format. Ussage: `Rscript shapeit2bgl.R fileroot inputfolder outputfolder`

--------------------------------------------------------------------------------------------------------
*Python scripst are written in Python 2.7. and needs [argparse-1.2.1](https://docs.python.org/dev/library/argparse.html)