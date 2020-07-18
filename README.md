# CYP2D6 Simulator

Simulate CYP2D6 genotypes and phenotypes based on existing star alleles.

Three input files are required, all within "data".  You can also specify the number of samples to generate and a prefix for the output label file.

```
python bin/cyp2d6_simulator.py -f data/gnomad.cyp2d6.pops.no_indel.vcf  -S data/star_definitions.tsv -a data/star_scores -n <COUNT>  --prefix <PREFIX> > <PREFIX>.vcf
```
