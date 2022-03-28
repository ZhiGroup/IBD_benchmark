## Create a VCF panel with switch errors:

A genotype/haplotype panel can be generated using msprime or any other simulator. stdpopsim uses msprime as its engine by default.
```sh
stdpopsim HomSap -c chr20 -o output_file -g HapMapII_GRCh37 -d OutOfAfrica_2T12 0 500000
```
The command line above generates 500000 haplotypes of European ancestry. You can convert the output of the simulator into VCF format using tskit:

```sh
tskit vcf output_file | gzip --stdout  > 500k_hapmap_chr20.vcf.gz
```
The sites with very low allele freuencies can be then filtered out using bcftools:

```sh
bcftools view -q 0.1:minor 500k_hapmap_chr20.vcf.gz
```

[SHAPEIT4](https://odelaneau.github.io/shapeit4/) or other phasing tools can be used to phase the genotype data.
