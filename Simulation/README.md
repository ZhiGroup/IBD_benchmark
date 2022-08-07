# IBD_benchmark
An example command of running simulation program:
```
python -u ./msprime_simulation.py 20 ./genetic_map_GRCh37_chr20.txt 1.38e-8 56 4000 17197 5000 1.0 ./output/
```


Modify line 280 to simulate mix or different population.

We used ```individuals_from_populations_to_sample_dictionary = {"CEU": 1334,"CHB":1333,"YRI":1333}``` for mixed populaiton.
