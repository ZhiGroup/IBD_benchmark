# IBD Detection Tool Benchmark Project

This open source project benchmarks most popular IBD detection tools with multiple measurements. 

For quick and simple demonstration precompiled executable versions and a set of test data are uploaded to ```IBD_Benchmark_ExperiencePackage``` directory.

The project does:
1. Simulate human genotype data (VCF file).
2. Add genotyping errors and phasing errors to VCF file.
3. Run IBD detection tools with the data, and collect the IBD results.
4. Evaluate the IBD results according to ground truth by multiple metrics.

Simulation code could be found in ```Simulation``` directory.

Error insertion code could be found in  ```Error_Insertion``` directory.

Evaluation software code could be found in  ```IBD_Benchmark``` directory.


