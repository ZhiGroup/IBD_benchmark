# IBD Detection Tool Benchmark Project

This open source project benchmarks most popular IBD detection tools with multiple measurements. 

For quick and simple demonstration precompiled executable versions and a set of test data are uploaded to ```IBD_Benchmark_ExperiencePackage``` directory.

The project does:
1. Simulate human genotype data (VCF file).
2. Create array data by downsampling sequencing data from step 1.
3. Add genotyping errors and phasing errors to VCF files.
4. Run IBD detection tools with the data, and collect the IBD results.
5. Evaluate the IBD results according to ground truth by multiple metrics.

Data simulation code could be found in ```Simulation``` directory.

Downsample method could be found in ```Downsample``` directory.

Error insertion code could be found in  ```Error_Insertion``` directory.

Evaluation software code could be found in  ```IBD_Benchmark``` directory.

