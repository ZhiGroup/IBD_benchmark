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


## Inputs:
1. IBD output from an IBD detection software.
2. VCF file that used as the input to the IBD detection software in 1.
3. Ground truth IBD file.

Details instructions could be found and configured in ```/IBD_Benchmark_ExperiencePackage/Config.txt```
The program will read the configure file to execute.

## Input Formats:
### IBD output from an IBD decetion software
We prebuilt 5 format pasers: FastSMC,hap-IBD,iLash,RaPID, and TPBWT. The choice of input type could be defined in ```/IBD_Benchmark_ExperiencePackage/Config.txt```.
For other formats of IBD output you may add additional parser in ```Loader.cs```, or you may also convert your IBD output into one of the supported format. Details of each format could be found in their GitHub page. 

### Ground truth IBD format
An example of ground truth file could be found in ```IBD_benchmark/IBD_Benchmark_ExperiencePackage/arr.gt.txt.gz``` as:

#individual_1_id,individual_1_haplotype_id,individual_2_id,individual_2_haplotype_id,chromosome_id,true_ibd_physical_position_start,true_ibd_physical_position_end,genetic_length

1710,1,44,1,20,69403,352433,1.012187

2197,0,44,1,20,69403,352433,1.012187

947,1,59,0,20,69403,352433,1.012187

Same as the input format, you may add additional parser in ```Loader.cs```, or convert your ground truth file into this format.


