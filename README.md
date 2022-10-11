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

### Formats:
#### IBD output from an IBD decetion software
We prebuilt 5 format pasers: FastSMC,hap-IBD,iLash,RaPID, and TPBWT. The choice of input type could be defined in ```/IBD_Benchmark_ExperiencePackage/Config.txt```.
For other formats of IBD output you may add additional parser in ```Loader.cs```, or you may also convert your IBD output into one of the supported format. Details of each format could be found in their GitHub page. 

#### Ground truth IBD format
An example of ground truth file could be found in ```IBD_benchmark/IBD_Benchmark_ExperiencePackage/arr.gt.txt.gz``` as:
```
#individual_1_id,individual_1_haplotype_id,individual_2_id,individual_2_haplotype_id,chromosome_id,true_ibd_physical_position_start,true_ibd_physical_position_end,genetic_length

1710,1,44,1,20,69403,352433,1.012187

2197,0,44,1,20,69403,352433,1.012187

947,1,59,0,20,69403,352433,1.012187
```
Same as the input format, you may add additional parser in ```Loader.cs```, or convert your ground truth file into this format.

## Output:
7 evaluation results binned by cM ranges:

### Each row: 
Accuracy,Length Accuracy,Length Accuracy,Length Discepancy,Recall,Power,Accumulative Recall,and Accumulative Power.

### Each column: 
Values for each bin. The binning could be modified in the ```/IBD_Benchmark_ExperiencePackage/Config.txt``` file by adjusting starting bin, ending bin and bin length. 

e.g, ```minBin=2; maxBin=7; binLen=1;``` This means the eveluation will only collect IBD length >=2cM, and the last bin will have results for all IBD that >= 7cM.

```
Bin_0	Bin_1	Bin_2	Bin_3	Bin_4	Bin_5	

Accuracy
ToolName	0.998434048366469	0.99970583368572	0.999488444322646	0.99975333873639	0.99988419894621	0.999829646394588	

Length Accuracy
ToolName	0.979041095412064	0.98614326906595	0.989468088410164	0.991318139999737	0.993142297654614	0.995252082425287	

Length Discepancy
ToolName	1.14525992959265	1.72141823296979	2.27453118780375	3.05301635367978	3.62253222715984	7.63100423669796	

Recall
ToolName	0.278648090397311	0.781735470693817	0.892241197257458	0.87602218496066	0.86575977604759	0.742540508230885	

Power
ToolName	0.248235241129757	0.64423443971452	0.738102138000221	0.752536732567245	0.750271981641258	0.67780075119114	

Accumulative Recall
ToolName	0.278648090397311	0.781747849229436	0.898703346520435	0.910589449245453	0.931108389467238	0.978142517891472	

Accumulative Power
ToolName	0.26727813509913	0.647393065974722	0.745047223562128	0.783646771365684	0.81164894727716	0.890121550075723	

```

