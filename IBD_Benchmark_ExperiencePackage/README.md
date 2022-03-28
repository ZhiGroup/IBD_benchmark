## IBD Benchmark Experience Package
A quick demostration for array data with 0.1% genotyping error.

Download all files in a directory.

Unzip all the Gzip files(.gz):

    May use command: ```gzip -d *.gz``` , if Gzip(available both in windows and Linux) is installed.
  
    Or anyother unzip software.
  
To Execute, run:

  ```IBD_BM_u18 Config.txt``` or ```IBD_BM_u18 Config.txt``` for ubuntu 18
  
  ```IBD_BM_u21 Config.txt``` or ```IBD_BM_u21 Config.txt``` for ubuntu 21
  
  ```IBD_BM.exe Config.txt``` or ```IBD_BM.exe Config.txt``` for windows 10 tested version
  

Results will be saved in ```IBD_BM_Result.txt``` by default.  

For configured executions, please read instructions in ```Config.txt```.


## List of Files:

```Config.txt``` the configure file contains paths and setting, as e input to the program. If this execute the program without given the configuration file name, the program will try to find Config.txt in the same directory.

```F_s9.arr.e0.001_2.1.1.FastSMC.ibd.gz``` Result from FastSMC for array data with 0.1% genotyping error.

```H_s9.arr.e0.001_2cM.ibd.gz``` Result from hap-IBD for array data with 0.1% genotyping error.

```IBD_BM.exe``` Windows compiled executable file, tested in Windows 10.

```IBD_BM_u18``` Ubuntu 18 compiled executable file.

```IBD_BM_u21``` Ubuntu 21 compiled executable file.

```I_s9.arr.e0.001_2.match.gz``` Rsult from iLash for array data with 0.1% genotyping error.

```R_s9.arr.e0.001.vcf_2w_3.results.max.gz``` Rsult from RaPID for array data with 0.1% genotyping error.

```T_s9.arr.e0.001.vcf_2cM.csv.gz``` Rsult from TPBWT for array data with 0.1% genotyping error.

```arr.gt.txt.gz``` Ground truth file for array data.

```genetic_map_GRCh37_chr20.txt``` Genetic map used for IBD detection tools and benchmarking result calculation.

```s9.arr.e0.001.vcf.gz``` Array data with 0.1% genotyping error as the input to IBD detection tools, also needed for benchmarking result calculation module.
