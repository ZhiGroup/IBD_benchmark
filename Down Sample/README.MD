## This method down-samples sequencing data into array data.

Given a target number of markers, the program will try to take as much marker as possible and as even as possible in term of genetic location.


The program does:
1. Scan and calculate the genetic location and minor allele count in the sequencing data.
2. Create a list to collect by using a certain window size.
3. Extract the sites from the sequencing data by the list of collection.


Note:
In Step 2, using smaller window size will produce more even markers, but less number of sites.
bigger window size will give more uneven markers, but more sites in array data.
Try multiple rounds to find ideal window size.
We used 5 as window size in our project.
