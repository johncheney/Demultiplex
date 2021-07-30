# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1|
| 1294_S1_L008_R3_001.fastq.gz | Index 2|
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
![R1_histogram](https://user-images.githubusercontent.com/71104613/127625815-8115ed02-cba1-492b-86d4-9677a66c0e95.png)
![R2_histogram](https://user-images.githubusercontent.com/71104613/127625832-230988a3-1b5b-4292-b97b-2d49ec131b21.png)
![R3_histogram](https://user-images.githubusercontent.com/71104613/127625844-3cd50cef-9edf-4bc1-b91c-673a63185b23.png)
![R4_histogram](https://user-images.githubusercontent.com/71104613/127625856-3924f40a-a223-4ccd-baa8-789b2a99697f.png)
    3. I think that a quality score cutoff of 25 is the most appropriate; qscore of 20 correspsonds to 1 error in 100 for basecall accuracy and 30 being 1 in 1000. 25 seems sufficently high enough as to remove poor quality data but also retain the most useable data as well
    4. I got 7304664 total index sequences with undetermined base calls. 
    
    My code was:

    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n "2~4p" | grep -E "[N]" | wc -l 

    zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n "2~4p" | grep -E "[N]" | wc -l 

    
## Part 2
1. Define the problem
Index swapping can lead to drawing incorrect conclusion from the sequencing data we receive. We need to look at all of the indexes for each read, see if they are matched, and then sort them accordingly by read number and index and whether or not they've index swapped. 

2. Describe output
Output will be 50 separate files, each sorted by the 24 indexes, and whether the file was read 1 or read 2. 2 files will be for index swapped read 1 and read 2 records. 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).


4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
