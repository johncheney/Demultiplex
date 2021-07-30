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
Output will be 52 separate files, 48 of which will be sorted by the 24 indexes, and whether the file was read 1 or read 2. 2 files will be for index swapped read 1 and read 2 records. 2 files will be for the unknown adapter parts 

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).


4. Pseudocode

make all 50 files with adapter names from the list (plus 2 swapped files)  

initalize a known_index.dict that holds all of the adapter seqs

gzip open all 4 fq files    
    store 4 lines of code in memory for each file in variables 
        R1-1 headR1
        R1-2 seqR1
        R1-3 plusR1
        R1-4 qscoreR1

        R2-1 headR2
        R2-2 seqR2
        R2-3 plusR2
        R2-4 qscoreR2

        R3-1 headR3
        R3-2 seqR3
        R3-3 plusR3
        R3-4 qscoreR3

        R4-1 headR4
        R4-2 seqR4
        R4-3 plusR4
        R4-4 qscoreR4
    compare R2 and R3 adapter sequences to known_index.dict (rev comp?)
        if True: 
            compare R2 and R3 adapter sequences to one another 
        if different:
            move R1 and R4 to swapped file (count this)
        if same:  
            read the adapter and
                move R1 and R4 to corresponding adapter file (count this)
        if indexes have N in them:  
            move to unknown (count this)


5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

def file_creator(int, list):
"""makes files for each adapter sequence with nested R1 and R4 files withn each """
Input: 24, [list of adapter names]
Expected output: 48, [50 different adapter files] 
return

def fastq_parser_sort (filename, filename)
"""parses fq index file, compares index sequences """
Input: 2 fq index files 
Expected output: sorted into R1 and R4 files 
return 

def reverse_complment(string)
# YOU HAVE A ROSALIND PYTHON CODE THAT DOES THIS
"""takes a sequence in string and gives the reverse complement"""
Input: [AGTCN.+]
Expected output: [TCAGN.+]

def counter(int)
"""counts the times a file has been opened (VERY computationally heavy and a record added to it, saves number into a text file after the end of the fq input files. """
Input: int 
Expected output: int 