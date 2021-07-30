          
            #DEPENDENCIES
import matplotlib as plt
import argparse  #module for parsing arguments, very convenient!
import math 
from typing import Union 

#Welcome to Bioinfo.py! John's Bioinfomatics Python Module! 
#Index
#line 10 - convert phred
#line 26   pythag  
#line 36   convert phred string      
#line 43   qual score
#line 53   init list
#line 61   populate list 
#line 75   get args (argparse)
#line 93   median calc
#line 110  kmerIZE
#line 128  knorm 
#line 158  fasta newline remover 
#line 175  fasta contigs parser 
#line 252  validate DNA seq 
#line 259  gc content 
#line 266  validate base seq 


def convert_phred(letter: str) -> str:   
    """Converts a single character into a phred score"""
    score = ord(letter) -33
    return score

def pythag(a: Union[int, float], b: Union[int, float]) -> float:   
    '''Computes the length of the hypotenuse of a right triangle given the lengths
    of the other two sides.'''
    c1= a**2 + b**2 
    c2=math.sqrt(c1)
    return float(c2)

def convert_phred_string(phred_score: str) -> str: 
    """This function converts a string of phred scores to their quality scores"""
    for char in range(0, len(phred_score)):
        print(char,': ',phred_score[char],' - ',convert_phred(phred_score[char]), sep="", end="\n") 
        return

def qual_score(phred_score: str) -> str: #
    """This function converts phred scores to their quality scores"""
    list1=[]
    for char in range(0, len(phred_score)):
        list1.append(convert_phred(phred_score[char]))
            #for item in list1: 
              #  list1.append())
    average=sum(list1)/len(list1)
    return(average)

def init_list(lst: list, value: float=0.0) -> list: #
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    lst = [value]*101
    return lst

def populate_list(file: str) -> str: #
    """This function initializes a list and populates it from a file"""
    lane1_list = []
    lane1_list1 = init_list(lane1_list)
    i=0
    with open(file, "r") as l1: 
        for line in l1:
            line=line.strip('\n')
            i+=1
            if i%4 == 0: 
                for y, x in enumerate(line):
                    lane1_list1[y]+=(convert_phred(x))
        return(lane1_list1, i)
     
def get_args():             
    """This function passes arguements to a python script to"""
    parser = argparse.ArgumentParser(description="A program to normalize kmerspec")
    parser.add_argument("-f", "--filename", help="Your filename", required=True)               
    parser.add_argument("-l", "--length", help="read length", required=True, type=int)
    parser.add_argument("-k", "--kmer", help="what is your kmer length?", type=int)
    parser.add_argument("-o", "--output_filename", help="what is your output file?", required=True)
    parser.add_argument("-c", "--coverage", help="what is your coverage?", type=int)
    return parser.parse_args()
	
    args = get_args()    
    filename = (args.filename)
    readlength = (args.length)
    kmersize = (args.kmer)
    output = (args.output_filename)
    coverage = (args.coverage)


def median_calc(lst: list) -> list:               
    """This function calculates the median given a list of integers"""
    length=len(lst) 
    lst.sort()   #sorts line line but doesnt return the variable 
    #sorted(list) sorting this way returns a full, sorted list 
    if (length%2==0)==True:
        prime_slice=int(length/2 -1)
        #print(prime_slice)
        first=lst[prime_slice]
        second=lst[prime_slice+1]
        median=((second+first)/2)
    else:
        prime_odd=length//2
        median=(lst[prime_odd])
    return(median)  


def kmerIZE(DNA: str) -> str:     
    """This function kmerizes data from a fasta file"""
    j=0
    k=kmersize
    kmerdict={}
    seq_range=range(0, (args.length-kmersize+1))
    for _ in seq_range:
        slice = DNA[j:k]
        if slice not in kmerdict:
            kmerdict[slice]=1
        else:
            kmerdict[slice]+=1
        j+=1
        k+=1
    # for _ in kmerdict.keys():
    #     print(str(_)+"\t"+str(kmerdict[_]))
    return(kmerdict)

def knorm(filename: str,output: str) -> str:
    """This function normalizes kspectra from fastq files"""
    knorm_dict={}
    with open(filename, "r") as fh, open(output, "w") as op:      #opening the fq file to read, output file to output to 
        while True: 
            header=fh.readline().strip()              #storing the 4 lines in a fq file in a variable 
            if header=="":                          #readlines fxn parses and ends @ a no space character: break before that
                break 
            seqline=fh.readline().strip()
            plus=fh.readline().strip()
            qscore=fh.readline().strip()                                 
            kmerdict=kmerIZE(seqline)                                    #kmerize line
            for kmer in kmerdict:                               
                if kmer not in knorm_dict:
                    knorm_dict[kmer]=1                                      #store it w/ a 1 count 
                else: 
                    knorm_dict[kmer]+=kmerdict[kmer]                         #increment it 
        #print("knormdict is:",knorm_dict)
            kmer_storage_list=[]                                      # init list 
            kmerdict=kmerIZE(seqline)
            for kmer in kmerdict:                                          #add values to list 
                kmer_storage_list.append(knorm_dict[kmer])
            kmer_storage_list.sort()
            mediancoverage=median_calc(kmer_storage_list)               #calcs median of a list 
            if mediancoverage < (args.coverage): 
                op.write(header+'\n'+seqline+'\n'+plus+'\n'+qscore+'\n')   #saves the desired info in an output file 
    return

#  Last edited on Friday July 16th, 2021 

def fasta_newline_remover(filename: str) -> str: 
    """This function removes newlines from fasta files"""
    with open(filename, "r") as fh:
        firstline = True			            #opening the fastq file to read it 
        for line in fh:                         
            if line[0] == '>':                  
                header = line.strip()           #removing end of line whitespace characters
                if firstline == False:
                    print('\n', end="")
                print(header)
                firstline = False
            else:
                seq = line.strip()
                print(seq, end="")
        print('\n', end="")
    return 

def fasta_contigs_parser(filename: str) -> str:
	"""This function parses a contigs.fasta file and returns statistics, distribution and a .pdf barplot of distribution, requires filename input and output"""
	with open(filename, "r") as fh:     
		#opening the fasta file to read it
		#innitilizing lists to hold contig lengths, physical lengths
		physlengths=[]
		lengths=[]
		coverages=[]
		weighted_coverages=[]
		for line in fh:
			#splitting and stripping contig lines 
			if line.startswith(">")==True:
				line=line.strip('\n')
				line=line.split('_')
			#assigning contig lengths and coverage values from split header, calculating physical length #ASSUMING kmer length 49 
				length=int(line[3])
				physlength=int(line[3]) + 49-1
				coverage=float(line[5])
				weighted_coverage= (int(physlength)*int(coverage)/int(length))*int(physlength)
			#appending variables to lists 
				physlengths.append(physlength)
				lengths.append(length)
				coverages.append(coverage)
				weighted_coverages.append(weighted_coverage)

	#assigning variable names for downstream calculations 
	total_contig_number=len(physlengths)
	genome_length=sum(physlengths)
	mean_contig_length=genome_length/total_contig_number
	print("mean contig length:",mean_contig_length)

	physlengths.sort(reverse=True)
	max_contig_length=physlengths[0]
	print("max contig length:",max_contig_length)

	sumphys=0
	n50=None
	half_genome_length=sum(physlengths)/2
	for physlength in physlengths:
		sumphys=sumphys+physlength
		if sumphys >= half_genome_length:
			n50=physlength
			print("N50 value:",n50)
			break

	mean_coverage_depth=sum(weighted_coverages)/genome_length
	print('mean coverage depth:', mean_coverage_depth)
	#visual partiton 
	print("__________________________")

	bucket_dict={} 
	frequency_dict={}
#creating dictionaries to store contig bins and contig frequencies 
	for length in lengths:
		bin_number=(length//100)*100
		if bin_number not in bucket_dict:
			bucket_dict[bin_number]=[]
		bucket_dict[bin_number].append(length)
	for key in bucket_dict:
		frequency_dict[key]=len(bucket_dict[key])

	sorted_frequency_dict={}

	for key in sorted(frequency_dict):
			sorted_frequency_dict[key]=frequency_dict[key]
#Plotting data in a barplot 
	plt.bar(sorted_frequency_dict.keys(), sorted_frequency_dict.values(), width=100.0)
	plt.xlabel('contig bins')
	plt.yscale('linear')
	plt.ylabel('number of contigs')
	plt.savefig(output.pdf)
#returning sorted data as well as other stats into standard output 
	print("The Distribution")
	for k,v in sorted_frequency_dict.items():
		print(k,v)
	return 

def validate_DNA_seq(seq: str) -> str:                                                  
    '''This function takes a string. Returns True if string is composed
    of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
    #here is a comment
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")

def gc_content(DNA: str) -> str:
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

def validate_base_seq(seq, RNAflag=False) -> str:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

