#!/usr/bin/env python
#Sbatch directeives 

from os import write
import Bioinfo
#import matplotlib.pyplot as plt 
import gzip 

#make a 2 dictionaries with all of the index file names fw and reverse and make a directory 
#with all of them that can be referenced. 

# def reverse(string):
#     for i in string:
#         string_new=string[::-1]
#     return string_new

# def complement(n):
#     """Creates the complement of a strand of DNA"""
#     rev_comp=""
#     for letter in n:
#         if letter == 'A':
#             rev_comp+='T'
#         elif letter == 'C': 
#             rev_comp+='G'
#         elif letter == 'T': 
#             rev_comp+='A'
#         elif letter == 'G': 
#             rev_comp+='C'
#         elif letter =='N':
#             rev_comp=+'N'
#     return(rev_comp)


index_forward_dict={}
index_reverse_dict={}

with open("indexes.txt", "r") as fh: 
    fh.readline()
    for line in fh:
        line=line.strip('\n')
        index_forward_dict[line.split('\t')[3]]=line.split('\t')[4]

index_forward_dict["UN"]=" "
index_forward_dict["SW"]=" "

for key in index_forward_dict.keys():
    index_reverse_dict[key]=Bioinfo.complement(Bioinfo.reverse(index_forward_dict[key]))

temp = "R2_"
index_reverse_dict = {temp + str(key): val for key, val in index_reverse_dict.items()}

# print(index_forward_dict)
# print(index_reverse_dict)

# for key in index_forward_dict.keys():
#     open(f"data_dir/{key}.fq", 'w')

# for key in index_reverse_dict.keys():
#     open(f"data_dir/{key}.fq", 'w')


#make a dictionary with indexes as keys, and their file handles as values 
full_index_dict={}
for key, val in index_forward_dict.items():
    full_index_dict[val]=open(f"data_dir/{key}.fq", 'w')        

for key, val in index_reverse_dict.items():
    full_index_dict[val]=open(f"data_dir/{key}.fq", 'w') 

#print(full_index_dict)
#full_index_dict['CTTAGGAC'].write('blah')
unR1=open("data_dir/UN.fq", 'w')
unR2=open("data_dir/R2_UN.fq", 'w')

swR1=open("data_dir/SW.fq", 'w')
swR2=open("data_dir/R2_SW.fq", 'w')

bdR1=open("data_dir/BD.fq", 'w')
bdR2=open("data_dir/R2_BD.fq", 'w')


with gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz","rt") as R1, \
    gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz","rt") as R2, \
        gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz","rt") as R3, \
            gzip.open("/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz","rt") as R4: 


# with open("./TEST-input_FASTQ/test_input_R1.fq","r") as R1, \
#     open("./TEST-input_FASTQ/test_input_R2.fq","r") as R2, \
#         open("./TEST-input_FASTQ/test_input_R3.fq","r") as R3, \
#             open("./TEST-input_FASTQ/test_input_R4.fq", "r") as R4:
    i=0 
    unknown_count=0
    swapped_count=0
    matched_index_count=0
    bad_quality_count=0
    while True:
        headR1=R1.readline().strip()
        if headR1=="":
            break 
        seqR1=R1.readline().strip() 
        plusR1=R1.readline().strip() 
        qscoreR1=R1.readline().strip() 
        #print(headR1,'\n',seqR1,'\n',plusR1,'\n',qscoreR1,'\n')
        headR2=R2.readline().strip() 
        seqR2=R2.readline().strip() 
        plusR2=R2.readline().strip() 
        qscoreR2=R2.readline().strip() 
        #print(headR2,'\n',seqR2,'\n',plusR2,'\n',qscoreR2,'\n')
        headR3=R3.readline().strip() 
        seqR3=R3.readline().strip() 
        plusR3=R3.readline().strip() 
        qscoreR3=R3.readline().strip() 
        #print(headR3,'\n',seqR3,'\n',plusR3,'\n',qscoreR3,'\n')
        headR4=R4.readline().strip() 
        seqR4=R4.readline().strip() 
        plusR4=R4.readline().strip() 
        qscoreR4=R4.readline().strip() 
        #print(headR4,'\n',seqR4,'\n',plusR4,'\n',qscoreR4,'\n')
        rev_seqR3=Bioinfo.complement(Bioinfo.reverse(seqR3))
        # print(seqR2, seqR3, rev_seqR3)
        # #print(index_forward_dict.keys(), index_reverse_dict.keys())
        # print("seqR2 in index_forward_dict.values()?", seqR2 in index_forward_dict.values())
        # print("N not in seq2?", "N" not in seqR2)
        # print("N not in seq3?", "N" not in seqR3)
        # print("rev_seqR3 in index_reverse_dict.values()?", rev_seqR3 in index_reverse_dict.values())
        # print("passes all tests", all(("N" not in seqR2 , "N" not in seqR3 , seqR2 in index_forward_dict.values() , rev_seqR3 in index_reverse_dict.values())))
        # qscore mean 
        seqR1_converted_qscore=Bioinfo.qual_score(qscoreR1)
        seqR2_converted_qscore=Bioinfo.qual_score(qscoreR4)
        if all(("N" not in seqR2 , "N" not in seqR3 , seqR2 in full_index_dict , rev_seqR3 in full_index_dict)):
            # print("passed first test")
            if (seqR1_converted_qscore <= 30) or (seqR2_converted_qscore  <= 30):
                bdR1.write(headR1+'\n'+seqR1+'\n'+plusR1+'\n'+qscoreR1+'\n')
                bdR2.write(headR4+'\n'+seqR4+'\n'+plusR4+'\n'+qscoreR4+'\n')
                bad_quality_count+=1
            elif (seqR2 != rev_seqR3) == True:
                #write to swapped reads file 
                swR1.write(headR1+'\n'+seqR1+'\n'+plusR1+'\n'+qscoreR1+'\n')
                swR2.write(headR4+'\n'+seqR4+'\n'+plusR4+'\n'+qscoreR4+'\n')
                #print("swapped")
                swapped_count+=1
            elif seqR2 == rev_seqR3: #the reverse comp, yeah? 
                full_index_dict[seqR2].write(headR1+" "+seqR2+" "+seqR3+'\n'+seqR1+'\n'+plusR1+'\n'+qscoreR1+'\n')
                full_index_dict[rev_seqR3].write(headR4+" "+seqR2+" "+seqR3+'\n'+seqR4+'\n'+plusR4+'\n'+qscoreR4+'\n')
                matched_index_count+=1
                #print("matched")
        else:       
            #print("in else")
            unR1.write(headR1+'\n'+seqR1+'\n'+plusR1+'\n'+qscoreR1+'\n')
            unR2.write(headR4+'\n'+seqR4+'\n'+plusR4+'\n'+qscoreR4+'\n')
            unknown_count+=1 
            #print("unknown")


    #print(index_forward_dict, index_reverse_dict)
    Total_read_count=unknown_count+swapped_count+matched_index_count+bad_quality_count
    print(unknown_count,"unknown reads")
    print(swapped_count,"index swapped reads")
    print(matched_index_count,"matched index reads")
    print(bad_quality_count, "bad quality reads")
    print(Total_read_count,"total reads")
    print((swapped_count/(Total_read_count))*100, "percent swapping occured") # fix this 
    
    #print(str(seqR3),index_reverse_dict)

#run some stats here****

unR1.close()
unR2.close()
swR1.close()
swR2.close()
bdR1.close()
bdR2.close()
for key in full_index_dict:
    full_index_dict[key].close()
