#!/usr/bin/env python3.9
'''
Created on 8 May 2024

@author: sdantu
'''
from multiprocessing import Pool
import tqdm
import kmerio.fileIO as fileio
from kmertools.extractkmers import Kmer
import kmertools.comparekmers as comp_kmers 
import pandas as pd
import re,timeit,os
from kmertools.comparekmers import dict_all_kmers
kmer_length=12
obj_kmer=Kmer()
#datafol="/Users/sdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D Sequence Evolution/data"
num_procs=56
datafol="/mnt/disk_c/scdantu/k-mers"
db_themo_fname="%s/thermomut/thermomutdb.json"%(datafol)
df_thermomut=pd.DataFrame()
cols_df=["PDB_wild","dtm","effect","mutation_type","mutation_code","ddg"]
data_out=""
dict_kmer_counts={}
dict_all_kmers={}
list_kmer_data=[]
def main():
    get_all_kmers()
    
    start=timeit.default_timer()
    
    get_kmers_list()
    process_kmers_list()
    end=timeit.default_timer()
    print("took : %12.3f"%(end-start))

def get_all_kmers():
    global dict_all_kmers,dict_kmer_counts
    dict_all_kmers=comp_kmers.get_all_kmers()
    
    total_num=0
    for key,value in dict_all_kmers.items():
        #print(key,len(value))
        num=len(value)
        dict_kmer_counts[key]=num
        total_num+=num
    print("Total No. of kmers : %12d"%(total_num))
    
def get_kmer_counts(kmer_aa):
    if(kmer_aa in dict_kmer_counts):
        return dict_kmer_counts[kmer_aa]
    else:
        return -1
def get_kmers_list():
    global list_kmer_data
    fname="/mnt/disk_c/scdantu/k-mers/thermomut/TMut_Out.txt"
    fdata=open(fname,'r')
    for line in fdata:
        line=line.strip().split(',')
        list_kmer_data.append(line)
    '''
    for i in range(0,10):
        print(list_kmer_data[i])
        process_kmer(list_kmer_data[i])
    '''
def process_kmers_list():
    global data_out
    kmer_results=[]
    pool=Pool(num_procs)
    test_len=96
    #for i in tqdm.tqdm(pool.imap_unordered(process_kmer,list_kmer_data[:test_len]), desc="Compare Kmers @ ",total=len(list_kmer_data[:test_len])):
    for i in tqdm.tqdm(pool.imap_unordered(process_kmer,list_kmer_data), desc="Compare Kmers @ ",total=len(list_kmer_data)):
       kmer_results.append(i)
    pool.close()
    pool.join()
    data_out+="%s,%s,%s,%s,%s/n"%("kmer","type","effect","levdist","seqratio")
    for result in kmer_results:
        data_out+="%s,%s,%s,%d,%.3f\n"%(result[0],result[1],result[2],result[4],result[5])
    fileio.save_data("/mnt/disk_c/scdantu/k-mers/thermomut/TMut_Kmer_Comparison.csv", data_out)

def process_kmer(kmer_data):
    global dict_all_kmers
    kmer_stats=comp_kmers.compare_kmer(kmer_data[0],dict_all_kmers[kmer_data[0][0]])
    return(kmer_data+kmer_stats)
    #return kmer_stats
main()
    
