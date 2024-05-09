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
datafol="/Users/sdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D Sequence Evolution/data"
db_themo_fname="%s/experimental/thermomut/thermomutdb.json"%(datafol)
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
    fname="TMut_Out.txt"
    fdata=open(fname,'r')
    for line in fdata:
        line=line.strip().split(',')
        list_kmer_data.append(line)
    '''
    for i in range(0,10):
        print(list_kmer_data[i])
        process_kmer(list_kmer_data[i])
    '''
    out=[]
    
    pool=Pool(2)
    for i in tqdm.tqdm(pool.imap_unordered(process_kmer,list_kmer_data[:50]), desc="Compare Kmers @ ",total=len(list_kmer_data[:50])):
        out.append(i)
    #r = pool.map(process_kmer, tqdm.tqdm(list_kmer_data[:100]))
    pool.close()
    pool.join()
    print(out)
    
    #print('Processing %12d kmers'%(len(list_kmer_data)))
    #exit()
    #data_out+="%15s%3s%4s%6d%6.3f\n"%(wt_kmer,"WT","NAT",wt_lev_d[0],wt_lev_d[1])
    #fileio.save_data("TMut_Out.txt", data_out)

def process_kmer(kmer_data):
    global dict_all_kmers
    kmer_stats=comp_kmers.compare_kmer(kmer_data[0],dict_all_kmers[kmer_data[0][0]])
    return(kmer_data+kmer_stats)
    #return kmer_stats
main()
    