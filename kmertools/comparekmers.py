#!/usr/bin/python3.9
'''
Created on 9 May 2024

@author: sdantu
'''
import Levenshtein as levenshtein
import kmerio.fileIO as fileio
import kmerio.hashmaps as hmaps
data_fol="/Users/sdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D Sequence Evolution/data/kmer"
data_fol="/mnt/disk_c/scdantu/k-mers/kmerdb"
fname=""
kmer_data=[]
dict_all_kmers={}
def get_all_kmers():
    global dict_all_kmers
    dict_all_kmers={}
    for aa in hmaps.get_aa_list():
        dbname="%s/%s-kmers.txt"%(data_fol,aa)
        kmer_natural_set=fileio.read_kmer_db(dbname)
        dict_all_kmers[aa]=kmer_natural_set
    return dict_all_kmers

def compare_kmer(target_kmer,kmer_natural_set):
    #dbname="%s/%s-kmers.txt"%(data_fol,target_kmer[0])
    
    #kmer_natural_set=fileio.read_kmer_db(dbname)
    
    #print(len(kmer_natural_set))
    lev_d=[]
    for kmer in kmer_natural_set:
        kmer=kmer.strip()
        #if(kmer==target_kmer):
            #print("Kmer match found")
        lev_d.append([levenshtein.distance(target_kmer,kmer),levenshtein.seqratio(target_kmer,kmer)])
        #out+="%15s %4d\n"%(k,distance(kmer,k))
    return(min(lev_d))
    
