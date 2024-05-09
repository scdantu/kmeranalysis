'''
Created on 8 May 2024

@author: sdantu
'''
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

def main():
    get_all_kmers()
    start=timeit.default_timer()
    process_thermomut_db()
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
def process_thermomut_db():
    global df_thermomut
    df_thermomut=fileio.read_thermomutdb_as_df(db_themo_fname)
    df_sub=df_thermomut[["PDB_wild","dtm","effect","mutation_type","mutation_code","ddg"]]
    for index,row in df_sub.iterrows():
        if(index>-1):
            print("Processing row number %d"%(index))
            pdb_fname   =   row["PDB_wild"]
            effect      =   row["effect"]
            mut_code    =   row["mutation_code"]
            dict_kmers={}
            
            if(check_pdb_fname(pdb_fname)==True):
                dict_kmers=get_kmers_for_pdb(pdb_fname)
                if(len(dict_kmers)>1):
                    process_mutations(mut_code,dict_kmers,effect)
        if(index==-1):
            break
    fileio.save_data("TMut_Out.txt", data_out)
def check_pdb_fname(pdb_fname):
    is_pdb=False
    pdb_fpath="%s/experimental/thermomut/pdbs/%s.pdb"%(datafol,pdb_fname)
    if(pdb_fname is not None) and(os.path.exists(pdb_fpath)==True):
        if(len(pdb_fname)==4):
            is_pdb=True
    elif(pdb_fname is None):
        print(pdb_fname)
    return is_pdb
def process_mutations(mut_code,dict_kmers,effect):
    global data_out
    mut_code=mut_code.split(',')
    
    print(mut_code)
    for mutation in mut_code:
        if(len(mutation)>=3 and (len(mutation)<=6)):
            mut_split=re.split('(\d+)',mutation)
            if(len(mut_split)==3):
                
                loc=int(mut_split[1])
                wt_v=mut_split[0].strip()
                mt_v=mut_split[2].strip()
                eout="STA"
                if(len(wt_v)<2 and len(mt_v)<2 and len(mt_v)>0 and len(wt_v)>0) and (loc in dict_kmers.keys()):
                    
                    wt_kmer=dict_kmers[loc][:kmer_length]
                    mt_kmer="%s%s"%(mt_v,wt_kmer[1:])
                    
                    data_out+="%s,%s,%s,%d\n"%(wt_kmer,"WT",eout,get_kmer_counts(wt_v))
                    
                    if(effect=="destabilizing"):
                        eout="DES"
                    data_out+="%s,%s,%s,%d\n"%(mt_kmer,"MT",eout,get_kmer_counts(mt_v))
                    '''
                    wt_lev_d=comp_kmers.compare_kmer(wt_kmer)
                    mt_lev_d=comp_kmers.compare_kmer(mt_kmer)
                    eout="STA"
                    if(effect=="destabilizing"):
                        eout="DES"
                    data_out+="%15s%3s%4s%6d%6.3f\n"%(wt_kmer,"WT","NAT",wt_lev_d[0],wt_lev_d[1])
                    data_out+="%15s%3s%4s%6d%6.3f\n"%(mt_kmer,"MT",eout,mt_lev_d[0],mt_lev_d[1])
                    print(data_out)
                    '''
    
    #exit()
def get_kmers_for_pdb(pdbfname):    
    #pdbfname="%s/test/4LYZ.pdb"%(datafol)
    pdb_fpath="%s/experimental/thermomut/pdbs/%s.pdb"%(datafol,pdbfname)
    dict_kmer=obj_kmer.get_kmers(pdb_fpath)
    return(dict_kmer)
    
main()
    