'''
Created on 8 May 2024

@author: sdantu
'''
import pandas as pd

def save_data(filename,data):
    f=open(filename,'w')
    f.write(data)
    f.close()
def read_thermomutdb_as_df(filename):
    df=pd.read_json(filename)
    return df

def read_kmer_db(fname):
    kmer_data=[]
    kmer_data=open(fname,'r').readlines()
    return(kmer_data)