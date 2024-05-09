'''
Created on 8 May 2024

@author: sdantu
'''
dict_aa={"ALA":"A",
          "ARG":"R",
          "ASN":"N",
          "ASP":"D",
          "CYS":"C",
          "GLN":"Q",
          "GLU":"E",
          "GLY":"G",
          "HIS":"H",
          "HIP":"H",
          "HIE":"H",
          "ILE":"I",
          "LEU":"L",
          "LYS":"K",
          "MET":"M",
          "PHE":"F",
          "PRO":"P",
          "SER":"S",
          "THR":"T",
          "TYR":"Y",
          "TRP":"W",
          "VAL":"V",

         }
def get_aa_code(aa_name):
    if(aa_name in dict_aa):
        return dict_aa[aa_name]
    else:
        return "-"
def get_aa_list():
    list_aa=list(set(list(dict_aa.values())))
    list_aa.sort()
    return list_aa