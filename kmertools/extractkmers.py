'''
Created on 8 May 2024

@author: sdantu
'''
from sklearn.neighbors import KDTree
import MDAnalysis as md
import kmerio.hashmaps as hmaps
import kmerio.fileIO as fileio

class Kmer(object):
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self._data_fol="/Users/sdantu/Library/CloudStorage/OneDrive-BrunelUniversityLondon/Work/Research/Manuscript/Me/3D Sequence Evolution/data/kmer"
        fol_natural_set=""
        self.dict_kmer={}
    
    def get_kmers(self,filename_pdb):
        self.process_pdb(filename_pdb)
        return(self.dict_kmer)
    
    def process_pdb(self,filename_pdb):
        pdb_universe = md.Universe(filename_pdb)
        ca_atoms_selection = pdb_universe.select_atoms("name CA and protein and not (resname MSE or resname HYP)")
        reslist=[]
        for res in ca_atoms_selection.resnames:
            reslist.append(hmaps.dict_aa[res])
        #print(ca_atoms_selection.resnums)
        if(len(ca_atoms_selection)>0):
            indices, distances = self._nearest_neighbours(reslist, ca_atoms_selection.positions)
            _kmer_list = [''.join([reslist[i] for i in ind]) for ind in indices]
            count=0
            for resid in ca_atoms_selection.resnums:
                self.dict_kmer[resid]=_kmer_list[count]
                count+=1
        
    def _nearest_neighbours(self,_residues, coordinates, search_radius=15):
        """ 
        Returns 2 lists, each containing lists of indices of residues and distances.
        For every residue, the list will be a sorted list of all residues within the search radius.
        default=15 angstroms
        """
        tree = KDTree(coordinates)
        indices, distances = tree.query_radius(coordinates, r=search_radius, return_distance=True, sort_results=True)

        return indices, distances