#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
from pathlib import Path
import shutil
from Bio.PDB import *
from scipy.spatial.distance import *
import argparse
import itertools
import numpy as np
from collections import Counter

polar_amks = ['ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'HIS', 'LYS', 'SER', 'THR', 'TYR']

def arg_parser():
    parser = argparse. ArgumentParser()
    parser.add_argument("pdb", type=str)
    parser.add_argument('-f', '--from_file', action='store_true')
    return parser.parse_args()

def get_structure():
    args = arg_parser()
    print(args)
    parser = PDBParser()
    return parser.get_structure(args.pdb, args.pdb)
    if not args.from_file:
        file_name = args.pdb + ".pdb"
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(sys.argv[1])
        shutil.copy(file_name, 'pdbs')
        Path(file_name).unlink(missing_ok = True)

    return parser.get_structure(args.pdb, "pdbs/"+file_name)

def get_info(structure):
        
    print("name of the pdb:", sys.argv[1])
    print("number of models:", len(structure))
    chains = 0
    residues = 0
    atoms = 0
    for model in structure:
        chains += len(model)
        for chain in model:
            residues += len(chain)
            for residue in chain:
                atoms += len(residue)
    print("number of chains:", chains)
    print("number of residues:", residues)
    print("number of atoms:", atoms)


def get_widht(structure):
    atoms = [a.get_coord() for a in structure.get_atoms()]
    return (max(pdist(atoms, 'euclidean')))


def get_nearest_substructures(hetatm, distance:float, type_of_substr :str, atom_list):
    ns = NeighborSearch(atom_list)
    return ns.search(hetatm.get_coord(), distance, level = type_of_substr)
                     
                     
def get_nearest_atoms(hetatm, distance:float):
    return get_nearest_substructures(hetatm, distance, 'A', Selection.unfold_entities(structure, 'A'))


def get_nearest_residues(hetatm, distance:float, structure):
    return get_nearest_substructures(hetatm, distance, 'R', Selection.unfold_entities(structure, 'A'))


def get_diameter(structure):
    atoms = structure.get_atoms()
    #print(len(list(atoms)))
    diam = 0
    i = 0
    for atom_pair in itertools.combinations(atoms, 2):
        if i % 100000 == 0:
            print(i)
        i += 1
        distance = np.linalg.norm(atom_pair[0] - atom_pair[1])
        diam = distance if distance > diam else diam
    return diam

def get_buried_ratio(structure):
    buried =  len(list(i[1]  for i in HSExposure.ExposureCN(structure[0], radius=10, offset=0) if i[1] >=20))
    exposed =  len(list(i[1]  for i in HSExposure.ExposureCN(structure[0], radius=10, offset=0) if i[1] <20))
    return {"buried_ratio":round(buried/(buried+exposed),3),"exposed_ratio":round(exposed/(buried+exposed),3)}


def get_buried_exposed_by_amk(structure):
    buried = Counter(list(i[0].resname  for i in HSExposure.ExposureCN(structure[0], radius=10, offset=0) if i[1]>=20))
    exposed = Counter(list(i[0].resname  for i in HSExposure.ExposureCN(structure[0], radius=10, offset=0) if i[1]<20))
    return buried,exposed
    
def get_histogram(structure):
    buried, exposed = get_buried_exposed_by_amk(structure)
    amas = dict.fromkeys(set(exposed.keys()).union(set(buried.keys())),0)
    for ama in amas.copy():
        if ama in buried.keys() and ama in exposed.keys():
            amas[ama] = round(buried[ama]/(buried[ama] + exposed[ama]),3)
        elif ama in buried.keys():
            amas[ama] = 1
    
    plt.bar(amas.keys(),amas.values())
    plt.title("Ratio of buried to all aminoacids")
    plt.show()
    
    
def get_ratio(amk_type):
    print(amk_type)
    amk_pol = sum(num for amk,num in amk_type.items() if amk in polar_amks)
    print(amk_pol)
    print(sum(num for num in amk_type.values()))
    return amk_pol/sum(num for num in amk_type.values())

def get_polarity(structure):
    buried, exposed = get_buried_exposed_by_amk(structure)
    buried_pol_ratio = get_ratio(buried)
    exposed_pol_ratio = get_ratio(exposed)
    print(buried_pol_ratio)
    print(exposed_pol_ratio)

structure = get_structure()
print(get_buried_ratio(structure))
#print(get_diameter(structure))
