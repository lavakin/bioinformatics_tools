import sys
from pathlib import Path
import shutil
from Bio.PDB import *
from scipy.spatial.distance import *
import argparse


def arg_parser():
    parser = argparse. ArgumentParser()
    parser.add_argument("pdb", type=str)
    parser.add_argument('-f', '--from_file', action='store_true')
    return parser.parse_args()

def get_structure():
    args = arg_parser()
    print(args)
    file_name = args.pdb + ".pdb"
    if not args.from_file:
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(sys.argv[1])
        shutil.copy(file_name, 'pdbs')
        Path(file_name).unlink(missing_ok = True)
    parser = PDBParser()
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


structure = get_structure()
get_info(structure)
