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


class PDB:
    def __init__(self, file_name, from_file=True):
        """
        :param file_name: name of the pdb file or pdb if it should be downloaded
        :param from_file: from file stored on the machine or should be downloaded
        """
        self.structure = PDB.get_structure(file_name, from_file)

    def get_info(self):
        """
        :return: information about the pdb
        """
        print("number of models:", len(self))
        chains = 0
        residues = 0
        atoms = 0
        for model in self:
            chains += len(model)
            for chain in model:
                residues += len(chain)
                for residue in chain:
                    atoms += len(residue)
        return {"number of chains": chains, "number of residues" : residues, "number of atoms": atoms}


    def get_structure(file_name, from_file=True):
        """
        :param from_file: from file stored on the machine or should be downloaded
        :return: parsed structure from pdb
        """
        parser = PDBParser()
        if from_file:
            return parser.get_structure(file_name, file_name)
        if not from_file:
            comp_file_name = file_name + ".pdb"
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(file_name)
            shutil.copy(comp_file_name, 'pdbs')
            Path(comp_file_name).unlink(missing_ok=True)
        return parser.get_structure(file_name, "pdbs/" + file_name)


    def get_diameter(self):
        """
        :return: get the distance between two most distant atoms
        """
        atoms = [a.get_coord() for a in self.structure.get_atoms()]
        return (max(pdist(atoms, 'euclidean')))


    def get_nearest_substructures(hetatm, distance: float, type_of_substr: str, atom_list):
        """
        :param distance: max distance from query atom
        :param type_of_substr: the type of substructure (atom/residue)
        :param atom_list: list of all atoms
        :return: list of nearest substructures
        """
        ns = NeighborSearch(atom_list)
        return ns.search(hetatm.get_coord(), distance, level=type_of_substr)


    def get_nearest_residues(self, hetatm, distance: float):
        """
        :param hetatm: query atom
        :param distance: max distance
        :return: nearest residues of a query atom
        """
        # print(structure['A'][111]['CA'])
        return PDB.get_nearest_substructures(hetatm, distance, 'R', Selection.unfold_entities(self.structure, 'A'))


    def get_nearest_atoms(self, hetatm, distance: float):
        """
        :param hetatm: query atom
        :param distance: may distance
        :return: nearest atoms from a query atom
        """
        return PDB.get_nearest_substructures(hetatm, distance, 'A', Selection.unfold_entities(self.structure, 'A'))


    def get_buried_ratio(self):
        """
        :return: nubmer of buried and exposed CA in the structure
        """
        buried = len(list(i[1] for i in HSExposure.ExposureCN(self.structure[0], radius=10, offset=0) if i[1] >= 20))
        exposed = len(list(i[1] for i in HSExposure.ExposureCN(self.structure[0], radius=10, offset=0) if i[1] < 20))
        return {"buried_ratio": round(buried / (buried + exposed), 3),
                "exposed_ratio": round(exposed / (buried + exposed), 3)}

    def get_buried_exposed_by_amk(self):
        """
        :return: number of buried and exposed CA for each aminoacid
        """
        buried = Counter(
            list(i[0].resname for i in HSExposure.ExposureCN(self.structure[0], radius=10, offset=0) if i[1] >= 20))
        exposed = Counter(
            list(i[0].resname for i in HSExposure.ExposureCN(self.structure[0], radius=10, offset=0) if i[1] < 20))
        return buried, exposed


    def get_histogram(self):
        """
        :return: histogram of the ratio of buried and exposed CA for each aminoacid
        """
        amas = self.get_buried_exposed_ratio_by_amk()
        plt.bar(amas.keys(), amas.values())
        plt.title("Ratio of buried to all aminoacids")
        return plt
    

    def get_buried_exposed_ratio_by_amk(self):
        """
        :return: dict of aminoacids with their respective ratios of buried and exposed CAs
        """
        buried, exposed = self.get_buried_exposed_by_amk()
        amas = dict.fromkeys(set(exposed.keys()).union(set(buried.keys())), 0)
        for ama in amas.copy():
            if ama in buried.keys() and ama in exposed.keys():
                amas[ama] = round(buried[ama] / (buried[ama] + exposed[ama]), 3)
            elif ama in buried.keys():
                amas[ama] = 1
        return amas

    def get_polarity(self):
        """
        :return: ratio of polar to nonpolar aminoacids that are buried or exposed
        """
        def get_ratio(amk_type):
            amk_pol = sum(num for amk, num in amk_type.items() if amk in polar_amks)
            return amk_pol / sum(num for num in amk_type.values())

        buried, exposed = self.get_buried_exposed_by_amk()
        buried_pol_ratio = get_ratio(buried)
        exposed_pol_ratio = get_ratio(exposed)
        return {"buried ratio": buried_pol_ratio, "exposed ratio": exposed_pol_ratio}



