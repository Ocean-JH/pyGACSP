#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import List
from ase import Atoms
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher


class Pool:
    def __init__(self):
        self.structures = []
        self.atoms = []

    @property
    def size(self) -> int:
        return len(self.structures)

    @property
    def struct_list(self):
        return self.structures

    @property
    def atoms_list(self):
        return self.atoms

    def add(self, structure):
        if isinstance(structure, Structure):
            self.structures.append(structure)
            self.atoms.append(AseAtomsAdaptor.get_atoms(structure))
        elif isinstance(structure, Atoms):
            self.structures.append(AseAtomsAdaptor.get_structure(structure))
            self.atoms.append(structure)

    # def _delete(self, num):
    #     del self.structures[:-num]
    #     del self.atoms[:-num]


def check_duplication(structure, pool_list):
    if len(pool_list) == 0:
        return False

    matcher = StructureMatcher(scale=False, comparator=ElementComparator())
    for struct in pool_list:
        is_equivalent = matcher.fit(structure, struct)
        if is_equivalent:
            return True
        else:
            continue

    return False


if __name__ == "__main__":
    pool = Pool()
    indiv1_path = r'CONTCAR'
    indiv2_path = r'CONTCAR'
    # indiv3_path = r'POSCAR'
    indiv1 = Structure.from_file(indiv1_path)
    indiv2 = Structure.from_file(indiv2_path)
    # indiv3 = Structure.from_file(indiv3_path)
    pool.add(indiv1)
    # pool.add(indiv2)
    dup = check_duplication(indiv1, pool.struct_list)
    print(dup)
