#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

import spglib
from ase.io import read
from pymatgen.core import Structure

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# class MyStructure:
#     def __init__(self, file):
#         with open(file, 'r') as f:
#             content = f.readlines()
#             f.close()
#
#         content = pd.DataFrame(content)
#         content.replace(r'\n', '', regex=True, inplace=True)
#
#         self.file = file
#         self.poscar = content
#
#     @property
#     def species(self):
#         """
#                 Element species of the structure.
#
#         @rtype: list of string
#         """
#         species = self.poscar.iloc[5, 0].split()
#
#         return species
#
#     @property
#     def composition(self):
#         """
#                 The numbers of atoms in the structure.
#
#         @rtype: list of int
#         """
#         composition = self.poscar.iloc[6, 0].split()
#         composition = [int(num) for num in composition]
#
#         return composition
#
#     @property
#     def lattice(self):
#         """
#                 Lattice vector of the structure.
#
#         @rtype: matrix
#         """
#         vector = self.poscar.iloc[2:5, 0].str.split(expand=True)
#         vector = vector.values.tolist()
#
#         vector = [[float(num) for num in row] for row in vector]
#         vector = np.array(vector)
#
#         return vector
#
#     @property
#     def lattice_parameters(self):
#         """
#                 Lattice parameters of the structure.
#
#         @rtype: tuple of floats
#         """
#         lattice = self.lattice
#         a = lattice[0]
#         b = lattice[1]
#         c = lattice[2]
#
#         _a = np.linalg.norm(a)
#         _b = np.linalg.norm(b)
#         _c = np.linalg.norm(c)
#
#         alpha = np.rad2deg(np.arccos(np.dot(b, c) / (_b * _c)))
#         beta = np.rad2deg(np.arccos(np.dot(a, c) / (_a * _c)))
#         gamma = np.rad2deg(np.arccos(np.dot(a, b) / (_a * _b)))
#
#         return _a, _b, _c, alpha, beta, gamma
#
#     @property
#     def _is_cartesian(self):
#         """
#                 Check whether the type of atomic coordinates is Cartesian.
#
#         @rtype: bool
#         """
#         coordinate_form = self.poscar.iloc[7, 0].split()[0].lower()
#
#         if coordinate_form == 'Cartesian'.lower() or coordinate_form == 'C'.lower():
#             is_cartesian = True
#         elif coordinate_form == 'Direct'.lower() or coordinate_form == 'D'.lower():
#             is_cartesian = False
#         else:
#             raise ValueError('The coordinate form cannot be recognized: {}'.format(coordinate_form))
#
#         return is_cartesian
#
#     @property
#     def positions(self):
#         """
#                 The Cartesian coordinates of all atoms in the structure.
#
#         @rtype: matrix
#         """
#         is_cartesian = self._is_cartesian
#
#         atomic_coordinates = self.poscar.iloc[8:, 0].str.split(expand=True)
#         atomic_coordinates = atomic_coordinates.drop([3], axis=1).values.tolist()
#         atomic_coordinates = [[float(num) for num in row] for row in atomic_coordinates]
#
#         if is_cartesian:
#             atomic_coordinates = np.array(atomic_coordinates)
#
#         else:
#             vector = self.poscar.iloc[2:5, 0].str.split(expand=True)
#             vector = vector.values.tolist()
#             vector = [[float(num) for num in row] for row in vector]
#             vector = np.array(vector)
#
#             atomic_coordinates = np.array(atomic_coordinates)
#
#             atomic_coordinates = np.dot(atomic_coordinates, vector)
#
#         return atomic_coordinates


# class Symmetry:
#     def __init__(self, file):
#         crystal = read(file)
#         sym_info = spglib.get_symmetry_dataset(crystal)
#
#         self.crystal = crystal
#         self.sym_info = sym_info
#
#     @property
#     def space_group(self):
#         """
#                 The space group symbol and number for the given structure.
#
#         @rtype: tuple
#         """
#         spg = self.sym_info['international']
#         spg_num = self.sym_info['number']
#
#         return spg, spg_num
#
#     @property
#     def refine(self):
#         """
#             Return refined cell.
#
#         :return: a tuple of refined cell: (lattice, positions, numbers)
#         @rtype: tuple
#         """
#         lattice, scaled_positions, numbers = spglib.refine_cell(self.crystal)
#
#         return lattice, scaled_positions, numbers
#
#     @property
#     def primitive_cell(self):
#         """
#                 Find primitive cell of the structure.
#
#         :return: a tuple of primitive cell: (lattice, positions, numbers)
#         @rtype: tuple
#         """
#         lattice, positions, numbers = spglib.find_primitive(self.crystal)
#
#         return lattice, positions, numbers
#
#     @property
#     def transformation(self):
#         """
#                 Transformation matrix from input lattice to standardized lattice:
#                     L^original = L^standardized * Tmat.
#
#         :return: Transformation matrix
#         @rtype: ndarray
#         """
#         transformation_matrix = self.sym_info['transformation_matrix']
#
#         return transformation_matrix
#
#     @property
#     def rotations(self):
#         """
#                 Matrix (rotation) parts of space group operations.
#
#         :return: Rotation matrix
#         @rtype: ndarray
#         """
#         rotation_matrix = self.sym_info['rotations']
#
#         return rotation_matrix
#
#     @property
#     def translations(self):
#         """
#                 Vector (translation) parts of space group operations.
#         :return: Translation matrix
#         @rtype: ndarray
#         """
#         translation_matrix = self.sym_info['translations']
#
#         return translation_matrix
#
#     @property
#     def apply_random_symmetry_operations(self):
#         """
#                 Apply random symmetry operations of space groups.
#
#         :return: The atom coordinates matrix after transformation
#         @rtype: ndarray
#         """
#         v_old = self.crystal.positions          # Cartesian coordinates
#         rotation_matrix = self.sym_info['rotations']
#         translation_matrix = self.sym_info['translations']
#
#         index = np.random.randint(len(rotation_matrix))
#         v_new = np.dot(v_old, rotation_matrix[index]) + translation_matrix[index]
#
#         return v_new, index


def get_space_group_symbol(struc_path):
    structure = Structure.from_file(struc_path)
    spga = SpacegroupAnalyzer(structure)
    spg_symbol = spga.get_space_group_symbol()

    return spg_symbol


def get_space_group_number(struc_path):
    structure = Structure.from_file(struc_path)
    spga = SpacegroupAnalyzer(structure)
    spg_number = spga.get_space_group_number()

    return spg_number


if __name__ == '__main__':
    file = r'F:\GA_CSP\File\Structures\ini_population_relaxed\0-69_from_mp-976122-POSCAR'
    # structure = MyStructure(file)
    # sym = Symmetry(file)
    # print(structure.positions)
    # print(sym.primitive_cell[1])
