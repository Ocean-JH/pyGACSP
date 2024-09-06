#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random
import itertools
import pandas as pd

from pyxtal import pyxtal
from pyxtal import msg


def get_random_crystal_system(crystal_systems=None, systems_prob=None):
    """
    @type crystal_systems: list of string
    @type systems_prob: list of number
    :param crystal_systems: list of candidate crystal systems
    :param systems_prob: list of weights for crystal systems
    :return: Random crystal system and space group selected by weight
    """

    if crystal_systems is None:
        crystal_systems = ['Triclinic', 'Monoclinic', 'Orthorhombic', 'Tetragonal', 'Trigonal', 'Hexagonal', 'Cubic']
    if systems_prob is None:
        systems_prob = [2, 13, 59, 68, 7, 45, 36]       # Default weights

    crystal_system = random.choices(crystal_systems, weights=systems_prob, k=1)[0]

    if crystal_system == 'Triclinic':
        spg = random.randint(1, 2)
    elif crystal_system == 'Monoclinic':
        spg = random.randint(3, 15)
    elif crystal_system == 'Orthorhombic':
        spg = random.randint(16, 74)
    elif crystal_system == 'Tetragonal':
        spg = random.randint(75, 142)
    elif crystal_system == 'Trigonal':
        spg = random.randint(143, 167)
    elif crystal_system == 'Hexagonal':
        spg = random.randint(168, 194)
    else:
        spg = random.randint(195, 230)

    random_crystal_system = [crystal_system, spg]
    return random_crystal_system


def check_input_element(species):
    """
            Returns the element symbol based on atomic number.
    @type species: List of str or int or both
    :param species: List of Element symbol or corresponding atomic number
    :return: List of atomic number of the element
    @rtype: List of str
    """

    element_symbol = ['H', 'He',
                      'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                      'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                      'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
                      'Br', 'Kr',
                      'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
                      'I', 'Xe',
                      'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                      'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                      'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', ' Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
                      'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', ' Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
                      ]
    atomic_number = list(range(1, 119))
    element_dict = dict(zip(atomic_number, element_symbol))

    output_ele = []
    for element in species:
        if type(element) is int:
            if element not in element_dict:
                raise ValueError('Atomic number must be an integer between 1 and 119, not {}'.format(element))
            output_ele.append(element_dict[element])
        elif type(element) is str:
            if element not in element_symbol:
                raise ValueError('No {} element symbol in the periodic table, please check your input'.format(element))
            output_ele.append(element)
        else:
            raise ValueError('The input must be element symbol or corresponding atomic number, not {}'.format(element))

    return output_ele


def get_random_structure(comp_type, species, min_atoms, max_atoms):
    """
            Get random structure from pyXtal.

    @type comp_type: int
    @type species: list of str
    @type min_atoms: int
    @type max_atoms: int
    :param comp_type: Type of individual - (unary, binary, ternary)
    :param species: Species of the system
    :param min_atoms: Minimum number of atoms of structure
    :param max_atoms: Maximum number of atoms of structure
    :return: Random structure and corresponding information
    """

    # Generate species combinations
    comp_unary = list(itertools.combinations(species, 1))
    comp_binary = list(itertools.combinations(species, 2))
    comp_ternary = list(itertools.combinations(species, 3))

    crystal = pyxtal()
    rand_crystal_system = get_random_crystal_system()
    crystal_system = rand_crystal_system[0]
    space_group_number = rand_crystal_system[1]

    if comp_type == 1:

        while True:
            try:
                comp = random.choice(comp_unary)

                num_atoms = [random.randint(min_atoms, max_atoms)]

                crystal.from_random(3, space_group_number, comp, num_atoms)
            except msg.Comp_CompatibilityError:
                rand_crystal_system = get_random_crystal_system()
                crystal_system = rand_crystal_system[0]
                space_group_number = rand_crystal_system[1]
                continue
            else:
                break

        composition = '{}'.format(comp[0]) + '{}'.format(num_atoms[0])

        structure_data = {'nelements': [comp_type],
                          'crystal_system': [crystal_system],
                          'space_group_number': [space_group_number],
                          'from': 'Random',
                          'composition': [composition],
                          species[0]: [0],
                          species[1]: [0],
                          species[2]: [0]}
        structure_data = pd.DataFrame(structure_data)
        structure_data.loc[0, '{}'.format(comp[0])] = num_atoms

    elif comp_type == 2:

        while True:
            try:
                comp = random.choices(comp_binary, weights=[0.1, 0.45, 0.45], k=1)[0]

                a_atoms = random.randint(1, max_atoms - 1)
                if a_atoms < min_atoms:
                    b_atoms = random.randint(min_atoms - a_atoms, max_atoms - a_atoms)
                else:
                    b_atoms = random.randint(1, max_atoms - a_atoms)
                num_atoms = [a_atoms, b_atoms]

                crystal.from_random(3, space_group_number, comp, num_atoms)
            except msg.Comp_CompatibilityError:
                rand_crystal_system = get_random_crystal_system()
                crystal_system = rand_crystal_system[0]
                space_group_number = rand_crystal_system[1]
                continue
            else:
                break

        composition = ('{}'.format(comp[0]) + '{}'.format(num_atoms[0]) +
                       '{}'.format(comp[1]) + '{}'.format(num_atoms[1]))

        structure_data = {'nelements': [comp_type],
                          'crystal_system': [crystal_system],
                          'space_group_number': [space_group_number],
                          'from': 'Random',
                          'composition': [composition],
                          species[0]: [0],
                          species[1]: [0],
                          species[2]: [0]}
        structure_data = pd.DataFrame(structure_data)
        structure_data.loc[0, '{}'.format(comp[0])] = num_atoms[0]
        structure_data.loc[0, '{}'.format(comp[1])] = num_atoms[1]

    elif comp_type == 3:
        comp = comp_ternary[0]

        while True:
            try:
                a_atoms = random.randint(1, max_atoms - 2)
                b_atoms = random.randint(1, max_atoms - a_atoms - 1)
                if a_atoms + b_atoms < min_atoms:
                    c_atoms = random.randint(min_atoms - (a_atoms + b_atoms), max_atoms - (a_atoms + b_atoms))
                else:
                    c_atoms = random.randint(1, max_atoms - (a_atoms + b_atoms))
                num_atoms = [a_atoms, b_atoms, c_atoms]

                crystal.from_random(3, space_group_number, comp, num_atoms)
            except Exception:
                rand_crystal_system = get_random_crystal_system()
                crystal_system = rand_crystal_system[0]
                space_group_number = rand_crystal_system[1]
                continue
            else:
                break

        composition = ('{}'.format(comp[0]) + '{}'.format(num_atoms[0]) +
                       '{}'.format(comp[1]) + '{}'.format(num_atoms[1]) +
                       '{}'.format(comp[2]) + '{}'.format(num_atoms[2]))

        structure_data = {'nelements': [comp_type],
                          'crystal_system': [crystal_system],
                          'space_group_number': [space_group_number],
                          'from': 'Random',
                          'composition': [composition],
                          species[0]: [0],
                          species[1]: [0],
                          species[2]: [0]}
        structure_data = pd.DataFrame(structure_data)
        structure_data.loc[0, '{}'.format(comp[0])] = num_atoms[0]
        structure_data.loc[0, '{}'.format(comp[1])] = num_atoms[1]
        structure_data.loc[0, '{}'.format(comp[2])] = num_atoms[2]

    else:
        raise ValueError('nelements should be an integer in (1, 2, 3), not {}'.format(type))

    num_1 = structure_data.loc[0, species[0]]
    num_2 = structure_data.loc[0, species[1]]
    num_3 = structure_data.loc[0, species[2]]

    print(
        'Structure generated with symmetry {}-{}. Composition: {}   {}   {}.'
        .format(crystal_system, space_group_number, num_1, num_2, num_3))

    structure = crystal.to_pymatgen()
    return structure, structure_data, comp_type, crystal_system, space_group_number, composition, num_1, num_2, num_3


if __name__ == '__main__':
    obj = get_random_structure(3, ['Ge', 'Sb', 'Te'], 6, 20)
    print(type(obj[0]))
