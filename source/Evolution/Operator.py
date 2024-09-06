#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas as pd

from ase import Atoms
from ase.ga.utilities import CellBounds
from ase.ga.standardmutations import PermutationMutation, StrainMutation
from ase.ga.soft_mutation import SoftMutation
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from Initialization.StructureGenerator import get_random_structure
from Utility.Deduplication import check_duplication


# def lattice_matching(original_poscar, new_poscar):
#     """
#             Keep the lattice vectors of parents consistent to facilitate subsequent crossover operation.
#     :param original_poscar: Structure need to change lattice vectors
#     :param new_poscar: Structure with the target lattice vectors
#     :return: Structure after changing the lattice vector
#     @:rtype: Structure
#     """
#     lattice1 = original_poscar.lattice
#     lattice2 = new_poscar.lattice
#
#     matrix1 = np.array(original_poscar.lattice.matrix).T
#     matrix2 = np.array(new_poscar.lattice.matrix).T
#
#     transform_matrix = np.dot(np.linalg.inv(matrix1), matrix2)
#
#     sites1 = original_poscar.sites
#
#     transformed_structure = Structure(new_poscar.lattice, [], [])
#     for site in sites1:
#         transformed_coords = np.dot(site.coords, transform_matrix)
#         transformed_structure.append(site.specie, transformed_coords, properties=site.properties)
#
#     new_poscar = Poscar(transformed_structure)
#
#     return new_poscar


# def crossover(parent1, parent2):
#     """
#             Crossover to generate offsprings.
#     @type parent1: str
#     @type parent2: str
#     :param parent1: structure after selection
#     :param parent2:
#     :return: structure after crossover
#     :rtype: POSCAR
#     """
#     p1 = Structure.from_file(parent1)
#     p2 = Structure.from_file(parent2)
#
#     lattice1 = p1.lattice
#     lattice2 = p2.lattice
#
#     site1 = p1.sites
#     site2 = p2.sites
#
#     if lattice1.matrix != lattice2.matrix:


def crossover(father, mother):
    f = Structure.from_file(father)
    m = Structure.from_file(mother)

    f = AseAtomsAdaptor.get_atoms(f)
    m = AseAtomsAdaptor.get_atoms(m)
    f.info['confid'] = '{}'.format(father)
    m.info['confid'] = '{}'.format(mother)

    atom_numbers_to_optimize = f.get_atomic_numbers()
    n_to_optimize = len(atom_numbers_to_optimize)
    blmin = closest_distances_generator(atom_numbers_to_optimize, 0.5)

    slab = Atoms('', pbc=True)
    pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)

    child = pairing.cross(f, m)
    child = AseAtomsAdaptor.get_structure(child)

    return child


def _strainmutate(parent):
    f = Structure.from_file(parent)
    f = AseAtomsAdaptor.get_atoms(f)

    atom_numbers_to_optimize = f.get_atomic_numbers()
    blmin = closest_distances_generator(atom_numbers_to_optimize, 0.5)
    cellbounds = CellBounds(bounds={'alpha': [45, 135], 'beta': [45, 135],
                                    'gamma': [45, 135], 'a': [2, 40],
                                    'b': [2, 40], 'c': [2, 40]})

    strainmut = StrainMutation(blmin, cellbounds)

    offspring = strainmut.mutate(f)
    offspring = AseAtomsAdaptor.get_structure(offspring)

    print(
        'Structure generated from Strain Mutation.'
    )

    return offspring


def _softmutate(parent):
    f = Structure.from_file(parent)
    f = AseAtomsAdaptor.get_atoms(f)
    f.info['confid'] = '{}'.format(parent)

    atom_numbers_to_optimize = f.get_atomic_numbers()
    blmin_soft = closest_distances_generator(atom_numbers_to_optimize, 0.1)

    softmut = SoftMutation(blmin_soft, bounds=[1, 3], used_modes_file=None)

    offspring = softmut.mutate(f)
    offspring = AseAtomsAdaptor.get_structure(offspring)

    print(
        'Structure generated from Soft Mutation.'
    )

    return offspring


def _permutation(parent):
    f = Structure.from_file(parent)
    f = AseAtomsAdaptor.get_atoms(f)

    atom_numbers_to_optimize = f.get_atomic_numbers()
    blmin_soft = closest_distances_generator(atom_numbers_to_optimize, 0.1)

    permut = PermutationMutation(None, blmin=blmin_soft)

    offspring = permut.mutate(f)
    offspring = AseAtomsAdaptor.get_structure(offspring)

    print(
        'Structure generated from Permutation.'
    )

    return offspring


def gen_mutation(cur_gen, id, species, operator, candidates, last_struct_dir, cur_struct_dir, pop_data, pool):
    for indiv in candidates:
        if operator == 'strain':
            offspring = _strainmutate(os.path.join(last_struct_dir, indiv))
            while check_duplication(offspring, pool.struct_list):
                offspring = _strainmutate(os.path.join(last_struct_dir, indiv))
        elif operator == 'soft':
            offspring = _softmutate(os.path.join(last_struct_dir, indiv))
            while check_duplication(offspring, pool.struct_list):
                offspring = _softmutate(os.path.join(last_struct_dir, indiv))
        elif operator == 'permutation':
            offspring = _permutation(os.path.join(last_struct_dir, indiv))
            while check_duplication(offspring, pool.struct_list):
                offspring = _permutation(os.path.join(last_struct_dir, indiv))
        else:
            raise Exception('Unknown operator: {}'.format(operator))

        offspring.to(os.path.join(os.path.join(cur_struct_dir, '{}_{}-from_mutation-POSCAR'.format(cur_gen, id))))
        pool.add(offspring)
        analyzer = SpacegroupAnalyzer(offspring)
        n_elems = offspring.n_elems
        spg_num = offspring.get_space_group_info()[1]
        composition = str(offspring.composition)
        composition = composition.replace(" ", "")
        crystal_system = analyzer.get_crystal_system()

        comp_dict = offspring.composition.get_el_amt_dict()

        offspring_data = {'nelements': [n_elems],
                          'crystal_system': [crystal_system],
                          'space_group_number': [spg_num],
                          'from': 'mutation({})-'.format(operator) + indiv.split('_')[0],
                          'composition': [composition],
                          species[0]: [0],
                          species[1]: [0],
                          species[2]: [0]}
        offspring_data = pd.DataFrame(offspring_data)
        for i in range(n_elems):
            offspring_data.loc[0, '{}'.format(list(offspring.composition.keys())[i].symbol)] = int(
                comp_dict[list(offspring.composition.keys())[i].symbol])
        offspring_data.insert(loc=0, column='ID', value='{}_{}'.format(cur_gen, id))
        pop_data = pd.concat([pop_data, offspring_data], ignore_index=True)
        id += 1

    return id, pop_data, pool


def gen_random(gen, id, comp_type, species, min_atoms, max_atoms, struct_dir, pop_data, pool):
    struct, individual_data, *other = get_random_structure(comp_type, species, min_atoms, max_atoms)
    while check_duplication(struct, pool.struct_list):
        struct, individual_data = get_random_structure(comp_type, species, min_atoms, max_atoms)
    struct.to_file(os.path.join(struct_dir, '{}_{}-from_random-POSCAR'.format(gen, id)))
    pool.add(struct)
    individual_data.insert(loc=0, column='ID', value='{}_{}'.format(gen, id))
    pop_data = pd.concat([pop_data, individual_data], ignore_index=True)
    id += 1

    return id, pop_data, pool


if __name__ == '__main__':
    father = r'C:\Users\Wangjq\Downloads\0_86-from_random-POSCAR'
    mother = r'C:\Users\Wangjq\Downloads\0_74-from_mp-1197020-POSCAR'
    son = _permutation(father)
    comp_dict = son.composition.get_el_amt_dict()
    print(list(son.composition.keys())[1].symbol)
    print(comp_dict[list(son.composition.keys())[1].symbol])
    # son.to_file(r'C:\Users\Wangjq\Downloads\offspring-POSCAR')
