#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import pickle
import random

import pandas as pd

from Selection.EliteSelection import elite_selection
from Initialization.StructureGenerator import check_input_element
from Initialization.SeedsOperation import get_seeds_from_mp, get_structure_from_seed
from Evolution.Operator import gen_random, gen_mutation
from Model.CHGNet.Relaxation import relax
from Utility.Deduplication import Pool, check_duplication


def split_candidates(candidates, ratios=None):
    if ratios is None:
        ratios = [0.5, 0.3, 0.2]
    if sum(ratios) != 1.0:
        raise ValueError('Ratios must sum to 1.0')

    sublists = []
    total_num = len(candidates)
    random.shuffle(candidates)

    nums = [int(total_num * i) for i in ratios]
    idx = 0
    for num in nums:
        start, end = idx, idx + num
        sublists.append(candidates[start:end])
        idx += num

    return tuple(sublists)


def get_restart_gen(pop_dir=None):
    if pop_dir is None:
        pop_dir = r'F:\GA_CSP\File\Structures'

    pop_list = os.listdir(pop_dir)
    pop_relaxed_list = []

    for pop in pop_list:
        if pop.endswith('relaxed'):
            gen = int(pop.split('_')[0])
            pop_relaxed_list.append(gen)

    return max(pop_relaxed_list)


def get_restart_id(pop_dir=None):
    if pop_dir is None:
        pop_dir = r'F:\GA_CSP\File\Structures'

    restart_gen = get_restart_gen(pop_dir)
    cur_gen_dir = os.path.join(pop_dir, '{}_population'.format(restart_gen + 1))
    cur_list = os.listdir(cur_gen_dir)
    n = len(cur_list)

    return n


def pop_initialization(species, pop_size, min_atoms, max_atoms, fmax=0.05, steps=1000, type_weight=None, seeds=False):
    """
            Generate crystal structures from seeds or random.

    @type species: list of string
    @type pop_size: int
    @type min_atoms: int
    @type max_atoms: int
    @type type_weight: list of number
    @type seeds: bool

    :param species: Species of the system
    :param pop_size: Size of the initial population
    :param min_atoms: Minimum number of atoms of structure
    :param max_atoms: Maximum number of atoms of structure
    :param type_weight: The proportion of compound type - (unary, binary, ternary)
    :return: random crystal structures generated from random or seeds
    """
    # Check the rationality of the input element
    species = check_input_element(species)
    # Create DataFrame of population information
    column_title = ['ID', 'nelements', 'crystal_system', 'space_group_number', 'from', 'composition', species[0],
                    species[1], species[2]]
    pop_data = pd.DataFrame(columns=column_title)

    pop_dir = r'F:\GA_CSP\File\Structures'
    if not os.path.exists(pop_dir):
        os.makedirs(pop_dir)

    ini_pop_dir = os.path.join(pop_dir, '0_population')
    if not os.path.exists(ini_pop_dir):
        os.makedirs(ini_pop_dir)

    if type_weight is None:
        type_scale = [1, 4, 5]      # Default compound proportion of compound type[unary, binary, ternary]
    elif len(type_weight) == 3:
        type_scale = type_weight
    else:
        raise ValueError('Type weight should be a list with three numbers, not {}'.format(type_weight))

    pool = Pool()

    n = 1  # ID of individual

    if seeds is True:
        seeds_path = r'F:\GA_CSP\File\Seeds'
        if not os.path.exists(seeds_path):
            os.makedirs(seeds_path)
        if not os.listdir(seeds_path):
            get_seeds_from_mp(species)

        # Sort seed information by formation energy per atom
        seeds_info = pd.read_csv(r'F:\GA_CSP\File\seeds_info.csv')
        seeds_sorted = seeds_info.sort_values(by='formation_energy_per_atom', ignore_index=True)

        seeds_chosen = seeds_sorted[:]  # Add all seed structures.

        # Ensure that the initial population scale is at least 10 more than the number of seeds
        if len(seeds_sorted) > pop_size - 10:
            pop_size = len(seeds_chosen) + 10
        else:
            pop_size = pop_size - len(seeds_sorted)

        # Generate structures from seeds
        for mp_id in seeds_chosen['material_id']:
            struct, individual_data = get_structure_from_seed(species, mp_id)
            if check_duplication(struct, pool.struct_list):
                continue
            struct.to_file(os.path.join(ini_pop_dir, '0_{}-from_{}-POSCAR'.format(n, mp_id)))
            pool.add(struct)
            individual_data.insert(loc=0, column='ID', value='0_{}'.format(n))
            pop_data = pd.concat([pop_data, individual_data], ignore_index=True)
            n += 1

    unary_num = (type_scale[0] / sum(type_scale)) * pop_size
    binary_num = (type_scale[1] / sum(type_scale)) * pop_size
    # ternary_num = (type_scale[2] / sum(type_scale)) * pop_size

    # Generate structures from random
    while len(pop_data) < unary_num:
        n, pop_data, pool = gen_random(0, n, 1, species, min_atoms, max_atoms, ini_pop_dir, pop_data, pool)

    while len(pop_data) < sum([unary_num, binary_num]):
        n, pop_data, pool = gen_random(0, n, 2, species, min_atoms, max_atoms, ini_pop_dir, pop_data, pool)

    while len(pop_data) < pop_size:
        n, pop_data, pool = gen_random(0, n, 3, species, min_atoms, max_atoms, ini_pop_dir, pop_data, pool)

    # Write population information
    pop_info_dir = r'F:\GA_CSP\File\Population_info'
    if not os.path.exists(pop_info_dir):
        os.mkdir(pop_info_dir)

    pop_data.to_csv(os.path.join(pop_info_dir, '0_pop-info.csv'), index=False)

    pool_path = r'F:\GA_CSP\File\pool.pkl'
    with open(pool_path, 'wb') as f:
        pickle.dump(pool, f)

    relax(fmax, steps, species, ini_pop_dir, 0)


def pop_iteration(total_gen_num: int, species, pop_size, min_atoms, max_atoms, crossover_rate, mutation_rate, fmax=0.05, steps=1000, restart=False, type_weight=None):

    pop_info_dir = r'F:\GA_CSP\File\Population_info'
    if not os.path.exists(pop_info_dir):
        os.mkdir(pop_info_dir)

    pop_dir = r'F:\GA_CSP\File\Structures'
    # if not os.path.exists(pop_dir):
    #     os.makedirs(pop_dir)

    if type_weight is None:
        type_scale = [1, 4, 5]      # Default compound proportion of compound type[unary, binary, ternary]
    elif len(type_weight) == 3:
        type_scale = type_weight
    else:
        raise ValueError('Type weight should be a list with three numbers, not {}'.format(type_weight))

    if restart is True:
        gen = get_restart_gen(pop_dir) + 1
    else:
        gen = 1

    for cur_gen in range(gen, total_gen_num):
        pool_path = r'F:\GA_CSP\File\pool.pkl'
        if os.path.isfile(pool_path):
            with open(pool_path, 'rb') as f:
                pool = pickle.load(f)
        else:
            pool = Pool()

        print("""
            Generation-{}
        """.format(cur_gen))

        last_pop_file = '{}_relaxed_pop-info.csv'.format(cur_gen - 1)
        last_pop_dir = os.path.join(pop_info_dir, last_pop_file)
        last_pop_info = pd.read_csv(last_pop_dir)

        cur_struct_dir = os.path.join(pop_dir, '{}_population'.format(cur_gen))
        if not os.path.exists(cur_struct_dir):
            os.makedirs(cur_struct_dir)
        else:
            shutil.rmtree(cur_struct_dir)
            os.makedirs(cur_struct_dir)

        last_struct_dir = os.path.join(pop_dir, '{}_population_relaxed'.format(cur_gen - 1))
        last_struct_file = os.listdir(last_struct_dir)
        last_struct_id = {}
        for name in last_struct_file:
            name_id = name.split('_')[0]
            last_struct_id[name_id] = name

        selected_indi = elite_selection(last_pop_info, scale=0.2)
        selected_struct_file = []
        for id in selected_indi['ID']:
            if id in last_struct_id.keys():
                selected_struct_file.append(last_struct_id[id])

        # Create DataFrame of population information
        column_title = ['ID', 'nelements', 'crystal_system', 'space_group_number', 'from', 'composition', species[0],
                        species[1], species[2]]
        pop_data = pd.DataFrame(columns=column_title)

        # Start heredity and mutation on the structure

        n = 1       # ID of individual

        # crossover_comb = []
        # for _ in range(crossover_num):
        #     while True:
        #         combination = sorted(random.sample(selected_struct_file, 2))
        #
        #         if combination not in crossover_comb:
        #             crossover_comb.append(combination)
        #             break
        #
        # for parents in crossover_combi:
        #     crossover(parents)

        mutation_candidates = random.sample(selected_struct_file, k=int(mutation_rate * len(selected_struct_file)))
        strain_candi, soft_candi, permutation_candi = split_candidates(mutation_candidates)
        n, pop_data, pool = gen_mutation(cur_gen, n, species, 'strain', strain_candi, last_struct_dir, cur_struct_dir, pop_data, pool)
        n, pop_data, pool = gen_mutation(cur_gen, n, species, 'soft', soft_candi, last_struct_dir, cur_struct_dir, pop_data, pool)
        n, pop_data, pool = gen_mutation(cur_gen, n, species, 'permutation', permutation_candi, last_struct_dir, cur_struct_dir, pop_data, pool)

        unary_num = (type_scale[0] / sum(type_scale)) * pop_size
        binary_num = (type_scale[1] / sum(type_scale)) * pop_size
        # ternary_num = (type_scale[2] / sum(type_scale)) * pop_size

        # Generate structures from random
        while len(pop_data) < unary_num:
            n, pop_data, pool = gen_random(cur_gen, n, 1, species, min_atoms, max_atoms, cur_struct_dir, pop_data, pool)

        while len(pop_data) < sum([unary_num, binary_num]):
            n, pop_data, pool = gen_random(cur_gen, n, 2, species, min_atoms, max_atoms, cur_struct_dir, pop_data, pool)

        while len(pop_data) < pop_size:
            n, pop_data, pool = gen_random(cur_gen, n, 3, species, min_atoms, max_atoms, cur_struct_dir, pop_data, pool)

        # Write population information
        pop_data.to_csv(os.path.join(pop_info_dir, '{}_pop-info.csv'.format(cur_gen)), index=False)

        pool_path = r'F:\GA_CSP\File\pool.pkl'
        with open(pool_path, 'wb') as f:
            pickle.dump(pool, f)

        relax(fmax, steps, species, cur_struct_dir, cur_gen)


if __name__ == '__main__':
    # pop_initialization(['Ge', 'Sb', 'Te'], 30, 6, 30)
    pop_iteration(3, ['Ge', 'Sb', 'Te'], 30, 6, 30, 1, 5)
