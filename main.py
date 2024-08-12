#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
        This is the main part of the program
        Author: Jianghai Wang@BUAA
"""

import time

from Utility import logo
from Initialization.PopulationGenerator import pop_iteration, pop_initialization, get_restart_gen


def main(species, pop_size, generations, min_atom, max_atom, crossover_rate, mutation_rate, restart=False):
    logo.main_logo()

    if restart is False:
        logo.pop_initialization_logo()
        pop_initialization(species, pop_size, min_atom, max_atom)

        logo.pop_iteration_logo()
        pop_iteration(generations, species, pop_size, min_atom, max_atom, crossover_rate, mutation_rate)
    else:
        cur_gen = get_restart_gen()
        logo.restart_logo(cur_gen)
        pop_iteration(generations, species, pop_size, min_atom, max_atom, crossover_rate, mutation_rate, restart=True)

    logo.finish_logo(generations)


if __name__ == '__main__':
    species = ['Ge', 'Sb', 'Te']
    pop_size = 50
    generations = 100
    min_atom = 4
    max_atom = 40
    crossover_rate = 0
    mutation_rate = 1

    main(species, pop_size, generations, min_atom, max_atom, crossover_rate, mutation_rate, restart=True)
