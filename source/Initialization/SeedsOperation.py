#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import itertools

import pandas as pd
from pyxtal import pyxtal
from mp_api.client import MPRester
from Utility.SplitComposition import splitcomposition


def get_seeds_from_mp(species, species_types=None, api=None):
    """
            Fetch seed data from Materials Project via API.

    @type species: list of string
    @type species_types: list of integer
    @type api: string
    :param species: Species of the system
    :param species_types: Type of individual - (unary, binary, ternary)
    :param api: API key of Materials Project
    """
    if species_types is None:
        species_types = [i for i in range(1, len(species) + 1)]     # Default type of species
    else:
        species_types = species_types

    if api is None:
        if os.path.exists(r'F:\GA_CSP\File\API_KEY'):
            with open(r'F:\GA_CSP\File\API_KEY') as f:
                api_key = f.read()
        else:
            raise ValueError('Please provide your API key of Materials Project.')
    else:
        api_key = api

    seeds_path = r'F:\GA_CSP\File\Seeds'
    if not os.path.exists(seeds_path):
        os.makedirs(seeds_path)

    if max(species_types) > len(species):
        raise ValueError('Impossible to form {}-component compound with {} species.'
                         .format(max(species_types), len(species)))

    # Generate species combinations
    for species_type in species_types:
        species_combinations = list(itertools.combinations(species, species_type))
        unit_style = '-{}'
        for combination in species_combinations:
            style = '{}' + unit_style * (species_type - 1)
            chemistry_system = style.format(*combination)

            with MPRester(api_key) as mpr:
                # list_of_available_fields = mpr.materials.summary.available_fields
                '''
                list_of_available_fields = ['builder_meta', 'nsites', 'elements', 'nelements', 'composition',
                                            'composition_reduced', 'formula_pretty', 'formula_anonymous', 'chemsys',
                                            'volume', 'density', 'density_atomic', 'symmetry', 'property_name',
                                            'material_id', 'deprecated', 'deprecation_reasons', 'last_updated',
                                            'origins', 'warnings', 'structure', 'task_ids',
                                            'uncorrected_energy_per_atom', 'energy_per_atom',
                                            'formation_energy_per_atom', 'energy_above_hull', 'is_stable',
                                            'equilibrium_reaction_energy_per_atom', 'decomposes_to', 'xas',
                                            'grain_boundaries', 'band_gap', 'cbm', 'vbm', 'efermi', 'is_gap_direct',
                                            'is_metal', 'es_source_calc_id', 'bandstructure', 'dos', 'dos_energy_up',
                                            'dos_energy_down', 'is_magnetic', 'ordering', 'total_magnetization',
                                            'total_magnetization_normalized_vol',
                                            'total_magnetization_normalized_formula_units', 'num_magnetic_sites',
                                            'num_unique_magnetic_sites', 'types_of_magnetic_species', 'k_voigt',
                                            'k_reuss', 'k_vrh', 'g_voigt', 'g_reuss', 'g_vrh', 'universal_anisotropy',
                                            'homogeneous_poisson', 'e_total', 'e_ionic', 'e_electronic', 'n',
                                            'e_ij_max', 'weighted_surface_energy_EV_PER_ANG2',
                                            'weighted_surface_energy', 'weighted_work_function', 'surface_anisotropy',
                                            'shape_factor', 'has_reconstructed', 'possible_species', 'has_props',
                                            'theoretical', 'database_IDs']
                '''

                docs = mpr.materials.summary.search(chemsys=chemistry_system,
                                                    fields=['material_id', 'database_IDs', 'nelements', 'chemsys',
                                                            'formula_anonymous', 'formula_pretty', 'composition',
                                                            'symmetry', 'energy_per_atom', 'formation_energy_per_atom',
                                                            'energy_above_hull', 'is_stable', 'theoretical', 'volume',
                                                            'density', 'possible_species', 'decomposes_to',
                                                            'deprecated', 'origins', 'warnings', 'last_updated',
                                                            'structure'])

                # Collect seeds information
                seeds_info = [
                    [doc.material_id, doc.database_IDs, doc.nelements, doc.chemsys, doc.formula_anonymous,
                     doc.formula_pretty, doc.composition, doc.symmetry, doc.energy_per_atom,
                     doc.formation_energy_per_atom, doc.energy_above_hull, doc.is_stable, doc.theoretical, doc.volume,
                     doc.density, doc.possible_species, doc.decomposes_to, doc.deprecated, doc.origins, doc.warnings,
                     doc.last_updated] for doc in docs]

            seeds_info = pd.DataFrame(seeds_info,
                                      columns=['material_id', 'database_IDs', 'nelements', 'chemsys',
                                               'formula_anonymous', 'formula_pretty', 'composition', 'symmetry',
                                               'energy_per_atom', 'formation_energy_per_atom', 'energy_above_hull',
                                               'is_stable', 'theoretical', 'volume', 'density', 'possible_species',
                                               'decomposes_to', 'deprecated', 'origins', 'warnings', 'last_updated'])

            if not os.path.exists(os.path.join(r'F:\GA_CSP\File', r'seeds_info.csv')):
                seeds_info.to_csv(os.path.join(r'F:\GA_CSP\File', r'seeds_info.csv'),
                                  mode='a', index=False)
            else:
                seeds_info.to_csv(os.path.join(r'F:\GA_CSP\File', r'seeds_info.csv'),
                                  mode='a', index=False, header=False)

            # Write structure information
            for doc in docs:
                if not os.path.exists(os.path.join(seeds_path, doc.material_id + "_POSCAR")):
                    doc.structure.to(os.path.join(seeds_path, doc.material_id + "_POSCAR"))

    # Retrieval symmetry information
    info = pd.read_csv(r'F:\GA_CSP\File\seeds_info.csv')
    info_symmetry = info['symmetry'].str.split('\'', expand=True)[[1, 3, 4, 5]]
    info_spg = info_symmetry[4].str.split('=', expand=True)[1].str.split(expand=True)[0]
    info_sym = pd.concat([info_symmetry, info_spg], axis=1)
    info = pd.concat([info, info_sym], axis=1)
    info = info.rename(columns={0: 'space_group_number', 1: 'crystal_system', 3: 'space_group', 5: 'point_group'})
    info.drop([4, 'symmetry'], axis=1, inplace=True)

    # Featurization - split composition
    info = splitcomposition(species, info)
    info = info.reindex(
        columns=['material_id', 'database_IDs', 'nelements', 'chemsys', 'formula_anonymous', 'formula_pretty',
                 'composition', '{}'.format(species[0]), '{}'.format(species[1]), '{}'.format(species[2]),
                 'crystal_system', 'space_group', 'space_group_number', 'point_group', 'energy_per_atom',
                 'formation_energy_per_atom', 'energy_above_hull', 'is_stable', 'theoretical', 'volume', 'density',
                 'possible_species', 'decomposes_to', 'deprecated', 'origins', 'warnings', 'last_updated'])

    info['composition'] = info['composition'].replace(' ', '', regex=True)

    info.to_csv(os.path.join(r'F:\GA_CSP\File', r'seeds_info.csv'), mode='w', index=False)


def get_structures_from_seeds(seeds_path=None):
    """
            Apply Perturbation to all seeds to generate new structures.

    @type seeds_path: str
    :param seeds_path: Path to the seeds file
    """
    if seeds_path is None:
        seeds_path = r'F:\GA_CSP\File\Seeds'

    if not os.path.exists(seeds_path):
        raise FileNotFoundError('No such path exists, please check it.   {}'.format(seeds_path))

    if not os.listdir(seeds_path):
        raise FileNotFoundError('No seed file exists, please check it.   {}'.format(seeds_path))

    seeds = os.listdir(seeds_path)

    structure_path = r'F:\GA_CSP\File\Structures\from_seeds'
    if not os.path.exists(structure_path):
        os.mkdir(structure_path)

    for seed in seeds:
        mp_id = seed.split('_')[0]
        crystal = pyxtal()
        crystal.from_seed(os.path.join(seeds_path, seed))
        crystal.apply_perturbation()
        crystal.to_ase().write(os.path.join(structure_path, 'from_{}'
                                            .format(mp_id + '_POSCAR')), format='vasp', vasp5=True)


def get_structure_from_seed(species, seed, seeds_path=None):
    """
            Apply Perturbation to one specific seed to generate a new structure.

    :param species: Species of the system
    :param seed: mp-id of the seed
    :param seeds_path: Path to the seeds file
    :return: The new structure after the perturbation, the corresponding structural information
    """
    if seeds_path is None:
        seeds_path = r'F:\GA_CSP\File\Seeds'

    if not os.path.exists(seeds_path):
        raise FileNotFoundError('No such path exists, please check it.   {}'.format(seeds_path))

    seed_name = seed + '_POSCAR'
    seeds = os.listdir(seeds_path)
    if seed_name not in seeds:
        raise FileNotFoundError('No such seed exists, please check it.   {}'.format(seeds_path))

    # structure_path = r'F:\GA_CSP\File\Structures\from_seed'
    # if not os.path.exists(structure_path):
    #     os.mkdir(structure_path)

    seeds_data = pd.read_csv(os.path.join(r'F:\GA_CSP\File', r'seeds_info.csv'))
    seed_data = seeds_data.query("material_id == '{}'".format(seed))
    crystal = pyxtal()
    crystal.from_seed(os.path.join(seeds_path, seed_name))
    # crystal.apply_perturbation()
    structure = crystal.to_pymatgen()

    structure_data = {'nelements': seed_data['nelements'],
                      'crystal_system': seed_data['crystal_system'],
                      'space_group_number': seed_data['space_group_number'],
                      'from': 'Seed-{}-{}'.format(seed, crystal.source),
                      'composition': seed_data['composition'],
                      species[0]: seed_data['{}'.format(species[0])],
                      species[1]: seed_data['{}'.format(species[1])],
                      species[2]: seed_data['{}'.format(species[2])]}
    structure_data = pd.DataFrame(structure_data)
    print('Structure generated with seed: {}'.format(seed))

    return structure, structure_data


if __name__ == '__main__':
    # get_seeds_from_mp(['Ge', 'Sb', 'Te'])
    # get_structures_from_seeds()
    obj = get_structure_from_seed(['Ge', 'Sb', 'Te'], 'mp-32')
    print(obj[1])
