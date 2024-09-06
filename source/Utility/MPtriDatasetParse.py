#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ijson
'''
        Extract structure information from large json files.
        Author: JHWang@BHU 2024
        License: CC BY-NC-SA 4.0 International
'''


def extract_json_data(path):

    mp_id_set = []
    frame_id = []

    """
    labels = ['structure', 'uncorrected_total_energy', 'corrected_total_energy', 'energy_per_atom', 'ef_per_atom',
              'e_per_atom_relaxed', 'ef_per_atom_relaxed', 'force', 'stress', 'magmom', 'bandgap']
    structure_labels = ['@module', '@class', 'charge', 'lattice', 'sites']
    lattice_labels = ['matrix', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'volume']
    sites_labels = ['species', 'abc', 'xyz', 'label', 'properties']
    """

    with open(path, 'r') as file:
        parser = ijson.parse(file)
        for prefix, event, value in parser:
            if prefix == '':
                mp_id_set.append(value)

            if mp_id_set:
                if prefix == mp_id_set[-1] and event == 'map_key':
                    frame_id.append(value)

                if frame_id:
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.@module':
                        module = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.@class':
                        class_type = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.charge':
                        charge = value

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.matrix' and event == 'start_array':
                        matrix = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.matrix.item' and event == 'start_array':
                        vector = []
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.matrix.item.item':
                        vector.append(value)
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.matrix.item' and event == 'end_array':
                        matrix.append(vector)

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.a':
                        a = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.b':
                        b = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.c':
                        c = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.alpha':
                        alpha = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.beta':
                        beta = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.gamma':
                        gamma = value
                    elif prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.lattice.volume':
                        volume = value

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites' and event == 'start_array':
                        elements = []
                        label = []
                        occu = []
                        abcs = []
                        xyzs = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.species.item.element':
                        elements.append(value)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.species.item.occu':
                        occu.append(value)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.abc' and event == 'start_array':
                        abc = []
                        xyz = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.abc.item':
                        abc.append(value)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.abc' and event == 'end_array':
                        abcs.append(abc)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.xyz.item':
                        xyz.append(value)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.xyz' and event == 'end_array':
                        xyzs.append(xyz)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.structure.sites.item.label':
                        label.append(value)

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.uncorrected_total_energy':
                        uncorrected_total_energy = value
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.corrected_total_energy':
                        corrected_total_energy = value
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.energy_per_atom':
                        energy_per_atom = value
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.ef_per_atom':
                        ef_per_atom = value
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.e_per_atom_relaxed':
                        e_per_atom_relaxed = value
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.ef_per_atom_relaxed':
                        ef_per_atom_relaxed = value

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.force' and event == 'start_array':
                        forces = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.force.item' and event == 'start_array':
                        force = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.force.item.item':
                        force.append(value)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.force.item' and event == 'end_array':
                        forces.append(value)

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.stress' and event == 'start_array':
                        stresses = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.stress.item' and event == 'start_array':
                        stress = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.stress.item.item':
                        stress.append(value)
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.stress.item' and event == 'end_array':
                        stresses.append(stress)

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.magmom' and event == 'start_array':
                        magmom = []
                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.magmom.item':
                        magmom.append(value)

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.bandgap':
                        bandgap = value

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] + '.mp_id':
                        mp_id = value

                    if prefix == mp_id_set[-1] + '.' + frame_id[-1] and event == 'end_map':
                        info = {
                            "mp_id": mp_id,
                            "frame_id": frame_id[-1],
                            "structure": {
                                "@module": module,
                                "@class": class_type,
                                "charge": charge,
                                "lattice": {
                                    "matrix":matrix,
                                    "a": a,
                                    "b": b,
                                    "c": c,
                                    "alpha": alpha,
                                    "beta": beta,
                                    "gamma": gamma,
                                    "volume": volume
                                },
                                "sites": {
                                    "species": elements,
                                    "occu": occu,
                                    "abc": abcs,
                                    "xyz": xyzs
                                }
                            },
                            "uncorrected_total_energy": uncorrected_total_energy,
                            "corrected_total_energy": corrected_total_energy,
                            "energy_per_atom": energy_per_atom,
                            "ef_per_atom": ef_per_atom,
                            "e_per_atom_relaxed": e_per_atom_relaxed,
                            "ef_per_atom_relaxed": ef_per_atom_relaxed,
                            "force": forces,
                            "stress": stresses,
                            "magmom": magmom,
                            "bandgap": bandgap
                        }

                        yield info


if __name__ == '__main__':
    file_path = r'F:\Dataset\MPtrj_2022.9_full.json'
    for i in extract_json_data(file_path):
        print(i)
