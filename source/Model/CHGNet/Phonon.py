#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from jobflow import SETTINGS
from jobflow import run_locally

from atomate2.forcefields.jobs import CHGNetRelaxMaker, CHGNetStaticMaker
from atomate2.forcefields.flows.phonons import PhononMaker

from pymatgen.core.structure import Structure
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos
from pymatgen.phonon.plotter import PhononBSPlotter, PhononDosPlotter


def plot_phonon_bs(candidates_path, phonon_path):
    if not os.path.exists(phonon_path):
        os.makedirs(phonon_path)

    candidates = os.listdir(candidates_path)
    phonon_calculator = PhononMaker(min_length=15.0, store_force_constants=False,
                                    bulk_relax_maker=CHGNetRelaxMaker(steps=1000,
                                                                      relax_kwargs={"fmax": 0.001, "verbose": False}),
                                    phonon_displacement_maker=CHGNetStaticMaker(),
                                    static_energy_maker=None)
    for candidate in candidates:
        indiv_name = candidate.rstrip('-POSCAR')
        cur_list = os.listdir(phonon_path)
        cur_name = f'BS-{indiv_name}.png'
        if cur_name in cur_list:
            continue
        structure = Structure.from_file(os.path.join(candidates_path, candidate))
        phonon_flow = phonon_calculator.make(structure=structure)
        run_locally(phonon_flow, log=False, create_folders=False, root_dir=phonon_path, store=SETTINGS.JOB_STORE)
        store = SETTINGS.JOB_STORE
        store.connect()

        result = store.query_one(
            {"name": "generate_frequencies_eigenvectors"},
            properties=[
                # "output.phonon_dos",
                "output.phonon_bandstructure",
            ],
            load=True,
            sort={"completed_at": -1}  # to get the latest computation
        )

        ph_bs = PhononBandStructureSymmLine.from_dict(result['output']['phonon_bandstructure'])
        # ph_dos = PhononDos.from_dict(result['output']['phonon_dos'])

        bs_plot = PhononBSPlotter(bs=ph_bs)
        bs_plot.save_plot(os.path.join(phonon_path, f'BS-{indiv_name}.png'), units='cm^-1')

        # dos_plot = PhononDosPlotter()
        # dos_plot.add_dos(label='a', dos=ph_dos)
        # dos_plot.save_plot(f'DOS-{indiv_name}.png', units='cm^-1')

        store = None
        phonon_flow = None


if __name__ == '__main__':
    in_path = r'Candidates\optimized_structure\USPEX'
    out_path = r'Candidates\optimized_structure\USPEX_Phonon'
    plot_phonon_bs(in_path, out_path)
