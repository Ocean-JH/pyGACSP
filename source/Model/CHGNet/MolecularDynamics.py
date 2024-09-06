#!/usr/bin/env python
# -*- coding: utf-8 -*-
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import MolecularDynamics
from pymatgen.core import Structure

from ase.io.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor
from chgnet.utils import solve_charge_by_mag

import warnings


warnings.filterwarnings("ignore", module="pymatgen")
warnings.filterwarnings("ignore", module="ase")

structure = Structure.from_file("examples/mp-18767-LiMnO2.cif")
chgnet = CHGNet.load()

md = MolecularDynamics(
    atoms=structure,
    model=chgnet,
    ensemble="nvt",
    temperature=1000,  # in K
    timestep=2,  # in femto-seconds
    trajectory="md_out.traj",
    logfile="md_out.log",
    loginterval=100,
)
md.run(50)  # run a 0.1 ps MD simulation

# visualization
traj = Trajectory("md_out.traj")
mag = traj[-1].get_magnetic_moments()

# get the non-charge-decorated structure
structure = AseAtomsAdaptor.get_structure(traj[-1])
print(structure)

# get the charge-decorated structure
struct_with_chg = solve_charge_by_mag(structure)
print(struct_with_chg)
