#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pickle
from chgnet.model.model import CHGNet
from pymatgen.io.vasp.outputs import Oszicar, Vasprun
from chgnet.data.dataset import StructureData, get_train_val_test_loader
from chgnet.trainer import Trainer


"""
        Fine-tune CHGNet model to adapt to Specific(Ge-Sb-Te) system.
"""


def check_convergence(calc_dir):
    file = os.path.join(calc_dir, 'OUTCAR')
    with open(file, 'r') as f:
        content = f.read()
        return 'reached required accuracy' in content


def get_calc_dir(root_dir):
    """
            Get the calculation directory for each structure.
    :param root_dir: Path to dataset
    :return: List of directory names of structures.
    """
    if not os.path.exists(root_dir):
        raise FileNotFoundError("No such directory: " + root_dir)

    # Initialize the calculation directory
    calc_dirs = []
    for root, _, _ in os.walk(root_dir):
        if os.path.exists(os.path.join(root, 'OUTCAR')):
            calc_dirs.append(root)
            # if check_convergence(root):
            #     calc_dirs.append(root)
            # else:
            #     warnings.warn(f"The calculations in {root} did NOT converge!")
        else:
            print(f"Can't find the required files('OSZICAR' or 'vasprun.xml') at {root}, please check it.")

    return calc_dirs


def prepare_data(calc_dirs):
    """
            Prepare training data from ab initio calculation.
    @type calc_dirs: list
    :param dataset_dir: Structure set directory for model fine-tuning
    :return: Structure and the corresponding energy per atom, force and stress
    """
    # Create lists for data storage
    structures, energies_per_atom, forces, stresses = [], [], [], []
    # magmoms = []

    xml = 'vasprun.xml'
    oszicar = 'OSZICAR'
    for calc_dir in calc_dirs:
        try:
            vasprun_file = os.path.join(calc_dir, xml)
            dataset = Vasprun(vasprun_file)

            structures.append(dataset.ionic_steps[-1]['structure'])

            n_atoms = len(dataset.ionic_steps[0]["structure"])
            energy = Oszicar(os.path.join(calc_dir, oszicar)).final_energy
            energy_per_atom = energy / n_atoms
            energies_per_atom.append(energy_per_atom)

            forces.append(dataset.ionic_steps[-1]['forces'])
            stresses.append(dataset.ionic_steps[-1]['stress'])
        except Exception as error:
            print(f"Error occurred in {calc_dir}: {error}")
            continue

    return structures, energies_per_atom, forces, stresses


def dataset_split(structures, energies_per_atom, forces, stresses):
    """
            Split the dataset into training sets, test sets and validation sets.

    :param structures: List of structures from ab initio calculation
    :param energies_per_atom: List of energies per atom from ab initio calculation
    :param forces: List of forces from ab initio calculation
    :param stresses: List of stresses from ab initio calculation
    :return: train_loader, val_loader, test_loader
    """
    dataset = StructureData(
        structures=structures,
        energies=energies_per_atom,
        forces=forces,
        stresses=stresses,
        # magmoms=magmoms,
    )
    train_loader, val_loader, test_loader = get_train_val_test_loader(
        dataset, batch_size=16, train_ratio=0.9, val_ratio=0.05
    )

    return train_loader, val_loader, test_loader


if __name__ == '__main__':
    """
            Prepare training data from ab initio calculation.
    """
    # dataset_dir = r"F:\Research_Files\Structure_Search\USPEX\30_30_pop\DataAnalysis\30_30_SCF"
    # dataset_dir = r"F:\Dataset\ElementSubstitution"
    # calc_dirs = get_calc_dir(dataset_dir)
    # with open(r'F:\GA_CSP\File\calc_dir.pkl', 'wb') as f:
    #     pickle.dump(calc_dirs, f)

    # with open(r'F:\GA_CSP\File\calc_dir.pkl', 'rb') as f:
    #     calc_dirs = pickle.load(f)
    #
    # structures, energies_per_atom, forces, stresses = prepare_data(calc_dirs)
    # dataset = (structures, energies_per_atom, forces, stresses)
    # with open(r'F:\GA_CSP\File\elesub_dataset.pkl', 'wb') as f:
    #     pickle.dump(dataset, f)

    with open(r'F:\GA_CSP\File\elesub_dataset.pkl', 'rb') as f:
        dataset = pickle.load(f)

    structures, energies_per_atom, forces, stresses = dataset

    """
            Split the dataset into training sets, test sets and validation sets.
    """
    train_loader, val_loader, test_loader = dataset_split(structures, energies_per_atom, forces, stresses)

    """
            Specify the trainer.
    """
    # Load initial CHGNet model
    chgnet = CHGNet.load()

    # Optionally fix the weights of some layers
    for layer in [
        chgnet.atom_embedding,
        chgnet.bond_embedding,
        chgnet.angle_embedding,
        chgnet.bond_basis_expansion,
        chgnet.angle_basis_expansion,
        chgnet.atom_conv_layers[:-1],
        chgnet.bond_conv_layers,
        chgnet.angle_layers,
    ]:
        for param in layer.parameters():
            param.requires_grad = False

    # Specify Trainer
    trainer = Trainer(
        model=chgnet,
        targets="efs",
        optimizer="Adam",
        scheduler="CosLR",
        criterion="MSE",
        epochs=50,
        learning_rate=1e-3,
        use_device="cpu",
        print_freq=100,
    )

    """
            Start training.
    """
    trainer.train(train_loader, val_loader, test_loader, save_dir=r'F:\GA_CSP\File\model')

    model = trainer.model
    best_model = trainer.best_model  # best model based on validation energy MAE
