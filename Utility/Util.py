#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from ase.io import read, write


def convert_format(path_in: str, path_out: str, out_format: str):
    os.listdir(path_in)

    for filename in os.listdir(path_in):
        new_filename = filename.rstrip('-POSCAR') + '.' + out_format
        structure = read(os.path.join(path_in, filename))
        write(os.path.join(path_out, new_filename), structure, format=out_format)


if __name__ == '__main__':
    path_in = r'F:\GA_CSP\File\Candidates'
    path_out = r''
    convert_format(path_in, path_out, 'cif')
