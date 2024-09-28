#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd


def splitcomposition(species, info):
    data = pd.DataFrame(info)
    data_comp = data['composition']
    composition = pd.DataFrame(0, index=range(len(data)), columns=species)
    data = pd.concat([data, composition], axis=1)
    comp_elements = data_comp.str.findall(r'([A-Za-z]+)')
    comp_numbers = data_comp.str.findall(r'(\d+)')
    for index, row in data.iterrows():
        for ele, num in zip(comp_elements[index], comp_numbers[index]):
            row['{}'.format(ele)] = int(num)
            data.iloc[index] = row
    return data


if __name__ == '__main__':
    info = pd.read_csv(r'..\..\file\seeds_info.csv')
    splitcomp = splitcomposition(['Ge', 'Sb', 'Te'], info)
    print(splitcomp)
