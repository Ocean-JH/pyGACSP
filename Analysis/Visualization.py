#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import ConvexHull


def comp_triangle(species, df):
    """
            Convert the composition to ternary coordinates.
    :param species: Species of the system
    :param df: Dataframe for population information
    :return: coordinates (x, y)
    """
    composition = df[[species[0], species[1], species[2]]]
    composition['sum'] = composition.sum(axis=1)
    comp_frac_0 = composition[species[0]] / composition['sum']
    comp_frac_1 = composition[species[1]] / composition['sum']
    df['X1'] = comp_frac_0 + comp_frac_1 / 2
    df['X2'] = comp_frac_1 * np.sqrt(3) / 2
    df['Y'] = df['formation_energy_per_atom']
    return df


def plot_convex_hull(df):
    data = df[['X1', 'X2', 'Y']]
    # data = data[~(data[['Y']] > 0.2).any(axis=1)]
    x = data['X1']
    y = data['X2']
    z = data['Y']

    hull_data = df[['X1', 'X2', 'Y']]
    # hull_data = hull_data[~(hull_data[['Y']] > 0).any(axis=1)]
    ele = [0, 0, 0, 1, 0, 0, 0.5, 0.866, 0]
    ele = pd.DataFrame(np.array(ele).reshape(3, 3), columns=['X1', 'X2', 'Y'])
    hull_data = pd.concat([hull_data, ele])
    hull_data = np.array(hull_data)
    hull = ConvexHull(hull_data)

    simplices = hull.simplices
    # simplices = np.delete(simplices, 1, axis=0)

    plt.figure(figsize=(16, 9))
    # print(plt.style.available)
    plt.style.use('seaborn-colorblind')
    ax = plt.axes(projection='3d')
    ax.scatter(x, y, z)

    for i in simplices:
        hull_data[i[1], 1] += 0.0001
        ax.plot_trisurf(hull_data[i, 0], hull_data[i, 1], hull_data[i, 2],
                        color=(0.2745, 0.5098, 0.7059, 0), alpha=0.5)

    ax.set_title('Convex Hull Diagram', fontsize=20)
    # ax.set_zlabel('Formation Energy', fontsize=16)

    ax.grid(None)
    ax.axis('off')
    ax.set(xticklabels=[],
           yticklabels=[],
           zticklabels=[])
    ax.view_init(elev=10, azim=-95)

    # plt.savefig('ConvexHullDiagram.png', dpi=600)
    plt.show()


if __name__ == '__main__':
    data = pd.read_csv(r'F:\GA_CSP\File\Population_info\0_relaxed_pop-info.csv')
    data = comp_triangle(['Ge', 'Sb', 'Te'], data)
    plot_convex_hull(data)
    # print(data)
