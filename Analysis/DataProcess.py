#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil

import numpy as np
import pandas as pd

from scipy.spatial import ConvexHull


def merge_pop_info(path=None):
    info_dir = r'F:\GA_CSP\File\Population_info'
    filename = 'Population_info.csv'
    if path is None:
        path = r'F:\GA_CSP\File'

    file_list = os.listdir(info_dir)
    file_relaxed_list = []

    for file in file_list:
        if file.endswith('relaxed_pop-info.csv'):
            file_relaxed_list.append(file)

    df = pd.read_csv(os.path.join(info_dir, file_relaxed_list[0]))

    df.to_csv(os.path.join(path, filename), index=False)

    for i in range(1, len(file_relaxed_list)):
        df = pd.read_csv(os.path.join(info_dir, file_relaxed_list[i]))
        df.to_csv(os.path.join(path, filename), index=False, header=False, mode='a+')

    df = pd.read_csv(os.path.join(path, filename))
    df = df.sort_values(by=['fitness'], ascending=False, ignore_index=True)
    df.to_csv(os.path.join(path, filename), index=False)


def comp_triangle(composition):
    comp_frac = composition / sum(composition)
    return comp_frac[0] + comp_frac[1] / 2, comp_frac[1] * np.sqrt(3) / 2


def calc_convex_hull(pop_info=None):
    if pop_info is None:
        pop_info = pd.read_csv(r'F:\GA_CSP\File\Population_info.csv')

    pts = np.array(pop_info)

    '''
                Coordinate transformation
    '''
    pts_Tri = []
    for pt in pts:
        x, y = comp_triangle(pt[6:9])
        pts_Tri.append([x, y, pt[10]])
    pts_Tri = np.array(pts_Tri)

    '''
                Convex hull determination
    '''
    hull = ConvexHull(pts_Tri)

    '''
                Set of convex hull vertices
    '''
    hull_vts_all = []
    for pt in hull.vertices:
        vertex = pts[pt]
        hull_vts_all.append(vertex)
    hull_vts_all = pd.DataFrame(hull_vts_all)

    hull_vts = hull_vts_all[~(hull_vts_all[[10]] > 0).any(axis=1)]

    # '''
    #             Removal of elemental components
    # '''
    # hull_vts = hull_vts[hull_vts[0] != 1]
    # hull_vts = hull_vts[hull_vts[1] != 1]
    # hull_vts = hull_vts[hull_vts[2] != 1]

    '''
                Addition of stable elemental components
    '''
    ele = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    ele = pd.DataFrame(np.array(ele).reshape(3, 11))
    hull_vts = pd.concat([hull_vts, ele])
    hull_vts = hull_vts.to_numpy()

    '''
                Vertex coordinate transformation
    '''
    hull_vts_Tri = []
    for pt in hull_vts:
        x, y = comp_triangle(pt[6:9])
        hull_vts_Tri.append([x, y, pt[10]])
    hull_vts_Tri = np.array(hull_vts_Tri)
    '''
                Hyperplane equation determination
    '''
    HULL = ConvexHull(hull_vts_Tri)
    HULL_Eqs = pd.DataFrame(HULL.equations)
    HULL_Eqs = np.array(HULL_Eqs)
    # HULL_Eqs = np.unique(HULL_Eqs, axis=0)
    # HULL_Eqs = np.delete(HULL_Eqs, [0], axis=0)  # Remember to delete plane c=0 !!!
    mask = np.abs(HULL_Eqs[:, 2]) < 0.000001
    HULL_Eqs = np.delete(HULL_Eqs, np.where(mask)[0], axis=0)

    '''
                Calculation of E_hull
    '''
    a, b, c, d, x, y = [], [], [], [], [], []

    for i in range(len(HULL_Eqs)):
        a_i = HULL_Eqs[i][0]
        a.append(a_i)
        b_i = HULL_Eqs[i][1]
        b.append(b_i)
        c_i = HULL_Eqs[i][2]
        c.append(c_i)
        d_i = HULL_Eqs[i][3]
        d.append(d_i)

    for j in range(len(pts_Tri)):
        x_j = pts_Tri[j][0]
        x.append(x_j)
        y_j = pts_Tri[j][1]
        y.append(y_j)

    E_on_hull_set = []
    for i in range(len(pts_Tri)):
        for j in range(len(HULL_Eqs)):
            e_on_hull = -(x[i] * a[j] + y[i] * b[j] + d[j]) / c[j]
            E_on_hull_set.append(e_on_hull)
    E_on_hull_set = pd.DataFrame(np.array(E_on_hull_set).reshape(len(pts_Tri), len(HULL_Eqs)))

    E_on_hull_set[E_on_hull_set >= 0] = -100  # Ignore point where convex hull energy >= 0
    E_on_hull = E_on_hull_set.max(axis=1)  # Select the maximum value as convex hull energy

    Ef = pop_info['formation_energy_per_atom']
    Ehull_set = []
    for i in range(len(Ef)):
        Ehull = Ef[i] - E_on_hull[i]  # Calculation of E_hull
        Ehull_set.append(Ehull)

    '''
                File output
    '''
    pop_info['Ehull'] = Ehull_set
    pop_info = pd.DataFrame(pop_info)
    pop_info.to_csv(r'F:\GA_CSP\File\pop_info_Ehull.csv', index=False)


def get_candidates(pop_info=None, criterion='formation_energy_per_atom', cutoff=None):
    candidates_dir = r'F:\GA_CSP\File\Candidates'
    if not os.path.exists(candidates_dir):
        os.mkdir(candidates_dir)

    if pop_info is None:
        pop_info = pd.read_csv(r'F:\GA_CSP\File\pop_info_Ehull.csv')

    if criterion == 'formation_energy_per_atom':
        if cutoff is None:
            cutoff = 0
        if cutoff > 0:
            raise ValueError(
                'The criterion of formation energy per atom must < 0, but current cutoff is {}'.format(cutoff))
        else:
            pop_info = pop_info.sort_values(by='formation_energy_per_atom', ascending=False, ignore_index=True)
    elif criterion == 'fitness':
        if cutoff is None:
            cutoff = 0.1
        if cutoff < 0:
            raise ValueError('The criterion of fitness must > 0, but current cutoff is {}'.format(cutoff))
        else:
            pop_info = pop_info.sort_values(by='fitness', ascending=False, ignore_index=True)
    else:
        raise ValueError(
            "The criterion of candidates must be either 'formation_energy_per_atom' or 'fitness'! Not {}".format(
                criterion))

    pop_info = pop_info[pop_info['{}'.format(criterion)] < cutoff]
    for _, row in pop_info.iterrows():
        gen_num = row['ID'].split('-')[0]
        pop_num = row['ID'].split('-')[1]
        dir_list = os.listdir(os.path.join(r'F:\GA_CSP\File\Structures', f'{gen_num}_population_relaxed'))
        for dir in dir_list:
            if dir.startswith('{}-{}_'.format(gen_num, pop_num)):
                shutil.copyfile(os.path.join(r'F:\GA_CSP\File\Structures\{}_population_relaxed'.format(gen_num), dir),
                                os.path.join(candidates_dir, dir))


if __name__ == '__main__':
    # merge_pop_info()
    # calc_convex_hull()
    get_candidates()
