#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd


def elite_selection(candidates, scale):
    """
            Select offspring individuals by elite selection.
    @type candidates: DataFrame
    @type scale: float | int
    :param candidates: DataFrame containing individual composition, formation energy and fitness
    :param scale: float: Fraction of candidates to next generation
                  int:   number of candidates to next generation
    :return: Candidates for next generation
    @rtype: DataFrame
    """
    candidates = pd.DataFrame(candidates)

    if type(scale) is int:
        num_of_elites = scale
    else:
        num_of_elites = int(scale * len(candidates))

    elites = candidates.nlargest(num_of_elites, 'fitness')
    return elites


if __name__ == '__main__':
    file_name = r'..\..\file\output.csv'
    df = pd.read_csv(file_name)
    selected_candidates = elite_selection(df, scale=10)
    selected_candidates.to_csv(r'..\..\file\next_gen_by_elite.csv', index=False)
