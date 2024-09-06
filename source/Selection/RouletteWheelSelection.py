#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd


def roulette_selection(candidates, scale):
    """
            Select offspring individuals by roulette method.
    @type candidates: DataFrame
    @type scale: float | int
    :param candidates: DataFrame containing individual composition, formation energy and fitness
    :param scale: float: Fraction of candidates to next generation
                  int:   number of candidates to next generation
    :return: Candidates for next generation
    @rtype: DataFrame
    """
    candidates = pd.DataFrame(candidates)

    total_fitness = candidates['fitness'].sum()
    probabilities = candidates['fitness'] / total_fitness

    if type(scale) is float:
        df_selected = candidates.sample(frac=scale, weights=probabilities)
    else:
        df_selected = candidates.sample(n=scale, weights=probabilities)

    return df_selected


if __name__ == '__main__':
    file_name = r'F:\GA_CSP\File\output.csv'
    df = pd.read_csv(file_name)
    selected_candidates = roulette_selection(df, scale=10)
    selected_candidates.to_csv(r'F:\GA_CSP\File\next_gen_by_roulette.csv', index=False)
