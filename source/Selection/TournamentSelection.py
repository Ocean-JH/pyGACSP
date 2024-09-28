#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd


def tournament_selection(candidates, scale, tournament_size=5, tournament_rounds=5):
    """
            Select offspring individuals by tournament method.

    @type candidates: DataFrame
    @type scale: float | int
    @type tournament_size: int
    @type tournament_rounds: int
    :param candidates: DataFrame containing individual composition, formation energy and fitness
    :param scale: float: Fraction of candidates to next generation
                  int:   number of candidates to next generation
    :param tournament_size: Size of tournament
    :param tournament_rounds: number of tournament rounds
    :return: Candidates for next generation
    @rtype: DataFrame
    """
    candidates = pd.DataFrame(candidates)

    if type(scale) == float:
        num_of_winner = int(scale * len(candidates))
    else:
        num_of_winner = scale

    selected_candidates = pd.DataFrame()

    for _ in range(num_of_winner):
        candidate = None

        for _ in range(tournament_rounds):
            tournament_individuals = candidates.sample(n=tournament_size)
            tournament_round_winner = tournament_individuals.loc[tournament_individuals['fitness'].idxmax()]

            if candidate is None or tournament_round_winner['fitness'] > candidate['fitness']:
                candidate = tournament_round_winner
                tournament_winner = tournament_individuals.loc[[tournament_individuals['fitness'].idxmax()]]
                candidates = candidates.drop(tournament_winner['ID'])

        selected_candidates = pd.concat([selected_candidates, tournament_winner])

    selected_candidates = selected_candidates.sort_values(by='fitness', ascending=False)

    return selected_candidates


if __name__ == '__main__':
    file_name = r'..\..\file\output.csv'
    df = pd.read_csv(file_name)
    selected_candidates = tournament_selection(df, scale=10)
    selected_candidates.to_csv(r'..\..\file\next_gen_by_tournament.csv', index=False)
