#!/usr/bin/python3
"""This file provides a method to generate sample data of time series based on a reference occurrence
dataset. Both files must be 
"""

import pandas as pd
import argparse

PARSER = argparse.ArgumentParser()
PARSER.add_argument('--color',
                    nargs='*',
                    type=str,
                    metavar=['red','green','blue','ir','swir1','swir2'],
                    help='Data band of the time series to extract samples from.')
PARSER.add_argument('--ts_path',
                    nargs=1,
                    type=str,
                    default='./',
                    help='Time series directory path.')
PARSER.add_argument('--ref_path',
                    nargs=1,
                    type=str,
                    default='./Presences_only_train_sample.csv',
                    help='Reference CSV path.')

ARGS = PARSER.parse_args()

def gen_samples(colors, ts_path, ref_path):
    for color in colors:
        gen_sample(color, ts_path, ref_path)

def gen_sample(color:str, ts_path:str = './', ref_path:str = 'Presences_only_train_sample.csv'):
    """Generate data sample of given data band based of a refernce CSV.

    Compares columns 'timeSerieID' in both a reference occurrences CSV and a
    time series data CSV (specific to a databand) to create a samplestr
    Args:
        color (str): dsata band of the time series to extract samples from.
        ts_path (str): path to the directory containing the time series CSV.
        ref_path (str, optional): file path to the ref. CSV. Defaults to 'Presences_only_train_sample.csv'.
    """
    ts_ref = pd.read_csv(ref_path, sep=';')
    ts = pd.read_csv(f'{ts_path}/time_series_{color}.csv', sep=';')
    inter_id = []
    for tsid in ts_ref['timeSerieID']:
        inter_id.append(ts[ts['timeSerieID']==tsid].index[0])
    ts_sample = ts.loc[inter_id]
    ts_sample.to_csv(f'time_series_{color}_sample.csv', sep=';', index=False)
    print(f'Generated sample for time_series_{color} in {ts_path}')

if __name__=='__main__':
    gen_samples(ARGS.color, ARGS.ts_path[0], ARGS.ref_path[0])