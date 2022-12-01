"""
This module provides a function to iterate over multiple inputs to the
Lattice Estimator, and graph the resulting bits of security of the associated
lattice instances.

Required arguments:
* nrange: [minimum n, inclusive]:[maximum n, exclusive]:[increment interval]
  Inputs must all be integers.
* xerange: [minimum xe, inclusive]:[maximum xe, exclusive]:[increment interval]
  Inputs can be floats or integers, they will be parsed as floats.
  Xe is the log_2 of the standard deviation of the error.

Example command for a simple sweep from n=600 (inclusive) to 700 (exclusive),
by increments of 20, and a log_2 of standard deviation or error from xe=6
(inclusive) to 7 (exclusive), by increments of 0.2:
  sage --python param_sweep.py 600:700:20 6:7:0.2

Optional flags that can be passed in:
* --security_cutoff: the number of bits of security of interest - values at or
  above this cutoff will be colored differently in the cutoff graph.
* --output_dir: the directory to output the resulting graphs and pickle files to.
  Defaults to the lattice-esimator directory if not specified.
* --num_proc: the number of parallel processes to launch while performing
  computations using the lattice-estimator. Defaults to 8.
* --make_pickle: append this flag if you want a pickle file of the intermediate
  computation results. This is useful if you want to generate a different graph
  from the same data, without having to recompute the data.
* --load_pickle: generate a graph from an existing pickle file, instead of
  computing the data from scratch.

Example command using many of the flags, generating a pickle file:
  sage --python param_sweep.py 600:1140:20 6:20:1 --security_cutoff=120 \
    --output_dir=param_sweep/ --num_proc=4 --make_pickle

Example command to generate graphs from a pre-existing pickle file:
  sage --python param_sweep.py 600:1140:20 6:20:1 --load_pickle=path/to/file

Note: running the parameter sweep can take minutes to hours, depending on the
number of parameter combinations (the larger the ranges and the smaller the 
increment intervals, the more parameter combinations will be computed).
"""

import argparse
from estimator import LWE
from estimator.io import Logging
from estimator.lwe_parameters import LWEParameters
from estimator.nd import NoiseDistribution
from estimator.schemes import TFHE630
from matplotlib import pyplot
import pickle
import time
import math
import multiprocessing
import numpy as np
import os
import re

attack_keys = [
    'arora-gb', 'bkw', 'usvp', 'bdd', 'bdd_hybrid', 'bdd_mitm_hybrid', 'dual',
    'dual_hybrid', 'dual_mitm_hybrid'
]


def PerformCalculation(input_params: tuple[int, float],
                       results: dict,
                       log_level=None) -> None:
"""
Call the lattice-estimator for a given set of input parameters, and append the
output to the `results` dict. 
The `input_params` field looks like (n, xe).
The q and Xs parameters are hardcoded for common LWE settings, but can be
modified to fit other settings.
"""
    n, xe = input_params
    Logging.log('sweep', log_level, f'Running for {n=}, {xe=}')
    LWEParams = LWEParameters(
        n=n,
        q=2**32,
        Xs=NoiseDistribution.UniformMod(2),
        Xe=NoiseDistribution.DiscreteGaussian(stddev=2**xe),
        tag=f'N{n}_Xe{xe}',
    )
    attacks = LWE.estimate(LWEParams)
    results[(n, xe)] = min(
        [math.log(attacks[i].rop, 2) for i in attack_keys if i in attacks])
    Logging.log('sweep', log_level, f'Finished running {n=}, {xe=}')


def MakeGraph(result_dict: dict,
              file_name: str,
              bit_cutoff: int,
              log_level=None) -> None:
    x = list(result_dict.items())
    x = sorted(x)
    x, values = zip(*x)

    x, y = zip(*x)
    mat = np.flip(np.array(values).reshape(
        (len(set(x)), len(set(y)))).transpose(),
                  axis=0)
    binary_mat = mat >= bit_cutoff

    fig, ax = pyplot.subplots(figsize=(20, 20), dpi=80)

    ax.imshow(mat)
    ax.set_xticks(np.arange(0, len(set(x)), 1))
    ax.set_yticks(np.arange(0, len(set(y)), 1))
    ax.set_xticklabels(sorted(list(set(x))))
    ax.set_yticklabels(sorted(list(set(y)), reverse=True))

    pyplot.title('Minimum bit length')

    fig2, ax2 = pyplot.subplots(figsize=(20, 20), dpi=80)
    ax2.imshow(binary_mat)
    ax2.set_xticks(np.arange(0, len(set(x)), 1))
    ax2.set_yticks(np.arange(0, len(set(y)), 1))
    ax2.set_xticklabels(sorted(list(set(x))))
    ax2.set_yticklabels(sorted(list(set(y)), reverse=True))

    for (j, i), label in np.ndenumerate(mat):
        ax.text(i, j, round(label, 1), ha='center', va='center', color='white')
        ax2.text(i,
                 j,
                 round(label, 1),
                 ha='center',
                 va='center',
                 color='white' if label < bit_cutoff else 'black')

    pyplot.title(f'Minimum bit length with cutoff at {bit_cutoff}')

    ax.set_xlabel('LWE length')
    ax2.set_xlabel('LWE length')

    ax.set_ylabel('log2(sigma) of error')
    ax2.set_ylabel('log2(sigma) of error')

    fig.savefig(file_name + '_gradient.png')
    Logging.log('sweep', log_level, 'Saved the gradient graph to: %s',
                file_name + '_gradient.png')

    fig2.savefig(file_name + '_cutoff.png')
    Logging.log('sweep', log_level, 'Saved the cutoff graph to: %s',
                file_name + '_cutoff.png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Iterate over possible options.')
    parser.add_argument('nrange',
                        type=str,
                        help='range for N, like 500:600:1',
                        default="500:600:1")
    parser.add_argument(
        'xerange',
        type=str,
        help='range for Xe: log(sigma), like 10:20:1. Logarithmic powers',
        default="10:20:1")
    parser.add_argument(
        '--security_cutoff',
        type=int,
        help='Mark values with higher bits of security than this with a'
        ' different color in the graph',
        default=128,
        required=False)
    parser.add_argument(
        '--output_dir',
        type=str,
        help='The directory to output the resulting pickle and graphs',
        default='',
        required=False)
    parser.add_argument('--num_proc',
                        type=int,
                        help='The number of parallel processes to launch',
                        default=8,
                        required=False)
    parser.add_argument(
        '--make_pickle',
        help=
        'Whether to dump a pickle file of the intermediate computation results',
        default=False,
        action='store_true',
        required=False)
    parser.add_argument(
        '--load_pickle',
        type=str,
        help='Path to a pickle to load. Overrides any computation',
        default='',
        required=False)
    args = parser.parse_args()
    log_level = 1

    # Validate the input arguments
    nrange = [int(i) for i in args.nrange.split(':')]
    xerange = [float(i) for i in args.xerange.split(':')]
    assert len(
        nrange) == 3, 'n range needs to be 3 integers separated by a colon'
    assert len(
        xerange) == 3, 'Xe range needs to be 3 floats separated by a colon'
    output_dir = args.output_dir
    if not output_dir:
        output_dir = os.path.dirname(os.path.realpath(__file__))
    assert args.num_proc >= 1, 'need at least one process to execute'
    assert args.security_cutoff >= 1, 'needs at least one bit of security'

    file_name = time.strftime('%d-%m-%Y_%H-%M-%S')
    file_name = file_name + '_n' + re.sub(
        ':', '-', args.nrange) + '_xe' + re.sub(':', '-', args.xerange)
    file_name = os.path.join(output_dir, file_name)

    if not args.load_pickle:
        # Create the permutations of n and Xe
        result_dict = multiprocessing.Manager().dict()
        work = []
        for n in range(*nrange):
            for xe in np.arange(*xerange):
                work.append((n, xe))

        # Parallel process the calculations
        pool = multiprocessing.Pool(processes=args.num_proc)
        for i in range(len(work)):
            pool.apply_async(PerformCalculation,
                             (work[i], result_dict, log_level))
        pool.close()
        pool.join()

        # Possibly pickle the intermediate computation results
        if args.make_pickle:
            pickle.dump(dict(result_dict),
                        open(file_name + '_result.pickle', 'wb'))
            Logging.log('sweep', log_level,
                        'Pickled the intermediate computations to: %s',
                        file_name + '_result.pickle')
    else:
        result_dict = pickle.load(open(args.load_pickle, 'rb'))

    # Make the graphs of the permutations
    MakeGraph(result_dict, file_name, args.security_cutoff, log_level)
