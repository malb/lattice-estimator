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

from typing import Iterable, Union
import argparse
from estimator import LWE
from estimator.io import Logging
from estimator.lwe_parameters import LWEParameters
from estimator.nd import NoiseDistribution as ND
import time
import math
import multiprocessing
import numpy as np
import os
from matplotlib import pyplot as plt


def parse_arg(arg_val: str) -> Iterable:
    if ':' in arg_val:
        # value is a range
        arg_val = [float(i.strip()) for i in arg_val.split(':')]
        if len(arg_val) <= 3:
            return [round(i, 4) for i in np.arange(*arg_val)]
        else:
            raise ValueError('argument of type range (containing colons :) '
                'cannot be longer than 3 values.')
    elif ',' in arg_val:
        # value is a list
        arg_val = [float(val) if '.' in val else int(val) 
            for val in arg_val.split(',')]
        return arg_val
    else:
        # value is a single value
        return [float(arg_val) if '.' in arg_val else int(arg_val)]


def GraphResults(result_dict: dict,
                params: dict[str, list[Union[int, float]]],
                file_name: str,
                log_level: int):
    """
    Creates one or more graphs, depending on how many variables are being
    iterated over.

    In the case of a single variable being a range, the output is a simple plot.
    In the case of two variables, a heatmap plot is created.
    In the case of three variables, a collection of heatmap plots are created.
    For four or more variables, an error is raised, as it is not practical to
    convey that information in a usable plot.

    Args:
        result_dict: a dictionary of results, like {(1, 45, 7, 4): 5.7}. The
                     order of the params is: n, q, Xe, Xs
        params: the range of params, like {'n': [1, 43, 67], 'q': [4] ...}

    Returns:
        None
    """
    params_to_idx = {'n': 0, 'q': 1, 'Xe': 2, 'Xs': 3}
    # First, figure out how many variables we need to plot.
    axis_vars = {}

    for param in params:
        if len(params[param]) > 1:
            axis_vars[param] = len(params[param])

    if len(axis_vars) == 0:
        raise ValueError('Cannot plot when there are no variables. '
                         'Call the lattice estimator directly for security.')

    elif len(axis_vars) > 2:
        raise ValueError('Cannot plot more than two variables. '
                         'Try freezing one or more of them.')

    elif len(axis_vars) == 1:
        variable_param = list(axis_vars.keys())[0]  # a string like 'Xe'
        x = sorted(params[variable_param])
        y = sorted(result_dict.items(), key=lambda x: x[params_to_idx[variable_param]])
        y = [_ for _ in zip(*y)][1]
        fig, ax = plt.subplots(figsize=(20, 20), dpi=80)

        ax.plot(x, y)
        ax.set_xlabel(f'Parameter: {variable_param}')
        ax.set_ylabel('Security')
        desc_params = {i:v[0] for i, v in params.items() if len(v) == 1}
        plt.title(f'Security with parameters {desc_params}')

        Logging.log('sweep', log_level, 'Saved the line plot graph to: %s',
                    file_name + '_plot.png')
        fig.savefig(file_name + '_plot.png')

    elif len(axis_vars) == 2:
        axis_vars = sorted(axis_vars.items(), key=lambda x: params_to_idx[x[0]])
        x_param = axis_vars[1][0]
        x = sorted(params[x_param])
        y_param = axis_vars[0][0]
        y = sorted(params[y_param])

        values = [i[1] for i in sorted(result_dict.items(), key=lambda x: x[0])]
        mat = np.flip(np.array(values).reshape(
            (len(set(x)), len(set(y)))).transpose(),
                      axis=0)

        fig, ax = plt.subplots(figsize=(20, 20), dpi=80)

        ax.imshow(mat)
        ax.set_xticks(np.arange(0, len(set(x)), 1))
        ax.set_yticks(np.arange(0, len(set(y)), 1))
        ax.set_xticklabels(sorted(list(set(x))))
        ax.set_yticklabels(sorted(list(set(y)), reverse=True))

        desc_params = {i:v[0] for i, v in params.items() if len(v) == 1}
        plt.title(f'Security with parameters {desc_params}')

        for (j, i), label in np.ndenumerate(mat):
            ax.text(i, j, round(label, 1), ha='center', va='center', color='white')

        ax.set_xlabel(f'Parameter: {x_param}')
        ax.set_ylabel(f'Parameter: {y_param}')

        fig.savefig(file_name + '_gradient.png')
        Logging.log('sweep', log_level, 'Saved the gradient graph to: %s',
                    file_name + '_gradient.png')


def PerformCalculation(input_params: tuple[int, float],
                       results: dict,
                       uniform_s: bool,
                       log_level=None) -> None:
    """
    This function calls the lattice-estimator for a given set of input
    parameters, and appends the output to the `results` dict. 
    """
    Logging.log('sweep', log_level, f'Running for n, q, Xe, Xs = {input_params}')
    Xs_ = ND.DiscreteGaussian if uniform_s else ND.UniformMod
    params = LWE.Parameters(
        n=int(input_params[0]),
        q=int(input_params[1]),
        Xe=ND.DiscreteGaussian(input_params[2]),
        Xs=Xs_(input_params[3]),
        m=input_params[0],
    )
    results[input_params] = security_level(params, f=LWE.estimate.rough)


def security_level(params, f=LWE.estimate, *args, **kwargs):
    res = f(params, *args, **kwargs)
    return min([math.log(res_.rop, 2) for res_ in res.values()])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Iterate over possible options.')
    parser.add_argument('-n',
                        type=str,
                        help='Value(s) for the LWE dimension n. This can be a\n'
                        '* list: looks like -n=25,70,80,100\n'
                        '* range: looks like -n=25:125:25\n'
                        '* single value: looks like -n=200\n',
                        required=True)
    parser.add_argument('-e', '--sigma_e',
                        type=str,
                        help='Value(s) for sigma of error e, as inputs '
                        'to a discrete gaussian distribution (unless specified '
                        'otherwise.) This can be a\n'
                        '* list: looks like --sigma_e=25,70,80,100\n'
                        '* range: looks like --sigma_e=25:125:25\n'
                        '* single value: looks like --sigma_e=200\n',
                        required=True)
    parser.add_argument('-q',
                        type=str,
                        help='Value(s) for the modulus q. This can be a\n'
                        '* list: looks like -q=25,70,80,100\n'
                        '* range: looks like -q=25:125:25\n'
                        '* single value: looks like -q=200\n',
                        default='4294967296', # q=2**32
                        )
    parser.add_argument('-s', '--sigma_s',
                        type=str,
                        help='Value(s) for sigma of secret key s, as inputs '
                        'to a discrete gaussian distribution (unless specified '
                        'otherwise.) This can be a\n'
                        '* list: looks like --sigma_s=25,70,80,100\n'
                        '* range: looks like --sigma_s=25:125:25\n'
                        '* single value: looks like --sigma_s=200.0\n',
                        default='2.0')
    parser.add_argument('-u', '--uniform_s',
                        help='Whether to use a uniform distribution for secret '
                        'key s. If unset, defaults to using a discrete '
                        'gaussian distribution.',
                        action='store_true')
    parser.add_argument('-f', '--fast',
                        help='Use LWE.Estimate.Rough, for a faster estimate. '
                        'If unset, defaults to using LWE.Estimate is slower.',
                        action='store_true')
    parser.add_argument(
        '--output_dir',
        type=str,
        help='The directory to output the resulting pickle and graphs',
        default='')
    parser.add_argument('--num_proc',
                        type=int,
                        help='The number of parallel processes to launch',
                        default=8)
    parser.add_argument(
        '--make_pickle',
        help=
        'Whether to dump a pickle file of the intermediate computation results',
        default=False,
        action='store_true')
    parser.add_argument(
        '--load_pickle',
        type=str,
        help='Path to a pickle to load. Overrides any computation',
        default='')
    args = parser.parse_args()
    log_level = 0

    output_dir = args.output_dir
    if not output_dir:
        output_dir = os.path.dirname(os.path.realpath(__file__))
    assert args.num_proc >= 1, 'need at least one process to execute'
    file_name = time.strftime('%d-%m-%Y_%H-%M-%S')
    file_name = os.path.join(output_dir, file_name)

    if not args.load_pickle:
        result_dict = multiprocessing.Manager().dict()
        work = []
        n = parse_arg(args.n)
        q = parse_arg(args.q)
        sigma_e = parse_arg(args.sigma_e)
        sigma_s = parse_arg(args.sigma_s)

        for n_ in n:
            for q_ in q:
                for sigma_e_ in sigma_e:
                    for sigma_s_ in sigma_s:
                        work.append((n_, q_, sigma_e_, sigma_s_))

        # Parallel process the calculations
        pool = multiprocessing.Pool(processes=min(args.num_proc, len(work)))
        for i in range(len(work)):
            pool.apply_async(PerformCalculation,
                             (work[i], result_dict, args.uniform_s, log_level))
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

    GraphResults(
        result_dict,
        {'n': n, 'q': q, 'Xe': sigma_e, 'Xs': sigma_s},
        file_name,
        log_level)

    print(result_dict)
