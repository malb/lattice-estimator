from typing import Iterable
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

def parse_arg(arg_val: str) -> Iterable:
    if ':' in arg_val:
        # value is a range
        arg_val = [int(i.strip()) for i in arg_val.split(':')]
        if len(arg_val) <= 3:
            return range(*arg_val)
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
        n=input_params[0],
        q=input_params[1],
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
                        help='Value(s) for sigma of error e. This can be\n'
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
                        help='Value(s) for sigma of secret key s. This can be\n'
                        '* list: looks like --sigma_s=25,70,80,100\n'
                        '* range: looks like --sigma_s=25:125:25\n'
                        '* single value: looks like --sigma_s=200.0\n',
                        default='2.0')
    parser.add_argument('-u', '--uniform_s',
                        help='Whether to use a uniform distribution for secret '
                        'key s. If unset, defaults to using a gaussian '
                        'distribution.',
                        action='store_true')
    parser.add_argument(
        '--security_cutoff',
        type=int,
        help='Mark values with higher bits of security than this with a'
        ' different color in the graph',
        default=128)
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
    log_level = 1

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
        print('work: ', work)

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

    print(result_dict)
