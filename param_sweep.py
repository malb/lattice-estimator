"""
This module provides a function to iterate over multiple inputs to the
Lattice Estimator, and graph the resulting bits of security of the associated
lattice instances.

Note: running the parameter sweep can take minutes to hours, depending on the
number of parameter combinations (the larger the ranges and the smaller the 
increment intervals, the more parameter combinations will be computed). For
faster and rougher results, use the `LWE.estimate.rough` function.
"""

from estimator import LWE
from estimator.io import Logging
from estimator.lwe_parameters import LWEParameters
from estimator.nd import NoiseDistribution as ND
from matplotlib import pyplot as plt
from typing import Iterable, Union
import argparse
import pickle
import time
import math
import multiprocessing
import numpy as np
import os


class ParameterSweep:
    def parameter_sweep(
        n: int,
        q: int,
        e: float,
        s: float,
        m=None,
        Xe=ND.DiscreteGaussian,
        e_log=True,
        Xs=ND.DiscreteGaussian,
        s_log=True,
        tag=None,
        f=LWE.estimate,
        num_proc=8,
        log_level=0,
    ) -> dict:
        try:
            iter(n)
        except TypeError:
            n = [n]
        try:
            iter(q)
        except TypeError:
            q = [q]
        try:
            iter(m)
        except TypeError:
            m = [m]
        try:
            iter(e)
        except TypeError:
            e = [e]
        try:
            iter(s)
        except TypeError:
            s = [s]

        result_dict = multiprocessing.Manager().dict()
        work = []

        for n_ in n:
            for q_ in q:
                for e_ in e:
                    for s_ in s:
                        for m_ in m:
                            if m_ is None:
                                m_ = n_
                            work.append((n_, q_, e_, s_, m_))

        # Parallel process the calculations
        pool = multiprocessing.Pool(processes=min(num_proc, len(work)))
        for i in range(len(work)):
            pool.apply_async(ParameterSweep.security_level,
                (work[i], result_dict, Xe, e_log, Xs, s_log, tag, f, log_level))
        pool.close()
        pool.join()

        return result_dict


    def graph_parameter_sweep(
        n: int,
        q: int,
        e: float,
        s: float,
        m=None,
        Xe=ND.DiscreteGaussian,
        e_log=True,
        Xs=ND.DiscreteGaussian,
        s_log=True,
        tag=None,
        f=LWE.estimate,
        num_proc=8,
        log_level=0,
        make_pickle=False,
        load_pickle=False,
        directory=None,
        security_cutoff=None,
    ):
        if not directory:
            directory = os.path.dirname(os.path.realpath(__file__))
        file_name = time.strftime('%d-%m-%Y_%H-%M-%S')
        file_name = os.path.join(directory, file_name)

        assert num_proc >= 1, 'need at least one process to execute'

        if not load_pickle:
            result_dict = ParameterSweep.parameter_sweep(
                n, q, e, s, m,
                Xe, e_log, Xs, s_log,
                tag, f, num_proc, log_level
            )
            if make_pickle:
                # Pickle the intermediate computation results
                pickle.dump(result_dict, open(file_name + '.pickle', 'wb'))
                Logging.log('sweep', log_level,
                            'Pickled the intermediate computations to: %s',
                            file_name + '.pickle')
        else:
            result_dict = pickle.load(open(load_pickle + '.pickle', 'rb'))

        Xe_string = 'log_2(Xe)' if e_log else 'Xe'
        Xs_string = 'log_2(Xs)' if s_log else 'Xs'
        params = {
            'n': (n, 0),
            'q': (q, 1),
            Xe_string: (e, 2, 2),
            Xs_string: (s, 2, 3),
            'm': (m, 4),
        }

        ParameterSweep.graph_results(
            result_dict,
            params,
            file_name,
            security_cutoff,
            log_level
        )


    def security_level(
        input_params: tuple[int, float],
        result_dict: dict,
        Xe=ND.DiscreteGaussian,
        e_log=True,
        Xs=ND.DiscreteGaussian,
        s_log=True,
        tag=None,
        f=LWE.estimate,
        log_level=0,
    ) -> None:
        """
        This function calls the lattice-estimator for a given set of input
        parameters, and appends the output to the `result_dict` dict. 
        """
        n_ = int(input_params[0])
        q_ = int(input_params[1])
        e_ = 2**float(input_params[2]) if e_log else float(input_params[2])
        s_ = 2**float(input_params[3]) if s_log else float(input_params[3])
        m_ = int(input_params[4])

        lwe_params = LWE.Parameters(
            n=n_,
            q=q_,
            Xe=Xe(e_),
            Xs=Xs(s_),
            m=m_,
            tag=tag,
        )
        estimator_result = f(lwe_params)
        security = min([math.log(res_.rop, 2) for res_ in estimator_result.values()])
        result_dict[(n_, q_, float(input_params[2]), float(input_params[3]), m_, tag)] = security
        Logging.log('sweep', log_level, f'Parameters = {lwe_params}; security = {security}')


    def graph_results(result_dict: dict,
                    params: dict[str, (list[Union[int, float]], int)],
                    file_name: str,
                    security_cutoff: int,
                    log_level: int):
        """
        Creates a two graphs: one with a heatmap plot, and one with a security cutoff.

        Args:
            result_dict: a dictionary of results, like {(1, 45, 7, 4, 1): 5.7}. The
                order of the params is: n, q, Xe, Xs, m.
            params: a dictionary of the parameters that were passed in,
                and associated titles, to assist with printing axes and titles for the graph.

        Returns:
            None
        """
        
        # Convert parameter iterators to lists.
        # Also keep track of the axis variables and fixed variables, for labeling.
        axis_vars = {}
        fixed_vars = {}
        for p in params:
            try:
                # The variable is an axis variable
                params[p] = (list(params[p][0]), params[p][1])
                axis_vars[p] = len(params[p][0])
            except TypeError:
                # The variable is a fixed variable
                fixed_vars[p] = params[p][0]
                params[p] = ([params[p][0]], params[p][1])

        match len(axis_vars):
            case 0:
                raise ValueError('Cannot plot when there are no variables. '
                                 'Call the lattice estimator directly for security.')
            case 1:
                variable_param = list(axis_vars.keys())[0]  # a string like 'Xe'
                x = sorted(params[variable_param][0])
                y = sorted(result_dict.items(), key=lambda x: x[params[variable_param][1]])
                y = [_ for _ in zip(*y)][1]
                fig, ax = plt.subplots(figsize=(20, 20), dpi=80)

                ax.plot(x, y)
                ax.set_xlabel(f'Parameter: {variable_param}')
                ax.set_ylabel('Security')
                
                plt.title(f'Security with parameters {fixed_vars}')

                Logging.log('sweep', log_level, 'Saved the line plot graph to: %s',
                            file_name + '_plot.png')
                fig.savefig(file_name + '_plot.png')
            case 2:
                axis_vars = sorted(axis_vars.items(), key=lambda x: params[x[0]][1])
                x_param = axis_vars[1][0]
                x = sorted(params[x_param][0])
                y_param = axis_vars[0][0]
                y = sorted(params[y_param][0])

                values = [i[1] for i in sorted(result_dict.items(), key=lambda x: x[0])]
                mat = np.flip(np.array(values).reshape(
                    (len(set(y)), len(set(x)))),
                              axis=0)

                fig, ax = plt.subplots(figsize=(20, 20), dpi=80)

                ax.imshow(mat)
                ax.set_xticks(np.arange(0, len(set(x)), 1))
                ax.set_yticks(np.arange(0, len(set(y)), 1))
                ax.set_xticklabels(sorted(list(set(x))))
                ax.set_yticklabels(sorted(list(set(y)), reverse=True))

                plt.title(f'Security with parameters {fixed_vars}')

                for (j, i), label in np.ndenumerate(mat):
                    ax.text(i, j, round(label, 1), ha='center', va='center', color='white')

                ax.set_xlabel(f'Parameter: {x_param}')
                ax.set_ylabel(f'Parameter: {y_param}')

                fig.savefig(file_name + '_gradient.png')
                Logging.log('sweep', log_level, 'Saved the gradient graph to: %s',
                            file_name + '_gradient.png')

                if security_cutoff:
                    fig2, ax2 = plt.subplots(figsize=(20, 20), dpi=80)
                    binary_mat = mat >= security_cutoff
                    ax2.imshow(binary_mat)
                    ax2.set_xticks(np.arange(0, len(set(x)), 1))
                    ax2.set_yticks(np.arange(0, len(set(y)), 1))
                    ax2.set_xticklabels(sorted(list(set(x))))
                    ax2.set_yticklabels(sorted(list(set(y)), reverse=True))

                    for (j, i), label in np.ndenumerate(mat):
                        ax2.text(i,
                                 j,
                                 round(label, 1),
                                 ha='center',
                                 va='center',
                                 color='white' if label < security_cutoff else 'black')

                    for (j, i), label in np.ndenumerate(mat):
                        ax.text(i, j, round(label, 1), ha='center', va='center', color='white')
                    plt.title(f'Security with parameters {fixed_vars}'
                              f' and cutoff at {security_cutoff}')

                    ax2.set_xlabel(f'Parameter: {x_param}')
                    ax2.set_ylabel(f'Parameter: {y_param}')

                    fig2.savefig(file_name + '_cutoff.png')
                    Logging.log('sweep', log_level, 'Saved the cutoff graph to: %s',
                                file_name + '_cutoff.png')
            case _:
                raise ValueError('Cannot plot more than two variables. '
                                 'Try freezing one or more of them.')
