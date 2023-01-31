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
from estimator.nd import NoiseDistribution as ND
from matplotlib import pyplot as plt
from typing import Iterable, Union, Optional, Callable
import pickle
import time
import math
import multiprocessing
import numpy as np
import os


class ParameterSweep:
    """
    A class that provides utilities for performing and graphing the results
    from performing parameter sweeps using the lattice estimator. To import in `sage`:
    ```
    from param_sweep import ParameterSweep as PS
    ```

    Example of a parameter sweep (without graphing), with:
    - `n` from 600 to 1140 incrementing by 20,
    - `log_2(Xe)` from 6 to 20 incrementing by 1,
    - `Xs` is a uniform binary distribution,
    - `q` is 2**32,
    - `m`==`n` (the default if `m` is not specified),
    - `f` as `LWE.estimate.rough` for faster computation
    From `sage`, run:
    ```
    from estimator import LWE, nd
    results = PS.parameter_sweep(n=range(600,1140,20), q=2**32, e=(6,20), s=2, s_log=False, Xs=nd.NoiseDistribution.UniformMod, f=LWE.estimate.rough)
    ```

    Example of a parameter sweep to generate a heatmap graph, using:
    - the same LWE parameters as the previous example,
    - changing the function to LWE.estimate for a more precise estimation,
    - pickling the intermediate results
    From `sage`, run:
    ```
    from estimator import nd
    PS.graph_parameter_sweep(n=range(600,1140,20), q=2**32, e=range(6,20), s=2, s_log=False, Xs=nd.NoiseDistribution.UniformMod, make_pickle=True, file_name='demo')
    ```

    Load the pickled results from the previous run, and additionally generate
    a security cutoff graph. From `sage`, run:
    ```
    PS.graph_parameter_sweep(n=range(600,1140,20), q=2**32, e=range(6,20), s=2, load_pickle=True, file_name='demo', security_cutoff=128)
    ```
    """

    def parameter_sweep(
        n: Union[int, Iterable],
        q: Union[int, Iterable],
        e: Union[float, Iterable],
        s: Union[float, Iterable],
        m: Optional[Union[int, Iterable]] = None,
        Xe: Callable = ND.DiscreteGaussian,
        e_log: bool = True,
        Xs: Callable = ND.DiscreteGaussian,
        s_log: bool = True,
        tag: str = None,
        f: Callable = LWE.estimate,
        num_proc: int = 8,
        log_level: int = 0,
    ) -> dict:
        """
        Performs a sweep over the parameters specified.

        :param n: the dimension of the LWE sample vector (Z/qZ)^n.
        :param q: the modulus of the space Z/qZ of integers for the LWE samples.
        :param e: the input to the distribution function for the error term.
        :param s: the input to the distribution function for the secret term.
        :param m: the number of LWE samples allowed to an attacker.
        :param Xe: the distribution function for the error term.
        :param e_log: whether to plot the error on a logarithmic scale.
        :param Xs: the distribution function for the secret term.
        :param s_log: whether to plot the secret on a logarithmic scale.
        :param tag: a name for the patameter set
        :param f: the estimation function. Use `LWE.estimate.rough` for speed.
        :param num_proc: the number of parallel processes for computation.
        :param log_level: the logging level.
        :returns: a dictionary mapping from a set of parameters, to the
            estimated security level for those parameters. The ordering of
            the parameters in the dict key is: (n, q, e, s, m).

        EXAMPLE ::

            >>> from estimator import LWE, nd
            >>> from param_sweep import ParameterSweep as PS
            >>> n_list = [600, 900]
            >>> e_list = [7, 9]
            >>> results = PS.parameter_sweep(\
                    n=n_list,\
                    q=2**32,\
                    e=e_list,\
                    s=2,\
                    s_log=False,\
                    Xs=nd.NoiseDistribution.UniformMod,\
                    f=LWE.estimate.rough,\
                    tag='test',\
                    log_level=2,\
                    num_proc=1,\
                )
            usvp                 :: rop: ≈2^45.6, red: ≈2^45.6, δ: 1.007290, β: 156, d: 1120, tag: usvp
            dual_hybrid          :: rop: ≈2^46.6, mem: ≈2^42.7, m: 579, β: 159, d: 1169, ↻: 1, ζ: 10, tag: dual_hybrid
            usvp                 :: rop: ≈2^51.7, red: ≈2^51.7, δ: 1.006767, β: 177, d: 1124, tag: usvp
            dual_hybrid          :: rop: ≈2^52.7, mem: ≈2^48.3, m: 571, β: 180, d: 1160, ↻: 1, ζ: 11, tag: dual_hybrid
            usvp                 :: rop: ≈2^82.9, red: ≈2^82.9, δ: 1.005021, β: 284, d: 1661, tag: usvp
            dual_hybrid          :: rop: ≈2^83.0, mem: ≈2^77.8, m: 830, β: 284, d: 1711, ↻: 1, ζ: 19, tag: dual_hybrid
            usvp                 :: rop: ≈2^92.6, red: ≈2^92.6, δ: 1.004667, β: 317, d: 1650, tag: usvp
            dual_hybrid          :: rop: ≈2^92.4, mem: ≈2^87.4, m: 816, β: 316, d: 1694, ↻: 1, ζ: 22, tag: dual_hybrid
            >>> results[(600, 4294967296, 9.0, 2.0, 600, 'test')]
            51.684000000000005
            >>> results[(600, 4294967296, 7.0, 2.0, 600, 'test')]
            45.552
            >>> results[(900, 4294967296, 7.0, 2.0, 900, 'test')]
            82.928
            >>> results[(900, 4294967296, 9.0, 2.0, 900, 'test')]
            92.36860677483823
        """
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

        if num_proc <= 1:
            result_dict = {}
            for i in range(len(work)):
                ParameterSweep.security_level(work[i], result_dict, Xe, e_log,
                                              Xs, s_log, tag, f, log_level)
        else:
            # Parallel process the calculations
            pool = multiprocessing.Pool(processes=min(num_proc, len(work)))
            for i in range(len(work)):
                pool.apply_async(ParameterSweep.security_level,
                                 (work[i], result_dict, Xe, e_log, Xs, s_log,
                                  tag, f, log_level))
            pool.close()
            pool.join()

        return dict(result_dict)

    def security_level(
        input_params: tuple[int, float],
        result_dict: dict,
        Xe: Callable = ND.DiscreteGaussian,
        e_log: bool = True,
        Xs: Callable = ND.DiscreteGaussian,
        s_log: bool = True,
        tag: str = None,
        f: Callable = LWE.estimate,
        log_level: int = 0,
    ) -> None:
        """
        Calls the lattice-estimator for a given set of input
        parameters, and appends the output to the `result_dict` dict.

        :param input_params: a tuple of params to estimate security for.
            Parameter ordering is: (n: int, q: int, e: float, s: float, m: int).
        :param result_dict: the dictionary to append the output to. The key is
            the tuple of parameters from input_params.
        :param Xe: the distribution function for the error term.
        :param e_log: whether to plot the error on a logarithmic scale.
        :param Xs: the distribution function for the secret term.
        :param s_log: whether to plot the secret on a logarithmic scale.
        :param tag: a name for the patameter set
        :param f: the estimation function. Use `LWE.estimate.rough` for speed.
        :param log_level: the logging level.
        """
        n_ = int(input_params[0])
        q_ = int(input_params[1])
        e_ = 2**input_params[2] if e_log else input_params[2]
        s_ = 2**input_params[3] if s_log else input_params[3]
        # If m = infinity, pass infinity to the estimator (since infinity can't be cast to an int).
        m_ = float('inf') if input_params[4] == float('inf') else int(
            input_params[4])

        lwe_params = LWE.Parameters(
            n=n_,
            q=q_,
            Xe=Xe(e_),
            Xs=Xs(s_),
            m=m_,
            tag=tag,
        )
        estimator_result = f(lwe_params)
        security = min(
            [math.log(res_.get('rop', 0), 2) for res_ in estimator_result.values()])
        if not security:
            raise ValueError('ROP for a estimator result was 0, estimator failed')
        result_dict[(n_, q_, float(input_params[2]), float(input_params[3]),
                     m_, tag)] = security
        Logging.log('sweep', log_level,
                    f'Parameters = {lwe_params}; security = {security}')

    def graph_parameter_sweep(
        n: Union[int, Iterable],
        q: Union[int, Iterable],
        e: Union[float, Iterable],
        s: Union[float, Iterable],
        m: Optional[Union[int, Iterable]] = None,
        Xe: Callable = ND.DiscreteGaussian,
        e_log: bool = True,
        Xs: Callable = ND.DiscreteGaussian,
        s_log: bool = True,
        tag: str = None,
        f: Callable = LWE.estimate,
        num_proc: int = 8,
        log_level: int = 0,
        make_pickle: bool = False,
        load_pickle: bool = False,
        security_cutoff: int = None,
        directory: str = None,
        file_name: str = None,
        extension: str = '.png',
    ) -> None:
        """
        Gets the results of a parameter sweep, and creates graph visualizations
        of the data. The type of graph depends on the number of variables.

        :param n: the dimension of the LWE sample vector (Z/qZ)^n.
        :param q: the modulus of the space Z/qZ of integers for the LWE samples.
        :param e: the input to the distribution function for the error term.
        :param s: the input to the distribution function for the secret term.
        :param m: the number of LWE samples allowed to an attacker.
        :param Xe: the distribution function for the error term.
        :param e_log: whether to plot the error on a logarithmic scale.
        :param Xs: the distribution function for the secret term.
        :param s_log: whether to plot the secret on a logarithmic scale.
        :param tag: a name for the patameter set
        :param f: the estimation function. Use `LWE.estimate.rough` for speed.
        :param num_proc: the number of parallel processes for computation.
        :param log_level: the logging level.
        :param make_pickle: whether to make a pickle file of the results dict.
        :param load_pickle: whether to load a pickle file of the results dict.
        :param security_cutoff: makes a separate graph with a security cutoff.
        :param directory: the directory to load files from and/or save files to.
        :param file_name: the file name to load files from and/or save files to.
        :param extension: the extension of the graph(s). Ex: .png, .pdf, .svg.

        EXAMPLE ::

            >>> from estimator import LWE, nd
            >>> from param_sweep import ParameterSweep as PS
            >>> import uuid
            >>> from pathlib import Path
            >>> e_range = range(7, 10, 2)
            >>> s_range = range(2, 5, 2)
            >>> file_name = 'test_file_' + str(uuid.uuid4())
            >>> _ = PS.graph_parameter_sweep(\
                    n=700,\
                    q=2**32,\
                    e=e_range,\
                    s=s_range,\
                    f=LWE.estimate.rough,\
                    tag='test',\
                    directory='/tmp',\
                    make_pickle=True,\
                    security_cutoff=128,\
                    file_name=file_name,\
                    num_proc=1,\
                )
            usvp                 :: rop: ≈2^69.2, red: ≈2^69.2, δ: 1.005647, β: 237, d: 1396, tag: usvp
            dual_hybrid          :: rop: ≈2^78.5, mem: ≈2^75.8, m: 700, β: 267, d: 1392, ↻: 1, ζ: 8, tag: dual_hybrid
            usvp                 :: rop: ≈2^78.8, red: ≈2^78.8, δ: 1.005191, β: 270, d: 1396, tag: usvp
            dual_hybrid          :: rop: ≈2^92.5, mem: ≈2^90.3, m: 700, β: 314, d: 1392, ↻: 1, ζ: 8, tag: dual_hybrid
            usvp                 :: rop: ≈2^78.8, red: ≈2^78.8, δ: 1.005191, β: 270, d: 1391, tag: usvp
            dual_hybrid          :: rop: ≈2^87.2, mem: ≈2^84.2, m: 700, β: 297, d: 1392, ↻: 1, ζ: 8, tag: dual_hybrid
            usvp                 :: rop: ≈2^90.5, red: ≈2^90.5, δ: 1.004738, β: 310, d: 1394, tag: usvp
            dual_hybrid          :: rop: ≈2^102.5, mem: ≈2^100.0, m: 700, β: 349, d: 1392, ↻: 1, ζ: 8, tag: dual_hybrid
            >>> Path(f'/tmp/{file_name}.pickle').exists()
            True
            >>> Path(f'/tmp/{file_name}_gradient.png').exists()
            True
            >>> Path(f'/tmp/{file_name}_cutoff.png').exists()
            True
            >>> results = pickle.load(open(f'/tmp/{file_name}.pickle', 'rb'))
            >>> results[(700, 4294967296, 9.0, 2.0, 700, 'test')]
            78.83999999999999
            >>> results[(700, 4294967296, 9.0, 4.0, 700, 'test')]
            90.52000000000001
            >>> results[(700, 4294967296, 7.0, 4.0, 700, 'test')]
            78.83999999999999
            >>> results[(700, 4294967296, 7.0, 2.0, 700, 'test')]
            69.204
        """
        if not directory:
            directory = os.path.dirname(os.path.realpath(__file__))
        if not file_name:
            file_name = time.strftime('%d-%m-%Y_%H-%M-%S')
        file_name = os.path.join(directory, file_name)
        assert num_proc >= 1, 'need at least one process to execute'

        if not load_pickle:
            result_dict = ParameterSweep.parameter_sweep(
                n, q, e, s, m, Xe, e_log, Xs, s_log, tag, f, num_proc,
                log_level)
            if make_pickle:
                # Pickle the intermediate computation results
                pickle.dump(result_dict, open(file_name + '.pickle', 'wb'))
                Logging.log('sweep', log_level,
                            'Pickled the intermediate computations to: %s',
                            file_name + '.pickle')
        else:
            result_dict = pickle.load(open(file_name + '.pickle', 'rb'))

        Xe_string = 'log_2(Xe)' if e_log else 'Xe'
        Xs_string = 'log_2(Xs)' if s_log else 'Xs'
        params = {
            'n': (n, 0),
            'q': (q, 1),
            Xe_string: (e, 2),
            Xs_string: (s, 3),
            'm': (m, 4),
        }

        ParameterSweep.graph_results(
            result_dict,
            params,
            file_name,
            security_cutoff,
            log_level,
            extension,
        )

    def graph_results(
        result_dict: dict,
        params: dict[str, (list[Union[int, float]], int)],
        file_name: str,
        security_cutoff: int = None,
        log_level: int = 0,
        extension: str = '.png',
    ):
        """
        Graph the security estimate results in `result_dict`, using the
        parameters in `params` for labeling the plot axes and titles.
        - For 1 variable: creates a line plot.
        - For 2 variables: creates a heatmap plot, and an optional cutoff plot.

        :param result_dict: a mapping from a set of parameters, to security.
            Parameter ordering: (n: int, q: int, e: float, s: float, m: int).
        :param params: a mapping from the string representation of a parameter,
            to a tuple of the parameter and its associated order in the
            `result_dict` key. Example: {'n': (600, 0), 'q'': (4294967296, 1)}
        :param file_name: the file name to write the output graphs to.
        :param security_cutoff: makes a separate graph with a security cutoff.
        :param log_level: the logging level
        :param extension: the extension of the graph(s). Ex: .png, .pdf, .svg.
        """

        # Convert parameter iterators to lists.
        # Also keep track of the axis and fixed variables, for labeling.
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

        if len(axis_vars) == 0:
            raise ValueError(
                'Cannot plot when there are no variables. '
                'Call the lattice estimator directly for security.')
        elif len(axis_vars) == 1:
            variable_param = list(axis_vars.keys())[0]  # a string like 'Xe'
            x = sorted(params[variable_param][0])
            y = sorted(result_dict.items(),
                       key=lambda x: x[params[variable_param][1]])
            y = [_ for _ in zip(*y)][1]
            fig, ax = plt.subplots(figsize=(20, 20), dpi=80)

            ax.plot(x, y)
            ax.set_xlabel(f'Parameter: {variable_param}')
            ax.set_ylabel('Security')

            plt.title(f'Security with parameters {fixed_vars}')

            Logging.log('sweep', log_level, 'Saved the line plot graph to: %s',
                        file_name + '_plot' + extension)
            fig.savefig(file_name + '_plot' + extension)
        elif len(axis_vars) == 2:
            axis_vars = sorted(axis_vars.items(),
                               key=lambda x: params[x[0]][1])
            x_param = axis_vars[0][0]
            x = sorted(params[x_param][0])
            y_param = axis_vars[1][0]
            y = sorted(params[y_param][0])

            values = [
                i[1] for i in sorted(result_dict.items(), key=lambda x: x[0])
            ]
            mat = np.flip(np.array(values).reshape((len(set(x)), len(set(y)))),
                          axis=1).transpose()

            fig, ax = plt.subplots(figsize=(20, 20), dpi=80)

            ax.imshow(mat)
            ax.set_xticks(np.arange(0, len(set(x)), 1))
            ax.set_yticks(np.arange(0, len(set(y)), 1))
            ax.set_xticklabels(sorted(list(set(x))))
            ax.set_yticklabels(sorted(list(set(y)), reverse=True))

            plt.title(f'Security with fixed parameters {fixed_vars}')

            for (j, i), label in np.ndenumerate(mat):
                ax.text(i,
                        j,
                        round(label, 1),
                        ha='center',
                        va='center',
                        color='white')

            ax.set_xlabel(f'Parameter: {x_param}')
            ax.set_ylabel(f'Parameter: {y_param}')

            fig.savefig(file_name + '_gradient' + extension)
            Logging.log('sweep', log_level, 'Saved the gradient graph to: %s',
                        file_name + '_gradient' + extension)

            if security_cutoff:
                fig2, ax2 = plt.subplots(figsize=(20, 20), dpi=80)
                binary_mat = mat >= security_cutoff
                ax2.imshow(binary_mat)
                ax2.set_xticks(np.arange(0, len(set(x)), 1))
                ax2.set_yticks(np.arange(0, len(set(y)), 1))
                ax2.set_xticklabels(sorted(list(set(x))))
                ax2.set_yticklabels(sorted(list(set(y)), reverse=True))

                for (j, i), label in np.ndenumerate(mat):
                    ax2.text(
                        i,
                        j,
                        round(label, 1),
                        ha='center',
                        va='center',
                        color='white' if label < security_cutoff else 'black')

                for (j, i), label in np.ndenumerate(mat):
                    ax.text(i,
                            j,
                            round(label, 1),
                            ha='center',
                            va='center',
                            color='white')
                plt.title(f'Security with fixed parameters {fixed_vars}'
                          f' and cutoff at {security_cutoff}')

                ax2.set_xlabel(f'Parameter: {x_param}')
                ax2.set_ylabel(f'Parameter: {y_param}')

                fig2.savefig(file_name + '_cutoff' + extension)
                Logging.log('sweep', log_level,
                            'Saved the cutoff graph to: %s',
                            file_name + '_cutoff' + extension)
        else:
            raise ValueError('Cannot plot more than two variables. '
                             'Try freezing one or more of them.')
