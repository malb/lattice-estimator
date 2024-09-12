# -*- coding: utf-8 -*-
"""
Cost estimates for lattice redution.
"""

from sage.all import ZZ, RR, pi, e, find_root, ceil, floor, log, oo, round, sqrt
from scipy.optimize import newton

from .cost import Cost


class ReductionCost:
    @staticmethod
    def _delta(beta):
        """
        Compute δ from block size β without enforcing β ∈ ZZ.

        δ for β ≤ 40 were computed as follows:

        ```
        # -*- coding: utf-8 -*-
        from fpylll import BKZ, IntegerMatrix

        from multiprocessing import Pool
        from sage.all import mean, sqrt, exp, log, cputime

        d, trials = 320, 32

        def f((A, beta)):

            par = BKZ.Param(block_size=beta, strategies=BKZ.DEFAULT_STRATEGY, flags=BKZ.AUTO_ABORT)
            q = A[-1, -1]
            d = A.nrows
            t = cputime()
            A = BKZ.reduction(A, par, float_type="dd")
            t = cputime(t)
            return t, exp(log(A[0].norm()/sqrt(q).n())/d)

        if __name__ == '__main__':
            for beta in (5, 10, 15, 20, 25, 28, 30, 35, 40):
                delta = []
                t = []
                i = 0
                  while i < trials:
                    threads = int(open("delta.nthreads").read()) # make sure this file exists
                    pool = Pool(threads)
                    A = [(IntegerMatrix.random(d, "qary", beta=d//2, bits=50), beta) for j in range(threads)]
                    for (t_, delta_) in pool.imap_unordered(f, A):
                        t.append(t_)
                        delta.append(delta_)
                    i += threads
                    print u"β: %2d, δ_0: %.5f, time: %5.1fs, (%2d,%2d)"%(beta, mean(delta), mean(t), i, threads)
                print
        ```

        """
        small = (
            (2, 1.02190),
            (5, 1.01862),
            (10, 1.01616),
            (15, 1.01485),
            (20, 1.01420),
            (25, 1.01342),
            (28, 1.01331),
            (40, 1.01295),
        )

        if beta <= 2:
            return RR(1.0219)
        elif beta < 40:
            for i in range(1, len(small)):
                if small[i][0] > beta:
                    return RR(small[i - 1][1])
        elif beta == 40:
            return RR(small[-1][1])
        else:
            return RR(beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))

    @staticmethod
    def delta(beta):
        """
        Compute root-Hermite factor δ from block size β.

        :param beta: Block size.
        """
        beta = ZZ(round(beta))
        return ReductionCost._delta(beta)

    @staticmethod
    def _beta_secant(delta):
        """
        Estimate required block size β for a given root-Hermite factor δ based on [PhD:Chen13]_.

        :param delta: Root-Hermite factor.

        EXAMPLE::

            >>> from estimator.reduction import ReductionCost
            >>> ReductionCost._beta_secant(1.0121)
            50
            >>> ReductionCost._beta_secant(1.0093)
            100
            >>> ReductionCost._beta_secant(1.0024) # Chen reports 800
            808

        """
        # newton() will produce a "warning", if two subsequent function values are
        # indistinguishable (i.e. equal in terms of machine precision). In this case
        # newton() will return the value beta in the middle between the two values
        # k1,k2 for which the function values were indistinguishable.
        # Since f approaches zero for beta->+Infinity, this may be the case for very
        # large inputs, like beta=1e16.
        # For now, these warnings just get printed and the value beta is used anyways.
        # This seems reasonable, since for such large inputs the exact value of beta
        # doesn't make such a big difference.
        try:
            beta = newton(
                lambda beta: RR(ReductionCost._delta(beta) - delta),
                100,
                fprime=None,
                args=(),
                tol=1.48e-08,
                maxiter=500,
            )
            beta = ceil(beta)
            if beta < 40:
                # newton may output beta < 40. The old beta method wouldn't do this. For
                # consistency, call the old beta method, i.e. consider this try as "failed".
                raise RuntimeError("β < 40")
            return beta
        except (RuntimeError, TypeError):
            # if something fails, use old beta method
            beta = ReductionCost._beta_simple(delta)
            return beta

    @staticmethod
    def _beta_find_root(delta):
        """
        Estimate required block size β for a given root-Hermite factor δ based on [PhD:Chen13]_.

        :param delta: Root-Hermite factor.

        TESTS::

            >>> from estimator.reduction import ReductionCost, RC
            >>> ReductionCost._beta_find_root(RC.delta(500))
            500

        """
        # handle beta < 40 separately
        beta = ZZ(40)
        if ReductionCost._delta(beta) < delta:
            return beta

        try:
            beta = find_root(
                lambda beta: RR(ReductionCost._delta(beta) - delta), 40, 2**16, maxiter=500
            )
            beta = ceil(beta - 1e-8)
        except RuntimeError:
            # finding root failed; reasons:
            # 1. maxiter not sufficient
            # 2. no root in given interval
            beta = ReductionCost._beta_simple(delta)
        return beta

    @staticmethod
    def _beta_simple(delta):
        """
        Estimate required block size β for a given root-Hermite factor δ based on [PhD:Chen13]_.

        :param delta: Root-Hermite factor.

        TESTS::

            >>> from estimator.reduction import ReductionCost, RC
            >>> ReductionCost._beta_simple(RC.delta(500))
            501

        """
        beta = ZZ(40)

        while ReductionCost._delta(2 * beta) > delta:
            beta *= 2
        while ReductionCost._delta(beta + 10) > delta:
            beta += 10
        while ReductionCost._delta(beta) >= delta:
            beta += 1

        return beta

    def beta(delta):
        """
        Estimate required block size β for a given root-hermite factor δ based on [PhD:Chen13]_.

        :param delta: Root-hermite factor.

        EXAMPLE::

            >>> from estimator.reduction import RC
            >>> 50 == RC.beta(1.0121)
            True
            >>> 100 == RC.beta(1.0093)
            True
            >>> RC.beta(1.0024) # Chen reports 800
            808

        """
        # TODO: decide for one strategy (secant, find_root, old) and its error handling
        beta = ReductionCost._beta_find_root(delta)
        return beta

    @classmethod
    def svp_repeat(cls, beta, d):
        """
        Return number of SVP calls in BKZ-β.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.

        .. note :: Loosely based on experiments in [PhD:Chen13].

        .. note :: When d ≤ β we return 1.

        """
        if beta < d:
            return 8 * d
        else:
            return 1

    @classmethod
    def LLL(cls, d, B=None):
        """
        Runtime estimation for LLL algorithm based on [AC:CheNgu11]_.

        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        """
        if B is None:
            return d**3  # ignoring B for backward compatibility
        else:
            return d**3 * B**2

    def short_vectors(self, beta, d, N=None, B=None, preprocess=True):
        """
        Cost of outputting many somewhat short vectors.

        The output of this function is a tuple of four values:

        - `ρ` is a scaling factor. The output vectors are expected to be longer than the shortest
          vector expected from an SVP oracle by this factor.
        - `c` is the cost of outputting `N` vectors
        - `N` the number of vectors output, which may be larger than the value put in for `N`.
        - `β'` the cost parameter associated with sampling, here: 2

        This baseline implementation uses rerandomize+LLL as in [EC:Albrecht17]_.

        :param beta: Cost parameter (≈ SVP dimension).
        :param d: Lattice dimension.
        :param N: Number of vectors requested.
        :param B: Bit-size of entries.
        :param preprocess: Include the cost of preprocessing the basis with BKZ-β.
               If ``False`` we assume the basis is already BKZ-β reduced.
        :return: ``(ρ, c, N, β')``

        EXAMPLES::

            >>> from estimator.reduction import RC
            >>> RC.CheNgu12.short_vectors(100, 500, N=1)
            (1.0, 1.67646...e17, 1, 2)
            >>> RC.CheNgu12.short_vectors(100, 500, N=1, preprocess=False)
            (1.0, 1, 1, 2)
            >>> RC.CheNgu12.short_vectors(100, 500)
            (2.0, 1.67646...e17, 1000, 2)
            >>> RC.CheNgu12.short_vectors(100, 500, preprocess=False)
            (2.0, 125000000000, 1000, 2)
            >>> RC.CheNgu12.short_vectors(100, 500, N=1000)
            (2.0, 1.67646...e17, 1000, 2)
            >>> RC.CheNgu12.short_vectors(100, 500, N=1000, preprocess=False)
            (2.0, 125000000000, 1000, 2)

        """

        if preprocess:
            cost = self(beta, d, B=B)
        else:
            cost = 0

        if N == 1:  # just call SVP
            return 1.0, cost + 1, 1, 2
        elif N is None:
            N = 1000  # pick something

        return 2.0, cost + N * RC.LLL(d), N, 2

    def short_vectors_simple(self, beta, d, N=None, B=None, preprocess=True):
        """
        Cost of outputting many somewhat short vectors.

        The output of this function is a tuple of four values:

        - `ρ` is a scaling factor. The output vectors are expected to be longer than the shortest
          vector expected from an SVP oracle by this factor.
        - `c` is the cost of outputting `N` vectors
        - `N` the number of vectors output, which may be larger than the value put in for `N`.
        - `β'` the cost parameter associated with sampling, here: `β`

        This naive baseline implementation uses rerandomize+BKZ.

        :param beta: Cost parameter (≈ SVP dimension).
        :param d: Lattice dimension.
        :param N: Number of vectors requested.
        :param B: Bit-size of entries.
        :param preprocess: This option is ignore.
        :return: ``(ρ, c, N, β')``

        EXAMPLES::

            >>> from estimator.reduction import RC
            >>> RC.CheNgu12.short_vectors_simple(100, 500, 1)
            (1.0, 1.67646160799173e17, 1, 100)
            >>> RC.CheNgu12.short_vectors_simple(100, 500)
            (1.0, 1.67646160799173e20, 1000, 100)
            >>> RC.CheNgu12.short_vectors_simple(100, 500, 1000)
            (1.0, 1.67646160799173e20, 1000, 100)

        """
        if N == 1:
            if preprocess:
                return 1.0, self(beta, d, B=B), 1, beta
            else:
                return 1.0, 1, 1, beta
        elif N is None:
            N = 1000  # pick something
        return 1.0, N * self(beta, d, B=B), N, beta

    def _short_vectors_sieve(self, beta, d, N=None, B=None, preprocess=True, sieve_dim=None):
        """
        Cost of outputting many somewhat short vectors.

        The output of this function is a tuple of four values:

        - `ρ` is a scaling factor. The output vectors are expected to be longer than the shortest
          vector expected from an SVP oracle by this factor.
        - `c` is the cost of outputting `N` vectors
        - `N` the number of vectors output, which may be larger than the value put in for `N`.
        - `β'` the cost parameter associated with sampling, here: `β` or ``sieve_dim``

        This implementation uses that a sieve outputs many somehwat short vectors [Kyber17]_.

        :param beta: Cost parameter (≈ SVP dimension).
        :param d: Lattice dimension.
        :param N: Number of vectors requested.
        :param B: Bit-size of entries.
        :param preprocess: Include the cost of preprocessing the basis with BKZ-β.
               If ``False`` we assume the basis is already BKZ-β reduced.
        :param sieve_dim: Explicit sieving dimension.
        :return: ``(ρ, c, N, β')``

        EXAMPLES::

            >>> from estimator.reduction import RC
            >>> RC.ADPS16.short_vectors(100, 500, 1)
            (1.0, 6.16702733460158e8, 1, 100)
            >>> RC.ADPS16.short_vectors(100, 500)
            (1.1547..., 6.16702733460158e8, 1763487, 100)
            >>> RC.ADPS16.short_vectors(100, 500, 1000)
            (1.1547..., 6.16702733460158e8, 1763487, 100)


        """

        if sieve_dim is None:
            sieve_dim = beta

        if N == 1:
            if preprocess:
                return 1.0, self(beta, d, B=B), 1, sieve_dim
            else:
                return 1.0, 1, 1, sieve_dim
        elif N is None:
            N = floor(2 ** (0.2075 * beta))  # pick something

        c0 = RR(N)
        c1 = RR(2 ** (RR(0.2075 * beta)))
        c = c0 / c1

        rho = sqrt(4 / 3.0) * RR(
            self.delta(sieve_dim) ** (sieve_dim - 1) * self.delta(beta) ** (1 - sieve_dim)
        )

        # arbitrary choice
        if c > 2**1000:
            # set c = oo
            return (
                rho,
                oo,
                oo,
                sieve_dim,
            )

        return (
            rho,
            ceil(c) * self(beta, d),
            ceil(c) * floor(c1),
            sieve_dim,
        )


class BDGL16(ReductionCost):
    __name__ = "BDGL16"
    short_vectors = ReductionCost._short_vectors_sieve

    @classmethod
    def _small(cls, beta, d, B=None):
        """
        Runtime estimation given β and assuming sieving is used to realise the SVP oracle for small
        dimensions following [SODA:BDGL16]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        TESTS::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.BDGL16._small(500, 1024), 2.0)
            222.9

        """
        return cls.LLL(d, B) + ZZ(2) ** RR(0.387 * beta + 16.4 + log(cls.svp_repeat(beta, d), 2))

    @classmethod
    def _asymptotic(cls, beta, d, B=None):
        """
        Runtime estimation given `β` and assuming sieving is used to realise the SVP oracle following [SODA:BDGL16]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        TESTS::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.BDGL16._asymptotic(500, 1024), 2.0)
            175.4
        """
        # TODO we simply pick the same additive constant 16.4 as for the experimental result in [SODA:BDGL16]_
        return cls.LLL(d, B) + ZZ(2) ** RR(0.292 * beta + 16.4 + log(cls.svp_repeat(beta, d), 2))

    def __call__(self, beta, d, B=None):
        """
        Runtime estimation given `β` and assuming sieving is used to realise the SVP oracle following [SODA:BDGL16]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.BDGL16(500, 1024), 2.0)
            175.4

        """
        # TODO this is somewhat arbitrary
        if beta <= 90:
            return self._small(beta, d, B)
        else:
            return self._asymptotic(beta, d, B)


class LaaMosPol14(ReductionCost):
    __name__ = "LaaMosPol14"
    short_vectors = ReductionCost._short_vectors_sieve

    def __call__(self, beta, d, B=None):
        """
        Runtime estimation for quantum sieving following [EPRINT:LaaMosPol14]_ and [PhD:Laarhoven15]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.LaaMosPol14(500, 1024), 2.0)
            161.9

        """
        return self.LLL(d, B) + ZZ(2) ** RR(
            (0.265 * beta + 16.4 + log(self.svp_repeat(beta, d), 2))
        )


class CheNgu12(ReductionCost):
    __name__ = "CheNgu12"

    def __call__(self, beta, d, B=None):
        """
        Runtime estimation given β and assuming [CheNgu12]_ estimates are correct.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        The constants in this function were derived as follows based on Table 4 in
        [CheNgu12]_::

            >>> from sage.all import var, find_fit
            >>> dim = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250]
            >>> nodes = [39.0, 44.0, 49.0, 54.0, 60.0, 66.0, 72.0, 78.0, 84.0, 96.0, \
                         99.0, 105.0, 111.0, 120.0, 127.0, 134.0]
            >>> times = [c + log(200,2).n() for c in nodes]
            >>> T = list(zip(dim, nodes))
            >>> var("a,b,c,beta")
            (a, b, c, beta)
            >>> f = a*beta*log(beta, 2.0) + b*beta + c
            >>> f = f.function(beta)
            >>> f.subs(find_fit(T, f, solution_dict=True))
            beta |--> 0.2701...*beta*log(beta) - 1.0192...*beta + 16.10...

        The estimation 2^(0.18728 β⋅log_2(β) - 1.019⋅β + 16.10) is of the number of enumeration
        nodes, hence we need to multiply by the number of cycles to process one node. This cost per
        node is typically estimated as 64.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.CheNgu12(500, 1024), 2.0)
            365.70...

        """
        repeat = self.svp_repeat(beta, d)
        cost = RR(
            0.270188776350190 * beta * log(beta)
            - 1.0192050451318417 * beta
            + 16.10253135200765
            + log(100, 2)
        )
        return self.LLL(d, B) + repeat * ZZ(2) ** cost


class ABFKSW20(ReductionCost):
    __name__ = "ABFKSW20"

    def __call__(self, beta, d, B=None):
        """
        Enumeration cost according to [C:ABFKSW20]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.ABFKSW20(500, 1024), 2.0)
            316.26...

        """
        if 1.5 * beta >= d or beta <= 92:  # 1.5β is a bit arbitrary, β≤92 is the crossover point
            cost = RR(0.1839 * beta * log(beta, 2) - 0.995 * beta + 16.25 + log(64, 2))
        else:
            cost = RR(0.125 * beta * log(beta, 2) - 0.547 * beta + 10.4 + log(64, 2))

        repeat = self.svp_repeat(beta, d)

        return self.LLL(d, B) + repeat * ZZ(2) ** cost


class ABLR21(ReductionCost):
    __name__ = "ABLR21"

    def __call__(self, beta, d, B=None):
        """
        Enumeration cost according to [C:ABLR21]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC
            >>> log(RC.ABLR21(500, 1024), 2.0)
            278.20...

        """
        if 1.5 * beta >= d or beta <= 97:  # 1.5β is a bit arbitrary, 97 is the crossover
            cost = RR(0.1839 * beta * log(beta, 2) - 1.077 * beta + 29.12 + log(64, 2))
        else:
            cost = RR(0.1250 * beta * log(beta, 2) - 0.654 * beta + 25.84 + log(64, 2))

        repeat = self.svp_repeat(beta, d)

        return self.LLL(d, B) + repeat * ZZ(2) ** cost


class ADPS16(ReductionCost):
    __name__ = "ADPS16"
    short_vectors = ReductionCost._short_vectors_sieve

    def __init__(self, mode="classical"):
        if mode not in ("classical", "quantum", "paranoid"):
            raise ValueError(f"Mode {mode} not understood.")

        self.mode = mode

    def __call__(self, beta, d, B=None):
        """
        Runtime estimation from [USENIX:ADPS16]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC, ADPS16
            >>> log(RC.ADPS16(500, 1024), 2.0)
            146.0
            >>> log(ADPS16(mode="quantum")(500, 1024), 2.0)
            132.5
            >>> log(ADPS16(mode="paranoid")(500, 1024), 2.0)
            103.75

        """

        c = {
            "classical": 0.2920,
            "quantum": 0.2650,  # paper writes 0.262 but this isn't right, see above
            "paranoid": 0.2075,
        }

        c = c[self.mode]

        return ZZ(2) ** RR(c * beta)


class ChaLoy21(ReductionCost):

    __name__ = "ChaLoy21_CoreSVP"
    short_vectors = ReductionCost._short_vectors_sieve

    def __call__(self, beta, d, B=None):
        """
        Runtime estimation for solving the Core-SVP using quantum sieving, based on the work in [AC:ChaLoy21]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.
        """

        return ZZ(2) ** RR(0.2570 * beta)

class ChaLoy21_BKZ(ReductionCost):
    __name__ = "ChaLoy21_BKZ"
    short_vectors = ReductionCost._short_vectors_sieve

    def __call__(self, beta, d, B=None):
        """
        Total runtime estimation for the BKZ (Block Korkine-Zolotarev) algorithm, 
        including multiple calls for solving Core-SVP using quantum sieving, based on the work in [AC:ChaLoy21]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.
        """

        return self.LLL(d, B) + ZZ(2) ** RR(
            (0.2570 * beta + 16.4 + log(self.svp_repeat(beta, d), 2))
        )

class Kyber(ReductionCost):
    __name__ = "Kyber"

    # These are not asymptotic expressions but compress the data in [AC:AGPS20]_ which covers up to
    # β = 1024

    # import glob
    # import pandas
    # var("a,b,d")
    # f = (a*d + b).function(d)

    # for filename in sorted(glob.glob("data/cost-estimate-*.csv")):
    #     if "sieve_size" in filename:
    #         continue
    #     costs = list(pandas.read_csv(filename)[["d", "log_cost"]].itertuples(index=False, name=None))
    #     key = filename.replace("data/cost-estimate-","").replace(".csv", "")
    #     values = find_fit(costs, f, solution_dict=True)
    #     print(f"\"{key}\": {{\"a\": {values[a]}, \"b\": {values[b]}}},")

    NN_AGPS = {
        "all_pairs-classical": {"a": 0.4215069316613415, "b": 20.1669683097337},
        "all_pairs-dw": {"a": 0.3171724396445732, "b": 25.29828951733785},
        "all_pairs-g": {"a": 0.3155285835002801, "b": 22.478746811528048},
        "all_pairs-ge19": {"a": 0.3222895263943544, "b": 36.11746438609666},
        "all_pairs-naive_classical": {"a": 0.4186251294633655, "b": 9.899382654377058},
        "all_pairs-naive_quantum": {"a": 0.31401512556555794, "b": 7.694659515948326},
        "all_pairs-t_count": {"a": 0.31553282515234704, "b": 20.878594142502994},
        "list_decoding-classical": {"a": 0.2988026130564745, "b": 26.011121212891872},
        "list_decoding-dw": {"a": 0.26944796385592995, "b": 28.97237346443934},
        "list_decoding-g": {"a": 0.26937450988892553, "b": 26.925140365395972},
        "list_decoding-ge19": {"a": 0.2695210400018704, "b": 35.47132142280775},
        "list_decoding-naive_classical": {"a": 0.2973130399197453, "b": 21.142124058689426},
        "list_decoding-naive_quantum": {"a": 0.2674316807758961, "b": 18.720680589028465},
        "list_decoding-t_count": {"a": 0.26945736714156543, "b": 25.913746774011887},
        "random_buckets-classical": {"a": 0.35586144233444716, "b": 23.082527816636638},
        "random_buckets-dw": {"a": 0.30704199612690264, "b": 25.581968903639485},
        "random_buckets-g": {"a": 0.30610964725102385, "b": 22.928235564044563},
        "random_buckets-ge19": {"a": 0.31089687599538407, "b": 36.02129978813208},
        "random_buckets-naive_classical": {"a": 0.35448283789554513, "b": 15.28878540793908},
        "random_buckets-naive_quantum": {"a": 0.30211421791887644, "b": 11.151745013027089},
        "random_buckets-t_count": {"a": 0.30614770082829745, "b": 21.41830142853265},
    }

    def __init__(self, nn="classical"):
        """
        :param nn: Nearest neighbor cost model. We default to "ListDecoding" (i.e. BDGL16) and to
                   the "depth × width" metric. Kyber uses "AllPairs".

        """
        if nn == "classical":
            nn = "list_decoding-classical"
        elif nn == "quantum":
            nn = "list_decoding-dw"
        self.nn = nn

    @staticmethod
    def d4f(beta):
        """
        Dimensions "for free" following [EC:Ducas18]_.

        :param beta: Block size ≥ 2.

        If β' is output by this function then sieving is expected to be required up to dimension β-β'.

        EXAMPLE::

            >>> from estimator.reduction import RC
            >>> RC.Kyber.d4f(500)
            42.597...

        """
        return max(float(beta * log(4 / 3.0) / log(beta / (2 * pi * e))), 0.0)

    def __call__(self, beta, d, B=None, C=5.46):
        """
        Runtime estimation from [Kyber20]_ and [AC:AGPS20]_.

        :param beta: Block size ≥ 2.
        :param d: Lattice dimension.
        :param B: Bit-size of entries.
        :param C: Progressive overhead lim_{β → ∞} ∑_{i ≤ β} 2^{0.292 i + o(i)}/2^{0.292 β + o(β)}.

        EXAMPLE::

            >>> from math import log
            >>> from estimator.reduction import RC, Kyber
            >>> log(RC.Kyber(500, 1024), 2.0)
            176.61534319964488
            >>> log(Kyber(nn="list_decoding-ge19")(500, 1024), 2.0)
            172.68208507350872

        """

        if beta < 20:  # goes haywire
            return CheNgu12()(beta, d, B)

        # "The cost of progressive BKZ with sieving up to blocksize b is essentially C · (n − b) ≈
        # 3340 times the cost of sieving for SVP in dimension b." [Kyber20]_
        svp_calls = C * max(d - beta, 1)
        # we do not round to the nearest integer to ensure cost is continuously increasing with β which
        # rounding can violate.
        beta_ = beta - self.d4f(beta)
        # "The work in [5] is motivated by the quantum/classical speed-up, therefore it does not
        # consider the required number of calls to AllPairSearch. Naive sieving requires a
        # polynomial number of calls to this routine, however this number of calls appears rather
        # small in practice using progressive sieving [40, 64], and we will assume that it needs to
        # be called only once per dimension during progressive sieving, for a cost of C · 2^137.4
        # gates^8." [Kyber20]_

        gate_count = C * 2 ** (
            RR(self.NN_AGPS[self.nn]["a"]) * beta_ + RR(self.NN_AGPS[self.nn]["b"])
        )
        return self.LLL(d, B=B) + svp_calls * gate_count

    def short_vectors(self, beta, d, N=None, B=None, preprocess=True):
        """
        Cost of outputting many somewhat short vectors using BKZ-β.

        The output of this function is a tuple of four values:

        - `ρ` is a scaling factor. The output vectors are expected to be longer than the shortest
          vector expected from an SVP oracle by this factor.
        - `c` is the cost of outputting `N` vectors
        - `N` the number of vectors output, which may be larger than the value put in for `N`.
        - `β'` the cost parameter associated with sampling

        This is using an observation insprired by [AC:GuoJoh21]_ that we can run a sieve on the
        first block of the basis with negligible overhead.

        :param beta: Cost parameter (≈ SVP dimension).
        :param d: Lattice dimension.
        :param N: Number of vectors requested.
        :param preprocess: Include the cost of preprocessing the basis with BKZ-β.
               If ``False`` we assume the basis is already BKZ-β reduced.
        :return: ``(ρ, c, N, β')``

        EXAMPLES::

            >>> from estimator.reduction import RC
            >>> RC.Kyber.short_vectors(100, 500, 1)
            (1.0, 2.7367476128136...19, 100, 1)
            >>> RC.Kyber.short_vectors(100, 500)
            (1.1547, 2.7367476128136...19, 176584, 84)
            >>> RC.Kyber.short_vectors(100, 500, 1000)
            (1.1547, 2.7367476128136...19, 176584, 84)

        """
        beta_ = beta - floor(self.d4f(beta))

        if N == 1:
            if preprocess:
                return 1.0, self(beta, d, B=B), beta, 1
            else:
                return 1.0, 1, beta, 1
        elif N is None:
            N = floor(2 ** (0.2075 * beta_))  # pick something

        c = N / floor(2 ** (0.2075 * beta_))
        return 1.1547, ceil(c) * self(beta, d), ceil(c) * floor(2 ** (0.2075 * beta_)), beta_


class GJ21(Kyber):
    __name__ = "GJ21"

    def short_vectors(self, beta, d, N=None, preprocess=True, B=None, C=5.46, sieve_dim=None):
        """
        Cost of outputting many somewhat short vectors according to [AC:GuoJoh21]_.

        The output of this function is a tuple of four values:

        - `ρ` is a scaling factor. The output vectors are expected to be longer than the shortest
          vector expected from an SVP oracle by this factor.
        - `c` is the cost of outputting `N` vectors
        - `N` the number of vectors output, which may be larger than the value put in for `N`.
        - `β'` the cost parameter associated with sampling

        This runs a sieve on the first β_0 vectors of the basis after BKZ-β reduction
        to produce many short vectors, where β_0 is chosen such that BKZ-β reduction and the sieve
        run in approximately the same time. [AC:GuoJoh21]_

        :param beta: Cost parameter (≈ SVP dimension).
        :param d: Lattice dimension.
        :param N: Number of vectors requested.
        :param preprocess: Include the cost of preprocessing the basis with BKZ-β.
               If ``False`` we assume the basis is already BKZ-β reduced.
        :param B: Bit-size of entries.
        :param C: Progressive overhead lim_{β → ∞} ∑_{i ≤ β} 2^{0.292 i + o(i)}/2^{0.292 β + o(β)}.
        :param sieve_dim: Explicit sieving dimension.
        :return: ``(ρ, c, N, β')``

        EXAMPLES::

            >>> from estimator.reduction import RC
            >>> RC.GJ21.short_vectors(100, 500, 1)
            (1.0, 2.7367476128136...19, 1, 100)
            >>> RC.GJ21.short_vectors(100, 500)
            (1.04228014727497, 5.56224438...19, 36150192, 121)
            >>> RC.GJ21.short_vectors(100, 500, 1000)
            (1.04228014727497, 5.56224438...19, 36150192, 121)

        """
        beta_ = beta - floor(self.d4f(beta))
        if sieve_dim is None:
            sieve_dim = beta_
            if beta < d:
                # set beta_sieve such that complexity of 1 sieve in dim sieve_dim is approx
                # the same as the BKZ call
                sieve_dim = min(
                    d, floor(beta_ + log((d - beta) * C, 2) / self.NN_AGPS[self.nn]["a"])
                )

        # MATZOV, p.18 (they call the slope δ_β, we call it α_β)
        # gh(β) ≈ √(β/2πe)
        # α_β ≈ (β/2πe)^(1/(β-1))
        # λ_1' = gh(sieve_dim) ⋅ α_β^{(d-sieve_dim)/2} ⋅ vol(Λ)^{1/d}
        #      = α_{sieve_dim}^((sieve_dim-1)/2)  ⋅ α_β^{(d-sieve_dim)/2} ⋅ vol(Λ)^{1/d}
        # shortest vector in BKZ-β reduced basis
        # λ_1 =  α_β^{(d-1)/2} ⋅ vol(Λ)^{1/d}
        #      = α_β^((sieve_dim-1)/2)  ⋅ α_β^{(d-sieve_dim)/2} ⋅ vol(Λ)^{1/d}
        # λ_1'/λ_1  = α_{sieve_dim}^((sieve_dim-1)/2)  ⋅ α_β^{(d-sieve_dim)/2} ⋅ vol(Λ)^{1/d}
        #           / α_{β}^((sieve_dim-1)/2)          ⋅ α_β^{(d-sieve_dim)/2} ⋅ vol(Λ)^{1/d}
        #           = α_{sieve_dim}^((sieve_dim-1)/2) / α_{β}^((sieve_dim-1)/2)
        #           ≈ δ_{sieve_dim}^(sieve_dim-1) / δ_{β}^(sieve_dim-1)

        rho = sqrt(4 / 3.0) * RR(
            self.delta(sieve_dim) ** (sieve_dim - 1) * self.delta(beta) ** (1 - sieve_dim)
        )

        if N == 1:
            if preprocess:
                return 1.0, self(beta, d, B=B), 1, beta
            else:
                return 1.0, 1, 1, beta
        elif N is None:
            N = floor(2 ** (0.2075 * sieve_dim))  # pick something

        c0 = RR(N)
        c1 = RR(2 ** RR(0.2075 * sieve_dim))
        c = c0 / floor(c1)
        sieve_cost = C * 2 ** RR((self.NN_AGPS[self.nn]["a"] * sieve_dim + self.NN_AGPS[self.nn]["b"]))

        # arbitrary choice
        if c > 2**1000:
            # set c = oo
            return (
                rho,
                oo,
                oo,
                sieve_dim,
            )

        return (
            rho,
            ceil(c) * (self(beta, d) + sieve_cost),
            ceil(c) * floor(c1),
            sieve_dim,
        )


class MATZOV(GJ21):
    """
    Improved enumeration routine in list decoding from [MATZOV22]_.
    """

    __name__ = "MATZOV"

    # These are not asymptotic expressions but compress the data in [AC:AGPS20]_ with the fix and
    # improvement from [MATZOV22]_ applied which covers up to β = 1024
    NN_AGPS = {
        "all_pairs-classical": {"a": 0.4215069316732438, "b": 20.166968300536567},
        "all_pairs-dw": {"a": 0.3171724396445733, "b": 25.2982895173379},
        "all_pairs-g": {"a": 0.31552858350028, "b": 22.478746811528104},
        "all_pairs-ge19": {"a": 0.3222895263943547, "b": 36.11746438609664},
        "all_pairs-naive_classical": {"a": 0.41862512941897706, "b": 9.899382685790897},
        "all_pairs-naive_quantum": {"a": 0.31401512571180035, "b": 7.694659414353819},
        "all_pairs-t_count": {"a": 0.31553282513562797, "b": 20.87859415484879},
        "list_decoding-classical": {"a": 0.29613500308205365, "b": 20.387885985467914},
        "list_decoding-dw": {"a": 0.2663676536352464, "b": 25.299541499216627},
        "list_decoding-g": {"a": 0.26600114174341505, "b": 23.440974518186337},
        "list_decoding-ge19": {"a": 0.26799889622667994, "b": 30.839871638418543},
        "list_decoding-naive_classical": {"a": 0.29371310617068064, "b": 15.930690682515422},
        "list_decoding-naive_quantum": {"a": 0.2632557273632713, "b": 15.685687713591548},
        "list_decoding-t_count": {"a": 0.2660264010780807, "b": 22.432158856991474},
        "random_buckets-classical": {"a": 0.3558614423344473, "b": 23.08252781663665},
        "random_buckets-dw": {"a": 0.30704199602260734, "b": 25.58196897625173},
        "random_buckets-g": {"a": 0.30610964725102396, "b": 22.928235564044588},
        "random_buckets-ge19": {"a": 0.31089687605567917, "b": 36.02129974535213},
        "random_buckets-naive_classical": {"a": 0.35448283789554536, "b": 15.28878540793911},
        "random_buckets-naive_quantum": {"a": 0.3021142178390157, "b": 11.151745066682524},
        "random_buckets-t_count": {"a": 0.3061477007403873, "b": 21.418301489775203},
    }


def cost(cost_model, beta, d, B=None, predicate=True, **kwds):
    """
    Return cost dictionary for computing vector of norm` δ_0^{d-1} Vol(Λ)^{1/d}` using provided lattice
    reduction algorithm.

    :param cost_model:
    :param beta: Block size ≥ 2.
    :param d: Lattice dimension.
    :param B: Bit-size of entries.
    :param predicate: if ``False`` cost will be infinity.

    EXAMPLE::

        >>> from estimator.reduction import cost, RC
        >>> cost(RC.ABLR21, 120, 500)
        rop: ≈2^68.9, red: ≈2^68.9, δ: 1.008435, β: 120, d: 500
        >>> cost(RC.ABLR21, 120, 500, predicate=False)
        rop: ≈2^inf, red: ≈2^inf, δ: 1.008435, β: 120, d: 500

    """
    # convenience: instantiate static classes if needed
    if isinstance(cost_model, type):
        cost_model = cost_model()

    cost = cost_model(beta, d, B)
    delta_ = ReductionCost.delta(beta)
    cost = Cost(rop=cost, red=cost, delta=delta_, beta=beta, d=d, **kwds)
    cost.register_impermanent(rop=True, red=True, delta=False, beta=False, d=False)
    if predicate is False:
        cost["red"] = oo
        cost["rop"] = oo
    return cost


beta = ReductionCost.beta
delta = ReductionCost.delta


class RC:
    beta = ReductionCost.beta
    delta = ReductionCost.delta

    LLL = ReductionCost.LLL
    ABFKSW20 = ABFKSW20()
    ABLR21 = ABLR21()
    ADPS16 = ADPS16()
    BDGL16 = BDGL16()
    CheNgu12 = CheNgu12()
    Kyber = Kyber()
    MATZOV = MATZOV()
    GJ21 = GJ21()
    LaaMosPol14 = LaaMosPol14()
    ChaLoy21 = ChaLoy21()
    ChaLoy21_BKZ = ChaLoy21_BKZ()
