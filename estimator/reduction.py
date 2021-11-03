# -*- coding: utf-8 -*-
from sage.all import ZZ, RR, pi, e, find_root, ceil, log, oo, round
from scipy.optimize import newton


class BKZ:
    """
    Cost estimates for BKZ.
    """

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
                delta_0 = []
                t = []
                i = 0
                  while i < trials:
                    threads = int(open("delta_0.nthreads").read()) # make sure this file exists
                    pool = Pool(threads)
                    A = [(IntegerMatrix.random(d, "qary", beta=d//2, bits=50), beta) for j in range(threads)]
                    for (t_, delta_0_) in pool.imap_unordered(f, A):
                        t.append(t_)
                        delta_0.append(delta_0_)
                    i += threads
                    print u"β: %2d, δ_0: %.5f, time: %5.1fs, (%2d,%2d)"%(beta, mean(delta_0), mean(t), i, threads)
                print
        ```

        """
        small = (
            (2, 1.02190),  # noqa
            (5, 1.01862),  # noqa
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

    @classmethod
    def delta(cls, beta):
        """
        Compute root-Hermite factor δ from block size β.
        """
        beta = ZZ(round(beta))
        return cls._delta(beta)

    @staticmethod
    def _beta_secant(delta):
        """
        Estimate required block size β for a given root-Hermite factor δ based on [PhD:Chen13]_.

        :param delta: root-Hermite factor

        EXAMPLE::

            >>> from estimator.reduction import BKZ
            >>> 50 == BKZ._beta_secant(1.0121)
            True
            >>> 100 == BKZ._beta_secant(1.0093)
            True
            >>> BKZ._beta_secant(1.0024) # Chen reports 800
            808

        .. [PhD:Chen13] Yuanmi Chen. Réduction de réseau et sécurité concrète du chiffrement
                        complètement homomorphe. PhD thesis, Paris 7, 2013.
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
                lambda beta: RR(BKZ._delta(beta) - delta),
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
            beta = BKZ._beta_simple(delta)
            return beta

    @classmethod
    def _beta_find_root(cls, delta):
        """
        TESTS::

            >>> from estimator.reduction import BKZ
            >>> BKZ._beta_find_root(BKZ.delta(500))
            500

        """
        # handle beta < 40 separately
        beta = ZZ(40)
        if cls.delta(beta) < delta:
            return beta

        try:
            beta = find_root(lambda beta: RR(BKZ._delta(beta) - delta), 40, 2 ** 16, maxiter=500)
            beta = ceil(beta - 1e-8)
        except RuntimeError:
            # finding root failed; reasons:
            # 1. maxiter not sufficient
            # 2. no root in given interval
            beta = BKZ._beta_simple(delta)
        return beta

    @classmethod
    def _beta_simple(cls, delta):
        beta = ZZ(40)

        while cls._deltaf(2 * beta) > delta:
            beta *= 2
        while cls._deltaf(beta + 10) > delta:
            beta += 10
        while True:
            if cls._deltaf(beta) < delta:
                break
            beta += 1

        return beta

    @classmethod
    def beta(cls, delta):
        """
        Estimate required blocksize β for a given root-hermite factor δ based on [PhD:Chen13]_.

        :param delta: root-hermite factor

        EXAMPLE::

            >>> from estimator.reduction import BKZ
            >>> 50 == BKZ.beta(1.0121)
            True
            >>> 100 == BKZ.beta(1.0093)
            True
            >>> BKZ.beta(1.0024) # Chen reports 800
            808

        .. [PhD:Chen13] Yuanmi Chen. Réduction de réseau et sécurité concrète du chiffrement
                        complètement homomorphe. PhD thesis, Paris 7, 2013.
        """
        # TODO: decide for one strategy (secant, find_root, old) and its error handling
        beta = cls._beta_find_root(delta)
        return beta

    # BKZ Estimates

    @staticmethod
    def svp_repeat(beta, d):
        """
        Return number of SVP calls in BKZ-β.

        :param beta: block size
        :param d: dimension

        .. note :: loosely based on experiments in [PhD:Chen13]

        """
        if beta < d:
            return 8 * d
        else:
            return 1

    @staticmethod
    def LLL(d, B=None):
        """
        Runtime estimation for LLL algorithm.

        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [CheNgu11] Chen, Y., & Nguyen, P.  Q.  (2011).  BKZ 2.0: better lattice security
            estimates.  In D.  H.  Lee, & X.  Wang, ASIACRYPT~2011 (pp.  1–20).  : Springer,
            Heidelberg.
        """
        if B:
            return d ** 3 * B ** 2
        else:
            return d ** 3  # ignoring B for backward compatibility

    @staticmethod
    def _BDGL16_small(beta, d, B=None):
        """
        Runtime estimation given β and assuming sieving is used to realise the SVP oracle for small dimensions.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [BDGL16] Becker, A., Ducas, L., Gama, N., & Laarhoven, T.  (2016).  New directions in
        nearest neighbor searching with applications to lattice sieving.  In SODA 2016, (pp.
        10–24).

        """
        return BKZ.LLL(d, B) + ZZ(2) ** RR(0.387 * beta + 16.4 + log(BKZ.svp_repeat(beta, d), 2))

    @staticmethod
    def _BDGL16_asymptotic(beta, d, B=None):
        """
        Runtime estimation given `β` and assuming sieving is used to realise the SVP oracle.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [BDGL16] Becker, A., Ducas, L., Gama, N., & Laarhoven, T.  (2016).  New directions in
        nearest neighbor searching with applications to lattice sieving.  In SODA 2016, (pp.
        10–24).

        """
        # TODO we simply pick the same additive constant 16.4 as for the experimental result in [BDGL16]
        return BKZ.LLL(d, B) + ZZ(2) ** RR(0.292 * beta + 16.4 + log(BKZ.svp_repeat(beta, d), 2))

    @staticmethod
    def BDGL16(beta, d, B=None):
        """
        Runtime estimation given `β` and assuming sieving is used to realise the SVP oracle.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [BDGL16] Becker, A., Ducas, L., Gama, N., & Laarhoven, T.  (2016).  New directions in
            nearest neighbor searching with applications to lattice sieving.  In SODA 2016, (pp.
            10–24).
        """
        # TODO this is somewhat arbitrary
        if beta <= 90:
            return BKZ._BDGL16_small(beta, d, B)
        else:
            return BKZ._BDGL16_asymptotic(beta, d, B)

    @staticmethod
    def LaaMosPol14(beta, d, B=None):
        """
        Runtime estimation for quantum sieving.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [LaaMosPol14] Thijs Laarhoven, Michele Mosca, & Joop van de Pol.  Finding shortest
            lattice vectors faster using quantum search.  Cryptology ePrint Archive, Report
            2014/907, 2014.  https://eprint.iacr.org/2014/907.

        ..  [Laarhoven15] Laarhoven, T.  (2015).  Search problems in cryptography: from
            fingerprinting to lattice sieving (Doctoral dissertation).  Eindhoven University of
            Technology. http://repository.tue.nl/837539
        """
        return BKZ.LLL(d, B) + ZZ(2) ** RR((0.265 * beta + 16.4 + log(BKZ.svp_repeat(beta, d), 2)))

    @staticmethod
    def CheNgu12(beta, d, B=None):
        """
        Runtime estimation given β and assuming [CheNgu12]_ estimates are correct.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        The constants in this function were derived as follows based on Table 4 in
        [CheNgu12]_::

            >>> from sage.all import var, find_fit
            >>> dim = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250]
            >>> nodes = [39.0, 44.0, 49.0, 54.0, 60.0, 66.0, 72.0, 78.0, 84.0, 96.0, 99.0, 105.0, 111.0, 120.0, 127.0, 134.0]  # noqa
            >>> times = [c + log(200,2).n() for c in nodes]
            >>> T = list(zip(dim, nodes))
            >>> var("a,b,c,beta")
            (a, b, c, beta)
            >>> f = a*beta*log(beta, 2.0) + b*beta + c
            >>> f = f.function(beta)
            >>> f.subs(find_fit(T, f, solution_dict=True))
            beta |--> 0.2701...*beta*log(beta) - 1.0192...*beta + 16.10...

        The estimation

        2^(0.270188776350190*beta*log(beta) - 1.0192050451318417*beta + 16.10253135200765)

        is of the number of enumeration nodes, hence we need to multiply by the number of
        cycles to process one node.  This cost per node is typically estimated as 100 [FPLLL]_.

        ..  [CheNgu12] Yuanmi Chen and Phong Q.  Nguyen.  BKZ 2.0: Better lattice security
        estimates (Full Version). 2012. http://www.di.ens.fr/~ychen/research/Full_BKZ.pdf

        ..  [FPLLL] The FPLLL development team.  fplll, a lattice reduction library.  2016.
        Available at https://github.com/fplll/fplll
        """
        repeat = BKZ.svp_repeat(beta, d)
        cost = RR(
            0.270188776350190 * beta * log(beta)
            - 1.0192050451318417 * beta
            + 16.10253135200765
            + log(100, 2)
        )
        return BKZ.LLL(d, B) + repeat * ZZ(2) ** cost

    @staticmethod
    def ABFKSW20(beta, d, B=None):
        """
        Enumeration cost according to [ABFKSW20]_.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [ABFKSW20] Martin R. Albrecht, Shi Bai, Pierre-Alain Fouque, Paul Kirchner, Damien
        Stehlé and Weiqiang Wen.  Faster Enumeration-based Lattice Reduction: Root Hermite
        Factor $beta^{1/(2k)}$ in Time $beta^{beta/8 + o(beta)}$.  CRYPTO 2020.
        """
        if 1.5 * beta >= d or beta <= 92:  # 1.5β is a bit arbitrary, β≤92 is the crossover point
            cost = RR(0.1839 * beta * log(beta, 2) - 0.995 * beta + 16.25 + log(100, 2))
        else:
            cost = RR(0.125 * beta * log(beta, 2) - 0.547 * beta + 10.4 + log(100, 2))
        repeat = BKZ.svp_repeat(beta, d)

        return BKZ.LLL(d, B) + repeat * ZZ(2) ** cost

    @staticmethod
    def ADPS16(beta, d, B=None, mode="classical"):
        """
        Runtime estimation from [ADPS16]_.

        :param beta: block size
        :param d: lattice dimension
        :param B: bit-size of entries

        ..  [ADPS16] Edem Alkim, Léo Ducas, Thomas Pöppelmann, & Peter Schwabe (2016).
            Post-quantum key exchange - A New Hope.  In T. Holz, & S. Savage, 25th USENIX
            Security Symposium, USENIX Security 16 (pp. 327–343). USENIX Association.
        """

        if mode not in ("classical", "quantum", "paranoid"):
            raise ValueError(f"Mode {mode} not understood.")

        c = {
            "classical": 0.2920,
            "quantum": 0.2650,  # paper writes 0.262 but this isn't right, see above
            "paranoid": 0.2075,
        }

        c = c[mode]

        return ZZ(2) ** RR(c * beta)

    @classmethod
    def cost(cls, cost_model, beta, d, B=None, predicate=None, **kwds):
        """
        Return cost dictionary for returning vector of norm` δ_0^{d-1} Vol(Λ)^{1/d}` using provided lattice
        reduction algorithm.

        :param lattice_reduction_estimate:
        :param beta: block size ≥ 2
        :param d: lattice dimension
        :param B: bit-size of entries
        :param predicate: if ``False`` cost will be infinity.

        """
        from .cost import Cost

        cost = cost_model(beta, d, B)
        delta = cls.delta(beta)
        cost = Cost(rop=cost, red=cost, delta=delta, beta=beta, d=d, **kwds)
        Cost.register_impermanent(rop=True, red=True, delta=False, beta=False, d=False)
        if predicate is not None and not predicate:
            cost["red"] = oo
            cost["rop"] = oo
        return cost


BKZ.default = BKZ.BDGL16
BKZ.classical_poly_space = BKZ.ABFKSW20
