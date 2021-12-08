from .conf import red_cost_model_default, red_shape_model_default
from .lwe_dual import dual, dual_hybrid, dual_mitm_hybrid
from .lwe_primal import primal_usvp, primal_hybrid, primal_bdd
from .gb import arora_gb
from .lwe_brute_force import exhaustive_search, mitm


def estimate_lwe(params,
                 red_shape_model=red_shape_model_default,
                 red_cost_model=red_cost_model_default,
                 skip={arora_gb, exhaustive_search, mitm},
                 log_level=1):
    """
    A very basic top-level function which prints a variety of cost estimates
    :param params: the input LWE parameters
    :param red_shape_model: the reduction shape, i.e. GSA or simulator
    :param red_cost_model: the cost model for BKZ

    EXAMPLE:
    >>> from estimator import *
    >>> _ = estimate_lwe(LWEParameters(n=200, q=127, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3)),)
    >>> _["dual"]
    rop: ≈2^109.2, mem: ≈2^77.3, m: 234, red: ≈2^108.6, δ: 1.005166, β: 272, d: 433, ↻: ≈2^75.8, tag: dual
    >>> _["dual_hybrid"]
    rop: ≈2^102.0, mem: ≈2^99.0, m: 219, red: ≈2^101.6, δ: 1.005500, β: 247, d: 398, ↻: ≈2^65.7, ...
    >>> _["dual_mitm_hybrid"]
    rop: ≈2^228.0, mem: ≈2^228.1, m: 35, k: 139, ↻: 228, red: ≈2^46.5, δ: 1.010545, β: 75, d: 38, ...
    >>> _["primal_usvp"]
    rop: ≈2^90.1, red: ≈2^90.1, δ: 1.006187, β: 205, d: 378, tag: usvp
    >>> _["primal_bdd"]
    rop: ≈2^86.2, red: ≈2^84.6, svp: ≈2^85.7, β: 185, η: 216, d: 370, tag: bdd
    >>> _["primal_hybrid"]
    rop: ≈2^114.6, red: ≈2^113.6, svp: ≈2^113.5, β: 191, η: 2, ζ: 87, |S|: ≈2^137.9, d: 302, ...
    """

    # all implemented algorithms
    all_algs = [
        arora_gb,
        exhaustive_search,
        mitm,
        dual,
        dual_hybrid,
        dual_mitm_hybrid,
        primal_bdd,
        primal_usvp,
        primal_hybrid]

    # all not-skipped algorithms
    algs = [alg for alg in all_algs if alg not in skip]

    # a dictionary to store the costs
    costs = dict()

    for alg in algs:
        try:
            if alg.__name__[0:4] == "dual":
                # evaluate dual attacks
                est = alg(params, red_cost_model=red_cost_model)
                costs["{}".format(alg.__name__)] = est
            elif alg.__name__[0:6] == "primal":
                # evaluate the primal attacks
                est = alg(params, red_cost_model=red_cost_model, red_shape_model=red_shape_model)
                costs["{}".format(alg.__name__)] = est
            else:
                # evaluate all other attacks
                est = alg(params)
                costs["{}".format(alg.__name__)] = est
        except:  # noqa
            # skip any estimate that crashes
            pass

    return costs
