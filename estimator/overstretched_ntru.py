#   Instruction to generate the predictions in the Figures.
#   Fig 1.
#   Left: 
#       combined_attack_prob(q, 127, 2/3., ntru="matrix", fixed_tours=8, only="SKR")
#       for several q
#   Right: combined_attack_prob(q, 127, 2/3., ntru="matrix", fixed_tours=8)
#       for q in [601, 739, 1373]
#   Fig 3.
#       DSLI_vols(dsl_logvol, FL_shape)
#       where dsl_logvol is the dense sublattice log-volume and FL_shape is the log-profile from the actual experiments 
#   Fig 4.
#       DSL_logvol_matrix(n, 2/3.)
#       for n in primes(70, 200)
#   Fig 5.
#       DSL_logvol_circulant(n, 2/3.)
#       for n in primes(70, 200)
#   Fig 6.
#       combined_attack_prob(q, 127, 2/3., ntru="matrix", fixed_tours=8)
#       combined_attack_prob(q, 127, 2/3., ntru="matrix", fixed_tours=8, only="SKR")
#       combined_attack_prob(q, 127, 2/3., ntru="matrix", fixed_tours=8, only="DSD")
#       for several q
#   Fig 7.
#       combined_attack_prob(2003, 151, 2/3., ntru="circulant", fixed_tours=8)
#       for the bottom right plot we used the experimental values of the sublattice volumes and used
#       combined_attack_prob(2003, 151, 2/3., ntru="circulant", fixed_tours=8, dsl_logvol=...)
#   Fig 8.
#       find_fatigue(n, 2/3., ntru="matrix", fixed_tours=8, DSD_ratio=ratio)
#       for n in primes(90, 500) and ratio in [0.01, 0.5, .99]

from sage.all import cached_function, RealDistribution, minimize, log, exp, floor, pi, RR
from math import lgamma
from scipy.special import digamma, zeta, gamma
import numpy as np

from .conf import ntru_fatigue_lb, ntru_fatigue_ub

max_n_cache = 10000

@cached_function
def ball_log_vol(n):
    return RR((n/2.)) * RR(log(pi)) - RR(lgamma(n/2. + 1))

gh_constant = {1:0.00000,2:-0.50511,3:-0.46488,4:-0.39100,5:-0.29759,6:-0.24880,7:-0.21970,8:-0.15748,9:-0.14673,10:-0.07541,11:-0.04870,12:-0.01045,13:0.02298,14:0.04212,15:0.07014,16:0.09205,17:0.12004,18:0.14988,19:0.17351,20:0.18659,21:0.20971,22:0.22728,23:0.24951,24:0.26313,25:0.27662,26:0.29430,27:0.31399,28:0.32494,29:0.34796,30:0.36118,31:0.37531,32:0.39056,33:0.39958,34:0.41473,35:0.42560,36:0.44222,37:0.45396,38:0.46275,39:0.47550,40:0.48889,41:0.50009,42:0.51312,43:0.52463,44:0.52903,45:0.53930,46:0.55289,47:0.56343,48:0.57204,49:0.58184,50:0.58852}
def log_gh(d, logvol=0):
    if d < 49:
        return RR(gh_constant[d]) + RR(logvol)/d

    return RR(1./d) * RR(logvol - ball_log_vol(d))


def delta(k):
    assert(k>=60)
    delta = exp(log_gh(k)/(k-1))
    return RR(delta)

small_slope_t8 = {2:0.04473,3:0.04472,4:0.04402,5:0.04407,6:0.04334,7:0.04326,8:0.04218,9:0.04237,10:0.04144,11:0.04054,12:0.03961,13:0.03862,14:0.03745,15:0.03673,16:0.03585,17:0.03477,18:0.03378,19:0.03298,20:0.03222,21:0.03155,22:0.03088,23:0.03029,24:0.02999,25:0.02954,26:0.02922,27:0.02891,28:0.02878,29:0.02850,30:0.02827,31:0.02801,32:0.02786,33:0.02761,34:0.02768,35:0.02744,36:0.02728,37:0.02713,38:0.02689,39:0.02678,40:0.02671,41:0.02647,42:0.02634,43:0.02614,44:0.02595,45:0.02583,46:0.02559,47:0.02534,48:0.02514,49:0.02506,50:0.02493,51:0.02475,52:0.02454,53:0.02441,54:0.02427,55:0.02407,56:0.02393,57:0.02371,58:0.02366,59:0.02341,60:0.02332}

@cached_function
def slope(beta):
    if beta<=60:
        return RR(small_slope_t8[beta])
    if beta<=70:
        # interpolate between experimental and asymptotics
        ratio = RR((70-beta)/10.)
        return RR(ratio*small_slope_t8[60])+RR(1.-ratio)*2*RR(log(delta(70)))
    else:
        return 2 * RR(log(delta(beta)))

def chi2_CDF(n, x):
    if x > 100 * n:
        return 1.
    return RR(1.) - RR(gamma(n/2., x/2.))/RR(gamma(n/2.))

chisquared_table = {i: None for i in range(2*max_n_cache+1)}
for i in range(2*max_n_cache+1):
    chisquared_table[i] = RealDistribution('chisquared', i)


def conditional_chi_squared(d1, d2, lt, l2):
    """
    Probability that a gaussian sample (var=1) of dim d1+d2 has length at most
    lt knowing that the d2 first cordinates have length at most l2
    """
    D1 = chisquared_table[d1].cum_distribution_function
    D2 = chisquared_table[d2].cum_distribution_function
    l2 = RR(l2)

    PE2 = D2(l2)
    # In large dim, we can get underflow leading to NaN
    # When this happens, assume lifting is successfully (underestimating security)
    if PE2==0:
        raise ValueError("Numerical underflow in conditional_chi_squared")

    steps = 5 * (d1 + d2)

    # Numerical computation of the integral
    proba = 0.
    for i in range(steps)[::-1]:
        l2_min = i * l2 / steps
        l2_mid = (i + .5) * l2 / steps
        l2_max = (i + 1) * l2 / steps

        PC2 = (D2(l2_max) - D2(l2_min)) / PE2
        PE1 = D1(lt - l2_mid)

        proba += PC2 * PE1

    return proba

def binary_search_min_sucess(f, xmin, xmax):

    if xmax - xmin < 2: 
        return xmin+1, f(xmin+1)

    xmid = floor((xmax + xmin)/2)
    if not f(xmid):
        return binary_search_min_sucess(f, xmid, xmax)
    return binary_search_min_sucess(f, xmin, xmid)

def zshape(q, n, beta):
    logq = RR(log(q))
    L = n*[logq] + n * [0]
    slope_ = slope(beta)
    diff = slope(beta)/2.

    for i in range(n):
        if diff > logq/2.: break
        L[n-i-1] = logq/2. + diff
        L[n+i  ] = logq/2. - diff

        diff += slope_
    return L

# log loss of length when projecting out k dimension out of d 
@cached_function
def proj_logloss(d, k):
    return (RR(digamma((d-k)/2.))-RR(digamma(d/2.)))/2.

def DSL_logvol_matrix(n,sigmasq):
    total = n*(RR(log(sigmasq))+RR(log(2.))+RR(digamma(n)))/2.
    proj_loss = np.sum([(digamma((2*n-i)/2.)-digamma(n)) for i in range(n)])/2.
    return total+proj_loss

def DSL_logvol_circulant(n,sigmasq):
    lambda0 = RR((np.log(2)-np.euler_gamma+np.log(n)+np.log(sigmasq))/2.)
    lambdai = (n-1)*(1-np.euler_gamma+np.log(n)+np.log(sigmasq))/2.
    return lambda0+lambdai

def DSL_logvol_circulant_fixed(n,R):
    lambda0 = (-np.euler_gamma+np.log(R))/2.
    lambdai = (n-1)*(1-np.euler_gamma+np.log(R)-np.log(2))/2.
    return lambda0+lambdai

@cached_function
def DSL_logvol(n,sigmasq,ntru="circulant"):
    if ntru=="matrix":
        return DSL_logvol_matrix(n,sigmasq)
    if ntru=="circulant":
        return DSL_logvol_circulant(n,sigmasq)
    if ntru=="fixed":
        return DSL_logvol_circulant_fixed(n,sigmasq)
    print("non implemented ntru type")

@cached_function
def zeta_prime(x):
    h = 1e-5
    return (zeta(x+h,1) - zeta(x-h,1))/(2*h)

zeta_precomputed = [zeta(i) for i in range(max_n_cache+1)]
zeta_prime_precomputed = [zeta_prime(i) for i in range(max_n_cache+1)]

def DSLI_vols(dsl_logvol, FL_shape):
    n = len(FL_shape)//2
    vols = (2*n+1)*[None]

    dsl_dim = n
    vols[2*n] = dsl_logvol
    
    # Going to a intersection of dimension s
    for s in range(2*n-1, n, -1):
        # Negate cause it's a dual vector really
        x = - FL_shape[s]
        x += proj_logloss(s+1,n)
        x += zeta_prime_precomputed[dsl_dim]/zeta_precomputed[dsl_dim] # primitivity
        dsl_logvol += x
        vols[s] = dsl_logvol
        dsl_dim -= 1

    assert(dsl_dim==1)
    assert(s==n+1)
    
    return vols

def SKR_attack_prob(beta, q, n, sk_variance=2/3.):    
    logq = RR(log(q))
    slope_ = slope(beta)
    threshold = RR(log(beta * sk_variance)/2.)

    def score(m):
         return - (RR(m[0]*logq/(n+m[0])) - slope_ * RR((m[0]+n-1)/2. - beta))
    m = minimize(score, [n])[0]
    m = floor(m)
    m = min(m, n)
    m = max(m, 1)

    proba_one = 1.
    for b in range(beta, min(max(2*n, 300), 3*beta), beta-1):
        threshold = RR(log(b * sk_variance)/2.)
        log_len_gs = RR(m*logq/(n+m)) - slope_ * RR((m+n-1)/2. - b)
        bound_chi2 = b*exp(2*(log_len_gs - threshold))
        proba_one *= chi2_CDF(b, bound_chi2)
    prob_all_not = (1.-proba_one)**n

    prob_pos = np.zeros(2*n, dtype='double')
    prob_pos[2*n-beta-1] = 1.-prob_all_not
    return RR(1.-prob_all_not), prob_pos

def prob_add(x,y):
    return RR(1.-(1.-x)*(1.-y))

def DSD_attack_prob(beta, q, n, sk_variance=2/3., ntru="matrix", dsl_logvol=None, zshapef=None):    
    if dsl_logvol==None:
        dsl_logvol = DSL_logvol(n, sk_variance, ntru=ntru)
    
    # breakpoint()
    if zshapef is None:
        B_shape = zshape(q, n, beta)

    else:
        # zshapef comes from estimator API, need to take the sqrt and log.
        B_shape = [log(r_)/2 for r_ in zshapef(q=q, n=n, beta=beta)]

    dsli_vols = DSLI_vols(dsl_logvol, B_shape)

    prob_all_not = RR(1.)
    prob_pos = np.zeros(2*n, dtype='double')
    for i in range(1, n+1):
        s = n + i

        dslv_len = log_gh(i, dsli_vols[s])
        sigma_sq = exp(2*dslv_len)/s
        # if beta == 40:
            # print(sigma_sq)

        # print(sigma_sq)
        if sigma_sq > 10**10:
            prob_pos[s-beta] = 0.
            continue

        norm_threshold = exp(2*(B_shape[s-beta]))/sigma_sq
        proba_one = chisquared_table[beta].cum_distribution_function(norm_threshold)

        if proba_one <= 10e-8:
            continue

        # account for pulling back probability if beta small
        if beta <= 20:
            for j in range(2, int(s/beta+1)):
                if proba_one < 10**(-6):
                    proba_one = 0.
                    break
                ind = s - j*(beta-1)-1
                norm_bt = exp(2*B_shape[ind])/sigma_sq
                norm_b2 = exp(2*B_shape[ind+beta-1])/sigma_sq
                proba_one *= conditional_chi_squared(beta-1,s-ind-(beta-1), norm_bt, norm_b2)

        prob_pos[s-beta] = proba_one
        prob_all_not *= max(1.-proba_one, 0.)
    return RR(1.-prob_all_not), prob_pos

# To obtain an estimate for circulant NTRU with parameters q, n, and variance ss, using t tours run
# combined_attack_prob(q, n, ss, ntru="circulant", fixed_tours=t)
# It returns
#   - The expected successful beta
#   - The probability that the SKR event is triggered first
#   - The probability that the DSD event is triggered first
#   - The distribution of the event positions kappa.
#
# For matrix NTRU change to ntru="matrix"
# To only account for the SKR or DSD event use only="SKR" or only="DSD"
# To run for a specific sublattice volume supply dsl_logvol=...
def combined_attack_prob(q, n, sk_variance=2/3., ntru="matrix", fixed_tours=None, only=None, verbose=True, dsl_logvol=None, zshapef=None):
    if n > max_n_cache:
        print("Please increase the hardcoded value of max_n_cache to run the predictor for such large n")
        return

    remaining_proba = RR(1.)
    average_beta = RR(0.)
    total_SKR_prob = RR(0.)
    SKR_prob = RR(0.)
    total_DSD_prob = RR(0.)
    DSD_prob = RR(0.)
    prob_pos_total = np.zeros(2*n, dtype='double')

    for beta in range(2,n):
        if fixed_tours is None:
            tours = floor(n**2 / beta**2)+3
        else:
            tours = fixed_tours

        if only!="DSD":
            SKR_prob, SKR_prob_pos = SKR_attack_prob(beta, q, n, sk_variance)
        if only!="SKR":
            DSD_prob, DSD_prob_pos = DSD_attack_prob(beta, q, n, sk_variance, ntru, dsl_logvol=dsl_logvol, zshapef=zshapef)

        if SKR_prob > 10e-8 or DSD_prob > 10e-8:
            for t in range(tours):
                for i in range(2*n):

                    if only!="DSD":
                        # SKR
                        prob_pos = SKR_prob_pos[i]
                        average_beta += RR(beta) * remaining_proba * prob_pos
                        prob_pos_total[i] += remaining_proba * prob_pos
                        total_SKR_prob += remaining_proba * prob_pos
                        remaining_proba *= (1.-prob_pos)

                    if only!="SKR":
                        # DSD
                        prob_pos = DSD_prob_pos[i]
                        average_beta += RR(beta) * remaining_proba * prob_pos
                        prob_pos_total[i] += remaining_proba * prob_pos
                        total_DSD_prob += remaining_proba * prob_pos
                        remaining_proba *= (1.-prob_pos)
        if verbose:
            print("β= %d,\t pr=%.4e, \t rem-pr=%.4e"%(beta, prob_add(SKR_prob, DSD_prob), remaining_proba))
        if remaining_proba < 0.001:
            average_beta += beta * remaining_proba
            break

    return average_beta, total_SKR_prob, total_DSD_prob, prob_pos_total

def find_fatigue(n, sk_variance=2/3., ntru="circulant", fixed_tours=None, DSD_ratio=0.5):
    if n > max_n_cache:
        print("Please increase the hardcoded value of max_n_cache to run the predictor for such large n")
        return

    print(n)
    def dsl_wins(q):
        average_beta, p_SKR, p_DSD, pos = combined_attack_prob(q, n, sk_variance, ntru, fixed_tours, verbose=False)
        return p_DSD >= DSD_ratio
    return binary_search_min_sucess(dsl_wins, ntru_fatigue_lb(n), ntru_fatigue_ub(n))[0]
