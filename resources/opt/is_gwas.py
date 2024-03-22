from pyspark.sql.functions import udf
from pyspark.sql.types import StructType, StructField, FloatType
import math
import scipy.stats as stats
import numpy as np
from firth_correction import Ka, Kb, Kaa, Kab, Kbb
from typing import List, Tuple, Literal

def check_input_vars(N: float, Nstar: float, MAF: float, MAFstar: float) -> float:
    return (N / (N - Nstar)) * (MAF - Nstar / N * MAFstar)

def check_initial_values(alpha: float, beta: float, g: int) -> float:
    return math.exp(alpha) / (1 + math.exp(alpha)) * abs(math.exp(beta * g) - 1)

def pr_g_under_hwe(maf: float, g: int) -> float:
    '''
    Probability of observing genotype g under HWE
    '''
    g0 = (1 - maf) ** 2
    g1 = 2 * maf * (1 - maf)
    g2 = maf ** 2
    return [g0, g1, g2][g]

def pr_case(alpha: float, beta: float, g: int) -> float:
    '''
    Probability of observing a case, given alpha, beta, and g 
    '''
    return (math.exp(alpha + beta * g)) / (1 + math.exp(alpha + beta * g))

def e_int(alpha: float, beta: float, cexp: int, maf: float, g: int) -> float:
    pr_y = pr_case(alpha = alpha, beta = beta, g = g)
    return pr_y * (1 - pr_y) * g**cexp * pr_g_under_hwe(maf = maf, g = g)

def compute_e(alpha: float, beta: float, cexp: int, maf: float) -> float:
    '''
    Compute e^{c, k}, an entry in the Fisher Information matrix for the kth iteration
    (see lines 485-486 from isGWAS preprint)
    '''
    return sum((e_int(alpha = alpha, beta = beta, cexp = cexp, maf = maf, g = 0),
                e_int(alpha = alpha, beta = beta, cexp = cexp, maf = maf, g = 1),
                e_int(alpha = alpha, beta = beta, cexp = cexp, maf = maf, g = 2)))

def tilde_e_int(alpha: float, beta: float, cexp: int, maf: float, g: int) -> float:
    pr_y = pr_case(alpha = alpha, beta = beta, g = g)
    return pr_y * g**cexp * pr_g_under_hwe(maf, g)

def compute_tilde_e(alpha: float, beta: float, cexp: int, maf: float) -> float:
    '''
    Compute \\tilde{e}^{c, k}, a component of the Score function for the kth iteration
    (see lines 489-490 from isGWAS preprint)
    '''
    return sum((tilde_e_int(alpha = alpha, beta = beta, cexp = cexp, maf = maf, g = 0),
                tilde_e_int(alpha = alpha, beta = beta, cexp = cexp, maf = maf, g = 1),
                tilde_e_int(alpha = alpha, beta = beta, cexp = cexp, maf = maf, g = 2)))

def salnr(alpha: float, beta: float, MAF: float, MAFstar: float, N: float, Nstar: float, firth: bool) -> dict:
    '''
    SaLN-R update function. See lines 507-508 from isGWAS preprint
    '''
    e0 = compute_e(alpha = alpha, beta = beta, cexp = 0, maf = MAF)
    e1 = compute_e(alpha = alpha, beta = beta, cexp = 1, maf = MAF)
    e2 = compute_e(alpha = alpha, beta = beta, cexp = 2, maf = MAF)
    tilde_e0 = compute_tilde_e(alpha = alpha, beta = beta, cexp = 0, maf = MAF)
    tilde_e1 = compute_tilde_e(alpha = alpha, beta = beta, cexp = 1, maf = MAF)
    if firth:
        pr_y = Nstar / N
        K_a = Ka(alpha, beta, MAF)
        K_b = Kb(alpha, beta, MAF)
        K_a_wrt_a = Kaa(alpha, beta, MAF)
        K_a_wrt_b = Kab(alpha, beta, MAF)
        K_b_wrt_b = Kbb(alpha, beta, MAF)
        scalar_a = float(1 / (((e0 - (K_a_wrt_a / N)) * (e2 - (K_b_wrt_b / N)) - (e1 - (K_a_wrt_b / N))**2)))
        matrix_b = np.array([[e2 - (K_a_wrt_a / N),
                              -e1 + (K_a_wrt_b / N)],
                             [-e1 + (K_a_wrt_b / N),
                             e0 - (K_b_wrt_b / N)]], dtype = float)
        vector_c = np.array([pr_y - tilde_e0 + (K_a / N),
                             (2 * pr_y * MAFstar) - tilde_e1 + (K_b / N)], dtype = float)
    else:
        scalar_a = 1 / (N * (e0 * e2 - (e1)**2))
        matrix_b = np.array([[e2, -e1], [-e1, e0]])
        vector_c = np.array([Nstar - (N * tilde_e0), (2 * Nstar * MAFstar) - (N * tilde_e1)])
    sigma = np.array([alpha, beta]) + np.array(scalar_a * np.matmul(matrix_b, vector_c))
    sigma_se = np.array([math.sqrt(e2 / (N * (e0 * e2 - e1**2))),
                         math.sqrt(e0 / (N * (e0 * e2 - e1**2)))])
    return dict(estimates = sigma, se = sigma_se)

def lrt(N: int, Nstar: int, MAFstar: float, MAF: float, alpha: float, beta: float) -> Tuple[float, float]:
    '''
    Liklihood ratio test statistic and p-value
    '''
    alpha_null = math.log(Nstar / (N - Nstar))
    pr_case0 = pr_case(alpha, beta, 0)
    pr_case1 = pr_case(alpha, beta, 1)
    pr_case2 = pr_case(alpha, beta, 2)
    pr_g0 = pr_g_under_hwe(MAF, 0)
    pr_g1 = pr_g_under_hwe(MAF, 1)
    pr_g2 = pr_g_under_hwe(MAF, 2)
    pr_case_null = pr_case(alpha_null, 0, 0)
    summ_a = math.log(1 - pr_case0) * pr_g0 + math.log(1 - pr_case1) * pr_g1 + math.log(1 - pr_case2) * pr_g2
    lrtstat = 2 * (Nstar * (alpha - alpha_null + (beta * 2 * MAFstar))) + N * (summ_a - math.log(1 - pr_case_null))
    pval = stats.chi2.sf(lrtstat, df = 1)
    pval = pval if pval > 0 else np.nextafter(0, 1)
    return lrtstat, pval

def wald(beta: float, se: float) -> Tuple[float, float]:
    '''
    Wald test statistic and p-value
    '''
    w_stat = beta**2 / se**2
    pval = stats.chi2.sf(w_stat, df = 1)
    pval = pval if pval > 0 else np.nextafter(0, 1)
    return w_stat, pval

def isGWAS(N_cases: float,
           N_controls: float,
           AF_cases: float,
           AF_controls: float,
           maf_min: float,
           standardized: bool,
           firth: bool,
           test: Literal['wald', 'lrt'],
           e_max: float,
           tau: int) -> Tuple[float, float, float]:
    # Step 0: Derive inputs to algorithm
    N = N_cases + N_controls
    Nstar = N_cases
    AF = (Nstar / N) * AF_cases + ((1 - Nstar / N) * AF_controls)

    # For the effect allele to always be the minor allele (as in the isGWAS preprint):
    # if AF > 0.5:
    #     MAF = 1 - AF
    #     MAFstar = 1 - AF_cases
    # else:
    #     MAF = AF
    #     MAFstar = AF_cases

    # For the effect allele to always be the "Alternative" allele:
    MAF = AF
    MAFstar = AF_cases

    if (MAF < maf_min) or ((1 - MAF) < maf_min):
        return np.nan, np.nan, np.nan

    # Step 1: Get initial values
    alpha0 = math.log(Nstar / (N - Nstar))
    beta_tmp1 = (MAFstar - MAF) / (MAF * (1 - MAF)) * (N / (N - Nstar))
    beta_tmp2 = math.log(2 + math.exp(-alpha0))
    beta0 = beta_tmp1 if check_initial_values(alpha0, beta_tmp1, 1) <= 1 else beta_tmp2
    if (beta0 > 0) & (beta0 > beta_tmp2):
        beta0 = beta_tmp2

    # Step 3: Run iterations
    i = 0
    k = 0
    delta_e = 0.1
    while (k <= tau) & (delta_e > e_max):
        if (k == 0):
            alpha = alpha0
            beta = beta0
        else:
            alpha = sigma['estimates'][0]
            beta = sigma['estimates'][1]
        sigma = salnr(alpha = alpha,
                      beta = beta,
                      MAF = MAF,
                      MAFstar = MAFstar,
                      N = N,
                      Nstar = Nstar,
                      firth = firth)
        if any(np.isnan(val) for val in sigma['estimates']) | any(np.isinf(val) for val in sigma['estimates']):
            i = i + 1
            sigma['estimates'] = np.array([alpha0 / (2 * i), math.log(2 + math.exp(-alpha0)) / 2 ** i])
        delta_e = max(abs(sigma['estimates'] - np.array([alpha, beta])))
        k = k + 1
    if test == 'wald':
        test_stat, pval = wald(beta = sigma['estimates'][1], se = sigma['se'][1])
    else:
        test_stat, pval = lrt(N, Nstar, MAFstar, MAF, sigma['estimates'][0], sigma['estimates'][1])
    if standardized:
        beta = float(sigma['estimates'][1] * (math.sqrt(2 * MAF * (1 - MAF))))
        se = float(sigma['se'][1] * (math.sqrt(2 * MAF * (1 - MAF))))
    else:
        beta = float(sigma['estimates'][1])
        se = float(sigma['se'][1])
    return beta, se, float(-math.log(pval, 10))

isGWAS_udf = udf(isGWAS, StructType([
    StructField("Beta", FloatType(), False),
    StructField("SE", FloatType(), False),
    StructField("log10p", FloatType(), False)
]))