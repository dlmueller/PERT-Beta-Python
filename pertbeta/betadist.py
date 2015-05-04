# -*- coding: utf-8 -*-
"""
Numerische Approximatino der GAMMA- sowie BETA-Funktion,
numerische Interation, um die BETA-Dichtefunktion (pdf)
und BETA-Verteiltungsfunktion (CDF) zu berechnen.

@copyright: David Lukas Müller (2013, 2015), Lanczos approximation from Wikipedia
"""

#---
#--- Python
from cmath import *
import doctest
import random

#---
MAX_DOTS = 40

#---
def gamma_lanczos(z):
    """
    Gamma function with the Lanczos approximation (found on Wikipedia)
    https://en.wikipedia.org/wiki/Lanczos_approximation
    """
    g = 7
    lanczos_coef = [ \
         0.99999999999980993,
       676.5203681218851,
     -1259.1392167224028,
       771.32342877765313,
      -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
         9.9843695780195716e-6,
         1.5056327351493116e-7]

    z = complex(z)
    if z.real < 0.5:
        return pi / (sin(pi*z)*gamma(1-z))
    else:
        z -= 1
        x = lanczos_coef[0] + \
            sum(lanczos_coef[i]/(z+i)
                for i in range(1, g+2))
        t = z + g + 0.5
        return sqrt(2*pi) * t**(z+0.5) * exp(-t) * x



def gamma(x) :
    """
    Testing against a known value (Γ(10) = 9! = 362880)::
        >>> int(gamma(10.0)) # 362880.00000000047+0j)
        362880
    """
    z = gamma_lanczos(x)
    return z.real


def INTEGRAL(f, lo, hi, epsilon = 0.001) :
    """
    Doctests::
        >>> eps = 0.00001
        >>> INTEGRAL(lambda x : x, 0, 1, epsilon = eps) - 0.5 < eps
        True
    """
    summe = 0.0
    for (A, tr) in iterIntegral(f, lo, hi, epsilon = epsilon) :
        summe += A
    return summe

def iterIntegral(f, lo, hi, epsilon = 0.001) :
    if lo >= hi :
        return
    t = lo
    while t <= hi :
        tr = min(hi, t + epsilon)
        y = f(t)
        yr = f(tr)
        m = (y + yr)/2.0
        A = m*(tr-t)
        yield (A, tr)
        if tr >= hi :
            break
        t = tr
        continue
    return

def BETA_GAMMA(alpha, beta, a, b) :
    gamma_quot = (gamma(alpha)*gamma(beta)/gamma(alpha+beta))
    return gamma_quot #* (b-a)**(alpha + beta - 1)

def BETA_INTEGRAL(alpha, beta, a, b, epsilon = 0.001) :
    f = lambda t : beta_pdf_nominator(t, alpha, beta, a, b)
    return INTEGRAL(f, a, b, epsilon = epsilon)

def BETA(alpha, beta, a, b) :
    return BETA_GAMMA(alpha, beta, a, b)
    # return BETA_INTEGRAL(alpha, beta, a, b)

beta_cache = {}
def BETA_CACHED(alpha, beta, a, b) :
    key = (alpha, beta, a, b)
    if not beta_cache.has_key(key) :
        B = BETA(alpha, beta, a, b)
        beta_cache[key] = B
    return beta_cache[key]

def beta_pdf_nominator(x, alpha, beta, a, b, infinite = 10.0) :
    """
    Notes
    -----
    The probability density function for `beta` is::

         beta.pdf(x, a, b) = gamma(a+b)/(gamma(a)*gamma(b)) * x**(a-1) *
        (1-x)**(b-1),

    (taken from https://github.com/scipy/scipy/blob/v0.11.0/scipy/stats/distributions.py)
    @param a, b, x: lower and upper bounds with a <= x <= b
    @param alpha, beta: shape parameters with alpha, beta > 0
    """
    NULL = 0.0
    if x < a or x > b:
        return NULL
    if alpha > 0 and beta > 0 :
        try :
            num = (x-a)**(alpha-1) * (b-x)**(beta-1)
        except ZeroDivisionError as E :
            num = infinite
        return num
    return NULL

def beta_pdf(x, alpha, beta, a, b, infinite = 10.0) :
    """
    Notes
    -----
    The probability density function for `beta` is::

         beta.pdf(x, a, b) = gamma(a+b)/(gamma(a)*gamma(b)) * x**(a-1) *
        (1-x)**(b-1),

    (taken from https://github.com/scipy/scipy/blob/v0.11.0/scipy/stats/distributions.py)
    @param a, b, x: lower and upper bounds with a <= x <= b
    @param alpha, beta: shape parameters with alpha, beta > 0
    """
    NULL = 0.0
    if x < a or x > b:
        return NULL
    if alpha > 0 and beta > 0 :
        num = beta_pdf_nominator(x, alpha, beta, a, b, infinite = infinite)
        den = BETA_CACHED(alpha,beta,a,b)*(b-a)**(alpha+beta-1)
        return num/den
    return NULL

def beta_cdf(x, alpha, beta, a, b, epsilon = 0.001) :
    """
    BETA(3.0, 4.0, 20.0) alpha=0.8 beta=3.2
    mean1 = 6.5 var = 8.0 sigma1 = 2.8
    mean2 = 6.6 var = 10.1 sigma2 = 3.2 (N = 1000)
    q25 = 4.0 q50 = 5.8 q75 = 8.6 | qConf = 9.4 (conf=80.0%)

    @param a, b, x: lower and upper bounds with a <= x <= b
    @param alpha, beta: shape parameters with alpha, beta > 0
    """
    if x < a or x > b:
        raise ValueError("x outside support [a,b]")
    if alpha > 0 and beta > 0 :
        f = lambda t : beta_pdf(t, alpha, beta, a, b)
        return INTEGRAL(f, a, x, epsilon = epsilon)
    raise ValueError("precondition violated: alpha, beta > 0")

def beta_inv(conf, alpha, beta, a, b, epsilon = 0.001) :
    """
    BETA(3.0, 4.0, 20.0) alpha=0.8 beta=3.2
    mean1 = 6.5 var = 8.0 sigma1 = 2.8
    mean2 = 6.6 var = 10.1 sigma2 = 3.2 (N = 1000)
    q25 = 4.0 q50 = 5.8 q75 = 8.6 | qConf = 9.4 (conf=80.0%)

    @param a, b, x: lower and upper bounds with a <= x <= b
    @param alpha, beta: shape parameters with alpha, beta > 0
    """
    if conf < 0 or conf > 1:
        raise ValueError("conf outside support [0,1]")
    if alpha <= 0 or beta <= 0 :
        raise ValueError("precondition violated: alpha, beta > 0")
    summe = 0.0
    f = lambda t : beta_pdf(t, alpha, beta, a, b)
    for (A, tr) in iterIntegral(f, a, b, epsilon = epsilon) :
        summe += A
        if summe >= conf :
            return tr
    return summe

def alphaBetaFromAmB(a, m, b) :
    first_numer_alpha = 2.0 * (b + 4 * m - 5 * a)
    first_numer_beta = 2.0 * (5 * b - 4 * m - a)
    first_denom = 3.0 * (b - a)
    second_numer = (m - a) * (b - m)
    second_denom = (b - a) ** 2
    second = (1 + 4 * (second_numer / second_denom))
    alpha = (first_numer_alpha / first_denom) * second
    beta = (first_numer_beta / first_denom) * second
    return alpha, beta


#---
class BetaDistribution(object) :

    @classmethod
    def FromAmB(cls, a, m, b) :

        def _alphaBetaFromAmB(a, m, b) :
            first_numer_alpha = 2.0 * (b + 4 * m - 5 * a)
            first_numer_beta = 2.0 * (5 * b - 4 * m - a)
            first_denom = 3.0 * (b - a)
            second_numer = (m - a) * (b - m)
            second_denom = (b - a) ** 2
            second = (1 + 4 * (second_numer / second_denom))
            alpha = (first_numer_alpha / first_denom) * second
            beta = (first_numer_beta / first_denom) * second
            return alpha, beta

        alpha, beta = _alphaBetaFromAmB(a, m, b)
        dist = BetaDistribution(a, b, alpha, beta)
        dist.m = m
        return dist

    def __init__(self, a, b, alpha, beta) :
        self.a = a
        self.b = b
        self.alpha = alpha
        self.beta = beta
        self.m = None

    def random(self) :
        r = random.betavariate(self.alpha, self.beta)
        return self.a + r * (self.b - self.a)

    def mean(self) :
        if self.m is None :
            return None
        a = self.a
        m = self.m # Modus
        b = self.b
        return (a + 4 * m + b) / 6.0

    def sigma(self) :
        a = self.a
        b = self.b
        return (b - a) / 6.0


    def iterpdf(self, epsilon = 0.001):
        """@rtype: C{[(int, float)]}"""
        a = self.a
        b = self.b
        alpha = self.alpha
        beta = self.beta
        f = lambda t:beta_pdf(t, alpha, beta, a, b)
        for x in xrange(a, b + 1) :
            prop = INTEGRAL(f, max(a, x - 0.5), min(x + 0.5, b), epsilon = epsilon)
            yield (x, prop)

    def iterPDFasHistogram(self, maxDots = MAX_DOTS, epsilon = 0.001) :
        bucketProps = list(self.iterpdf(epsilon))
        maxProp = max([b[1] for b in bucketProps])
        header = "%5s %6s %s" % ("units", "CDF", "PDF")
        yield header
        cdf = 0.0
        for (x, prop) in bucketProps :
            dots = round(maxDots * prop / maxProp)
            cdf += prop
            cdfPercent = 100.0 * cdf
            theDots = "*"*dots
            line = "%(x)5i %(cdfPercent)6.2f %(theDots)s" % locals()
            yield line

#def iterPDFasHistogram(a, b, alpha, beta, maxDots = MAX_DOTS, epsilon = 0.001) :
#    f = lambda t : betadist.beta_pdf(t, alpha, beta, a, b)
#    bucketProps = [(x, betadist.INTEGRAL(f, max(a, x - 0.5), min(x + 0.5, b), epsilon = epsilon)) for x in xrange(a, b + 1)]
#    maxProp = max([b[1] for b in bucketProps])
#    header = "%5s %6s %s" % ("units", "CDF", "PDF")
#    yield header
#    cdf = 0.0
#    for (x, prop) in bucketProps :
#        dots = round(maxDots * prop / maxProp)
#        cdf += prop
#        cdfPercent = 100.0 * cdf
#        theDots = "*"*dots
#        line = "%(x)5i %(cdfPercent)6.2f %(theDots)s" % locals()
#        yield line

def round(real) :
    return int(real + 0.5)

#---
def test_uniform(a, b, **keywords) :
    """
    @keyword alpha, beta: shape parameter
    @keyword nominal: N aus O/N/P
    """
    try :
        alpha = keywords['alpha']
        beta = keywords['beta']
        nominal = None
    except KeyError :
        nominal = keywords['nominal']
        alpha, beta = alphaBetaFromAmB(a, nominal, b)
    print "a = %(a)r, b = %(b)r, alpha = %(alpha)r, beta = %(beta)r" % locals()
    if nominal is not None :
        print "nominal = %(nominal)r" % locals()
    step = (b - a) / 20.0
    x = a
    q80 = beta_inv(0.8, alpha, beta, a, b)
    print "q80 = %(q80)f" % locals()
    print "%8s %16s %16s %s" % ("x", "pdf", "cdf", "scale")
    area = BETA_GAMMA(alpha, beta, a, b)
    while x <= b :
        pdf = beta_pdf(x, alpha, beta, a, b)
        cdf = beta_cdf(x, alpha, beta, a, b)
        print "%(x)8.3f %(pdf)16f %(cdf)16f %(area)f" % locals()
        if x == b :
            break
        x = min(b, x + step)
        continue
    print

def main() :
    doctest.testmod()
    print INTEGRAL(lambda x : x, 0, 1) # ca. 0.5
    # [0, 1]
    test_uniform(0, 1, alpha = 1, beta = 1)
    test_uniform(0, 1, nominal = 0.5)
    test_uniform(0, 1, nominal = 0.25)
    test_uniform(0, 1, nominal = 0.75)
    # [0, b]
    test_uniform(0, 2, alpha = 1, beta = 1)
    test_uniform(0, 2, nominal = 0.5)
    test_uniform(0, 2, nominal = 0.25)
    test_uniform(0, 2, nominal = 0.75)
    # [a, b]
    test_uniform(3, 20, alpha = 0.8, beta = 3.2)
    test_uniform(9, 15, alpha = 2.5, beta = 1.6)
    return

if __name__ == "__main__" :
    main()
