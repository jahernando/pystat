import math,random
from psparam import *

#----------------------------------------
#  Common pdf functions
#----------------------------------------

def fbinomial(N,p):
    """ returns a binomial pmf with parameters: N-number of trials, 
    and p-probability of success
    """
    if (not isinstance(N,int)): 
        print '>>> binomial parameter N not an int! '
        return None
    if (p<=0 or p>1):
        print '>>> binomial parameter p (probability) not valid! ',p
        return None
    def fun(n):
        if (not isinstance(n,int)): 
            print ">>> binomial vaiable n not an int! "
            return None
        if (n<0 or n>N): return 0.
        norma = math.factorial(N)/(math.factorial(n)*math.factorial(N-n))
        return norma*math.pow(p,n)*math.pow(1-p,N-n)
    return fun

def fpoisson(nu):
    """ returns poisson pmf with nu (average) parameter
    """
    if (nu<0): 
        print ">>> Negative Poisson parameter nu! "
        return None
    def fun(n):
        if (not isinstance(n,int)): 
            print ">>> Poisson variable not an int! "
            return None
        if (n<0): return 0.
        return math.pow(nu,n)*math.exp(-nu)/math.factorial(n)
    return fun
        
def fgauss(mu,sigma):
    """ returns gaussian pdf with average (mu) and sigma
    """
    if (sigma<=0): 
        print ">>> Gaussian parameter sigma negative! "
        return None
    def fun(x):
        x = float(x)
        norma = 1./(math.sqrt(2*math.pi)*sigma)
        val = norma*math.exp(-(x-mu)*(x-mu)/(2*sigma*sigma))
        return val
    return fun

def fcdfgauss(mu,sigma):
    """ returns gaussian cdf with average (mu) and sigma
    """
    if (sigma<=0): 
        print ">>> Gaussian parameter sigma negative! "
        return None
    def fun(x):
        x = float(x)
        val = 0.5*(1+math.erf((x-mu)/(sigma*math.sqrt(2))))
        return val
    return fun

def fcdfuniform(x0,xf):
    """ returns uniform cdf with average (mu) and sigma
    """               
    def fun(x):
        if (x<x0 or x>xf): return 0.
        x = float(x)
        val = float(x)/(float(xf)-float(x0))
        return val
    return fun

def fexponential(tau):
    """ returns exponetial pdf with mean lifetime (tau)
    """
    if (tau<=0): 
        print ">>> Exponential tau parameter negative or null! "
        return None
    def fun(t):
        if (t<0): return 0.
        t = float(t)
        return (1./tau)*math.exp(-1*t/tau)
    return fun

def fcdfexponential(tau):
    """ returns exponetial cdf with mean lifetime (tau)
    """
    if (tau<=0): 
        print ">>> Exponential tau parameter negative or null! "
        return None
    def fun(t):
        if (t<0): return 0.
        t = float(t)
        return 1-math.exp(-1*t/tau)
    return fun

def fchi2(n):
    """ returns chi2 pdf with n degrees of freedom
    """
    if (not isinstance(n,int)):
        print ">>> n-degrees of freedom in chi2 must be int"
        return None
    def fun(x):
        if (x<0): 
            print ">>> chi2 argument should be possitive"
            return None
        val = (1./(math.pow(2.,n/2.)*math.gamma(n/2.)))*math.pow(x,(n/2.-1.))*math.exp(-x/2.)
        return val
    return fun

def binomial(n,N,p):
    """ returns binomial probability for n-success of N-trials and p-probability
    """
    return fbinomial(N,p)(n)

def poisson(n,nu):
    """ returns probability of n-events in a Poisson with average nu
    """
    return fpoisson(nu)(n)

def gauss(x,mu,sigma):
    """ returns gaussian probability density for average mu and sigma
    """
    return fgauss(mu,sigma)(x)

def exponential(t,tau):
    """ returns exponential probability density for mean lifetime tau
    """
    return fexponential(tau)(t)

def uniform(x,x0,xf):
    """ return the probability density for a uniform pdf in the intercal x0,xf
    """
    if (x<x0 or x>xf): return 0.
    return 1./(xf-x0)

def funiform(x0,xf):
    fun = lambda x: uniform(x,x0,xf)
    return fun

def triangle(x,x0,xc,xf):
    if (x0>xc or xc>xf):
        raise ArithmeticError('Impossible to create triangle function!')
    if (x<x0 or x>xf): return 0.
    if (x<=xc): return 2.*(x-x0)/((xf-x0)*(xc-x0))
    else: return (2./(xf-x0))*(1.-(x-xc)/(xf-xc))

def ftriangle(x0,xc,xf):
    fun = lambda x: triangle(x,x0,xc,xf)
    return fun

def chi2(x,n):
    """ return the probability density for a uniform chi2 of n degrees of freedom
    """
    return fchi2(n)(x)
