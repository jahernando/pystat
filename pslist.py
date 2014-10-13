import math

"""
statistic with list of values

@version 1.0
@author Jose A. Hernando
@date   25/10/13
"""

#----------------------------------------
# Statistic with list of numbers
#-----------------------------------------

def average(a):
    """ returns the average of list a
    """
    v = sum(a)/float(len(a))
    return v

def mean(a):
    """ returns the average of list a
    """
    return average(a)

def variance(a):
    """ returns the variance of list a
    """
    av = average(a)
    var = sum(map(lambda x: (x-av)*(x-av),a))/float(len(a))
    return var

def svariance(a):
    """ returns the estimation of the variance of list a.
    It has the factor 1/(n-1)
    """
    av = average(a)
    svar = sum(map(lambda x: (x-av)*(x-av),a))/(float(len(a))-1)
    return svar

def srms(a):
    """ returns the estimator of the sqrt(variance) of a list a.
    The estimator of the variance has the 1/(n-1) factor
    """
    svar = svariance(a)
    return math.sqrt(svar)

def rms(a):
    """ returns the standar desviation of list a
    """
    var = variance(a)
    xrms = math.sqrt(var)
    return xrms

def sigma(a):
    """ returns the standar desviation of list a
    """
    return rms(a)

def meanrms(x):
    """ returns the rms of the mean of a list
    """
    mu = mean(x)
    rms = srms(x)
    n = len(x)
    smu = rms/math.sqrt(1.*n) 
    return smu

def nmoment(x,n):
    nn = float(len(x))
    mu = mean(x)
    mo = sum(map(lambda x: math.pow((x-mu),n)/nn,x))
    return mo

def varvariance(x):
    """ returns the variance of the variance
    """
    n = float(len(x))
    mu2 = nmoment(x,2)
    mu4 = nmoment(x,4)
    vvar =(mu4/n-(mu2/n)*(mu2*(n-3)/(n-1)))
    return vvar

def rmsrms(x):
    """ returns the rms of the rms
    """
    rms = srms(x)
    varrms = variancerms(x)
    rrms = varrms/(2*rms)
    return rrms

def variancerms(x):
    """ returns the rms of the variance
    """
    vvar = varvariance(x)
    if (vvar<0): return 0.
    return math.sqrt(vvar)

def median(a):
    """ returns the median of list a
    """
    b = a[:]
    b.sort()
    n = len(a)
    ni = int(n/2.)
    return b[ni]

def contaiment(a,cont=0.64):
    """ returns the interval of contaiment of list a
    """
    b = a[:]
    b.sort()
    n = len(a)
    n0,n1 = int((0.5-cont/2.)*n),int((0.5+cont/2.)*n)
    return b[n0],b[n1]

def cov(a,b):
    """ returns the covariance between two list a,b
    """
    am = average(a)
    bm = average(b)
    c = map(lambda x,y: (x-am)*(y-bm),a,b)
    c = sum(c)/float(len(a))
    return c

def rho(a,b):
    """ returns the correlation coeffiencent between two list a,b
    """
    cc = cov(a,b)
    sa = rms(a)
    sb = rms(b)
    return cc/(sa*sb)

