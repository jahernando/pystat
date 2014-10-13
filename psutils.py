import math
import random
import bisect

from psparam import *

#----------------------------------------
#  Utilities for pystat
#----------------------------------------

#  Note: the functions that return a function start with a 'f',
#        they are called 'functors'

def checktype(x,typ):
    """ check that all the objects in the list (xs) ar of a given type (typ).
    Returns True or False.
    """
    if (isinstance(x,list)): 
        return (sum(map(lambda xi:checktype(xi,typ),x))==len(x))
    return isinstance(x,typ)

def finterval(x0,xf):
    """ Returns a boolean function f(x) that check is x is inside the interval [x0,xf]
    """
    def fun(x): return ((x>=x0) and (x<=xf))
    return fun

def range_partition(x0,xf,npoints):
    """ makes a partition of the interval [x0,xf] in npoints.
    Returns a list.
    """
    if (xf<x0 or (not isinstance(npoints,int)) or npoints<0):
        raise ArithmeticError('Not possible to make partition! ')
    x0,xf,npoints=float(x0),float(xf),int(npoints)
    step = (xf-x0)/(npoints+1)
    xs = map(lambda i: x0+i*step,range(npoints+1))
    xs.append(xf)
    return xs

def range_minfun(fun,x0,step,fmin=PMIN):
    """ returns the range (x0,xf) of the fun(x)>PMIN.
    It increments x on units of step. Note that step can be negative!
    """
    if (abs(step)<ZERO):
        raise ArithmeticError('No valid step in rage minfun!')
    xf = x0+step; p = fun(xf)
    while (p>fmin): 
        xf +=step; p = fun(xf)
    xs = [x0,xf]
    xs.sort()
    return xs
        
def ffuniparam(fun,xs,i):
    """ It takes a function f(x1,...,xn) of n parameters and returns a function f(xi) of the ith parameter, the other are taken as constants.
    xs are the values of (x1,xn) (except xi) taken as constant.
    """
    xxs = list(xs)
    def ifun(x):
        xxs[i] = x
        return fun(*xxs)
    return ifun
    
#---------------------------------------------------------------
# Funxy classes : functions created from a list of (x,y) points
#---------------------------------------------------------------

class Funxy:
    """ base class for functions creatd with (x,y) points """

    types = [int,float,str]

    def __init__(self,points):
        xs = map(lambda x: x[0],points)
        ys = map(lambda x: x[1],points)
        self._xtype, self._ytype = None, None
        for type in Funxy.types: # check for the type of x and y
            if checktype(xs,type): self._xtype = type
            if checktype(ys,type): self._ytype = type
        if (not self._xtype or not self._ytype):
            raise ArithmeticError('Not valid types for points in Funxy! ')
        def xsort(p1,p2): 
            if (p1[0] < p2[0]): return -1
            elif (p1[0] > p2[0]): return 1
            return 0
        points.sort(xsort) # sort the points using the x variable
        self._points = tuple(points) # me a copy
        self._npoints = len(self._points)
        self._xs = map(lambda x:x[0],self._points)
        x0,xf = self._xs[0],self._xs[-1]
        self._isin = finterval(x0,xf)
    
    def points(self):
        return self._points
    
    def domain(self):
        x0,xf = self._points[0][0],self._points[-1][0] # interval
        return x0,xf

    def isindomain(self,x):
        return self._isin(x)

    def npoints(self):
        return len(self._points)

    def xs(self):
        return self._xs

    def ys(self):
        ys = map(lambda x : x[1],self._points)
        return ys

    def point(self,k):
        return self._points[k]

    def __getitem__(self,k):
        return self._points[k]

    def __str__(self):
        return str(self._points)

    def __call__(self,x):
        if (not isinstance(x,self._xtype)):
            raise TypeError('Not valid input type for Funxy!')
        k = (bisect.bisect_right(self._xs,x)-1)
        if (not self.isindomain(x)):return  self.foutside(x,k)
        else: return self.finside(x,k)

    def foutside(self,x,k):
        """ return 0 if x is outside domain """
        return self._ytype(0)

    def finside(self,x,k):
        """ (virtual) return value for x given k-index. 
        xk is the greatest x-value of the points with xk <= x.
        """
        raise ArithmeticError('function not defined!')
        return None

class Funxy_xdiscrete(Funxy):
    """ Funxy function that returns the value of the y point with same x """
    def __init__(self,points):
        Funxy.__init__(self,points)

    def isindomain(self,x):
        """ check if x is in the domain """
        return (x in self._xs)

    def finside(self,x,k):
        """ return the value of the point of k-th index """
        return self._points[k][1]    

class Funxy_xinterpolate(Funxy):
    """ Funxy function that interpolates between two (x,y) points """

    def __init__(self,points):
        """ constructor with a list of (x,y) points """
        Funxy.__init__(self,points)

    def finside(self,x,k):
        """ return interpolation between points """
        if (k==self._npoints-1): k-=1
        x0,xf = self._points[k][0],self._points[k+1][0]
        y0,yf = self._points[k][1],self._points[k+1][1]
        dx,ddx = x-x0, xf-x0
        if (ddx<=0.): raise ArithmeticError('null dx, division by zero!')
        slope = (yf-y0)/ddx
        dy = y0+slope*dx
        return self._ytype(dy)
    
class Funxy_xright(Funxy):
    """ Funxy function that returns the value of the point (xi,yi) to the right of x  """
    def __init__(self,points):
        """ constructor with points (xi,yi) """
        Funxy.__init__(self,points)

    def finside(self,x,k):
        """ return the y value of the point with same of higher x """
        if (k==self._npoints-1): k-=1
        yf = self._points[k+1][1]
        return self._ytype(yf)

#------------------------------------------------------
# functions to do comulation, averange and variance
#------------------------------------------------------

def nsum(fun,ns):
    """ return the sum(f(k)) in the ns points """
    return sum(map(f,ns))

def riemannsum(fun,xs):
    """ compute the Riemann's sum of the function (fun) in a partition xs (a list of x-values)
    """
    ns = len(xs)-1
    dxs = map(lambda i: (xs[i+1]-xs[i]),range(ns))
    xcs = map(lambda i: 0.5*(xs[i+1]+xs[i]),range(ns))
    ys = map(fun,xcs)
    ars = map(lambda y,dx: y*abs(dx),ys,dxs)
    return sum(ars)

def integral(fun,x0,x,epsilon=EPSILON):
    """ returns the integral int_x0^x fun(x) with a precision epsilon.
    It uses Riemann's sums.
    """
    if (not checktype([x0,x],float)):
        raise TypeError('Not valid types for integral! ')
    if (abs(x-x0)<=ZERO): return 0.
    if (x<x0): return -1.*integral(fun,x,x0,epsilon)
    ar0,ar = -INFTY,INFTY
    niter, npoints = 0,1
    while ((abs(ar0-ar)>epsilon) and niter<NITER):
        xs = range_partition(x0,x,npoints)
        ar0,ar = ar,riemannsum(fun,xs)
        niter+=1; npoints *=10;
        if (DEBUG>2): print  'integral (iter) ',niter, ar
    if ((abs(ar0-ar)>epsilon)):
        print '>>> limit # iterations, integral precision (', abs(ar-ar0),')'
    if (DEBUG>1): print " integral [",x0,x,"]=",ar 
    return ar

def fintegral(fun,x0,epsilon=EPSILON):
    """ returns a function F(x) = int_x0^x f(x) dx.
    The integral is done numerically till the Riemann's sum covergence is smaller than an epsilon (default value 1e-3)
    """
    ff = lambda x : integral(fun,x0,x,epsilon)
    return ff

def sumaverage(pf,ns):
    """ Compute the average of a function, pf, in the ns list.
    it does the sum pf(n) with n in ns list.
    """
    return nsum(pf,ns)    

def sumvariance(pf,ns):
    """ Compute the variance of a given function, pf, in the ns list.
    """
    xave = sumaverage(pf,ns)
    dvar = lambda x : (x-xave)*(x-xave)*pf(x)
    var = nsum(dvar,ns)
    if (DEBUG>2): print ' average ',xave,' variance ',var
    return var

def intaverage(pf,x0,xf):
    """ Compute the average of a function, pf, in the [x0,xf] range.
    it does the int pf(n) with [x0,xf] interval
    """
    return nsum(pf,ns)    

def intvariance(pf,x0,xf):
    """ Compute the variance of a function, pf, in the [x0,xf] range.
    it does the int (pf(n)-mu)^2 with [x0,xf] interval
    """
    xave = intaverage(pf,x0,xf)
    dvar = lambda x : (x-xave)*(x-xave)*pf(x)
    var = integral(dvar,x0,xf)
    if (DEBUG>2): print ' average ',xave,' variance ',var
    return var

#--------------------------------------------------
# functions for comulative and inverse comulative
#-------------------------------------------------

def funxy_pmf(pf,ns):
    """ returns a (Funxy object) probability mass function (pmf) in ns (a list of rvs)
    """
    ys = map(pf,ns)
    norma = sum(map(pf,ns))
    ys = map(lambda x: x/norma,ys)
    zs = zip(ns,ys)
    if (DEBUG): print " (pmf) zs ",zs
    xfun = Funxy_xdiscrete(zs)
    return xfun

def funxy_cmf(pmf,ns):
    """ returns a (Funxy object) comulative function (cmf) of a pmf in the partition ns (a list of rvs)
    """
    ys = map(pmf,ns)
    zs = map(lambda i: min(1.,sum(ys[:i+1])),range(len(ys)))
    z = zip(ns,zs)
    if (DEBUG): print " (cmf) z ",z
    xfun = Funxy_xdiscrete(z)
    return xfun

def funxy_invcmf(cmf,ns):
    """ returns a (Funxy object) inverse of the comulative mass function (invcmf) of a cmf in ns partition (ints).
    """
    ys = map(cmf,ns)
    zs = zip(ys,ns)
    if (zs[0][0]>0.): zs.insert(0,(0.,ns[0]))
    if (zs[-1][0]<1.): zs.append((1.,ns[-1]))
    if (DEBUG): print " (icmf) z ",zs
    xfun = Funxy_xright(zs)
    return xfun

def funxy_pdf(pf,xs,epsilon=EPSILON):
    """ returns a (Funxy object) probability mass function (pdf) using a function and a partition xs.
    It uses integral function with precision epsilon.
    """
    if (not checktype(xs,float)): 
        raise TypeError('Not valid types for funxy_pdf! ')
    x0,xf = xs[0],xs[-1]
    norma = integral(pf,x0,xf,epsilon)
    ys = map(lambda x: pf(x)/norma,xs)
    zs = zip(xs,ys)
    if (DEBUG): print " (pdf) zs ",zs
    xfun = Funxy_xinterpolate(zs)
    return xfun

def funxy_cdf(pdf,xs,epsilon=EPSILON):
    """ returns a (Funxy object) comulative density function (cdf) of a pdf, using the points given in the partition xs (a list).
    It uses integral function with precision epsilon.
    """
    if (not checktype(xs,float)): 
        raise TypeError('Not valid types for funxy_cdf! ')
    dar = lambda i: integral(pdf,xs[i],xs[i+1],epsilon)
    ys = map(dar,range(len(xs)-1))
    ys.insert(0,0.)
    # print ys
    zs = map(lambda i: min(1.,sum(ys[:i+1])),range(len(ys)))
    z = zip(xs,zs)
    if (DEBUG): print " (cmf) z ",z
    xfun = Funxy_xinterpolate(z)
    return xfun
    
def funxy_invcdf(cdf,xs):
    """ returns a (Funxy object) inverse comulative density function (invcdf) of a cdf, using the points given in the partition xs (a list).
    """
    if (not checktype(xs,float)): 
        raise TypeError('Not valid types for funxy_invcdf! ')
    ys = map(cdf,xs)
    if (ys[-1]<1.): ys[-1]=1.
    zs = zip(ys,xs)
    if (DEBUG): print " (icdf) z ",zs
    xfun = Funxy_xinterpolate(zs)
    return xfun

