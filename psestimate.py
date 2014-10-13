"""

Module to Estimate parameters using the Log Likelihood

@author Jose A Hernando
@date 7/11/2013
"""

from psparam import *
import math, bisect
import psutils, pslib

#--------------------------------------------------
# Base classe for estimate and some guess functions
#--------------------------------------------------

class Estimate():
    """ class Estimate to estore the estimate theta of a parameter
    center value, x, and range of validity, [x0,xf]
    """
    def __init__(self,x,x0,xf):
        """ Set the estimate x, and the range x0,xf
        """
        if ((x0>x) or (x0>xf) or (x>xf)):
            print 'x0,x,xf',x0,x,xf
            raise ArithmeticError('Invalid estimate! ')
        self._x = x
        self._x0 = x0
        self._xf = xf

    def __str__(self):
        """ converts to a string """
        s = '['+str(self._x)+' ('+str(self._x0)+','+str(self._xf)+')] '
        return s

    def __repr__(self):
        """ converts to a string """
        return str(self)
    
    def value(self):
        """ return the value """
        return self._x

    def range(self):
        """ return the range """
        return (self._x0,self._xf)

    def width(self):
        """ return the width """
        return (self._xf-self._x0)

    def sigma(self):
        """ return the sigma = width/2"""
        return self.width()/2.

    def xestimate(self):
        """ return (value,sigma) """
        return (self.value(),self.sigma())

def guessparam_mean(data,n=1.):
    """ guess the mean of data """
    mu = st.average(data)
    smu = st.meanrms(data)
    muhat = Estimate(mu,mu-n*smu,mu+n*smu)
    return muhat

def guessparam_rms(data,n=1.):
    """ guess the rms of data """
    mu = st.srms(data)
    smu = st.rmsrms(data)
    muhat = Estimate(mu,max(ZERO,mu-n*smu),mu+n*smu)
    return muhat

def guessparam_gauss(data,n=1):
    """ guess mu,sigma of data if gaussian """
    muhat = guessparam_mean(data,n)
    sghat = guessparam_rms(data,n)
    return [muhat,sghat]

def guessparam_expo(data,n=1.):
    """ guess tau of data if exponental """
    return guessparam_mean(data,n)

def guessparam_poisson(data,n=1.):
    """ guess nu of data if poisson """
    return guessparam_mean(data,n)

def guessparam_chi2(data,n=1.):
    """ guess nodf of data if chi2 """
    return guessparam_mean(data,n)
    
#---------------------------------------
# Base classes for Log Likehood method
#---------------------------------------

#-------------------------
class LLEstimatorBase():
    """ Base Class for parmater estimator of a pdf(x;theta) """
    def __init__(self,fun):
        """ contructor with a functor, fun(parameter), that returns f(x)=f(x;parameter). xhat is an Estimate object of the parameter with the initial guess.
        """
        self._fun = fun # should the like function
        self._flike = None 
        self._vals = []
        self._llmax = -INFTY
        self._xhat = None
        
    def xhatestimate(self):
        """ returns the estimate as Estimate object"""
        return self._xhat

    def xestimate(self):
        """ returns the estimate as a tupe (value,sigma) """
        if (isinstance(self._xhat,list)):
            return map(lambda x: x.xestimate(),self._xhat)
        return self._xhat.xestimate()

    def value(self):
        """ returns the value of the estimate"""
        if (isinstance(self._xhat,list)):
            return map(lambda x: x.value(),self._xhat)
        return self._xhat.value()
    
    def flikelihood(self):
        """ return the likelihood function f(x;theta) with theta is the estimate (theta-hat) """
        return self._flike

    def likelihood(self,data):
        """ returns the likelihood of the data (a value or list) """
        if (isinstance(data,list)):
            vals = map(self.likelihood,data)
            val = reduce(lambda x,y: x*y, vals)
            return val
        return self._flike(data)
    
    def llikelihood(self,data):
        """ returns the log likelihood of the data (a value or list) """
        if (isinstance(data,list)):
            vals = map(self.llikelihood,data)
            xval = sum(vals)
            return xval
        return math.log(max(ZERO,self.likelihood(data)))

    def llmax(self):
        """ return the log likelihood of the estimate"""
        return self._llmax

    def points(self):
        """ return a list of points (x,llikelihood(x) """
        return self._vals

    def set_initial_estimate(self,xhat):
        """ set the initial parameter guess (xhat is an Estimate object) """
        self._xhat = xhat
        self._set_value(self.value())

    def estimate(self,data,xhat,epsilon=ELLEPSILON):
        """ estimate the parameter using data and initial guess xhat """
        self.set_initial_estimate(xhat)
        print ' LL guess ',self.xestimate()
        tol,iter = INFTY,0
        while (abs(tol)>epsilon and iter<ELLNITER):
            self._scan(data)
            tol = self._update_estimate()
            # print 'LLE estimate iter ',iter,' tol ',tol
            iter+=1
        self._set_estimate()
        if (iter == ELLNITER):
            print '>>> limit # iterations reached, tol (',tol,') '
        print ' LL estimate ',self.xestimate(),' llmax ',self.llmax()
        #print ' LL iter ',iter
        return self._xhat
        
    def _set_value(self,x):
        """ (private) set the value of the parameter in the likelihood function """
        self._flike = self._fun(x)
        return

    def _scan(self,data):
        """ (private) does an scan of parameter values (xi,llikelihood(xi)) """
        return -INFTY

    def _update_estimate(self):
        """ (private) update the parameter estimate, returns tolerance """
        return INFTY
    
    def _set_estimate(self):
        """ private """
        x = self.value()
        self._set_value(x)
        return

#-------------------------------------
class LLEstimatorBisect(LLEstimatorBase):
    """ Class for Estimator of a pdf(x;parameter), one parameter using bisect of the range """
            
    def __init__(self,fun):
        LLEstimatorBase.__init__(self,fun)
    
    def _scan(self,data):
        """ scan the log likelihood in the range of the Estimate with npoints """
        xc = self._xhat.value()
        x0,xf = self._xhat.range()
        xs1 = psutils.range_partition(x0,xc,ELLNPOINTS)
        xs2 = psutils.range_partition(xc,xf,ELLNPOINTS)
        xs3 = xs1[:-1]+xs2
        xs0 = map(lambda x: x[0],self._vals)
        xs = []
        for xi in xs3:
            if (not xi in xs0): xs.append(xi)
        def llx(x):
            self._set_value(x)
            return self.llikelihood(data)
        # print 'len (xs) ',len(xs),npoints
        lls =map(llx,xs)
        zs = zip(xs,lls)
        self._vals += zs
        xsort = lambda v1,v2: v1[0]<v2[0]
        self._vals.sort()
        llmax = max(map(lambda x:x[1],self._vals))
        self._llmax0 = self._llmax
        self._llmax = llmax
        #print ' scan vals ',self._vals
        #print ' scan zs ',zs
        #print ' scan llmax ',self._llmax
        return llmax

    def _update_estimate(self):
        """ (private)  """     
        ys = map(lambda x: x[1],self._vals)
        xs = map(lambda x: x[0],self._vals)
        ymax = max(ys); imax = ys.index(ymax)
        xmax = xs[imax]; vmax = self._vals[imax]
        lvals = self._vals[:imax] ;rvals = self._vals[imax+1:]
        lval,rval =  lvals[-1],rvals[0]
        xn = lval[0] 
        if (rval[1]>lval[1]): xn=rval[0]
        xc = 0.5*(xmax+xn)
        lvals.append(vmax)
        rvals.insert(0,vmax)
        lvals = map(lambda x: (x[0],ymax-x[1]),lvals)
        rvals = map(lambda x: (x[0],ymax-x[1]),rvals)
        if (abs(lvals[0][1])<0.5 or abs(rvals[-1][1])<0.5):
            print ' scan range ',lvals[0],rvals[-1]
            raise ArithmeticError('Uncoverage!')
        def vsort(v1,v2): 
            if (abs(v1[1]-0.5) < abs(v2[1]-0.5)): return -1
            elif (abs(v1[1]-0.5) < abs(v2[1]-0.5)): return +1
            return 0
        lvals.sort(vsort)
        rvals.sort(vsort)
        xl = 0.5*(lvals[0][0]+lvals[1][0])
        xr = 0.5*(rvals[0][0]+rvals[1][0])
        if (xc <= xl or xc>=xr): xc =xmax
        xhat0 = self._xhat
        self._xhat = Estimate(xc,xl,xr)
        def utol(x0,x1):
            u1 = abs(x0._x0-x1._x0)
            u2 = abs(x0._x-x1._x)
            u3 = abs(x0._xf-x1._xf)
            return max(u1,u2,u3)
        tol = utol(xhat0,self._xhat) 
        #print ' update val max ',vmax
        #print ' update lval, rval',lval,rval
        #print ' update lvals ',lvals
        #print ' update rvals ',rvals
        #print ' update lvals (2) ',lvals[-2:]
        #print ' update rvals (2) ',rvals[:2]
        #print ' update estimate ',self._xhat
        #print ' update tolerance ',tol
        # raw_input('enter key ')
        return tol

#-------------------------------------
class LLEstimatorVertex(LLEstimatorBase):
    """ LLEstimator of one parameter using the computatiom of the vortex of the llike parabola"""

    def __init__(self,fun):
        """ constructor """
        LLEstimatorBase.__init__(self,fun)
        
    def _scan(self,data):
        """ (private) scan the log likelihood in the range of the Estimate with npoints """
        xc = self._xhat.value()
        x0,xf = self._xhat.range()
        xs1 = psutils.range_partition(x0,xc,3*ELLNPOINTS)
        xs2 = psutils.range_partition(xc,xf,3*ELLNPOINTS)
        xs3 = xs1[:-1]+xs2
        xs0 = map(lambda x:x[0],self._vals)
        xs = []
        for xi in xs3:
            if (not xi in xs0): xs.append(xi)
        def llx(x):
            self._set_value(x)
            return self.llikelihood(data)
        lls =map(llx,xs)
        zs = zip(xs,lls)
        self._vals += zs
        def xsort(v1,v2): 
            if (v1[0]<v2[0]): return -1
            elif (v1[0]>v2[0]): return +1
            return 0
        self._vals.sort(xsort)
        llmax = max(map(lambda x:x[1],self._vals))
        self._llmax = llmax
        # print 'scan vals ',self._vals
        # print 'scan llmax ',self._llmax
        return llmax

    def _update_estimate(self):
        """ (private)  """     
        ys = map(lambda x: x[1],self._vals)
        xs = map(lambda x: x[0],self._vals)
        x0,xf = xs[0],xs[-1]
        ymax = max(ys); imax = ys.index(ymax); xmax = xs[imax]
        lvals = self._vals[:imax+1] ;rvals = self._vals[imax:]
        y0,yf = self._vals[0][1], self._vals[-1][1]
        #print self._vals
        if (abs(y0-ymax)<ELLDELTAMIN or (abs(yf-ymax)<ELLDELTAMIN)):
            print ' y0,yf,ymax ',y0,yf,ymax
            raise ArithmeticError('LE Uncoverage!')
        xlt,ylt = lvals[-2]
        xrt,yrt = rvals[1]
        k = math.sqrt((ymax-yrt)/(ymax-ylt))
        xc = (xrt*xrt-(1-k)*xmax*xmax-k*xlt*xlt)/(2*(xrt-(1-k)*xmax-k*xlt))
        xlvals = map(lambda x: (ymax-x[1],x[0]),lvals)
        xrvals = map(lambda x: (ymax-x[1],x[0]),rvals)
        xlvals.reverse()
        flt = psutils.Funxy_xinterpolate(xlvals)
        #xxlt = flt(ELLDELTAMAX)
        x0lt = flt(ELLDELTAMIN)
        if (x0lt>xc): x0lt = x0
        #xxlt = x0lt+(xc-x0lt)
        frt = psutils.Funxy_xinterpolate(xrvals)
        #xxrt = frt(ELLDELTAMAX)
        x0rt = frt(ELLDELTAMIN)
        if (x0rt<xc): x0rt = xf
        if (xc<x0lt or xc>x0rt): xc = xmax
        #xxrt = x0rt+(x0rt-xc)
        # self._xhat = Estimate(xc,xxlt,xxrt)
        self._xhat = Estimate(xc,x0lt,x0rt)
        sig = self._xhat.width()/2.
        tol = abs(xc-xmax)/sig
        #print ' update xlvals ',xlvals
        #print ' update xrvals ',xrvals
        #print ' update xlt,xmax,xrt',xlt,xmax,xrt
        #print ' update k xc',k,xc
        #print ' update (2sigma) ',xxlt,xxrt
        #print ' update (2sigma) ',self._xhat,self._xhat.xestimate()
        #print ' update (1sigma) ',x0lt,x0rt
        #print ' LLV update estimate ',self._xhat.xestimate()
        #print ' LLV update tolerance ',tol
        #raw_input('enter key ')
        return  tol

#LLEstimator = LLEstimatorBisect
LLEstimator = LLEstimatorVertex

#---------------------------------------
# Estimator of N parameters
#---------------------------------------

class LLNEstimator(LLEstimatorBase):

    def __init__(self,fun):
        LLEstimatorBase.__init__(self,fun)
        
    def set_initial_estimate(self,xhats):
        """ Constructor using a functor, fun, that f(parameters) returns the function f(x) = f(x;parameters). xhats is the list of Estimate parameters as initial guess.
        """
        if (not isinstance(xhats,list)):
            raise TypeError("LLNEstimator requires a list of estimates ")
        self._npars = len(xhats)
        self._xhat = list(xhats)
        self._set_value(self.value())
        self._xestimators = []
        for i in range(len(xhats)):
            ifun,xhat = self._ifuns[i],self._xhat[i]
            llestimator = LLEstimator(ifun)
            llestimator.set_initial_estimate(xhat)
            self._xestimators.append(llestimator)
        #print ' set initial estimate ',xhats
        return

    def _set_value(self,xs):
        """ (private) """
        self._flike = self._fun(*xs)
        self._ifuns= map(lambda i:psutils.ffuniparam(self._fun,xs,i),
                         range(self._npars))
        # print ' set value ',xs
        return

    def _scan(self,data):
        """ (private) scan the log likelihood in the range of the Estimate with npoints 
        """
        self.set_initial_estimate(self._xhat)
        xs = self.value()
        llmax = -INFTY
        for i in range(self._npars):
            xhat = self._xhat[i]
            self._xestimators[i]._fun = self._ifuns[i] # update function
            xhat = self._xestimators[i].estimate(data,xhat)
            ll = self._xestimators[i]._llmax            
            xs[i] = xhat.value() # update values
            self._vals.append((xs,ll))
            if (ll>llmax): llmax = ll
            self._set_value(xs) # update values
            # print ' LLN scan par',i,' llmax ',ll,' val ',xhat.xestimate()
        self._llmax = ll
        # print ' scan llmax ',llmax
        return self._llmax

    def _update_estimate(self):
        """ (private) """
        mtol = 0.
        xhats = []
        for i in range(self._npars):
            x0 = self._xhat[i].value()
            xhat = self._xestimators[i]._xhat
            x,s = xhat.xestimate()
            tol = abs(x-x0)/s
            #print ' LLN update estimate i ',x0,x,s,tol
            if (tol>mtol): mtol = tol
            xhats.append(xhat)
        xhat0 = self._xhat
        self._xhatbest = list(xhats)
        xhats = map(lambda x,y: Estimate(x.value(),*y.range()),xhats,self._xhat)
        self._xhat = xhats
        #print ' LLN update estimate - init ',xhat0
        #print ' LLN update estimate - next ',self._xhat
        #print ' LLN update estimate - best ',self._xhatbest
        #print ' LLN update estimate - next range ',self.xestimate()
        #print ' LLN update tol ',mtol
        #raw_input('enter key')
        return mtol
    
    def _set_estimate(self):
        #print ' vals ',self._vals
        self._xhat = self._xhatbest
        LLEstimatorBase._set_estimate(self)
        return

        

#--------------------------------------------------
# LLEstimator for the common pdfs
#--------------------------------------------------

def lle_gauss(data,n=ELLNSIG):
    """ LLEstimator with gauss and data """
    fun = pslib.fgauss
    spar = guessparam_gauss(data,n)
    llm = LLNEstimator(fun)
    llm.estimate(data,spar)
    return llm

def lle_expo(data,n=ELLNSIG):
    """ LLEstimator with exponental and data """
    fun = pslib.fexponential
    spar = guessparam_expo(data,n)
    llm = LLEstimator(fun)
    llm.estimate(data,spar)
    return llm

def lle_poisson(data,n=ELLNSIG):
    """ LLEstimator with poisson and data """
    fun = pslib.fpoisson
    spar = guessparam_expo(data,n)
    llm = LLEstimator(fun)
    llm.estimate(data,spar)
    return llm

def lle_chi2(data,n=ELLNSIG):
    """ LLEstimator with chi2 and data """
    fun = pslib.fchi2
    spar = guessparam_chi2(data,n)
    llm = LLEstimator(fun)
    llm.estimate(data,spar)
    return llm

#---------------------
# check
#---------------------

def check(tau=1.,n=100):
    tau = 1.
    pe = st.pdf_exponential(tau)
    xs = Estimate(1.5*tau,tau/100.,3.*tau)
    fun = pslib.fexponential
    data = pe.random(n)
    llm = LLEstimator(fun)
    hat = llm.estimate(data,xs)
    return llm,data

def check2(mu=0.,sigma=1.,n=100):
    pg = st.pdf_gauss(mu,sigma)
    data = pg.random(n)
    xm = guessparam_mean(data,ELLNSIG)
    xs = guessparam_rms(data,ELLNSIG)
    xm = Estimate(1.2*mu,mu-2.*sigma,mu+2.*sigma)
    xs = Estimate(1.2*sigma,0.2*sigma,mu+2.5*sigma)
    fun = pslib.fgauss
    llm = LLNEstimator(fun)
    hat = llm.estimate(data,[xm,xs])
    return llm,data

#    pdf = llm.makepdf()
#    pval = llm.pvalue(data)
#    return llm
    #tg = tgraph(vals)
    #return tg        

