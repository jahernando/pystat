import math
import random

from psparam import * #import parameters
import psutils, pslib

"""
A module to deal with statistic in Python

@version 1.0
@author Jose A. Hernando
@date   25/10/13
"""

#----------------------------------------
#   Probability Function class (base for pmf and pdf classes)
#----------------------------------------        

class PF():
    """ Probability Function Base class
    """

    types = [int,float,str]

    def __init__(self,pf,rvs,cf=None,icf=None,name ='unknown'):
        """ create a probability function object. 
        pf is the probabilty function.
        rvs is the complete range of the randon variable.
        (for continous variables, it must be a small enough partition)
        cf is the comulative function.
        icf is the inverse comulative function.
        name is the name of the function.
        """
        self._xtype = None
        for type in PF.types:
            if (psutils.checktype(rvs,type)): self._xtype = type
        if (not self._xtype):
            raise TypeError('Not valid types for PF!')
        self._isnormal = psutils.finterval(0.,1.) # normal interval [0.,1.]
        self._pf,self._cf,self._icf = pf,cf,icf
        self._name = str(name)
        self._parameters = {}
        self._rvs = list(rvs)
        self._rvs.sort()

    def name(self):
        """ return the name of the pf (i.e Gauss, Poisson, etc) """
        return self._name

    def domain(self):
        """ returns the domain of PF """
        return self._rvs[0],self._rvs[-1]
    
    def isindomain(self,x):
        """ returns true if x is inside the domain of PF """
        print ">>> PF isindomain not implemented!"
        return False

    def checktype(self,x):
        """ returns true if x is of type of the random variable """
        if (not isinstance(x,self._xtype)):
            raise TypeError(' Not valid type for rv!')
        return 

    def rvariables(self):
        """ return a list of ordered values inside the domain """ 
        return self._rvs

    def points(self):
        zs = zip(self._rvs,self.likelihood(self._rvs))
        return zs

    def __str__(self):
        return str(self.points())

    def parameters(self):
        """ returns a dictionary with the parameters"""
        return dict(self._parameters)

    def parameter(self,name):
        """ returns the parameter value associated to that name """
        if (self._parameters.has_key(name)): return self._parameters[name]
        return None
    
    def set_parameter(self,name,value):
        """ set a parameter """
        self._parameters[name] = value
        setattr(self,name,value)
        return
    
    def likelihood(self,x):
        """ returns the likelihood of x, pf(x) """
        if (isinstance(x,list)): return map(self.likelihood,x)
        self.checktype(x)
        if (not self.isindomain(x)): return 0.
        return self._pf(x)
    
    def llikelihood(self,x):
        if (isinstance(x,list)): return map(self.llikelihood,x)
        lk =  self.likelihood(x)
        llk = math.log(max(ZERO,lk))
        return llk

    def __call__(self,x):
        """ returns likelihood(x)  """
        return self.likelihood(x)
    
    def comulative(self,x):
        """ returns the value of cumulative function at x, cf(x).
        x can be a value or a list of values.
        """
        if (isinstance(x,list)): return map(self.comulative,x)
        self.checktype(x)
        if (not self.isindomain(x)): return 0.
        if (not self._cf): self.make_comulative()
        return self._cf(x) 

    def inverse_comulative(self,y):
        """ returns the value of inverse cumulative function, icf(x).
        y can be a value or a list of values.
        y must be inside the interval [0.,1.].
        """
        if (isinstance(y,list)): return map(self.inverse_comulative,y)
        if (not isinstance(y,float)):
            raise TypeError("Not valid type for argument in inv-comulative")
        if (not self._isnormal(y)): 
            raise ArithmeticError('argument of inv-cdf out of [0.,1.]!')
        if (not self._icf): self.make_inverse_comulative()
        return self._icf(y)

    def random(self,n=1):
        """ returns n randoms numbers. 
        If n=1 returns a value, if n>1 returns a list of values.
        """
        n = int(n)
        ys = map(lambda i: random.uniform(0.,1.),range(n))
        xs = map(self.inverse_comulative,ys)
        if (n==1): return xs[0]
        return xs
    
    def median(self):
        """ returns the median """
        return self.inverse_comulative(0.5)
    
    def pvalue(self,xs):
        """ returns the p-value of a xs random variables (a list)"""
        raise TypeError('p-value virtual method!')
        return None
        

#----------------
class PNF(PF):
    """ probability class (PF) for discrete random variables non additive """
    
    def __init__(self,pf,rvs,cf=None,icf=None,name='PNF'):
        """ constructor of a PF object with discrete random variables non additive (i.e string).
        pf is the function p(x).
        rvs are the list of all values of the random variables (or hypothesis).
        """
        PF.__init__(self,pf,rvs,cf=cf,icf=icf,name=name)
        return

    def isindomain(self,x):
        """ check x in in the rvs """
        return (x in self._rvs)
    
    def make_comulative(self):
        """ make comulative function """
        self._cf = psutils.funxy_cmf(self._pf,self._rvs)
        return

    def make_inverse_comulative(self):
        """ make inv comulative function """
        if (not self._cf): self.comulative(self._rvs[0])
        self._icf = psutils.funxy_invcmf(self._cf,self._rvs)
        return

    def pvalue(self,xs):
        if (not isinstance(xs,list)): return self.inverse_comulative(xs)
        n = len(xs)
        if (n<10):
            raise ArithmeticError('Not p value for few values of rvs')
        zs = zip(self.comulative(self._rvs),self.likelihood(self._rvs))
        pf = pmf_points(zs)
        ave = pf.average()
        rms = pf.rms()
        pu = pdf_gauss(ave,rms/math.sqrt(n))
        acds = sum(self.comulative(xs))/(1.*n)
        return pu.comulative(acds)

#---------------
class PMF(PF):
    """ probability class (PF) for discrete additive random variables (PMF) """
    
    def __init__(self,pf,rvs,cf=None,icf=None,name='PMF'):
        """ constructor of a PMF object with discrete additive random variables.
        pf is the function p(x) (pmf(x)).
        rvs are the list of all values of the random variables (or hypothesis).
        """
        PF.__init__(self,pf,rvs,cf=cf,icf=icf,name=name)
        self._ave, self._var = None, None
        return

    def isindomain(self,x):
        """ check x in in the domain """
        return (x in self._rvs)
    
    def make_comulative(self):
        """ make comulative function (CMF) """
        self._cf = psutils.funxy_cmf(self._pf,self._rvs)
        return

    def make_inverse_comulative(self):
        """ make inverse comulative function (invCMF) """
        if (not self._cf): self.comulative(self._rvs[0])
        self._icf = psutils.funxy_invcmf(self._cf,self._rvs)
        return

    def average(self):
        """ returns average """
        if (self._ave): return self._ave
        dave = lambda x: x*self._pf(x)
        ave = sum(map(dave,self._rvs))
        if (DEBUG>2): print ' average ',ave
        self._ave = ave
        return self._ave

    def variance(self):
        """  returns variance """
        if (self._var): return self._var
        xave = self.average()
        dvar = lambda x : (x-xave)*(x-xave)*self._pf(x)
        var = sum(map(dvar,self._rvs))
        if (DEBUG>2): print ' variance ',var
        self._var = var
        return self._var

    def rms(self):
        """ returns standar deviation  """ 
        var = self.variance()
        return math.sqrt(var)

    def convolute(self,pg):
        """ convolute this PMF with pg (a PMF type)
        """
        if (not isinstance(pg,PMF)): 
            raise TypeError('Can not convolute PMFs!')
        dfz = lambda z,x: pg(z-x)*self._pf(x)
        def conv(z):
            df = lambda x: dfz(z,x)
            val = sum(map(df,self._rvs))
            return val                      
        n0,n1 = self.domain()
        m0,m1 = pg.domain()
        rvs = range(n0+m0,n1+m1+1)
        ps = map(conv,rvs)
        xfun = psutils.Funxy_xdiscrete(zip(rvs,ps))
        name = self.name()+'*'+pg.name()
        xpmf = PMF(xfun,rvs,name=name)
        return xpmf

    def __add__(self,pg):
        """ convolute with other PMF """
        return self.convolute(pg)
    
    def pvalue(self,xs):
        if (not isinstance(xs,list)): return self.inverse_comulative(xs)
        n = len(xs)
        if (n<10):
            raise ArithmeticError('Not p value for few values of rvs')
        zs = zip(self.comulative(self._rvs),self.likelihood(self._rvs))
        pf = pmf_points(zs)
        ave = pf.average()
        rms = pf.rms()
        pu = pdf_gauss(ave,rms/math.sqrt(n))
        acds = sum(self.comulative(xs))/(1.*n)
        return pu.comulative(acds)

#---------------
class PDF(PF):
    """ probability class (PF) for continuous additive random variables (PDF)"""
    
    def __init__(self,pf,rvs,cf=None,icf=None,name='PDF'):
        """ constructor of a PDF object with continuous additive random variables.
        pf is the function pdf(x).
        rvs is a small enough partition of the x random variable interval.
        """
        if (not psutils.checktype(rvs,float)):
            raise TypeError('Not valid rvs for PDF!')
        PF.__init__(self,pf,rvs,cf=cf,icf=icf,name=name)
        self._ave,self._var = None,None
        x0,xf = self.domain()
        self._isindomain = psutils.finterval(x0,xf)
        return
    
    def isindomain(self,x):
        """ check is x is in the domain """
        return self._isindomain(x)
    
    def make_comulative(self):
        """ creates a comulative function if non provided """
        self._cf = psutils.funxy_cdf(self._pf,self._rvs)
        return

    def make_inverse_comulative(self):
        """ creates an inv comulative function if non provided """
        if (not self._cf): self.comulative(self._rvs[0])
        self._icf = psutils.funxy_invcdf(self._cf,self._rvs)
        return
    
    def average(self):
        """ returns average """
        if (self._ave): return self._ave
        dave = lambda x: x*self._pf(x)
        x0,xf = self.domain()
        xave = psutils.integral(dave,x0,xf)
        if (DEBUG>2): print ' average ',xave
        self._ave = xave
        return self._ave    

    def variance(self):
        """  returns variance """
        if (self._var): return self._var
        x0,xf = self.domain()
        xave = self.average()
        dvar = lambda x : (x-xave)*(x-xave)*self._pf(x)
        var = psutils.integral(dvar,x0,xf)
        if (DEBUG>2): print ' variance ',var
        self._var = var
        return self._var
    
    def rms(self):
        """ returns standar deviation  """ 
        var = self.variance()
        xrms = math.sqrt(var)
        if (DEBUG>2): print ' rms ',xrms
        return xrms

    def convolute(self,pg):
        """ convolute this PDF with pg (a PDF)
        """
        if (not isinstance(pg,PDF)): 
            raise TypeError('Can not convolute PDFs!')
        x0,xf = self.domain()
        y0,yf = pg.domain()
        dfz = lambda z,x: pg(z-x)*self._pf(x)
        def conv(z):
            df = lambda x: dfz(z,x)
            val = psutils.integral(df,x0,xf)
            return val             
        rvs = psutils.range_partition(x0+y0,xf+yf,NPOINTS)
        ps = map(conv,rvs)
        xfun = psutils.Funxy_xinterpolate(zip(rvs,ps))
        name = self.name()+'*'+pg.name()
        xpdf = PDF(xfun,rvs,name=name)
        return xpdf

    def __add__(self,pg):
        """ convolute with other PDF """
        return self.convolute(pg)

    def pvalue(self,x):
        """ retruns the average comulative of a list of rvs.
        x can also be a value!
        """
        if (not isinstance(x,list)): return self.pvalue([x])
        n = len(x)
        acds = sum(self.comulative(x))/(1.*n)
        fun = pdf_un(n)
        return fun.comulative(acds)

#-----------------------------------------
# User clear interface to 
# create specific pnf, pmf, pdf objects
#-----------------------------------------

def pnf_list(rvs,ps,name='PN'):
    """ returns a PNF constructed with a list of random variables (rvs) and their probabilities (ps). 
    rvs is a list of int, float or str. ps a list of float"""
    zs = zip(rvs,ps)
    xfun = psutils.Funxy_xdiscrete(zs)
    xpf = PNF(xfun,rvs,name=name)
    return xpf

def pnf_points(points,name='PNPoints'):
    """ returns a PNF constructed with a list of points (xi,p(xi)).
    where xi can be str,int or float. p(xi) is float.    
    """
    rvs = map(lambda x:x[0],points)    
    norma = sum(map(lambda x:x[1],points))
    ps = map(lambda x: x[1]/norma,points)
    return pnf_list(rvs,ps,name='PNFPoints')

def pnf_uniform(rvs):
    """ returns a PNF constructed with a list of random variables (rvs) with uniform probabilities.
    """
    rvs.sort()
    nn = len(rvs)
    ps = nn*[1./(1.*nn)]
    return pnf_list(rvs,ps,name='PNUniform')

def pmf_list(rvs,ps,name='PMF'):
    """ returns a PMF constructed with a list random variables (rvs) and their probabilities (ps). 
    rvs is a list of int or float. ps a list of float"""
    zs = zip(rvs,ps)
    pf = psutils.Funxy_xdiscrete(zs)
    xpmf = PMF(pf,rvs,name=name)
    return xpmf

def pmf_points(points):
    """ creates a PDF using a list of points (xi,yi) 
    """
    rvs = map(lambda x:x[0],points)    
    norma = sum(map(lambda x:x[1],points))
    ps = map(lambda x: x[1]/norma,points)
    return pmf_list(rvs,ps,name='PMFPoints')

def pmf_coin():
    """ Creates a PMF object with a coin (0,1) results
    """
    rvs,ps = [0,1],[0.5,0.5]
    return pmf_list(rvs,ps)

def pmf_binomial(N,p):
    """ Creates a PMF object with a binomial distribution with p-probability and N trials
    """
    if (not isinstance(N,int)):
        raise TypeError('N of binomial must be int!')
    if (p<0. or p>=1.):
        raise ArithmeticError('probability must be in [0.,1.]!')
    xf = pslib.fbinomial(N,p)
    rvs = range(0,N+1)
    xpf = PMF(xf,rvs,name='binomial')
    xpf.set_parameter('p',p)
    xpf.set_parameter('N',N)
    return xpf

def pmf_ncoins(N=2):
    """ Creates a PMF object of N coins
    """
    if (not isinstance(N,int)):
        raise TypeError('N coins have to be an int!')
    xpf = pmf_binomial(N,0.5)
    return xpf
    
def pmf_dice(nsides=6):
    """ Creates a PMF object for a dice of n sides
    """
    if (not isinstance(nsides,int)):
        raise TypeError('# sides of dice have to be an int! ')    
    rvs,ps = range(1,nsides+1),nsides*[1./float(nsides)]
    xpf = pmf_list(rvs,ps,name='Dice')
    xpf.set_parameter('nsides',6)
    return xpf

def pmf_ndices(n=2,nsides=6):
    """ Creates a PMF objext of n dices with nsides
    """
    if (not isinstance(n,int)):
        print ">>> Not possible to create a not int number of dices!"
        return None
    c = pmf_dice(nsides)
    for i in range(1,n):
        ci = pmf_dice(nsides); c = c+ci; 
    return c

def pmf_ndifferentdices(ns):
    """ ns is a list of the sides of the dices (i.e [6,4,2] 3 dices, one with 6 sides, one with 4 and one with 3).
    Creates a PMF object with len(ns) dices of different sides
    """
    c = pmf_dice(ns[0])
    for nsides in ns[1:]:
        ci = pmf_dice(nsides); c = c+ci; 
    return c

def pmf_uniform(n0,nf):
    """ returns a PMF object with a uniform distribution between n0,nf (inclusives). n0,nf are int.
    """
    if ((not psutils.checktype([n0,nf],int)) or (nf<=n0)):
        raise TypeError('Invalid range for a uniform pmf!')
    ns = range(n0,nf+1)
    nn = len(ns)
    ps = nn*[1./(1.*nn)]
    return pmf_list(ns,ps,name='PMFUniform')

def pmf_poisson(nu,pmin=PMIN):
    """ Creates a PMF object with Poisson pdf with mean mu and sigma.
    The range of validity of the pdf is for pdf(n)>pmin
    """
    xpois = pslib.fpoisson(nu)
    n0,step = int(nu),1
    n0,nf = psutils.range_minfun(xpois,n0,step,pmin)
    n0,nc = psutils.range_minfun(xpois,n0,-1*step,pmin)
    if (n0<0): n0=0
    rvs = range(n0,nf+1)
    xpmf = PMF(xpois,rvs,name='Poisson')
    xpmf.set_parameter('nu',nu)
    return xpmf    

def pdf_partition(pf,x0,xf,npoints=NPOINTS,name='PDF'):
    """ creates a PDF object with a pdf function and an intercel [x0,xf].
    It creates a partition of npoints as random variables (rvs) of the PDF.
    """
    xs = psutils.range_partition(x0,xf,npoints)
    xpdf = PDF(pf,xs,name=name)
    return xpdf

def pdf_gauss(mu,sigma,pmin=PMIN,npoints=NPOINTS):
    """ creates a PDF object with a Gauss pdf with mean mu and sigma.
    The range of validity of the pdf is for pdf(x)>pmin
    """
    xgaus = pslib.fgauss(mu,sigma)
    # xgcdf = pslib.fcdfgauss(mu,sigma)
    x0,step = mu,0.1*sigma,
    x0,xf = psutils.range_minfun(xgaus,x0,step)
    x0,xc = psutils.range_minfun(xgaus,x0,-1.*step)
    xpdf = pdf_partition(xgaus,x0,xf,npoints=npoints,name='Gauss')
    xpdf.set_parameter('mu',mu)
    xpdf.set_parameter('sigma',sigma)
    return xpdf

def pdf_normal(pmin=PMIN,npoints=NPOINTS):
    """ creates a PDF object with a normal distribution, gaussian wit mu=0, sigma=1
    """
    return pdf_gauss(0.,1.,pmin=pmin,npoints=npoints)

def pdf_uniform(x0,xf,npoints=NPOINTS):
    """ creates a PDF object with a uniform pdf in the interval [x0,xf]
    """
    xuni = pslib.funiform(x0,xf)
    xpdf = pdf_partition(xuni,x0,xf,npoints,'Uniform')
    xpdf.set_parameter('x0',x0)
    xpdf.set_parameter('xf',xf)
    return xpdf

def pdf_triangle(x0,xc,xf,npoints=NPOINTS):
    xfun = pslib.ftriangle(x0,xc,xf)
    xpdf= pdf_partition(xfun,x0,xf,npoints,'Triangle')    
    xpdf.set_parameter('x0',x0)
    xpdf.set_parameter('xc',xc)
    xpdf.set_parameter('xf',xf)
    return xpdf

def pdf_un(n,npoints=NPOINTS):
    if (not isinstance(n,int) or n<=0):
        raise ArithmeticError('Not valid n for un-function!')
    if (n==1): return pdf_uniform(0.,1.,npoints)
    elif (n==2): return pdf_triangle(0.,0.5,1.,npoints)
    elif (n==3):
        u1,u2 = pdf_un(1,npoints),pdf_un(2,npoints)
        uu = u1+u2
        zs = uu.points()
        zs = map(lambda x: (x[0]/2.,2.*x[1]),zs)
        return pdf_points(zs)
    elif (n==4):
        u1,u2 = pdf_un(2,npoints),pdf_un(2,npoints)
        uu = u1+u2
        zs = uu.points()
        zs = map(lambda x: (x[0]/2.,2.*x[1]),zs)
        return pdf_points(zs)
    else:
        return pdf_gauss(0.5,1./math.sqrt(12.*n),npoints)
    return
    

def pdf_exponential(tau,pmin=PMIN,npoints=NPOINTS):
    """ creates a PDF object with an Exponential pdf with mean lifetime tau.
    The range of validity of the pdf is for pdf(x)>pmin
    """
    xexp = pslib.fexponential(tau)
    x0,step=0.,0.1*tau
    x0,xf = psutils.range_minfun(xexp,x0,step)
    xpdf = pdf_partition(xexp,x0,xf,npoints,'Exponential')
    xpdf.set_parameter('tau',tau)
    return xpdf

def pdf_chi2(n,pmin=PMIN,npoints=NPOINTS):
    """ creates a PDF object with a chiw with n dof
    """ 
    xchi = pslib.fchi2(n)
    x0,step = 0.,n*0.1
    x0,xf = psutils.range_minfun(xchi,x0,step)
    xpdf = pdf_partition(xchi,x0,xf,npoints,'Chi2')
    xpdf.set_parameter('ndof',n)    
    return xpdf    

def pdf_points(points):
    """ creates a PDF using a list of points (xi,yi) 
    """
    xfun = psutils.Funxy_xinterpolate(points)
    x0,xf = xfun.domain()
    norma = psutils.integral(xfun,x0,xf)
    xpoints = map(lambda x: (x[0],x[1]/norma),points)
    xfun = psutils.Funxy_xinterpolate(xpoints)
    xs = map(lambda x:x[0],xpoints)
    xpdf = PDF(xfun,xs,name='PDFpoints')
    return xpdf

def pfnconvolute(pfcreator,n):
    """ creates a pdf (pmf) convoluting a pdf (pmf) n times
    """
    if (n<=0):
        raise ArithmeticError('Can not nconvolute n times!')
    print 'be patient, it can be slow!'
    p1 = pfcreator()
    for i in range(1,n):
        p2 = pfcreator()
        p1 = p1+p2
    return p1

