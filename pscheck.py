#
# examples using mypystat
#
from math import *
import random
import pslist, psutils, pspdf, psbayes, psestimate, pslib, psplot

#---------------------------------
# extra drawing utilities
#---------------------------------

def thistograph(xs,fun,name='',nbins=100):
    h = thisto(xs,name,nbins)
    x0,xf = min(xs),max(xs)
    xs = stutils.range_partition(x0,xf,nbins+1)
    nn = len(ys); step = (xf-x0)/(nbins+1)
    zs = map(lambda x: nn*fun(x)*step,xs)
    tg = tgraph(zip(xs,zs))
    return h,tg

def tpull(tval,xpars,xname='var'):
    """ with the true value, tval, the list of (value,sigma) it returns a list of histograms with the values, sigmas, and pulls.
    """
    xs = map(lambda x: x[0],xpars)
    ss = map(lambda x: x[1],xpars)
    pu = map(lambda x,sx: (x-tval)/sx,xs,ss)
    hvals = [xs,ss,pu]
    names = ['','rms','pull']
    hs = []
    for i in range(len(names)):
        name,hval = names[i],hvals[i]
        hs.append(thisto(hval,xname+name))
    return hs  

#----------------------------------
# check pystatlist
#----------------------------------

def check_statlist(n=200,nn=2000):
    """ creates a poisson pmf and gaus pdf and compute the average, rms and pull of a list of n rvs in nn samples.
    """
    mu,sigma = 0.,1.
    p1 = pspdf.pdf_gauss(mu,sigma)
    hs2 = check_statlist_pf(p1,mu,sigma,n,nn)
    c = psplot.tncanvas(hs2,name='statlist')
    return c,hs2

def check_statlist_pf(pf,tmu,tsig,n=100,nn=1000):
    """ takes a pdf, atrue mean and a true sigma.
    generates nn samples with n events, and estimates mean, sigma, man rms and sigma rms.
    return the histograms of mean, sigma, mean rms, sigma rms, and the pull of the mean and the sigma.
    """
    funs = [pslist.mean,pslist.meanrms,pslist.srms,pslist.rmsrms]
    def sample(i):
        data = pf.random(n)
        return map(lambda f:f(data),funs)
    vals = map(sample,range(nn))
    mus = map(lambda x: (x[0],x[1]),vals)
    sis = map(lambda x: (x[2],x[3]),vals)
    hs1 = psplot.tpull(tmu,mus,'mu')
    hs2 = psplot.tpull(tsig,sis,'sigma')
    return hs1+hs2

#-------------------------------------------
# check pystatpdf
#-------------------------------------------

def check_pdfs():
    """ check the pmfs and pdfs defined in pystat
    """
    pfs = [pspdf.pmf_coin(),
           pspdf.pmf_ncoins(100),
           pspdf.pmf_binomial(10,0.2),
           pspdf.pmf_dice(), 
           pspdf.pmf_ndices(2),
           pspdf.pmf_poisson(2.),
           pspdf.pdf_gauss(0.,1.),
           pspdf.pdf_uniform(0.,1.),
           pspdf.pdf_exponential(1.), 
           pspdf.pdf_chi2(2)
           ]
    p1 = pspdf.pdf_gauss(0.,1.) 
    p2 = pspdf.pdf_gauss(0.,1.) 
    pfs.append(p1+p2) # add 2-gaus convolution 
    p3 = pspdf.pdf_gauss(0.,1.)
    p4 = pspdf.pdf_exponential(2.)
    pfs.append(p3+p4) # add gaus+expo convolution
    hs = map(lambda p: check_pf(p),pfs)
    return hs

def check_pf(pf,n=1000):
    """ check a pmf or pdf
    """
    print ">>> NEW CHECK"
    print '\t name ',pf.name()
    print '\t parameters ',pf.parameters()
    print '\t domain ',pf.domain()
    xs = pf.rvariables()
    zs = pf.points()
    tg1 = psplot.tgraph(zs,title='probability')
    zs = zip(xs,pf.comulative(xs))
    tg2 = psplot.tgraph(zs,title='cumulative')
    ys = psutils.range_partition(0.,1.,101)
    zs = zip(ys,pf.inverse_comulative(ys))
    tg3 = psplot.tgraph(zs,title='inverse cumulative')
    print '\t average ',pf.average()
    print '\t variance ',pf.variance()
    print '\t rms ',pf.rms()
    print '\t median ',pf.median()
    xs = pf.random(n)
    h1 = psplot.thisto(xs,'random')
    hs = [tg1,tg2,tg3,h1]
    c = psplot.tcanvas(hs,name=pf.name())
    raw_input('return to next check')
    return c,hs

#----------------------------
# check bayes
#----------------------------

def check_bayes_dices(nn=5):
    """ One player takes a dice from a bowl, there is one dice of 4 sides, one of 6, one of 8, one of 12 and one of 20. 
    What is the posterior probability if he gets a 6.
    And what is if he gets the following data: [6,6,,8,7,7,5,4]
    """
    hypos = [4,6,8,12,20] # the dices of Dragons and Dungeons!
    def dlike(d,h):
        # Likelihood of the data given an hypothesis, P(D|H), for the dices. 
        if (d>h): return 0.
        return 1./(1.*h)
    pnf = pspdf.pnf_uniform(hypos) # prior probability (uniform for all dices)
    pb = psbayes.PBayes(pnf,hypos,dlike) # create PBayes object!
    dice = random.sample(hypos,1)[0] # select a dice
    print ' selected dice ',dice
    pmf = pspdf.pmf_dice(dice) # generate nn events with the dices
    data = pmf.random(nn) 
    print ' data ',data
    pb.evidence(data) # update the bayes probability
    zs = zip(hypos,pb.posterior(hypos)) # points with the posterior
    print 'posterior ',zs
    tg = psplot.tgraph(zs,title='dice'+str(dice))
    c = psplot.tcanvas(tg)
    return c,tg

def check_bayes_expo():
    """ Consider N=10,100,1000 measurements of a exponential rv. with tau=1.
    And a prior probability for tau flat in [0.1,2.], use the data to compute the posterior probability.
    """
    cc = []
    NS = [10,100,1000,10000] # simple sizes
    RS = [(0.1,3.),(0.4,1.6),(0.8,1.2),(0.94,1.06)] # ranges
    for i in range(len(NS)):
        N = NS[i]
        u0,uf = RS[i]
        fu = pspdf.pdf_uniform(u0,uf) # pdf prior P(H)
        hypos = fu.rvariables() # all possible hypothesis! 
        fe = pspdf.pdf_exponential(1.) # 
        dlike = pslib.exponential # likelihood P(D|H)
        # bayes object for continuous hypothesis:
        pb = psbayes.PDBayes(fu,hypos,dlike) 
        data = fe.random(N) # generate data
        pb.evidence(data) # update bayes
        zs = zip(hypos,pb.posterior(hypos)) # points of the posterior
        xpf = pspdf.pdf_points(zs) # generate pdf with the posterior points
        print ' N points ',N
        print ' average ',xpf.average(),' rms ',xpf.rms()
        c1 = xpf.inverse_comulative(0.5-0.341)
        c2 = xpf.inverse_comulative(0.5+0.341) 
        print ' sigma-contaiment ',c2-c1
        tg = psplot.tgraph(zs,title='posterior_N'+str(N))
        hg = psplot.thisto(data,'d')
        c = psplot.tncanvas([tg,hg])
        cc.append((c,(tg,hg)))
    return cc

#------------------------------------------
# check pystatestimate
#------------------------------------------

def check_lle_gauss_scan(mu=0,sigma=1.,n=1000):
    pg = pslib.fgauss
    llm = psestimate.LLNEstimator(pg)
    g1 = pspdf.pdf_gauss(mu,sigma)
    data = g1.random(n)
    xm = psestimate.Estimate(mu,mu-sigma,mu+sigma)
    xs = psestimate.Estimate(0.8*sigma,0.2*sigma,2.5*sigma)
    ll = llm.estimate(data,[xm,xs])
    h = psplot.thisto(data)
    flike = llm.flikelihood()
    x0,xf = min(data),max(data)
    xs = psutils.range_partition(x0,xf,100)
    dx = (xf-x0)/102.
    ys = map(lambda x: n*dx*flike(x),xs)
    zs = zip(xs,ys)
    tg = psplot.tgraph(zs,'fit')
    c = psplot.tcanvas([tg,h])
    return c,(tg,h),flike

def check_lle_expo(tau=1.,n=200,nn=1000):
    """ check LL Estimation with an exponential, return histograms with mean, sigma and pull of tau
    """
    pe = pspdf.pdf_exponential(tau)
    vals = []
    for i in range(nn):
        data = pe.random(n)
        llm = psestimate.lle_expo(data)
        vals.append(llm.xestimate())
    hs = psplot.tpull(tau,vals,'tau')
    c = psplot.tncanvas(hs,name='lle_expo')
    return c,hs

def check_lle_gauss(mu=0.,sigma=1.,n=200,nn=1000):
    """ check LL Estimation with gauss, return histograms with mean, sigma and pull of mu, sigma
    """
    pg = pspdf.pdf_gauss(mu,sigma)
    vals = []
    for i in range(nn):
        data = pg.random(n)
        llm = psestimate.lle_gauss(data)
        vals.append(llm.xestimate())
    theta = [mu,sigma]
    names = ['mu','sigma']
    hs = []
    for i in range(2):
        xvals = map(lambda x:x[i],vals)
        hs1 = psplot.tpull(theta[i],xvals,names[i])
        hs+=hs1
    c = psplot.tncanvas(hs,name='lle_gauss')
    return c,hs

def check_lle_poisson(nu=1.,n=200,nn=1000):
    """ check LL Estimation with a poisson, return histograms with mean, sigma and pull of nu.
    """
    pg = pspdf.pmf_poisson(nu)
    vals = []
    for i in range(nn):
        data = pg.random(n)
        llm = psestimate.lle_poisson(data)
        vals.append(llm.xestimate())
    hs = psplot.tpull(nu,vals,'nu')
    c = psplot.tncanvas(hs,name='lle_poisson')
    return hs,c

def check_lle_custom(n=400):
    """ check LL Estimation with a gaus+expo custom function.
    It can be quite slow... be aware!
    return histogram with one sample and results of the LL fit.
    """
    tau,mu,sigma=2.,0.,1.
    pe = pspdf.pdf_exponential(tau)
    pg = pspdf.pdf_gauss(mu,sigma)
    print ' this check can be slow (be-patient!)!!'
    print ' building exp+gauss'
    pp = pg+pe
    def fun(t,sig): 
        p1 = pspdf.pdf_exponential(t)
        p2 = pspdf.pdf_gauss(mu,sig)
        p3 = p1+p2
        return p3
    print ' random '
    data = pp.random(n)
    xhats = [psestimate.Estimate(1.2*tau,0.2*tau,2.5*tau),
             psestimate.Estimate(1.2*sigma,0.2*sigma,2.5*sigma)]
    llm = psestimate.LLNEstimator(fun)
    print ' estimate! '
    llm.estimate(data,xhats)
    vals = llm.xestimate()
    flike = llm.flikelihood()
    h = psplot.thisto(data)
    x0,xf = min(data),max(data)
    xs = psutils.range_partition(x0,xf,100)
    dx = (xf-x0)/102.
    ys = map(lambda x: n*dx*flike(x),xs)
    zs = zip(xs,ys)
    tg = psplot.tgraph(zs,'fit')
    c = psplot.tcanvas([tg,h])
    return (c,tg,h)

