#
# examples using pystat
#
from math import *
import random
import pslist, pspdf, psbayes, pslib, psutils, psplot

#-----------------------
#  basic concepts
#-----------------------

def exncoins():
    """ flop N coins, get the number of head in n trials, make m experiments
    compute average and variance. Return histograms of mu, sigma.
    Notice that the rms is bias!
    """
    N = 100
    pf = stpdf.pmf_ncoins(N)
    nn = [3,6,10,40,100]
    nm, ns, nss = [],[],[]
    m = 1000
    hs = []
    for n in nn:
        name = 'n'+str(n)
        h1 = TH1F('mu_'+name,'mu_'+name,21,40,60)
        h2 = TH1F('rms_'+name,'rms_'+name,13,0,12)
        h3 = TH1F('rmss_'+name,'rmss_'+name,13,0,12)
        yym,yys,yyss = [],[],[]
        for k in range(m):
            rs = pf.random(n)
            ym = pslist.average(rs); yym.append(ym)
            ys = pslist.variance(rs); yys.append(ys)
            yss = pslist.svariance(rs); yyss.append(yss)
            h1.Fill(ym,1.)
            h2.Fill(ys,1.)
            h2.Fill(yss,1.)
        hs.append((h1,h2))
        nm.append(pslist.average(yym))
        ns.append(pslist.average(yys))
        nss.append(pslist.average(yyss))
    z = zip(nn,nm)
    tg1 = psplot.tgraph(z,title='mean')
    z = zip(nn,ns)
    tg2 = psplot.tgraph(z,title='variance')
    z = zip(nn,nss)
    tg3 = psplot.tgraph(z,title='svariance')
    return ((tg1,tg2,tg3),hs)

def excov():
    """ generate x,y gaussin in L, L/2. Compute cov and rho.
    Rotate by theta, compute cov and rho
    """
    L = 2.
    xpf = pspdf.pdf_gauss(0.,L)
    ypf = pspdf.pdf_gauss(0.,L/2)
    N = 10000
    xs, ys = xpf.random(N), ypf.random(N)
    h1 = TH2F('xy1','xy1',40,-3.*L,3.*L,40,-3.*L,3.*L)
    cov,rho = pslist.cov(xs,ys),pslist.rho(xs,ys)
    print ' cov ', cov, 'rho ',rho
    map(lambda x,y: h1.Fill(x,y),xs,ys)
    theta = pi/4.
    xps = map(lambda x,y: +x*cos(theta)+y*sin(theta),xs,ys)
    yps = map(lambda x,y: -x*sin(theta)+y*cos(theta),xs,ys)
    h2 = TH2F('xy2','xy2',40,-3.*L,3.*L,40,-3.*L,3.*L)
    cov,rho = pslist.cov(xps,yps),pslist.rho(xps,yps)
    print ' cov ', cov, 'rho ',rho
    map(lambda x,y: h2.Fill(x,y),xps,yps)
    return (h1,h2)

def exring():
    """ generate x,y in a ring of radious ra,rb
    Compute cov and rho
    """
    ra, rb=2.,3.
    rf = pspdf.pdf_uniform(ra,rb)
    pf = pspdf.pdf_uniform(-pi,pi)
    N = 10000
    rs = rf.random(N); phis = pf.random(N) 
    xs = map(lambda r,phi: r*cos(phi),rs,phis)
    ys = map(lambda r,phi: r*sin(phi),rs,phis)
    h1 = TH2F('r','r',40,-1.2*rb,1.2*rb,40,-1.2*rb,1.2*rb)
    map(lambda x,y: h1.Fill(x,y),xs,ys)
    cov = pslist.cov(xs,ys); rho = pslist.rho(xs,ys)
    print ' cov ',cov,' rho ',rho
    return h1
    
#-----------------------
#  probability
#-----------------------

def extwodices():
    """ Roll two dices, check Bayes' Theorem.
    A = Get a 6 in one dice, B = get 8 points
    Returns P(A), P(B), PA(A|B), P(B|A)
    """
    N = 10000
    d1,d2= pspdf.pmf_dice(),pspdf.pmf_dice()
    z = zip(d1.random(N),d2.random(N))
    k,n = 6,8
    aevt = lambda x: (k in x)
    bevt = lambda x: (sum(x)==n)
    ya = filter(aevt,z)
    yb = filter(bevt,z)
    yab = filter(aevt,yb)
    yba = filter(bevt,ya)
    pa,pb = len(ya)/(1.*N),len(yb)/(1.*N)
    pab,pba = len(yab)/(1.*len(yb)),len(yba)/(1.*len(ya))
    print ' P(A)=',pa,', P(B|A)=',pba
    print ' p(B)=',pb,', P(A|B)=',pab
    print ' P(A) P(B|A)=',pa*pba,' P(B) P(A|B)=',pb*pab
    return pa,pb,pab,pba

def exmontyhall(N=10):
    """ Monty Hall question.
    Behind nd=3 doors there is only one price. Select one door.
    Monty hall opens one of the door you did not select, and there is nothing behind it. Do you want to select other door?
    It ask you the question N=10 times, and compute the number of success.
    Consider two strategies: always keep your selection or always swith.
    What is the best strategy?
    """
    nd = 3
    n,p = range(nd),nd*[0]
    p[0] = 1 
    ys = []
    for i in range(N):
        y = nd*[0]
        y[random.sample(n,1)[0]]=1
        kks = filter(lambda h: y[h]==0,range(nd))
        #print '>>> There is a price behind a number! ',y
        print '>>> There is a price behind a number! '
        ii = raw_input('\t select a number (0=return) '+str(n)+' ')
        if (ii==''): ii = 0
        else: ii = int(ii)
        ops = filter(lambda h: h!=ii,range(nd))
        kops = filter(lambda h: y[h]==0,ops)
        kk = random.sample(kops,1)[0]
        # print '\t opts ',ops,' kk ',kops
        print '\t there is nothing at ',kk
        ops = filter(lambda h: h != kk, ops)
        # ok = raw_input('\t do you want to swith? (yes=return) '+str(ops)+' ')
        # if (ok == ''): ii = ops[0]
        ok = raw_input('\t do you want to swith? (yes=return) '+str(ops)+' ')
        if (ok!=''):
            if (len(ops)==1): ii = ops[0]
            else: ii = int(ok)
        print '\t you select ',ii, ' and got ',y[ii]
        if (y[ii]==1): print '\t congratulations!! '
        else: print '\t we are sooo sorry! '
        print '\t these were the prizes ',y
        ys.append(y[ii])
    print ' You got right ',sum(ys),'times, probability ',sum(ys)/(1.*N)
    return
    
def exbayesdices():
    """ One player takes a dice from a bowl, there is one dice of 4 sides, one of 6, one of 8, one of 12 and one of 20. 
    What is the posterior probability if he gets a 6.
    And what is if he gets the following data: [6,6,,8,7,7,5,4]
    """
    hypos = [4,6,8,12,20]
    def dlike(d,h):
        if (d>h): return 0.
        return 1./(1.*h)
    ys = [1./len(hypos)]*len(hypos)
    zs = zip(hypos,ys)
    xfun = psutils.Funxy_xdiscrete(zs)
    pb = psbayes.PBayes(xfun,hypos,dlike)
    datas = [6,6,8,7,7,5,4]
    pb.evidence(datas)
    print 'posterior ',zip(hypos,pb.posterior(hypos))
    return

def exbayesexpo():
    """ Consider N=10,100,1000 measurements of a exponential rv. with tau=1.
    And a prior probability for tau flat in [0.1,2.], use the data to compute the posterior probability.
    """
    tgs = []
    NS = [10,100,1000,10000]
    RS = [(0.1,3.),(0.4,1.6),(0.8,1.2),(0.94,1.06)]
    for i in range(len(NS)):
        N = NS[i]
        u0,uf = RS[i]
        fu = pspdf.pdf_uniform(u0,uf)
        fe = pspdf.pdf_exponential(1.)
        hypos = fu.rvariables()
        dlike = pslib.exponential
        pb = psbayes.PBayes(fu,hypos,dlike)
        data = fe.random(N)
        x0,xf = fe.domain()
        pb.evidence(data)
        xs = pb.hypothesis()
        zs = zip(xs,pb.posterior(xs))
        xpf = pspdf.pdf_points(zs)
        print ' N points ',N
        print ' average ',xpf.average(),' rms ',xpf.rms()
        print ' contaiment ',xpf.inverse_comulative(0.5-0.341),xpf.inverse_comulative(0.5+0.341) 
        tg = psplot.tgraph(zs,title='posterior_N'+str(N))
        h = psplot.thisto(data)
        tgs.append((h,tg))
    return tgs

#---------------------------------------------------
# Basic pdfs
#---------------------------------------------------

def ex2poisson(nu1=1.,nu2=1.,n=1000):
    p1 = pspdf.pmf_poisson(nu1)
    p2 = pspdf.pmf_poisson(nu2)
    x1 = p1.random(n)
    x2 = p2.random(n)
    xs = map(lambda x,y:x+y,x1,x2)
    h = psplot.thisto(xs)
    p = pspdf.pmf_poisson(nu1+nu2)
    ys = map(lambda x: n*p(x),p.rvariables())
    zs = zip(p.rvariables(),ys)
    tg = psplot.tgraph(zs)
    return (h,tg)

def ex2gauss(mu1=1.,sig1=1.,mu2=1.,sig2=1.,n=1000):
    g1 = pspdf.pdf_gauss(mu1,sig1)
    g2 = pspdf.pdf_gauss(mu2,sig2)
    sig = sqrt(sig1*sig1+sig2*sig2)
    x1 = g1.random(n)
    x2 = g1.random(n)
    xs = map(lambda x,y:x+y,x1,x2)
    g3 = pspdf.pdf_gauss(mu1+mu2,sig)
    hs = psplot.thistograph(xs,g3)
    return hs
    
def exgaussline(mu=1.,sig=1.,a=0.5,b=2.,n=1000):
    g = pspdf.pdf_gauss(mu,sig)
    xs = g.random(n)
    ys = map(lambda x: x*a+b,xs)
    gg = pspdf.pdf_gauss(a*mu+b,a*sig)
    hs = psplot.thistograph(ys,gg)
    return hs

def exbinomial(n=1000):
    p,n =0.5,100
    pf = pspdf.pmf_binomial(n,p)
    xs = pf.random(n)
    h1 = psplot.thisto(xs)
    p,n =0.01,100
    pf = pspdf.pmf_binomial(n,p)
    xs = pf.random(n)
    h2 = psplot.thisto(xs)
    return h1,h2
    


#-----------------------------------------------------
# Estimation of paramters using Log Likelihood
#-----------------------------------------------------

  

    
        
