import math
from psparam import *
import psutils

"""
Module to deal with Bayesian Statistics

@version 1.0
@author Jose A. Hernando
@date 05/11/13 
"""

class PBayes:
    """ Base class for Bayes Probabilities.
    Default implementation for dicrete probabilities (PNF and PMF).
    """
    def __init__(self,fprior,hypos,flike,name='PBayes'):
        """ constructor with the likelihood function pf(D|H), the prior pf(H), 
        and the list of hypothesis (H),
        that correspond to the possible random variables of the fprior.
        if hypothesis are continous provide a small enough partition.     
        """
        if (not isinstance(hypos,list)):
            raise TypeError('hypothesis should be a list! ')
        self._hypos = hypos
        zs = zip(hypos,map(fprior,hypos))
        self._fprior = self.normalize(zs)
        self._flike = flike
        self._name = name
        self._fpost = None
        
    def hypothesis(self):
        """ returns a list with the hypothesis """
        return self._hypos

    def prior(self,hypo):
        """ returns the probability of the hipothesis fprior(hip)
        """
        if (isinstance(hypo,list)):
            return map(self._fprior,hypo)
        val = self._fprior(hypo)
        if (DEBUG): print " prior ",hipo,val
        return self._fprior(hypo)

    def posterior(self,hypo):
        """ returns the pf of the hipothesis fposterior(hipo)
        """
        if (not self._fpost):
            raise ArithmeticError('Not valid posterior yet!')
        if (isinstance(hypo,list)):
            return map(self._fpost,hypo)
        val = self._fpost(hypo)
        if (DEBUG): print ' posterior ',hipo,val
        return self._fpost(hypo)
        
    def likelihood(self,data,hypo):
        """ returns the likelihood of the data for a given hipothesis
        """
        if (isinstance(data,list)): 
            vals = map(lambda d : self.likelihood(d,hypo),data)
            val = reduce(lambda x,y: x*y, vals)
            return val
        val = self._fprior(hypo)*self._flike(data,hypo)
        if (DEBUG): print ' P(H) P(D|H) ',data,hypo,' = ',val
        return val
    
    def llikelihood(self,data,hypo):
        """ returns the log-likelihood of the data for a given hipothesis
        """
        if (isinstance(data,list)): 
            vals = map(lambda d:self.llikelihood(d,hypo),data)
            val = sum(vals)
            return val
        val = math.log(max(ZERO,self.likelihood(data,hypo)))
        if (DEBUG): print ' log P(H) P(D|H) ',data,hipo,' = ',val
        return val
        
    def evidence(self,data):
        """ given a data compute the posterior
        """
        hs = self.hypothesis()
        if (self._fpost): self._prior = self._fpost
        lls = map(lambda h: self.llikelihood(data,h),hs)
        llmax = max(lls)
        llrats = map(lambda ll: ll-llmax,lls)
        numes = map(lambda lrat: math.exp(lrat),llrats)
        deno = sum(numes)
        ps = map(lambda nume: nume/deno,numes)
        zs = zip(hs,ps)
        self._fpost = self.normalize(zs)
        if (DEBUG): print ' evidence P(H|D) : ',zs
        return

    def normalize(self,zs):
        """ normalize for discrete hypothesis """
        return psutils.Funxy_xdiscrete(zs)
        
class PDBayes(PBayes):
    """ Bayes probability for continuous hypothesis """
    
    def __init__(self,fprior,hypos,flike,name='PDBayes'):
        """ constructor of Bayes probability for continuous hypothesis.
        fprior is the pdf of the prior.
        hypos is a small enough partition of all the possible hypothesis.
        flike is the likelihood p(data | hypo)
        """
        PBayes.__init__(self,fprior,hypos,flike)

    def normalize(self,zs):
        """ normalize pdf """
        return psutils.Funxy_xinterpolate(zs)
