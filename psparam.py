
#  PF parameters
NPOINTS = 100 # number of points to build (x,y) functions
NITER = 8 # maximun number of iteration
EPSILON = 5.e-4  # precision required to stop iterations
INFTY = 1.e24 # large real value associated to infinity
NIFTY = int(INFTY) # large int value associated to infinity
ZERO = 1.e-24 # zero, minimun difference between two numbers
PMIN = 1.e-7 # minimum probability

#  Log Likelihood estimation parameters
ELLNSIG = 5 # number of sigmas to get a guess of a parameter
ELLNITER = 12 # maximum number of iterations
ELLEPSILON = 1.e-2 # precision of the paramter estimate of LL
ELLNPOINTS = 1 # number of points in the LL scan 
ELLDELTAMAX = 1. # log likelihood delta range
ELLDELTAMIN = 0.5 # log likelihood delta range for one sigma

ELLNPOINTS = 5 # number of points in the LL scan 
ELLDELTAMAX = 1. # log likelihood delta range

# debug level
DEBUG = 0 # level of print out (0=none,1=info,2-verbose,3-very verbose)
