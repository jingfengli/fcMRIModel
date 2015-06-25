# -*- coding: utf-8 -*-
#!/usr/bin/env python
# rk methods

"""A variety of methods to solve first order ordinary differential equations.

AUTHOR:
    Jonathan Senning <jonathan.senning@gordon.edu>
    Gordon College
    Based Octave functions written in the spring of 1999
    Python version: March 2008, October 2008
"""
from boldodes import boldodes
import numpy
def rk2( f, x0, t, **kwargs):
    """Second-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        x = rk2a(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.

    NOTES:
        This version is based on the algorithm presented in "Numerical
        Analysis", 6th Edition, by Burden and Faires, Brooks-Cole, 1997.
    """
    
    n = len( t )
    x = numpy.array( [ x0 * 1.0 ] * n )
    
    for i in xrange( n - 1 ):
        h = t[i+1] - t[i]
        k1 = h * f( x[i], t[i], **kwargs) / 2.0

        x[i+1] = x[i] + h * f( x[i] + k1, t[i] + h / 2.0 ,**kwargs)

    return x

def rk4( f, x0, t, **kwargs):
    """Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        x = rk4(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a NumPy array.  In this
                case f must return a NumPy array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or NumPy array
                if a system of equations is being solved.
        t     - list or NumPy array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - NumPy array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.
    """
    
    n = len(t)

    x = numpy.array( [x0 * 1.0] * n)
    
    for i in xrange(n - 1):
        h = t[i+1] - t[i]
        
        k1 = h * f(x[i],t[i], **kwargs)
        k2 = h * f(x[i] + k1 / 2.0, t[i] + h / 2.0, **kwargs)
        k3 = h * f(x[i] + k2 / 2.0, t[i] + h / 2.0, **kwargs)
        k4 = h * f(x[i] + k3, t[i] + h, **kwargs)
        
        x[i+1] = x[i] + (k1 + 2.0 * k2 + 2.0 * k3 + k4)/6.0
    
    return x


        
        
        
        
        