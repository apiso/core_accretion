"""
general purpose tools that could be useful for a range of problems

TODO: should put zbrac and nextguess into `opttools.py` or something...
"""

from utils.zbrac import zbrac
from utils.myspace import mylogspace, loglinspace, linlogspace
from utils.constants import Me, Re, Myr

def nextguess(xnew, ys, xs, maxdy = 0):
    """
    predicts y(xnew) from linear extraplotion of two prior solutions passed as
    `xs` = [x0,x1] and `ys` = [y(x0),y(x1)].  Assumes ascending or descending
    order for x0,x1,xnew.  Limits step in dy if `maxdy` nonzero.  Step limit
    only works for scalar y values.
    """
    if type(ys[0]) == (list or tuple):
        ys = (np.array(ys[0]), np.array(ys[1]))
    dypred = (ys[1] - ys[0]) / (xs[1] - xs[0]) * (xnew - xs[1])
    if maxdy != 0:
        if abs(dypred) > maxdy*abs(ys[1] - ys[0]):
            print "maxdy step limited to", maxdy
            dypred = maxdy*abs(ys[1] - ys[0])
    return ys[1] + dypred

class nf(float):
    """
    Defines a class that forces representation of float to have only one digit
    after decimal.  Also removes trailing zero so '1.0' becomes '1'.
    Found here: http://matplotlib.org/examples/pylab_examples/contour_label_demo.html
    """
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1]=='0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
