"""
This script makes plot "IsoRoutMout.pdf" showing how isothermal atmosphere around higher mass cores cannot have large outer radii!
The disk model in these plots are the standard 10AU values.  The dashed curves are the masses corrected for the disk mass.  All atmospheres have a base at 0.1 RBc (Bondi radius as set by core mass).

Useage:
in python folder:
import userpath
cd isothermal
execfile("plot_rout_vs_m.py")

you may need a plt.show() if not in interactive mode...
"""

from isotherm import *
from isoout import *
import matplotlib.pyplot as plt

savethefig = 0

cmE = [1 , 3 , 5 , 10] #separate array in earth masses good for plot legends!
cmasses = np.array(cmE)*Me #range of "core masses" to consider, note that core is base of isothermal layer at 0.1 rBc
LPovermax = (4.39 , 4.5038 , 8 , 7) #log10 of max overpressure at base
nPr = 101 #number of base pressures to consider

mall = np.zeros( (len(cmasses) , nPr) )
rall = np.zeros( (len(cmasses) , nPr) )
mcorrall = np.zeros( (len(cmasses) , nPr) )

for i , mc in enumerate(cmasses):
    rall[ i , : ] , mall[ i , : ] , Pa = varyPbase( 1 , LPovermax[i] , nPr , mc)
    mcorrall[ i , : ] = mnorms(mall[ i , : ] , rall[ i , : ] , mc) #disk corrected masses


#begin plotting
f = plt.figure(figsize = (5,4))
colors = ['r' , 'g' , 'b' , 'k']
for i , mc in enumerate(cmasses):
    plt.loglog(mall[ i , : ] - 1 , rall[ i , : ] , c = colors[i] , label = str(cmE[i]) )
    plt.loglog(mcorrall[ i , : ] - 1 , rall[ i , : ] , '--' , c = colors[i])

plt.xlim(1e-4 , 6)
plt.ylim(ymax = 30)
plt.xlabel('$M_{\mathrm{atm}} / M_\mathrm{c}$')
plt.ylabel('$R_{\mathrm{out}} / R_\mathrm{B,c}$')
plt.legend(title = r'$M_\mathrm{c} [M_\oplus]$' , frameon = False , loc = 'upper left' , borderpad = 0)

plt.draw()
plt.tight_layout()
plt.draw()
if savethefig:
    plt.savefig("../../figs/ModelAtmospheres/isothermal/IsoRoutMout.pdf")
