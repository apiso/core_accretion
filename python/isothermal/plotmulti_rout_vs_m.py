"""
similar to plot_rout_vs_m.py but considers the the total Bondi (not core) and Hill radii.
Solid (dashed) curves use total (disk corrected) mass for both atmospheres mass and characteristic radii (Hill and Bondi)
"""

#or comment this line out if that script already run
#execfile("plot_rout_vs_m.py")

savethefigs = 1

#do all the radius normalizations (using isoout.rnorms)
rm_Bt , rm_Bds , rm_Hc , rm_Ht , rm_Hds = np.zeros( (len(cmasses) , nPr) ) , np.zeros( (len(cmasses) , nPr) ) , np.zeros( (len(cmasses) , nPr) ) ,  np.zeros( (len(cmasses) , nPr) ) ,  np.zeros( (len(cmasses) , nPr) )
for i , mc in enumerate(cmasses):
    rm_Bt[i,:] , rm_Bds[i,:] , rm_Hc[i,:] , rm_Ht[i,:] , rm_Hds[i,:] = rnorms(mall[ i , : ] , rall[ i , : ] , mc)

#start plotting
f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize= (5,7))

for i , mc in enumerate(cmasses):
    ax1.loglog(mall[ i , : ] - 1 , rm_Bt[ i , : ] , c = colors[i] , label = str(cmE[i]) )
    ax1.loglog(mcorrall[ i , : ] - 1 , rm_Bds[ i , : ] , '--' , c = colors[i])
    ax2.loglog(mall[ i , : ] - 1 , rm_Ht[ i , : ] , c = colors[i] , label = str(cmE[i]) )
    ax2.loglog(mcorrall[ i , : ] - 1 , rm_Hds[ i , : ] , '--' , c = colors[i])

plt.xlim(1e-4,10)
ax1.set_ylim(.01,30)
ax2.set_ylim(.02,5)
ax2.legend(title = r'$M_\mathrm{c} [M_\oplus]$' , frameon = False , loc = 'upper left' , borderpad = 0)
plt.xlabel('$M_{\mathrm{atm}} / M_\mathrm{c}$')
ax1.set_ylabel('$R_{\mathrm{out}} / R_\mathrm{B}$')
ax2.set_ylabel('$R_{\mathrm{out}} / R_\mathrm{H}$')
plt.draw()
plt.tight_layout()
f.subplots_adjust(hspace = 0)
plt.draw()
if savethefigs:
    plt.savefig("../../figs/ModelAtmospheres/isothermal/IsoRoutMoutBH.pdf")


#no-self-gravity data
LPmax_nsg = [4.328,4.298,4.268 , 4.194]
r_nsg , m_nsg , mcorr_nsg = np.zeros( (len(cmasses) , nPr) ) , np.zeros( (len(cmasses) , nPr) ) , np.zeros( (len(cmasses) , nPr) ) 
for i , mc in enumerate(cmasses):
    r_nsg[ i , : ] , m_nsg[ i , : ] , Pa = varyPbase( 1 , LPmax_nsg[i] , nPr , mc , sg = 0)
    mcorr_nsg[ i , : ] = mnorms(m_nsg[ i , : ] , r_nsg[ i , : ] , mc) #disk corrected masses

#no-self-gravity plot
f = plt.figure(figsize = (5,4))
colors = ['r' , 'g' , 'b' , 'k']
for i , mc in enumerate(cmasses):
    plt.loglog(m_nsg[ i , : ] - 1 , r_nsg[ i , : ] , c = colors[i] , label = str(cmE[i]) )
    plt.loglog(mcorr_nsg[ i , : ] - 1 , r_nsg[ i , : ] , '--' , c = colors[i])
    
plt.xlim(1e-4 , 10)
plt.ylim(ymax = 30)
plt.xlabel('$M_{\mathrm{atm}} / M_\mathrm{c}$')
plt.ylabel('$R_{\mathrm{out}} / R_\mathrm{B,c}$')
plt.legend(title = r'$M_\mathrm{c} [M_\oplus]$' , frameon = False , loc = 'upper left' , borderpad = 0)
plt.title('No Self-Gravity')
plt.draw()
plt.tight_layout()
plt.draw()

if savethefigs:
    plt.savefig("../../figs/ModelAtmospheres/isothermal/IsoRoutMoutNSG.pdf")
