getdat = 1
savefig = 1
figfolder = "/figs/ModelAtmospheres/isopoly/"
#id = "Hillmc5au10"
id = "HillLog200"
figfile = "IPcool_"+id+".pdf"

if getdat:
    ngrid = 100
    from utils.constants import Me, Re, G
    from IsoPoly import energetics as en 
    #file = "mc5au10_N1000.npz"
    file = id + ".npz"
    model , arr = en.modelseriesload(file)
    pfs = en.profiles(ngrid, file)
    LdtdM , LdtdM_cb = en.explicitLarr(pfs, model)
    Tcb , mco , R , Po , delad = model.Tcb[0] , model.mco[0] , model.R[0] , model.Po[0] , model.delad[0]
    Mcbav = 0.5*(arr.Mcb[:-1] + arr.Mcb[1:])
    Matm = Mcbav - mco    
    barP = 1e6

close('all')
f , ax = plt.subplots(1,figsize = (5,4),sharex = True, num = 1)

ax.loglog(Matm/Me, LdtdM, 'b', label = r'$+L / (dM/dt)$')
ax.loglog(Matm/Me, -LdtdM, 'b--')#, label = r'$-L / (dM/dt)$')
ax.loglog(Matm/Me, LdtdM_cb, 'g', label = r'$+L_\mathrm{CB} / (dM/dt) $')
ax.loglog(Matm/Me, -LdtdM_cb, 'g--')#, label = r'$-L_\mathrm{CB} / (dM/dt) $')
ax.set_xlabel(r'$M_\mathrm{atm} \; [M_\oplus]$')
ax.set_ylabel(r'$dE/dM \; [\mathrm{erg/g}]$')
ax.set_ylim(8e6,1e12)
ax.set_xlim(xmin = 3e-2)
ax.legend(frameon = False , loc = 'upper center',borderpad = 0 )

plt.draw()
plt.tight_layout()
#f.subplots_adjust(hspace=0)
plt.draw()

if savefig:
    from utils.userpath import userpath
    plt.savefig(userpath+figfolder+figfile)
