getdat = 1
savefig = 1
figfolder = "/figs/ModelAtmospheres/isopoly/"
#id = "Hillmc5au10"
id = "HillLog200"
figfile = "IPcoolglobal_"+id+".pdf"

if getdat:
    ngrid = 100
    from utils.constants import Me, Re, G
    #file = "mc5au10.npz"
    file = id+".npz"
    from IsoPoly import energetics as en    
    model , arr = en.modelseriesload(file)
    pfs = en.profiles(ngrid, file)
    vr = en.virialterms(file)
    LdtdM , LdtdM_cb = en.explicitLarr(pfs, model)
    dEdM_neg , eacc , work = en.globalLarr(pfs, model, vr)
    LdtdM_glob = dEdM_neg + eacc + work
    Tcb , mco , R , Po , delad = model.Tcb[0] , model.mco[0] , model.R[0] , model.Po[0] , model.delad[0]
    Mcbav = 0.5*(arr.Mcb[:-1] + arr.Mcb[1:])
    Matm = Mcbav - mco    
    barP = 1e6

close('all')
f , ax = plt.subplots(1,figsize = (5,6),sharex = True, num = 1)

ax.loglog(Matm/Me, LdtdM + LdtdM_cb, 'b', label = r'$+L_\mathrm{int} / (dM/dt)$')
ax.loglog(Matm/Me, -LdtdM -LdtdM_cb , 'b--')#, label = r'$-L / (dM/dt)$')
#ax.loglog(Matm/Me, LdtdM_cb, 'g', label = r'$+L_\mathrm{CB} / (dM/dt) $')
#ax.loglog(Matm/Me, -LdtdM_cb, 'g--')#, label = r'$-L_\mathrm{CB} / (dM/dt) $')
ax.loglog(Matm/Me, LdtdM_glob, 'k', label = r'$+L_\mathrm{glob} / (dM/dt)$')
ax.loglog(Matm/Me, -LdtdM_glob, 'k--')#, label = r'$-L / (dM/dt)$')
                                      #l1 = ax.legend(frameon = False , loc = 'upper center' , borderpad = 0)
ax.loglog(Matm/Me, dEdM_neg, 'r', label = r'$- dE/dM$')
ax.loglog(Matm/Me, -dEdM_neg, 'r--')#, label = r'$+ dE/dM$')
ax.loglog(Matm/Me, eacc, 'g', label = r'$+ e_\mathrm{acc}$')
ax.loglog(Matm/Me, -eacc, 'g--')#, label = r'$+ dE/dM$')
ax.loglog(Matm/Me, work, 'y', label = r'$ - P \Delta V/\Delta M$')
ax.loglog(Matm/Me, -work, 'y--')#, label = r'$+ dE/dM$')

ax.set_xlabel(r'$M_\mathrm{atm} \; [M_\oplus]$')
ax.set_ylabel(r'$dE/dM \; [\mathrm{erg/g}]$')
ax.set_ylim(3e7,2e11)
ax.set_xlim(xmin = 4e-2)
ax.legend(frameon = False , loc = 'upper right' , borderpad = 0, labelspacing = 0)
#h, l = ax.get_legend_handles_labels()
#ax.legend(h[2:],l[2:], frameon = False , loc = 'lower left', borderpad = 1)# , borderpad = 0)
#gca().add_artist(l1)

plt.draw()
plt.tight_layout()
#f.subplots_adjust(hspace=0)
plt.draw()

if savefig:
    from utils.userpath import userpath
    plt.savefig(userpath+figfolder+figfile)
