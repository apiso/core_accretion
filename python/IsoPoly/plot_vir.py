getdat = 1
savefig = 0
figfolder = "/figs/ModelAtmospheres/isopoly/"
#id = "Hillmc5au10"
id = "HillLog200"
figfile = "IPenergy_"+id+".pdf"

if getdat:
    from utils.constants import Me, Re, G
    file = id+".npz"
    #file = "mc10au10.npz"
    from IsoPoly import energetics as en
    vr = en.virialterms(file)
    model , arr = en.modelseriesload(file)
    Tcb , mco , R , Po , delad = model.Tcb[0] , model.mco[0] , model.R[0] , model.Po[0] , model.delad[0]
    Matm = arr.Mcb - mco
    Matm_vir = vr.Mvir - mco
    zeta = 3 * delad / (1 - delad)
    Ei = -vr.Iu/zeta
    barP = 1e6

close('all')
f , (ax1,ax2) = plt.subplots(2,figsize = (5,7),sharex = True, num = 1)
#f , ax1 = plt.subplots(1,figsize = (5,4),num = 1)
lw = 2
p1, = ax1.loglog(Matm_vir/Me, -vr.Eg, 'k', lw = lw, label = r'$-E_G$')
p2, = ax1.loglog(Matm_vir/Me, -vr.Iu, 'r--', lw = lw, label = r'$I_u$' )
l1 = ax1.legend(frameon = False , loc = 'upper left')# , borderpad = 0)
p3, = ax1.loglog(Matm_vir/Me, vr.top, lw = lw, label = r'$\mathrm{top}$')
p4, = ax1.loglog(Matm_vir/Me, -vr.bottom , '--', lw = lw,  label = r'$\mathrm{bottom}$')
ax1.set_ylabel(r'$\mathrm{Energy \; [erg]}$')
#ax1.set_xlabel(r'$M_\mathrm{atm} \; [M_\oplus]$')
h, l = ax1.get_legend_handles_labels()
ax1.legend(h[2:],l[2:], frameon = False , loc = 'lower right')# , borderpad = 0)
gca().add_artist(l1)
ax1.set_ylim(ymax = 1.2*max(-vr.Eg))

plt.draw()
plt.tight_layout()
#f.subplots_adjust(hspace=0)
plt.draw()

#f2 , ax2 = plt.subplots(1,figsize = (5,4),num = 2)
escale = 1e40
ax2.semilogx(Matm_vir/Me, vr.Etot/escale, lw = lw, label = r'$E_\mathrm{tot}$')
ax2.semilogx(Matm_vir/Me, vr.Eg/escale, lw = lw, label = r'$E_G$')
ax2.semilogx(Matm_vir/Me, Ei/escale, lw = lw, label = r'$E_i$')
ax2.set_ylabel(r'$\mathrm{Energy \; [10^{39} erg]}$')
ax2.set_xlabel(r'$M_\mathrm{atm} \; [M_\oplus]$')
ax2.legend(frameon = False , loc = 'upper left')
ax2.set_ylim(ymax = 9)
ax2.set_ylim(ymin = -11)
ax2.set_xlim(xmin = 2e-2)
plt.draw()
plt.tight_layout()
f.subplots_adjust(hspace=0)
plt.draw()

if savefig:
    plt.savefig(userpath+figfolder+figfile)

#semilogx(Matm_vir/Me, -vr.vsum/vr.Eg)
