getdat = 0
savefig = 1
figfolder = "/figs/ModelAtmospheres/isopoly/"
id = "Hillmc5au10"
figfile = "IPprofiles_"+id+".pdf"

if getdat:
    ngrid = 100
    from utils.constants import Me, Re, G
    #file = "mc5au10_N1000.npz"
    file = id+".npz"
    from IsoPoly import energetics as en    
    model , arr = en.modelseriesload(file)
    pfs = en.profiles(ngrid, file)
    Tcb , mco , R , Po , delad = model.Tcb[0] , model.mco[0] , model.R[0] , model.Po[0] , model.delad[0]
    Matm = arr.Mcb - mco    
    barP = 1e6

close('all')
f , (ax1,ax2) = plt.subplots(2,figsize = (5,7),sharex = True, num = 1)

step = 20
natm = arr.shape[0]
for i in range(natm-1,-1,-step):
#for i in range(0,natm,step):
    fi = float(i)/(natm -1)
    #print fi
    #ax1.loglog(pfs.P[i , :]/barP , pfs.TK[i , :] , c = cm.rainbow(fi))
    ax1.loglog(pfs.r[i , :]/Re , (pfs.m[i , :]- mco)/Me , '', c = cm.rainbow(fi))
    ax2.loglog(pfs.r[i , :]/Re , pfs.P[i , :]/barP , '' , c = cm.rainbow(fi))
    ax1.plot(pfs.r[i , -1]/Re, (pfs.m[i, -1] - mco)/Me, 'o', ms = 7, c = cm.rainbow(fi), mec = cm.rainbow(fi), mew = 0)
    ax2.plot(pfs.r[i , -1]/Re, pfs.P[i, -1]/ bar, 'o', ms = 7, c = cm.rainbow(fi), mec = cm.rainbow(fi), mew = 0)
     

ax1.set_ylabel(r'$M_\mathrm{atm}\; [M_\oplus]$')
ax2.set_ylabel(r'$P\;[\mathrm{bar}]$')
ax2.set_xlabel(r'$r \; [R_\oplus]$')
ax1.set_xlim(1.5e0, 4e4)
ax1.set_ylim(ymin = 4e-6 , ymax = 2e3)
    #ax1.invert_xaxis()

plt.draw()
plt.tight_layout()
f.subplots_adjust(hspace=0)
plt.draw()

if savefig:
    plt.savefig(userpath+figfolder+figfile)
