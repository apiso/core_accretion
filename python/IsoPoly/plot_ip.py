
getdat = 1
savefig = 0
figfolder = "/figs/ModelAtmospheres/isopoly/"
#id = "Hillmc5au10"
id = "HillLog200"
figfile = "IsoPoly_"+id+".pdf"

if getdat:
    from utils.constants import G, Me, Re
    from utils.userpath import userpath
    datfolder = userpath + "/dat/MODELS/isopoly/"
    #file = "mc5au10.npz"
    file = id+".npz"
    npzdat = np.load(datfolder + file)
    model = npzdat['model'].view(np.recarray)
    arr = npzdat['atmser'].view(np.recarray)
    Tcb , mco , R , Po = model.Tcb[0] , model.mco[0] , model.R[0] , model.Po[0]
    RBc = G * mco / (R *Tcb)
    RB = G * arr.Mcb / (R*Tcb)
    Matm = arr.Mcb - mco
    bar = 1e6

close(1)
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize= (5,10), num = 1)
#f = plt.figure(1, figsize = (5,4))
ax1.loglog(Matm/Me, arr.Pcb/bar)
ax1.set_ylabel(r'$P_\mathrm{CB}\;[\mathrm{bar}]$')
ax1.set_ylim(ymin= 0.8*Po/bar)
ax1.set_ylim(ymax= 1.2*max(arr.Pcb)/bar)

#f = plt.figure(2, figsize = (5,4))
ax2.loglog(Matm/Me, arr.Pc/bar)
#plt.xlabel(r'$M_\mathrm{atm} [M_\oplus]$')
ax2.set_ylabel(r'$P_\mathrm{c}\;[\mathrm{bar}]$')
ax2.set_ylim(ymin = 0.8*arr.Pc[0]/bar)
ax2.set_ylim(ymax= 1.2*max(arr.Pc)/bar)

#f = plt.figure(3, figsize = (5,4))
ax3.loglog(Matm/Me, arr.Rcb/RBc, label = r'$R_\mathrm{B,c}$')
ax3.loglog(Matm/Me, arr.Rcb/RB, label = r'$R_\mathrm{B,tot}$')
#plt.xlim(xmin = 4.)
ax3.set_ylim(5e-2, 2e1)
ax3.set_ylabel(r'$R_\mathrm{CB} / R_\mathrm{B}$')
ax3.legend(frameon = False , loc = 'upper left' )#, borderpad = 0)

ax3.set_xlabel(r'$M_\mathrm{atm} \; [M_\oplus]$')
plt.xlim(xmin=3e-2) # for mc = 5

ax1.set_title(r"IsoPoly "+id)

plt.draw()
plt.tight_layout()
f.subplots_adjust(hspace=0)
plt.draw()

if savefig:
    plt.savefig(userpath+figfolder+figfile)
