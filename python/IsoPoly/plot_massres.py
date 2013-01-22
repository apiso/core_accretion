getdat = 0
savefig = 1
figfolder = "/figs/ModelAtmospheres/isopoly/"
id1 = "Hillmc5au10"
id2 = "HillLog1000"
figfile = "IPmasses_Hill.pdf"
if getdat:
    from utils.constants import Me, Re, G
    from IsoPoly import energetics as en    
    file1, file2 = id1+".npz" , id2+".npz"
    model1 , arr1 = en.modelseriesload(file1)
    model2 , arr2 = en.modelseriesload(file2)

    Mcbav1 , Mcbav2 = 0.5*(arr1.Mcb[:-1] + arr1.Mcb[1:]) , 0.5*(arr2.Mcb[:-1] + arr2.Mcb[1:])
    Matm1, Matm2 = Mcbav1 - mco, Mcbav2 - mco
    dMonM1, dMonM2 = (arr1.Mcb[1:] - arr1.Mcb[:-1])/Matm1 , (arr2.Mcb[1:] - arr2.Mcb[:-1])/Matm2
    dPcb1 , dPcb2 = (arr1.Pcb[1:] - arr1.Pcb[:-1]) , (arr2.Pcb[1:] - arr2.Pcb[:-1])
    Pcbav1 , Pcbav2 = 0.5*(arr1.Pcb[:-1] + arr1.Pcb[1:]) , 0.5*(arr2.Pcb[:-1] + arr2.Pcb[1:])
    dlMdlP1 , dlMdlP2 = dMonM1 / dPcb1 * Pcbav1 , dMonM2 / dPcb2 * Pcbav2

close(4)
f , (ax, ax2) = plt.subplots(2,figsize = (5,7),sharex = False, num = 4)

ax.semilogy(dMonM1, 'b-', label = r"even log")
ax.semilogy(-dMonM1, 'b--')
ax.semilogy(dMonM2, 'r-', label = r"finer @ max")
ax.semilogy(-dMonM2, 'r--')
ax.set_xlabel(r"$i$")
ax.set_ylabel(r"$\Delta M_i/M_\mathrm{atm}$")
ax.legend(loc = 'lower right', title = r'$P_\mathrm{cb}$ spacing', frameon = False, borderpad = 0)

ax2.loglog(Matm2/Me , dlMdlP2, 'r-')
ax2.loglog(Matm2/Me , -dlMdlP2, 'r--')
ax2.loglog(Matm1/Me , dlMdlP1, 'b-')
ax2.loglog(Matm1/Me , -dlMdlP1, 'b--')

ax2.set_xlabel(r'$M_\mathrm{atm} \; [M_\oplus]$')
ax2.set_ylabel(r'$d [ln (M_\mathrm{atm})] /d [ln (P_\mathrm{cb})]$')
ax2.set_xlim(xmin = .03)
ax2.set_ylim(3e-3,2e4)

plt.draw()
plt.tight_layout()
#f.subplots_adjust(hspace=0)
plt.draw()

if savefig:
    from utils.userpath import userpath
    plt.savefig(userpath+figfolder+figfile)
