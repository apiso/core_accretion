"""
script to generate ErrMap2D.pdf, 2D errors in mass and luminsity for shooting problem
"""

getdat = 1
tspace = 'lin'

if getdat:
    from RadSGPoly import shoot2 as s2
    from utils.mynumsci import Me, Re, Myr
    Mold, Mc = 5.8*Me, s2.prms.mco
    prev = s2.makeprev(Mold)
    Mnew = Mc + 1.01*(Mold-Mc)
    lLg, dtg = s2.staticguess(Mnew, prev, dt_log=0)
    #lLm, dtm, errmap = s2.errormap(Mnew, .98*lLg, 1.02*lLg, .5*dtg, 2.5*dtg,
    #                               prev, tspace='lin', nL=5, ndt=5)
    #lLm, dtm, errmap = s2.errormap(Mnew, 24.6, 24.65, .6*dtg, 1.5*dtg,
    #                               prev, tspace=tspace, nL=8, ndt=8)
    lLm, dtm, errmap = s2.errormap(Mnew, 24.6, 24.65, .012*Myr, .026*Myr,
                                   prev, tspace=tspace, nL=30, ndt=30)

ext = (dtm[0]/Myr,dtm[-1]/Myr,lLm[0],lLm[-1])
DT, LL = np.meshgrid(dtm/Myr, lLm)
Lscale = 1#5e-3#factor by which to reduce luminosity errors
close(1)

f, (ax1, ax2, ax3) = plt.subplots(3, figsize = (5,10),sharex = True,num = 1)
im1 = ax1.imshow(errmap[:,:,0],origin = 'lower', extent = ext, 
                 vmin = -errmap[:,:,0].max())
ax1.set_aspect('auto')
if tspace == 'log':
    ax1.set_xscale('log') 
ax1.contour(DT, LL, errmap[:,:,0], colors = 'k', levels = (0,))
ax1.contour(DT, LL, errmap[:,:,1], colors = 'k', linestyles = 'dashed',
            levels = (0,))
ax1.set_ylabel(r'$\mathrm{log}_{10}(L)\;[\mathrm{erg/s}]$')
ax1.set_title(r'$M_\mathrm{c} \,\mathrm{relative \, error}$')
#ax1.set_title(r'mass error')
cbar = plt.colorbar(im1, ax = ax1,use_gridspec = True)

im2 = ax2.imshow(errmap[:,:,1]*Lscale,origin = 'lower', extent = ext,
                 vmin = -errmap[:,:,1].max()*Lscale)#, vmin=-10)
ax2.set_aspect('auto')
if tspace == 'log':
    ax2.set_xscale('log') 
ax2.contour(DT, LL, errmap[:,:,0], colors = 'k',linestyles = 'dashed', levels = (0,))
ax2.contour(DT, LL, errmap[:,:,1], colors = 'k', levels = (0,))
ax2.set_ylabel(r'$\mathrm{log}_{10}(L)\;[\mathrm{erg/s}]$')
ax2.set_title(r'$L_\mathrm{c}/L_\mathrm{out}\, \mathrm{error}$')
plt.colorbar(im2, ax = ax2, use_gridspec = True)

    #num = 3
    #close(num)
    #f, ax = plt.subplots(1,num = num)
im = ax3.imshow(np.log10((errmap[:,:,0]**2 + (errmap[:,:,1] * Lscale)**2)**0.5),origin = 'lower', extent = ext)#,vmax = 1e21)
cbar = f.colorbar(im, ax = ax3, use_gridspec = True)


ax3.set_aspect('auto')
ax3.set_title('$\mathrm{log_{10}(RMS\, error)}$')
if tspace == 'log':
    ax3.set_xscale('log')
    ax3.set_xlabel(r'$\mathrm{log}_{10}(\Delta t)\;[\mathrm{Myr}]$')
else:
    ax3.set_xlabel(r'$\Delta t\;[\mathrm{Myr}]$')
ax3.set_ylabel(r'$\mathrm{log}_{10}(L)\;[\mathrm{erg/s}]$')




ax3.contour(DT, LL, errmap[:,:,0], colors = 'k', levels = (0,))
ax3.contour(DT, LL, errmap[:,:,1], colors = 'k', linestyles = 'dashed',levels = (0,))
#cont = ax3.contour(DT, LL,np.log10((errmap[:,:,0]**2 + (errmap[:,:,1] * Lscale)**2)**0.5) , colors = 'k', levels = (-3,-2))
#ax3.clabel(cont, inline = 1, fmt = '%1.1f')

plt.draw()
plt.tight_layout()

savefig = 0
if savefig:
    from utils.userpath import userpath
    figfolder = userpath + "/figs/ModelAtmospheres/RadSelfGravPoly/Hill2/"  
    plt.savefig(figfolder+'ErrMap2D.pdf')
