"""
Basic plots of EOS tables.
To be used in pylab enviroment (since code assumes many things imported)

makes all plots with
>execfile('eosplotsscript.py') 
"""

#allow import from parent folder
import np, os, sys
lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from userpath import userpath
from EOSread import EOSread
from table_size import Nt, Np, Ncol

#read tables and say where the variables are
EOS , T , NpV = EOSread(Nt, Np, Ncol, userpath + 'dat/EOStables/H_HE_TAB_I_30_EXT.txt')
ent , delad , u = 4 , 10 , 5
P = EOS[-1 , : , 0]

#create subset of table for plotting
minp = 42 #if 42 then the extrapolation to low P (ugly since doesn't include dissociation) isn't shown
maxt, maxp = 45 , 78
eos = EOS[:maxt , minp : maxp , :]
t , p = T[:maxt], eos[-1 , : , 0]
pbar = 10**(p-6)
tk = 10**t
PP , TT = np.meshgrid(pbar,tk)


"""
plot entropy (log axes)
"""

arr = eos[: , : , 4]
fig = figure(figsize = (5,3.5))
ax = fig.add_subplot(111)
ext = [min(pbar),max(pbar),min(tk),max(tk)]
ax.set_yscale('log')
ax.set_xscale('log')
im = ax.imshow(arr, extent = ext, vmin = 8.3, vmax = 9.5, interpolation = 'nearest', origin = 'lower')
cbar = fig.colorbar(im)
cbar.set_label('log(S)')
ax.set_aspect('auto')
ax.set_ylabel('T [K]')
ax.set_xlabel('P [bar]')
cont = ax.contour(PP, TT, arr, colors = 'k', levels = np.arange(8.4,9.41,0.2))
plt.clabel(cont, inline = 1, fmt='%1.1f')

plt.draw()
plt.tight_layout()
fig.subplots_adjust(top = 0.90 , right = 0.74)

#savefig(userpath+'figs/EOS/S.pdf')


"""
plot entropy on delad contours
"""

from astropysics.constants import h, kb, mp

fig = figure(figsize = (5,3.5))
ax = fig.add_subplot(111)
ext = [min(pbar),max(pbar),min(tk),max(tk)]
ax.set_yscale('log')
ax.set_xscale('log')
arr = eos[: , : , delad]
im = ax.imshow(arr, extent = ext, vmin = np.min(arr[arr > 0]), interpolation = 'nearest', origin = 'lower')
cbar = fig.colorbar(im)
cbar.set_label(r'$\nabla_\mathrm{ad}$')
cont = ax.contour(PP, TT, eos[: , : , ent], colors = 'grey' , levels = arange(8.6,9.3,.1))
#cont = ax.contour(, extent = ext, origin = 'lower', colors = 'k', levels = arange(8.6,9.4,.05))
clabel(cont, inline = 1, fmt = '%1.1f')

ax.set_aspect('auto')
ax.set_ylabel('T [K]')
ax.set_xlabel('P [bar]')

Tdiss = 51800.
tem = 10**t
Pdiss = (kb * tem)**2.5 / 4. * (pi * mp / h**2)**1.5 * exp(-Tdiss / tem)
#ax.plot(Pdiss , tem, color = 'white')
#plt.xlim(PP.min(),PP.max())

#mask below a given entropy contour
minent = 8.5
cs = plt.contour(PP,TT,eos[:,:,ent],[minent])
v = cs.collections[0].get_paths()[0].vertices #gets contour path
x , y = v[:,0] , v[:,1]
plt.fill_between(x , y, TT.min() , color = 'white')
plt.draw()
plt.tight_layout()
fig.subplots_adjust(top = 0.90 , right = 0.74)
#savefig(userpath+'figs/EOS/delad_S.pdf')

"""
internal energy
"""

arr = eos[: , : , u]
fig = figure(figsize = (5,3.5))
ax = fig.add_subplot(111)
ext = [min(pbar),max(pbar),min(tk),max(tk)]
ax.set_yscale('log')
ax.set_xscale('log')
im = ax.imshow(arr, extent = ext, vmin = 8.5, interpolation = 'nearest', origin = 'lower')
cbar = fig.colorbar(im)
cbar.set_label('log(U)')
ax.set_aspect('auto')
ax.set_ylabel('T [K]')
ax.set_xlabel('P [bar]')

cont = ax.contour(PP, TT, arr, colors = 'k', levels = np.arange(8.0,13.01,0.5))
plt.clabel(cont, inline = 1, fmt='%1.1f')

plt.draw()
plt.tight_layout()
fig.subplots_adjust(top = 0.90 , right = 0.74)
#savefig(userpath+'figs/EOS/U.pdf')


"""
Plot Delad
"""

arr = eos[: , : , delad]
fig = figure(figsize = (5,3.5))
ax = fig.add_subplot(111)
ext = [min(pbar),max(pbar),min(tk),max(tk)]
ax.set_yscale('log')
ax.set_xscale('log')
im = ax.imshow(arr, extent = ext, vmin = np.min(arr[arr > 0]), interpolation = 'nearest', origin = 'lower')
cbar = fig.colorbar(im)
cbar.set_label(r'$\nabla_\mathrm{ad}$')
ax.contour(PP, TT, arr,[2./7, 2./5])
ax.set_ylabel('T [K]')
ax.set_xlabel('P [bar]')
ax.set_aspect('auto')
ax.set_ylabel('T [K]')
ax.set_xlabel('P [bar]')
plt.draw()
plt.tight_layout()
fig.subplots_adjust(top = 0.90 , right = 0.74)
#savefig(userpath+'figs/EOS/delad.pdf')





