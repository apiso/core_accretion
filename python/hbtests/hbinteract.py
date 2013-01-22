"""interactive script to implement code in module hbtests.py"""

import hbint as hb

#sample integration of adiabats
x , y , K = hb.goout(10**6 , .01 , 4000) #params are Pc, Pout & Tcore (cgs), w/ other optional keywords
x , y , K = hb.gooutrho(3. , 1e-10 , 1.) #params are rhoc, rhoout (cgs) & Tcore as fraction of virial estimate, w/ optional keywords

#check virial equilibrium
hb.vircheck(x,y)

#physical values
p , r , m = exp(x) , y[:,0] , y[:,1]
rr = r/hb.Re #and earth radii
mm = m/hb.Me #in earth masses

#plots
#...
