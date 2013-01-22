A = 6.02e23 #Avog #
ex = exp(1)
hf = 6.55e-27

mu = 200.6 # molecular weight of mercury
#Ri = A*kb #ideal gas const
Ri = 8.314e7
om = 1 #stat weight

temp = 630
press = 1e6
aa = (2*pi * mu)**1.5 * Ri**2.5 * om * ex**2.5 * hf**(-3) * A**(-4)
ent = Ri*(2.5*log(temp) - log(press) + log(aa))

answer = 191e7

print "agreement with Mercury result " + str(ent/answer)

mHe = 6.646442e-24
temp = 10**2.18
press = 10**4.2
L = h/sqrt(2*pi * mHe * kb * temp)
ent = kb / mHe * (log(kb * temp / press / L**3) + 2.5)
print log10(ent)
