from uncertainties import ufloat

m = ufloat(-0.620, 0.176)
r = ufloat(-2047234, 1343628)

print((m*r)/(m-r))