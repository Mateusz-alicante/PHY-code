from uncertainties import ufloat

def eq(v, i, r):
    return (v * r)/(i * r - v)

v = ufloat(6.437, 0.0005)
r = ufloat(463.6, 0.05)
i = ufloat(0.013883, 0.0000005)

print(eq(v, i, r))