from uncertainties import ufloat, unumpy

n = ufloat(6.40 * 10**(5), 0.37 * 10 ** (5))

print('Result = {:10.3f}'.format(n**2))