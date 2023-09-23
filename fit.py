import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *
import scipy.stats as stats

def chi2(observed_data, expected_data):
    chi_square_test_statistic1 = 0
    for i in range(len(observed_data)):
        chi_square_test_statistic1 = chi_square_test_statistic1 + \
            (np.square(observed_data[i]-expected_data[i]))/expected_data[i]
        
    return chi_square_test_statistic1, stats.chi2.ppf(1-0.05, df=6)

class fit():
    def __init__(self, x, x_err, y, y_err):

        # https://stackoverflow.com/questions/22670057/linear-fitting-in-python-with-uncertainty-in-both-x-and-y-coordinates

        self.y = np.array(y)
        self.x = np.array(x)

        def lin_func(p, x):
            m, c = p
            return m*x + c
        
        lin_model = Model(lin_func)

        data = RealData(x, y, sx=x_err, sy=y_err)

        odr = ODR(data, lin_model, beta0=[0., 1.])

        # Run the regression.
        odr.run()

        print(odr.output.beta[0])

        self.m = odr.output.beta[0]
        self.c = odr.output.beta[1]
        self.m_err = odr.output.sd_beta[0]
        self.c_err = odr.output.sd_beta[1]

    def print_equation(self, precision=None, units=""):
        '''Print the equation of the best fit line with uncertainties.'''
        if (precision == None):
            print(f"y = ({self.m} +- {self.m_err}) x + ({self.c} +- {self.c_err}) {units}")
        else:
            print(f"y = ({self.m:.{precision}f} +- {self.m_err:.{precision}f}) x + ({self.c:.{precision}f} +- {self.c_err:.{precision}f}) {units}")

    def preditc(self, x):
        '''Predict y value from x value.'''
        return self.m*x + self.c
    
    def plot(self, title, x_label, y_label, precision=3, units="", fit_font_size=12):
        '''Plot the data and the best fit line.'''
        plt.errorbar(x, y, xerr=x_err, yerr=y_err, ls = "None", color = "r", )
        plt.scatter(x, y, s = 5, marker = "h", color = "r", label='Experimental data')
        x_r = np.linspace(min(x), max(x), 200)
        y_r = [self.preditc(i) for i in x_r]
        plt.plot(x_r, y_r, label='Fit')
        plt.legend()
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': fit_font_size,
        }

        x_text_pos = min(x_r) + (max(x_r) - min(x_r)) * 0.0
        y_text_pos = min(y_r) + (max(y_r) - min(y_r)) * 0.0


        plt.text(x_text_pos, y_text_pos, f"y = ({self.m:.{precision}f} +- {self.m_err:.{precision}f}) x + ({self.c:.{precision}f} +- {self.c_err:.{precision}f}) {units}", fontdict=font)

        plt.show()

    def getIndicators(self):
        # R^2
        residuals = self.y- self.preditc(self.x)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((self.y-np.mean(self.y))**2)
        r_squared = 1 - (ss_res / ss_tot)

        

        # R^2, [chi pi value, chi critical value]
        return r_squared, chi2(self.y, self.preditc(self.x))
        




x = [6.443, 6.451, 6.475, 6.479]
x_err = [0.0005, 0.0005, 0.0005, 0.0005]

y = [0.06338, 0.029344, 0.013916, 0.002399]
y_err = [0.000005, 0.0000005, 0.0000005, 0.0000005]

res = fit(x, x_err, y, y_err)

print(res.getIndicators())

res.print_equation(3)
res.plot("Fit and experimental data", "Voltage (V)", "Current (A)", 3)

