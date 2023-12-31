import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *
import scipy.stats as stats


def chi2_online(observed_data, expected_data):
    # This is from stakc overflow
    chi_square_test_statistic1 = 0
    for i in range(len(observed_data)):
        chi_square_test_statistic1 = chi_square_test_statistic1 + \
            (np.square(observed_data[i]-expected_data[i]))/expected_data[i]

    return chi_square_test_statistic1, stats.chi2.ppf(1-0.05, df=6)


def chi2(observed_data, observed_data_error, expected_data):

    print("Observed: ", observed_data)
    print("Error: ", observed_data_error)
    print("Expected: ", expected_data)
    # This is from the handout
    chi_square = 0
    for i in range(len(observed_data)):
        chi_square += (expected_data[i] - observed_data[i]
                       )**2 / observed_data_error[i]**2

    m = 2  # Two parameters (m and c)
    N = len(observed_data)

    return chi_square, chi_square/(N-m)


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

        self.x_err = x_err
        self.y_err = y_err

    def print_equation(self, precision=None, units=""):
        '''Print the equation of the best fit line with uncertainties.'''
        if (precision == None):
            print(
                f"y = ({self.m} +- {self.m_err}) x + ({self.c} +- {self.c_err}) {units}")
        else:
            print(
                f"y = ({self.m:.{precision}f} +- {self.m_err:.{precision}f}) x + ({self.c:.{precision}f} +- {self.c_err:.{precision}f}) {units}")

    def preditc(self, x):
        '''Predict y value from x value.'''
        return self.m*x + self.c

    def plot(self, title, x_label, y_label, precision=3, units="", fit_font_size=12):
        '''Plot the data and the best fit line.'''

        plt.errorbar(self.x, self.y, xerr=self.x_err,
                     yerr=self.y_err, ls="None", color="r", )

        plt.scatter(self.x, self.y, s=5, marker="h",
                    color="r", label='Experimental data')

        x_r = np.linspace(min(x), max(x), 200)
        y_r = [self.preditc(i) for i in x_r]
        plt.plot(x_r, y_r, label='Fit')
        plt.legend()
        plt.title(title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        plt.ylim(0, 10)

        font = {'family': 'serif',
                'color':  'darkred',
                'weight': 'normal',
                'size': fit_font_size,
                }

        x_text_pos = min(x_r) + (max(x_r) - min(x_r)) * 0.0
        y_text_pos = min(y_r) + (max(y_r) - min(y_r)) * 0.0

        plt.text(x_text_pos, y_text_pos,
                 f"y = ({self.m:.{precision}f} +- {self.m_err:.{precision}f}) x + ({self.c:.{precision}f} +- {self.c_err:.{precision}f}) {units}", fontdict=font)

        # plt.show()

    def getIndicators(self):
        # R^2
        residuals = self.y - self.preditc(self.x)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((self.y-np.mean(self.y))**2)
        r_squared = 1 - (ss_res / ss_tot)

        # R^2, [chi pi value, chi critical value]
        return r_squared, chi2(self.preditc(self.x), self.y_err, self.y)


# First experiment:
y = [6.443, 6.451, 6.475, 6.479]
y_err = [0.0005, 0.0005, 0.0005, 0.0005]

x = [0.06338, 0.029344, 0.013916, 0.002399]
x_err = [0.000005, 0.0000005, 0.0000005, 0.0000005]

res = fit(x, x_err, y, y_err)

print(res.getIndicators())

res.print_equation(3)
res.plot("Fit and experimental data (first experiment)",
         "Current (A)", "Voltage (V)", 3, "V")


# Second experiment:

y = [6.399, 6.407, 6.437, 6.469]
y_err = [0.0005, 0.0005, 0.0005, 0.0005]

x = [0.06338, 0.029266, 0.013883, 0.002397]
x_err = [0.000005, 0.0000005, 0.0000005, 0.0000005]

res = fit(x, x_err, y, y_err)

print(res.getIndicators())

res.print_equation(3)
res.plot("Fit and experimental data (second experiment)",
         "Current (A)", "Voltage (V)", 3, "V")
