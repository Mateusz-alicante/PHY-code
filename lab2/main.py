import matplotlib.pyplot as plt
import numpy as np
import scipy
from uncertainties import ufloat, unumpy
from scipy.odr import *
import math
from decimal import Decimal
from sklearn.metrics import r2_score


u0 = math.pi * 4e-7


def chi2(observed_data, expected_data_error, expected_data):

    # This is from the handout
    chi_square = 0
    for i in range(len(observed_data)):
        chi_square += (expected_data[i] - observed_data[i]
                       )**2 / expected_data_error[i]**2

    m = 2  # Two parameters (m and c), for linear fit
    N = len(observed_data)

    return chi_square, chi_square/(N-m)


def r2(observed_data, expected_data):
    '''Compute the R^2 value for the fit.'''
    y_observed = unumpy.nominal_values(observed_data)
    y_real = unumpy.nominal_values(expected_data)
    return r2_score(y_real, y_observed)
    # residuals = y_observed - observed_data
    # ss_res = np.sum(residuals**2)
    # ss_tot = np.sum((expected_data_error-np.mean(expected_data_error))**2)
    # r_squared = 1 - (ss_res / ss_tot)
    # return r_squared


def plot_reisduals(x, y, m, c, title="Change title", xlabel="Change x label", ylabel="Change y label"):
    '''Plot the residuals of the fit.'''
    print("residuals: ", ylabel)
    plt.clf()
    plt.figure(constrained_layout=True)
    title = title + " residuals"
    plt.scatter(unumpy.nominal_values(x), unumpy.nominal_values(
        y) - (m * unumpy.nominal_values(x) + c), s=5, marker="h", color="r")
    plt.axhline(y=0, color='k', linestyle='--')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(
        f"/Users/mateuszkazimierczak/Downloads/{title}_residuals.png", dpi=300)


def get_goodness_indicators(observed_data, expected_data_error, expected_data):

    # R^2
    r_squared = r2(observed_data, expected_data)

    # Chi^2
    chi_square, reduced_chi_squared = chi2(
        observed_data, expected_data_error, expected_data)

    print("Chi^2: ", chi_square.nominal_value)
    print("Reduced Chi^2: ", reduced_chi_squared.nominal_value)
    print("R^2: ", r_squared)


def calc_r(dist):
    '''Calculate the radius of the path based on the number of increments to the left on the self-iluminating scale'''
    return (6 * 5 + 20 + 5 * dist)


class exp():
    def __init__(self, bulb_radius, coil_radius, coil_turns, coil_dist):
        self.bulb_radius = bulb_radius
        self.coil_radius = coil_radius
        self.coil_turns = coil_turns
        self.coil_dist = coil_dist

    def set_measurments(self, measurments, uncertainties, connst_values):

        # const I measurments
        self.const_I_value = ufloat(connst_values[0][0], connst_values[0][1])
        # The following has the format: [Voltage value, radius of electron path]
        self.const_I = \
            [[ufloat(v, uncertainties[0][0]), ufloat(calc_r(d), uncertainties[0][1]) / 1000]
             for v, d in measurments[0]]

        self.field_correction_factor_const_I = \
            [((1 - (r**4) / (self.coil_radius ** 4 * (0.6583 + 0.29 * r **
              2 / self.coil_radius ** 2) ** 2), r/self.coil_radius)) for v, r in self.const_I]

        self.const_I_field_c = [(4/5) ** (3/2) * u0 * self.coil_turns *
                                self.const_I_value / self.coil_radius * self.field_correction_factor_const_I[index][0] for index, const_i in enumerate(self.const_I)]

        print(np.array(self.const_I_field_c))
        # const V measurments
        self.const_V_value = ufloat(connst_values[1][0], connst_values[1][1])
        # The following has the format: [current_value, radius of electron path]
        self.const_V = \
            [[ufloat(v, uncertainties[1][0]), ufloat(calc_r(d), uncertainties[1][1]) / 1000]
             for v, d in measurments[1]]

        self.field_correction_factor_const_V = \
            [(1 - (r**4) / (self.coil_radius ** 4 * (0.6583 + 0.29 * r **
              2 / self.coil_radius ** 2) ** 2), r/self.coil_radius) for v, r in self.const_V]
        # For each mesurment, store [correction factor with uncertainty, ratio of electron radius to coil radius]

        self.const_V_field_c = [(4/5) ** (3/2) * u0 * self.coil_turns *
                                const_v[0] / self.coil_radius * self.field_correction_factor_const_V[index][0] for index, const_v in enumerate(self.const_V)]

        print("Constat V fields: ")
        for index, i in enumerate(self.const_V_field_c):
            i = i * self.field_correction_factor_const_V[index][1]
            print(f"{i.nominal_value / 10 ** (-5):.2f} +/- {i.std_dev/ 10 ** (-5):.2f}")
        print("Constat I fields: ")
        for i in self.const_I_field_c:
            i = i * self.field_correction_factor_const_V[index][1]
            print(f"{i.nominal_value / 10 ** (-5):.2f} +/- {i.std_dev/ 10 ** (-5):.2f}")

    def prepare_equation(self, m, c, m_err, c_err, units="", precision=None):
        '''Print the equation of the best fit line with uncertainties.'''
        if (precision == None):
            return (f"y = ({m} +- {m_err}) x + ({c} +- {c_err}) {units}")
        else:
            return (
                f"y = ({Decimal(m):.{precision}E} +- {Decimal(m_err):.{precision}E}) x + ({Decimal(c):.{precision}E} +- {Decimal(c_err):.{precision}E}) {units}")

    def compute_fit(self, x, y):

        o = scipy.stats.linregress(
            unumpy.nominal_values(x), unumpy.nominal_values(y))

        # def lin_func(p, x):
        #     m, c = p
        #     return m*x + c

        # lin_model = Model(lin_func)

        # data = RealData(unumpy.nominal_values(x), unumpy.nominal_values(y),
        #                 sx=unumpy.std_devs(x), sy=unumpy.std_devs(y))

        # odr = ODR(data, lin_model, beta0=[0., 1.])

        # # Run the regression.
        # odr.run()

        # return odr.output.beta[0], odr.output.beta[1], odr.output.sd_beta[0], odr.output.sd_beta[1], x, y

        return o.slope, o.intercept, o.stderr, o.intercept_stderr, x, y

    def predict(self, x, m, c):
        '''Predict y value from x value.'''
        return m*x + c

    def get_bc_r_graph(self, title="Change title", x_label="1/r", y_label="Bc", precision=3, units="T", fit_font_size=10):
        # if (type == "Const I"):
        #     m, c, m_err, c_err, x, y = self.compute_fit(
        #         self.const_I, self.const_I_field_c)
        # elif (type == "Const V"):
        #     m, c, m_err, c_err, x, y = self.compute_fit(
        #         self.const_V, self.const_V_field_c)
        # else:
        #     return print("Invalid type")

        # Get all points together on one graph
        plt.clf()
        m, c, m_err, c_err, x, y = self.compute_fit(
            1 / np.array([i[1] for i in self.const_V + self.const_I]),  np.array(self.const_V_field_c + self.const_I_field_c))

        self.b_e = ufloat(c, c_err)

        print(self.prepare_equation(m, c, m_err,
              c_err, units="T", precision=precision))

        plt.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), xerr=unumpy.std_devs(x),
                     yerr=unumpy.std_devs(y), ls="None", color="r", )
        plt.scatter(unumpy.nominal_values(x), unumpy.nominal_values(y), s=5, marker="h",
                    color="r", label='Experimental data')
        x_r = np.linspace(min(unumpy.nominal_values(x)),
                          max(unumpy.nominal_values(x)), 200)
        y_r = [self.predict(i, m, c) for i in x_r]
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

        x_text_pos = (min(unumpy.nominal_values(x)) +
                      max(unumpy.nominal_values(x))) / 2 * 0.7
        y_text_pos = min(unumpy.nominal_values(y))

        plt.text(x_text_pos, y_text_pos,
                 self.prepare_equation(m, c, m_err, c_err, units="T", precision=precision), fontdict=font)

        plt.savefig(
            f"/Users/mateuszkazimierczak/Downloads/{title}.png", dpi=300)

        plot_reisduals(x, y, m, c, title=title, xlabel=x_label, ylabel=y_label)

        print("Goodness indicators:", title, "=====================")
        get_goodness_indicators([self.predict(i, m, c) for i in x], unumpy.std_devs(
            y), unumpy.nominal_values(y))

    def get_ratio(self, title="Change title", x_label="1/r", y_label="Bc", precision=3, units="T", fit_font_size=10):
        plt.clf()
        k = 1/math.sqrt(2) * (4/5) ** (3/2) * u0 * \
            self.coil_turns / self.coil_radius
        I0 = self.b_e / k
        r = np.array([i[1] for i in self.const_I + self.const_V])
        Ivalues = np.array([i[0] for i in self.const_V] +
                           [self.const_I_value] * len(self.const_I))
        Vvalues = np.array([self.const_V_value] *
                           len(self.const_V) + [i[0] for i in self.const_I])
        x_raw = k * (Ivalues + 1/math.sqrt(2) * I0) / Vvalues
        y_raw = 1 / r
        m, c, m_err, c_err, x, y = self.compute_fit(x_raw, y_raw)

        print(self.prepare_equation(m, c, m_err,
              c_err, units="T", precision=precision))

        plt.errorbar(unumpy.nominal_values(x), unumpy.nominal_values(y), xerr=unumpy.std_devs(x),
                     yerr=unumpy.std_devs(y), ls="None", color="r", )
        plt.scatter(unumpy.nominal_values(x), unumpy.nominal_values(y), s=5, marker="h",
                    color="r", label='Experimental data')
        x_r = np.linspace(min(unumpy.nominal_values(x)),
                          max(unumpy.nominal_values(x)), 200)
        y_r = [self.predict(i, m, c) for i in x_r]
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

        x_text_pos = (min(unumpy.nominal_values(x)) +
                      max(unumpy.nominal_values(x))) / 2 * 0.45
        y_text_pos = min(unumpy.nominal_values(y))

        plt.text(x_text_pos, y_text_pos,
                 self.prepare_equation(m, c, m_err, c_err, units=units, precision=precision), fontdict=font)

        plt.savefig(
            f"/Users/mateuszkazimierczak/Downloads/{title}.png", dpi=300)

        plot_reisduals(x, y, m, c, title=title, xlabel=x_label, ylabel=y_label)

        print("Goodness indicators:", title, "=====================")
        get_goodness_indicators([self.predict(i, m, c) for i in x_raw], unumpy.std_devs(
            y_raw), unumpy.nominal_values(y_raw))


e = exp(
    ufloat(8.9, 0.05) / 1000,  # bulb radius in mm
    # coil radius (diameter of one coil from coil ends minus width of coil)
    ufloat(17.2 - 2.5, 0.05) / 100,
    ufloat(130, 0),  # coil turns
    # coil distance (distance from coil ends minus width of coil)
    ufloat(17.2 - 1.3, 0.05) / 100
)

e.set_measurments([
    [[113, 3.5], [119.4, 4], [127.8, 5], [142.7, 6], [155.7, 7], [
        173.5, 8], [191.1, 9], [210.7, 10]],  # measurments for const I

    [[1.850, 1], [1.735, 2], [1.575, 3], [1.457, 4], [1.349, 5], [1.269, 6], [
        1.206, 7], [1.184, 8], [1.097, 9], [0.960, 10]]  # measurments for const V
],
    # uncertainties, same for const I and const V
    [[0.001, 1], [0.001, 1]],

    [[1.200000, 0.001], [153, 0.1]]
    # constant values for const I and const V with uncertainties
)

# e.get_bc_r_graph(title="Corrected field value vs inverse radius", x_label="1/r (1/m)",
#                  y_label="corrected field value (T)", precision=3, units="T", fit_font_size=10)

# e.get_ratio(title="charge to mass ratio", y_label="1/r (1/m)",
#             x_label="k * (I + 1/sqrt(2)*I0)/V (1/sqrt(m*c))", precision=3, units="1/m",)
