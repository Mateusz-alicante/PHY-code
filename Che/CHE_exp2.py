from scipy import integrate, stats
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import csv
import statistics
from uncertainties import ufloat
import math

def get_denisty(pressure, temperature):
    R_air = 287
    return pressure / (R_air * temperature)


class exp():
    def __init__(self, heat_filename, mass_filename, custom_errors = None):
        Time = []
        mass_Time = []
        T = []
        P = []
        MassFlow = []
        P_mass_flow = []
        T_mass_flow = []
        EnergyTransfer = []
    
        with open(heat_filename) as f:
            reader = csv.reader(f, delimiter="\t")
            d = list(reader)
            for i in range(4, len(d)):
                Time.append(float(d[i][0]))
                # Add temperatures in K
                T.append(float(d[i][1]) + 273)
                # Add pressures in PSI
                P.append(float(d[i][5]) + 14.75034)

                # Get compund energy transfer in J
                EnergyTransfer.append(float(d[i][8]) * 1000)

            with open(mass_filename) as f:
                reader = csv.reader(f, delimiter="\t")
                d = list(reader)
                for i in range(4, len(d)):

                    mass_Time.append(float(d[i][0]))
                    # Get mass transfer per second (g / sec)
                    MassFlow.append(float(d[i][7]) / 60)

                    T_mass_flow.append(float(d[i][1]) + 273)
                    P_mass_flow.append(float(d[i][5]) + 14.75034)

        
        # Get standart deviation from the first 2 seconds
        if custom_errors:
            self.errors = custom_errors
        else: 
            self.errors = {"T": statistics.stdev(T[:200]), "P": statistics.stdev(P[:200]), "MassFlow": statistics.stdev(MassFlow[:200])}

        print("Errors: ", self.errors)

        self.data = {"Time": Time, "T": T, "mass_Time": mass_Time, "P": P, "MassFlow": MassFlow, "T_mass": T_mass_flow, "P_mass": P_mass_flow, "EnergyTransfer": EnergyTransfer}

    def get_fan_power(self):

        # get real and manufacturer gas denisty
        m_T = 25 + 273
        m_P = 14.75034
        r_T = ufloat(self.data["T"][0], self.errors["T"])
        r_P = ufloat(self.data["P"][0], self.errors["P"])

        # real and manufacturer fan speeds in rpm
        m_speed = 4200
        r_speed = 2000

        m_power = 0.001 * 745.7

        p = m_power * get_denisty(r_P, r_T) / get_denisty(m_P, m_T) * (r_speed / m_speed) ** 3

        print("Power: ", p)

    def calc_mass_flow(self, digitUncertainty = 0.05, digitUncertaintyMass = 0.005):
        mass_entered_value = integrate.simps(self.data["MassFlow"], self.data["mass_Time"])
        # Get the unceertainy value for the mass entered in the system by approximating the calculation as a single rectangle

        mean_flow_rate = statistics.mean(self.data["MassFlow"])

        start_flow_index = None
        end_flow_index = None
        for i, v in enumerate(self.data["MassFlow"]):
            if v > mean_flow_rate and start_flow_index == None:
                start_flow_index = i
            if v > mean_flow_rate:
                end_flow_index = i

        print("flow index: ", start_flow_index, end_flow_index)


        mass_entered_error = math.sqrt((statistics.mean(self.data["MassFlow"][start_flow_index:end_flow_index]) * max(self.errors["MassFlow"] * 60, digitUncertaintyMass * 60))**2 + ((end_flow_index - start_flow_index) * 0.1 * digitUncertainty)**2)
        self.mass_entered = ufloat(mass_entered_value, mass_entered_error)

        final_temp = ufloat(self.data["T_mass"][-1], self.errors["T"])
        initial_temp = ufloat(self.data["T_mass"][0], self.errors["T"])

        final_p = ufloat(self.data["P_mass"][-1], self.errors["P"])
        initial_p = ufloat(self.data["P_mass"][0], self.errors["P"])

        self.mass = self.mass_entered * (1  + 1/((final_p*initial_temp)/(initial_p*final_temp) - 1))

        print("mass: ", self.mass)

    def calc_heat_loss(self, title = "Change title", subtitle = None, start = 150, end = 320):
        x = self.data["Time"]
        y = self.data["EnergyTransfer"]

        fit = stats.linregress(x[start*10:end*10], y[start*10:end*10])

        # plt.rcParams["figure.figsize"] = (10, 5)
        plt.figure(constrained_layout=True)

        plt.scatter(x, y, s=2, alpha=0.5, color = "r", label='Measurments')
        plt.axline(xy1=(0, fit.intercept), slope= fit.slope, label=f'$y = ({fit.slope:.1f} +/- {fit.stderr:.1e})x + ({fit.intercept:.1f}+/-{fit.intercept_stderr:.1e})$')
        plt.legend()

        plt.suptitle(title)
        plt.title(subtitle)

        plt.xlabel("Time [s]")
        plt.ylabel("Cummulative heat transfer [J]")

        plt.savefig(f"/Users/mateuszkazimierczak/Downloads/{title}_{subtitle}.png", dpi=300)

        self.heat_loss_rate = ufloat(fit.slope, fit.stderr)

    def calc_cv(self):
        # Get temp values for heat series
        final_temp = ufloat(self.data["T"][-1], self.errors["T"])
        initial_temp = ufloat(self.data["T"][0], self.errors["T"])

        self.cv = self.heat_loss_rate / self.mass / (final_temp - initial_temp)

        print("cv: ", self.cv)


a = exp("/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - part 2a", "/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - part 1a")

a.calc_heat_loss("Heat loss rate", "Part a")

a.get_fan_power()

a.calc_mass_flow()

a.calc_cv()



b = exp("/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - part 2b", "/Users/mateuszkazimierczak/Downloads/lab-2/lab 2- part 1b")

b.calc_heat_loss("Heat loss rate", "Part b")

b.get_fan_power()

b.calc_mass_flow()

b.calc_cv()


c = exp("/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - part2c", "/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - 1c")

c.calc_heat_loss("Heat loss rate", "Part c")

c.get_fan_power()

c.calc_mass_flow()

c.calc_cv()


d = exp("/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - part 2d", "/Users/mateuszkazimierczak/Downloads/lab-2/lab 2 - part 1d")

d.calc_heat_loss("Heat loss rate", "Part d")

d.get_fan_power()

d.calc_mass_flow()


d.calc_cv()