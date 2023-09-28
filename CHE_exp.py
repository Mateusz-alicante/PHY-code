from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import csv
import statistics
from uncertainties import ufloat
import math

class exp():
    def __init__(self, filename, custom_errors = None):
        Time = []
        T1 = []
        T2 = []
        P1 = []
        P2 = []
        MassFlow = []
    
        with open(filename) as f:
            reader = csv.reader(f, delimiter="\t")
            d = list(reader)
            for i in range(4, len(d)):
                Time.append(float(d[i][0]))
                T1.append(float(d[i][1]) + 273)
                T2.append(float(d[i][2]) + 273)
                P1.append(float(d[i][3]) + 14.75034)
                P2.append(float(d[i][4]) + 14.75034)

                # Get mass transfer per second
                MassFlow.append(float(d[i][5]) / 60)
        
        # Get standart deviation from the first 2 seconds
        if custom_errors:
            self.errors = custom_errors
        else: 
            self.errors = {"T1": statistics.stdev(T1[:200]), "T2": statistics.stdev(T2[:200]), "P1": statistics.stdev(P1[:200]), "P2": statistics.stdev(P2[:200]), "MassFlow": statistics.stdev(MassFlow[:200])}

        self.data = {"Time": Time, "T1": T1, "T2": T2, "P1": P1, "P2": P2, "MassFlow": MassFlow}

    def calc_mass_flow(self, start = 0, end = 0, digitUncertainty = 0.05, digitUncertaintyMass = 0.005):
        mass_entered_value = integrate.simps(self.data["MassFlow"], self.data["Time"])
        mass_entered_error = math.sqrt((statistics.mean(self.data["MassFlow"][start:end]) * max(self.errors["MassFlow"] * 60, digitUncertaintyMass * 60))**2 + ((end - start) / 10 * digitUncertainty)**2)
        self.mass_entered = ufloat(mass_entered_value, mass_entered_error)
        print("Mass entered: ", self.mass_entered)

    def calc_vol(self, pressure_sampling_time = 100, temperature_sampling_time = 100, pressure_sampling_start = 0, temperature_sampling_start = 0):

        if not self.mass_entered:
            self.calc_mass_flow()

        Pi = ufloat(statistics.mean(self.data["P2"][pressure_sampling_start:(pressure_sampling_start + pressure_sampling_time)]), self.errors["P2"])
        Pf = ufloat(statistics.mean(self.data["P2"][-pressure_sampling_time:]), self.errors["P2"])
        Ti = ufloat(statistics.mean(self.data["T2"][temperature_sampling_start:(temperature_sampling_start + temperature_sampling_time)]), self.errors["T2"])
        Tf = ufloat(statistics.mean(self.data["T2"][-temperature_sampling_time:]), self.errors["T2"])

        print("Pi: ", Pi, "Pf: ", Pf, "Ti: ", Ti, "Tf: ", Tf, "mass entered: ", self.mass_entered)
        print("final values:    ", "Pi: ", Pi * 6894.76, "Pf: ", Pf * 6894.76, "Ti: ", Ti, "Tf: ", Tf, "mass entered: ", self.mass_entered / 1000)
        R_air = 287

        self.vol = (self.mass_entered / 1000) / ((Pf  * 6894.76 / (Tf * R_air)) - (Pi  * 6894.76 / (Ti * R_air)))
        print("Volue of container: ", self.vol)
        
        
    def plot_temp(self, title = "Modify title", subtitle = None):
        """Plot the temperatures over time."""

        plt.scatter(self.data["Time"], self.data["T1"], s=2, alpha=0.5, color = "r", label='Left tank')
        plt.errorbar(self.data["Time"], self.data["T1"], xerr=0,  alpha=0.5, yerr=self.errors["P1"], ls = "None", color = "salmon", zorder=-1, label='Left tank error')

        plt.scatter(self.data["Time"], self.data["T2"], s=2, alpha=0.5, color = "g", label='Right tank')
        plt.errorbar(self.data["Time"], self.data["T2"], xerr=0,  alpha=0.5, yerr=self.errors["P2"], ls = "None", color = "springgreen",  zorder=-1, label='Right tank error')

        plt.suptitle(title)
        plt.title(subtitle)
        plt.xlabel("Time [s]")
        plt.ylabel("Temperature [K]")

        plt.legend()
        plt.savefig(f"/Users/mateuszkazimierczak/Downloads/{title}_{subtitle}.png", dpi=200)

        plt.clf()

    def plot_p(self, title = "Modify title", subtitle = None):
        
        plt.scatter(self.data["Time"], self.data["P1"], s=2, alpha=0.5, color = "r", label='Left tank')
        plt.errorbar(self.data["Time"], self.data["P1"], xerr=0,  alpha=0.5, yerr=self.errors["P1"], ls = "None", color = "salmon", zorder=-1, label='Left tank error')

        plt.scatter(self.data["Time"], self.data["P2"], s=2, alpha=0.5, color = "g", label='Right tank')
        plt.errorbar(self.data["Time"], self.data["P2"], xerr=0,  alpha=0.5, yerr=self.errors["P2"], ls = "None", color = "springgreen",  zorder=-1, label='right tank error')

        plt.suptitle(title)
        plt.title(subtitle)
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Psi]")

        plt.legend()

        plt.savefig(f"/Users/mateuszkazimierczak/Downloads/{title}_{subtitle}.png", dpi=200)

        plt.clf()

    def plot_mass_flow(self, title = "Modify title", subtitle = None):
            
            plt.scatter(self.data["Time"], self.data["MassFlow"], s=2, alpha=0.5, color = "r", label='Left tank')
            plt.errorbar(self.data["Time"], self.data["MassFlow"], xerr=0,  alpha=0.5, yerr=self.errors["MassFlow"], ls = "None", color = "salmon", zorder=-1, label='Left tank error')
    
            plt.suptitle(title)
            plt.title(subtitle)
            plt.xlabel("Time [s]")
            plt.ylabel("mass flow rate [g/s]")
    
            plt.legend()
            plt.savefig(f"/Users/mateuszkazimierczak/Downloads/{title}_{subtitle}.png", dpi=200)

    def calc_vol_ratio(self, pressure_sampling_time = 100, temperature_sampling_time = 100,
                       pressure_sampling_start = 0, temperature_sampling_start = 0):
        # Compute pressure averages
        Pri = ufloat(statistics.mean(self.data["P2"][pressure_sampling_start:(pressure_sampling_start + pressure_sampling_time)]), self.errors["P2"]) * 6894.76
        Prf = ufloat(statistics.mean(self.data["P2"][-pressure_sampling_time:]), self.errors["P2"]) * 6894.76
        Pli = ufloat(statistics.mean(self.data["P1"][pressure_sampling_start:(pressure_sampling_start + pressure_sampling_time)]), self.errors["P1"]) * 6894.76
        Plf = ufloat(statistics.mean(self.data["P1"][-pressure_sampling_time:]), self.errors["P1"]) * 6894.76

        # Compute temp averages
        Tri = ufloat(statistics.mean(self.data["T1"][temperature_sampling_start:(temperature_sampling_start + temperature_sampling_time)]), self.errors["T1"])
        Trf = ufloat(statistics.mean(self.data["T1"][-temperature_sampling_time:]), self.errors["T1"])
        Tli = ufloat(statistics.mean(self.data["T2"][temperature_sampling_start:(temperature_sampling_start + temperature_sampling_time)]), self.errors["T2"])
        Tlf = ufloat(statistics.mean(self.data["T2"][-temperature_sampling_start:]), self.errors["T2"])

        print("Pri: ", Pri, "Prf: ", Prf, "Pli: ", Pli, "Plf: ", Plf)
        print("Tri: ", Tri, "Trf: ", Trf, "Tli: ", Tli, "Tlf: ", Tlf)
        self.vol_ratio = ((Prf / Trf) - (Pri / Tri)) / ((Pli / Tli) - (Plf / Tlf))
        print("Vol ratio: ", self.vol_ratio)
       


custom_errors = {'T1': 0.0740076731048306, 'T2': 0.06732319110308066, 'P1': 0.3585954594213208, 'P2': 0.37697553337560147, 'MassFlow': 1.662473787906861e-05}
first = exp("/Users/mateuszkazimierczak/Downloads/lab/Experiment 1 part 1.txt", custom_errors=custom_errors)
print(first.errors)
#first.plot_temp("Temperature over time")

first.calc_vol_ratio(200, 200, 2200, 2300)

first.plot_p("Pressure in both tanks over time", "first experiment first section")
first.plot_temp("Temperature in both tanks over time", "first experiment first section")

second = exp("/Users/mateuszkazimierczak/Downloads/lab/Experiment 1 part 2.txt", custom_errors=custom_errors)

second.calc_vol_ratio(200, 200, 1800, 2000)

second.plot_p("Pressure in both tanks over time", "first experiment second section")
second.plot_temp("Temperature in both tanks over time", "first experiment second section")

thrid = exp("/Users/mateuszkazimierczak/Downloads/lab/Experiment 2 v2.txt", custom_errors=custom_errors)

thrid.calc_mass_flow(20, 55)
thrid.calc_vol(100, 100, 0, 0)


thrid.plot_p("Pressure in both tanks over time", "second experiment")
thrid.plot_temp("Temperature in both tanks over time", "second experiment")
thrid.plot_mass_flow("Mass flow rate over time", "second experiment")
