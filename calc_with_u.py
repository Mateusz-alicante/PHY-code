from sympy import *
import math

class solver():
    def __init__(self, values):
        self.values = values
        self.symbols = [Symbol(el) for el in list(values.keys())]
        if (type(list(values.values())[0][0]) == int or type(list(values.values())[0][0]) == float):
            self.n_measurments = 1
        else:
            self.n_measurments = len(list(values.values())[0][0])


    def define_equation(self, equation):
        self.eq = equation(*self.symbols)
        self.compute_derivatives()

    def compute_derivatives(self):
        self.derivatives = []
        for el in self.symbols:
            self.derivatives.append(self.eq.diff(el))

    def solve(self):
        self.solved_values = []
        self.solved_errors = []
        for i in range(self.n_measurments):

            # For each measurment, substitute the values and solve the equation
            curr_exp = self.eq
            for el in self.symbols:
                curr_exp = curr_exp.subs(el, self.values[str(el)][0][i])
            self.solved_values.append(curr_exp.evalf())

            # For each measurment compute the error
            curr_err = 0
            for j in range(len(self.symbols)):
                curr_deriv = self.derivatives[j]
                for el in self.symbols:
                    curr_deriv = curr_deriv.subs(el, self.values[str(el)][0][i])
                curr_deriv = curr_deriv.evalf()

                curr_err += ((curr_deriv)*self.values[str(self.symbols[j])][1][i])**2

            self.solved_errors.append(math.sqrt(curr_err))
                
    def print_solution(self):
        for i in range(self.n_measurments):
            print(f"measurment {i}: {self.solved_values[i]} +- {self.solved_errors[i]}")

v = [6.443, 6.451, 6.475, 6.479]
v_err = [0.0005, 0.0005, 0.0005, 0.0005]

i = [0.06338, 0.029344, 0.013916, 0.002399]
i_err = [0.000005, 0.0000005, 0.0000005, 0.0000005]

r = [100.03, 218.89, 463.6, 2696.5]
r_err = [0.005, 0.005, 0.05, 0.05]

s = solver({"v": [v, v_err], "i": [i, i_err], "r": [r, r_err]})

def eq(v, i, r):
    return (v-i*r)/i

s.define_equation(eq)

s.solve()

s.print_solution()


