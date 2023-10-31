from sympy import *
import math
import statistics


class solver():
    def __init__(self, values):
        self.values = values
        self.symbols = [Symbol(el) for el in list(values.keys())]

    def define_equation(self, equation):
        self.eq = equation(*self.symbols)

    def get_steps(self, precision, func_name="f", units=""):
        """Return the steps to solve the equation."""
        res = ""
        res += f"We start with the eqation: \n $${func_name} = {latex(self.eq)}$$ \n \n"
        res += "First, we compute the value of the function at the values of the variables: \n"
        for symbol in self.symbols:
            res += f"$${symbol} = {self.values[str(symbol)][0]:.{precision}f}$$ \n"
        res += "\n \n"
        res += "Substituting the values: \n"
        final_func_val = self.eq.subs({key: self.values[str(key)][0] for key in self.symbols})
        res += f"$$ {func_name} = {final_func_val:.{precision}f}$$"
        res += "\n \n"
        res += "Now, we compute the derivative of the function with respect to each variable: \n \n"
        single_u_values = []
        for symbol in self.symbols:
            derivative = self.eq.diff(symbol)
            var_value = self.values[str(symbol)][1]
            res += f"The derivative of the function with respect to ${symbol}$ is : \n \
            ${latex(derivative)}$. \n"
            derivative_value = derivative.subs({key: self.values[str(key)][0] for key in self.symbols})
            single_u_values.append(derivative_value * var_value)
            res += f"Substituting the values: \n \
                $$({latex(derivative)}) \cdot \Delta {latex(symbol)} = {derivative_value:.{precision}f} \cdot {var_value:.{precision}f} = {(derivative_value * var_value):.{precision}f}$$ \n"
            res +=  "\n \n "
            
        final_res = math.sqrt(sum([el**2 for el in single_u_values]))

        res += "Computing the final uncertainty: \n $\sqrt{"
        for i, single_u_value in enumerate(single_u_values):
            res += f"({single_u_value:.{precision}f})^2"
            if i != len(single_u_values) - 1:
                res += " + "
        res += "}"
        res += f" = {final_res:.{precision}f}$ \n"

        res += "\n \n"

        res += "The final result is: \n"
        res += f"$$ {func_name} = {final_func_val:.{precision}f} \pm {final_res:.{precision}f} {units}$$"

        return res, final_res
            

# First part
# def eq(m, r):
#     return (m*r)/(r-m)

# s = solver({"m_1": [0.620, 0.176], "R_v": [-2047234, 1343628]})
# s.define_equation(eq)
# print(s.get_steps(3, "R_{01}", "Ohm")[0])

# Second part
def eq(m, r):
    return m-r

s = solver({"m_2": [1.059, 0.403], "R_a": [2.12, 0.2]})
s.define_equation(eq)
print(s.get_steps(3, "R_{02}", "Ohm")[0])

# val = {"v": [6.399, 0.0005], "i": [0.063381, 0.000005], "r": [100.03, 0.005]}
# sym = [Symbol(el) for el in list(val.keys())]
# eq = (sym[0] * sym[2])/(sym[1] * sym[2] - sym[0])

# print(eq.diff(sym[0]))
# print({key: val[str(key)][0] for key in sym})
# print(eq.diff(sym[1]).subs({key: val[str(key)][0] for key in sym}))