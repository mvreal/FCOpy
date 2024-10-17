import numpy as np
from scipy.optimize import fsolve

# Define the system of equations
def equations(vars):
    x, y, z = vars
    eq1 = x**2 + y**2 + z**2 - 10
    eq2 = x*y*z - 1
    eq3 = x + y + z - 6
    return [eq1, eq2, eq3]

# Initial guesses for the variables
initial_guess = [1.00, 1.00, 1.00]

# Solve the system of equations
solution = fsolve(equations, initial_guess)

print("Solution:")
print(f"x = {solution[0]:.4f}")
print(f"y = {solution[1]:.4f}")
print(f"z = {solution[2]:.4f}")

results = equations(solution)

print("Results:")
print(f"eq1 = {results[0]:.4f}") 
print(f"eq2 = {results[1]:.4f}") 
print(f"eq3 = {results[2]:.4f}")

