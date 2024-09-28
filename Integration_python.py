import numpy as np
import matplotlib.pyplot as plt
from scipy import constants, integrate


def trapezoidal_integration(func, a, b, n):
        h = (b - a) / n
    # Initialize the integral with the first and last terms
        integral = 0.5 * (func(a) + func(b))
    # Sum the middle terms
        for i in range(1, n):
            integral += func(a + i * h)
    # Multiply by the width of the trapezoids
        integral *= h
        return integral

def simpson_integration(func, a, b, n):
    if n % 2 == 1:
        raise ValueError("Number of intervals n must be even.")
        
    h = (b - a) / n
    integral = func(a) + func(b)

    for i in range(1, n, 2):
        integral += 4 * func(a + i * h)  # Odd indices are multiplied by 4

    for i in range(2, n-1, 2):
        integral += 2 * func(a + i * h)  # Even indices are multiplied by 2
    
    integral *= h / 3  # Final multiplication by h/3
    return integral

def riemann_integration(func, a, b, n):

    h = (b - a) / n
    integral = 0.0
    
    for i in range(n):
        integral += func(a + i * h)  # Sum the left endpoint values
    
    integral *= h  # Multiply by the width of the intervals
    return integral

def electric_field(r, R, Q): 
        
        epsilon_0 = constants.epsilon_0
        if r < R:        
               return (Q / (4 * np.pi * epsilon_0 * R**3)) * r 
        else:
                return (Q / (4 * np.pi * epsilon_0 * r**2)) 
def electric_potential_trapezoidal(electric_field, r,R,Q,n):   
        inf_approx = 1000*R          
        def field_to_integrate(x):
                return electric_field(x, R, Q) 
        if r <= R: 
                potential = trapezoidal_integration(field_to_integrate, inf_approx, R, n)     
        else:         
                potential = trapezoidal_integration(field_to_integrate, inf_approx, r,n)          
        return -potential

def electric_potential_riemann(electric_field,r,R,Q,n):     
        inf_approx = 1000*R         
        def field_to_integrate(x):         
                return electric_field(x, R, Q)          
        if r <= R:         
                potential = riemann_integration(field_to_integrate, inf_approx, R, n)     
        else:         
                potential = riemann_integration(field_to_integrate, inf_approx, r,n)          
        return -potential

def electric_potential_simpson(electric_field,r,R,Q,n):     
        inf_approx = 1000*R          
        def field_to_integrate(x):         
                return electric_field(x, R, Q)          
        if r <= R:         
                potential = simpson_integration(field_to_integrate, inf_approx, R, n)     
        else:         
                potential = simpson_integration(field_to_integrate, inf_approx, r,n)          
        return -potential  
def analytical_potential(r, R, Q):
        epsilon_0 = constants.epsilon_0 
        if r < R:
                return (Q / (4 * np.pi * epsilon_0 * R))
        else:         
                return Q / (4 * np.pi * epsilon_0 * r)
        
def electric_potential_from_field(electric_field, r, R, Q): 
    inf_approx = 1000 * R  # Approximation for infinity          
    def field_to_integrate(r):         
        return electric_field(r, R, Q)
    
    if r <= R:        
        potential, _ = integrate.quad(field_to_integrate, inf_approx,R)  
    else:      
        potential, _ = integrate.quad(field_to_integrate, inf_approx, r)  
        
    return -potential


def calculate_relative_error(numerical, analytical):
    return np.abs((numerical - analytical) / analytical)

R = 1.0  # Radius of the sphere in meters
Q = 1e-9  # Total charge in Coulombs
n_trapezoidal = 1000  # Number of intervals for trapezoidal integration
n_intervals = 1000

# Calculate and plot
r_values = np.linspace(0.01 * R, 3 * R, 1000)  # Start from 0.01*R to avoid division by zero
E_values = [electric_field(r, R, Q) for r in r_values]
V_values_numerical = [electric_potential_from_field(electric_field, r, R, Q) for r in r_values]
V_values_trapezoidal = [electric_potential_trapezoidal(electric_field, r, R, Q, n_trapezoidal) for r in r_values]
V_values_simpson = [electric_potential_simpson(electric_field, r, R, Q, n_trapezoidal) for r in r_values]
V_values_riemann = [electric_potential_riemann(electric_field, r, R, Q, n_trapezoidal) for r in r_values]
V_values_analytical = [analytical_potential(r, R, Q) for r in r_values]



# Plot the electric field and potential
plt.figure(figsize=(18, 6))

plt.plot()
plt.plot(r_values, E_values)
plt.title('Electric Field')
plt.xlabel('Distance from center (m)')
plt.ylabel('Electric Field (V/m)')
plt.axvline(x=R, color='r', linestyle='--', label='Sphere surface')
plt.legend()
plt.savefig("electric_field")

plt.figure()
plt.plot
plt.plot(r_values, V_values_numerical, label='Numerical (quad)')
plt.plot(r_values, V_values_trapezoidal, label='Numerical (trapezoidal)', linestyle=':')
plt.plot(r_values, V_values_simpson, label='Numerical (simpson)', linestyle=':')
plt.plot(r_values, V_values_riemann, label='Numerical (riemann)', linestyle=':')
plt.plot(r_values, V_values_analytical, label='Analytical', linestyle='--')
plt.title('Electric Potential')
plt.xlabel('Distance from center (m)')
plt.ylabel('Electric Potential (V)')
plt.axvline(x=R, color='r', linestyle='--', label='Sphere surface')
plt.legend()
plt.savefig("Potential.png")

plt.tight_layout()
plt.show()

