# Import the libraries needed
import numpy as np
import matplotlib.pyplot as plt


# Define physical constant
m = 9.10938356 * 10**(-31)  # Mass of the electron in kilograms
hbar =  1.0545718 * 10**(-34) # Angular Planck's constant in Js
e = 1.60217662 * 10**(-19)   # Charge of electron in Coulomb


def V(x):
	"""Finite square well potential function
	Input : x : position for x
	Output: Vx: Potential at position x"""

	##########################################################################
	### Because I want to plot the potential with input as an array          #
	### Python gives error when doing "if" to an array                       #
	### so I added this grand "if" to distinguish an array or a single number#
	##########################################################################

	### Different Output type (single float or array) depends on input type ###
	if np.size(x) == 1:
		### For a single number ###
		if np.abs(x) >= a/2:
			return V_0
		else:
			return 0
	else:          # For an array, comparing every element
		V = np.zeros(np.size(x))    # An zeroes array with the same size with x
		for i in range(np.size(x)):
			if np.abs(x[i]) >= a/2:  # Using absolute value for both sides
				V[i] = V_0          # Finite potential at two sides
			else:
				V[i] = 0
		return V
	print("Error, np.size(x)=", np.size(x))      # Debug if error
	return 0

def Schro_eq(r,x,E):
            """ One dimension Time independent Schrodinger Equation
            Input: r : Wavefunction corresponds to x position
                   x : Position of electron
                   e : Wavefunction energy           
            Output: psi_d,phi_d: First and second order derivative of the wave function"""
            psi=r[0]
            phi= r[1]
            psi_d = phi
            phi_d = 2*m*(V(x)-E)*psi/(hbar**2)
            return np.array([psi_d,phi_d])

            
def RungeKutta2d(r,xpoints,E):
	'''Fourth-order Runge-Kutta for two differential equations
	Inputs: r: initial guess of [wavefunction psi, first derivative phi]
			xpoints: array of x coordinates.
			E : Energy of the wave function
	Outputs: [psi, phi]: solutions for psi(x) and phi(x), numpy arrays one longer than xpoints'''
	psi = [] # initialise empty arrays
	phi = []
	for x in xpoints: # loops over all xpoints 
		psi.append(r[0])
		phi.append(r[1])
		# This RK function solves shcrodinger equation
		k1 = h*Schro_eq(r, x, E)
		k2 = h*Schro_eq(r+0.5*k1, x+0.5*h, E)
		k3 = h*Schro_eq(r+0.5*k2, x+0.5*h, E)
		k4 = h*Schro_eq(r+k3, x+h, E)
		r = r + (k1 + 2*k2 + 2*k3 + k4)/6
	# these next two lines calculate for the point at x = x_end (Adding up the last element)
	psi.append(r[0])
	phi.append(r[1])
	return np.array([psi, phi]) # convert output to numpy array with 2 rows and N+1 columns            



### Setting up position array ###
a = 10*10**(-11)       # half-width of the potential well
x_start = -a           # start of the well
x_end = a              # end of the well
N = 10000              # number of points for Runge-Kutta
h = (x_end - x_start)/N # size of Runge-Kutta steps
xpoints = np.arange(x_start, x_end, h)       # x position array of the well 
	
def secand(E1,E2):
		""" Find specific wavefunction solution using 2 initial guesses for each energy
		and normalise the result wavefunction
		Input: E1, E2: 2 intial random guess for the energy
		Output: psi :   wavefunction array (normalised)"""
		#Initial guess
		phi = 1
		wave1 = RungeKutta2d(np.array([0,phi]),xpoints,E1)[0,N]      # Solve equation with initial guess and extrate last wavefuction(psi) component
		wave2 = RungeKutta2d(np.array([0,phi]),xpoints,E2)[0,N]
		tolerance = e/1000                 # set the tolerance
		err = 1                          # set initial error 
		while err > tolerance:
			E3 = E2 - wave2*(E2-E1)/(wave2-wave1)
			err = abs(E2-E1) 
			# reset initial phi for the next iteration
			E1 = E2 
			E2 = E3 
			# obtain wavefunction at the end
			wave1 = RungeKutta2d(np.array([0, phi]),xpoints,E1)[0,N]
			wave2 = RungeKutta2d(np.array([0, phi]),xpoints,E2)[0,N]
		psi = RungeKutta2d(np.array([0, phi]),xpoints,E1)[0,]
		### Normalising ###
		I = h*(0.5*psi[N]**2+0.5*psi[0]**2+np.sum(psi[1:N-1]**2))   # Using trapezium rule to integrate wavefunction
		psi_n = psi/np.sqrt(I)        # Normalize the wavefunction
		print("Energy is",E1/e,"eV")  # Print energy 
		return psi_n     # return normalised wave function array
	
# Finite potential value 
V_0 = 50*e      # The step depth of finite potential is 500 eV
	
# Ploat a graph of potential 
plt.figure()
plt.plot(xpoints,V(xpoints)/e,label="Potential")
plt.ylabel("Potential Value (in eV)")
plt.xlabel("x (m)")
plt.title("Potential with respect to position")
plt.ylim(-50,500)
plt.savefig("plots/potential.pdf")
	
# Plot a graph of different energy level wavefunctions 
plt.figure(figsize=(9,4))
plt.plot(xpoints,secand(0*e,5*e)[0:N],label='A')
plt.plot(xpoints,secand(150*e,10*e)[0:N],label='B')# Higher initial guesses approch to higher quantum number
plt.plot(xpoints, secand(400*e, 30*e)[0:N],label='C')
plt.plot(xpoints, secand(650*e, 40*e)[0:N], label='D')
plt.xlabel("x (m)")
plt.ylabel(r"$\psi(x)$")
plt.title("Enrgy state wavefunctions for infinite square well")
plt.legend()
plt.savefig("plots/wave.pdf")

plt.show()
