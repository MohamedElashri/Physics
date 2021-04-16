#----------------------------------------------------------------------
## 2d Ising model Monte-Carlo Simulation
## Author: Mohamed Elashri 
## Email: elashrmr@mail.uc.edu
##  Algorithm 
##  1- Prepare some initial configrations of N spins. 
##  2- Flip spin of a lattice site chosen randomly 
##  3- Calculate the change in energy due to that 
##  4- If this change is negative, accept such move. If change is positive, accept it with probability exp^{-dE/kT}
##  5- repeat 2-4. 
##  6- calculate Other parameters and plot them 
#----------------------------------------------------------------------

'''
Lattice is a periodical structure of points that align one by one. 2D lattice can be plotted as: 

* * * * * * * *   
* * * * * * * * 
* * * * * * * *
* * * * * * * *
* * * * * * * *

The points in lattice are called lattice points, neareast lattice points of point ^ are those lattice points denoted by (*) shown in the graph below:

* * *(*)* * * *
* *(*)^(*)* * *
* * *(*)* * * *
* * * * * * * *

Each lattice point is denoted by a number i in the Harmitonian.

The expression for the Energy of the total system is (online latex formula)
http://melashri.net/url/a or (H = - J \sum_{ i = 0 }^{ N-1 } \sum_{ j = 0 }^{ N-1 } (s_{i,j}s_{i,j+1}+s_{i,j}s_{i+1,j}) )

* * * * * * * * 
* * * * * * * *
* * * * * * * * <-the i-th lattice point
* * * * * * * *
* * * * * * * *

Periodical strcture means that lattice point at(1,1) is the same as that at(1,9) if the lattice is 5 by 8. more e.g.(1,1)<=>(6,1),
(2,3)<=>(2,11). A 2D lattice can be any Nx by Ny. The location (x,y) here is another denotion of lattice point that 
is fundementally same as i-th lattice point denotation above.

* * * * * * * * 4
* * * * * * * * 3
* * * * * * * * 2
* * * * * * * * 1
1 2 3 4 5 6 7 8 



'''

#----------------------------------------------------------------------
##  Import needed python libraries 
#----------------------------------------------------------------------
# %matplotlib inline # if I'm working on jupyter notebook, comment if working from terminal.
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt



#----------------------------------------------------------------------
##  Define the main functions
#----------------------------------------------------------------------
def initialstate(N):   
    ''' generates a random spin configuration for a specific boundry condition'''
    state = 2*np.random.randint(2, size=(N,N))-1
    return state


# Calculate interaction energy between spins. Assume periodic boundaries
# Interaction energy will be the difference in energy due to flipping spin i,j 
# (Example: 2*spin_value*neighboring_spins)

def mcmove(config, beta):
    '''Monte Carlo implementation function using Metropolis algorithm '''
    for i in range(N):
        for j in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                s =  config[a, b] # spin value
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N] ## neighboring spins
                cost = 2*s*nb # interaction energy (difference in energy in Dr.Yashar's matlab implementation code)
                if cost < 0:
                    s *= -1
                elif rand() < np.exp(-cost*beta):
                    s *= -1
                config[a, b] = s
    return config


def calcEnergy(config):
    '''Calculate the Energy of a given configuration'''
    ''' The energy is the sum of interactions between spins divided by the total number of spins '''
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            S = config[i,j]
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S
    return energy/4.


def calcMag(config):
    '''Calculate the Magnetization of a given configuration'''
    ''' Magnetization is the sum of all spins divided by the total number of spins '''
    mag = np.sum(config)
    return mag




## Fine tuning paramters 
## change these parameters for a smaller (usually means faster) simulations
nt      = 32         #  number of temperature points
N       = 50         #  size of the lattice, N x N > we are working on L X L periodic lattice 
eqSteps = 1024       #  number of MC sweeps for equilibration
mcSteps = 1024       #  number of MC sweeps for calculation
B = 0                #  Magnetic field strength 
# The two numbers should be the same for our purposes 

T       = np.linspace(1.53, 3.28, nt);  # Calculate the temperature. 
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N) 
# divide by number of samples, and by system size to get intensive values.


# Main MC function definition
# Calculate interaction energy between spins. Assume periodic boundaries
# Interaction energy will be the difference in energy due to flipping spin i,j 

for tt in range(nt): # loop on Lattice. 
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N)
    iT=1.0/T[tt]; iT2=iT*iT;
    
    for i in range(eqSteps):         # equilibrate
        mcmove(config, iT)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(config, iT)           
        Ene = calcEnergy(config)     # calculate the energy
        Mag = calcMag(config)        # calculate the magnetisation

        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag 
        E2 = E2 + Ene*Ene

    E[tt] = n1*E1 # Energy equation
    M[tt] = n1*M1 # Magnetization equation
    C[tt] = (n1*E2 - n2*E1*E1)*iT2 # Specific heat equation
    X[tt] = (n1*M2 - n2*M1*M1)*iT # Susceptibility equation 


#----------------------------------------------------------------------
#  Plotting area
#----------------------------------------------------------------------

f = plt.figure(figsize=(18, 10)); # plot the calculated values  


sp =  f.add_subplot(2, 2, 1 );
plt.plot(T, E, marker='o', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');

sp =  f.add_subplot(2, 2, 2 );
plt.plot(T, abs(M), marker='o', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');

sp =  f.add_subplot(2, 2, 3 );
plt.plot(T, C, marker='o', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20);  
plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');   

sp =  f.add_subplot(2, 2, 4 );
plt.plot(T, X, marker='o', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');

plt.subplots_adjust(0.12, 0.11, 0.90, 0.81, 0.26, 0.56)
plt.suptitle("Simulation of 2D Ising Model by Metropolis Algorithm\n" + "Lattice Dimension:" + str(N) + "X" + str(
    N) + "\n" + "External Magnetic Field(B)=" + str(B) + "\n" + "Metropolis Step=" + str(mcSteps))


plt.show() # function to show the plots

  
