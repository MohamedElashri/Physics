clear;clc;close all;

tic
global MagMean;
global EnMean;
global Mag;
global En;
global En2Mean;
global Mag2Mean;
global Check;
L = 32; % The dimension of the LxL lattice. 
IC= 3; % positive (1), negative(2) or random(3).
Check = 0;
spin = Intialization(L, IC);
n=10^8;

spin = Metropolis(spin, 1, En, Mag, L, n, MagMean, EnMean);
n=10^7;
i = 1;
B = 0;

for T=1:0.002:4
    Check = 1;
    spin = Metropolis(spin, T, En, Mag, L, n);
    Energy(i) = EnMean;
    Magnetization (i) = MagMean;
    MagSus (i) = ((L^2)/T)*(Mag2Mean - (MagMean^2));
    SpesHeat (i) = ((L/T)^2)*(En2Mean - (EnMean^2));
    i = i + 1;
end 

Temp = 1:0.002:4;
figure(3)

subplot(2,2,1);
plot(Temp,Energy,'.')
title('Mean Energy')
xlabel('Temperature')
ylabel('Energy')

subplot(2,2,2);
plot(Temp,Magnetization,'.')
title('Mean Magnetization')
xlabel('Temperature')
ylabel('Magnetization')

subplot(2,2,3);
plot(Temp, SpesHeat, '.')
title('Spesific Heat')
xlabel('Temperature')
ylabel('Spesific Heat')

subplot(2,2,4);
plot(Temp,MagSus, '.')
title('Magnetic Susceptibility')
xlabel('Temperature')
ylabel('Magnetic Susceptibility')

sgtitle({
    ['Simulation of 2D Ising Model by Metropolis Algorithm' ] 
    ['Lattice Dimension = ' L ' and External Magnetic Field(B) = ' B ] 
    ['Metropolis Step  = ' n ]
    });


toc