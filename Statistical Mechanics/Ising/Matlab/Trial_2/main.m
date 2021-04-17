clc;
clear all;

numSpins = 5; %dimension
J = 1;%magnetic atomic momenta
h=0; %external magnetic field
numIters=5*10000;

% Temperatures to sample
numTemps = 2^8; 
kTc = 2*J / log(1+sqrt(2)); % Curie temperature
kT = linspace(0, 10, numTemps);
probSpinUp = 1;

% Preallocate to store results
Emean = zeros(size(kT));
CvVar=zeros(size(kT));
Mmean = zeros(size(kT));
MgSuspVar=zeros(size(kT));


parfor tempIndex = 1 : numTemps
    spin = GridBuilding(numSpins, probSpinUp);
    spin= metropolisEquilibrium(numIters,spin, kT(tempIndex), J);
    [spin,Energies,Magnetizations] = metropolisSampling(numIters,spin, kT(tempIndex), J);
    Emean(tempIndex) = mean(Energies);
    CvVar(tempIndex) =cvIsing(Energies,kT(tempIndex));
    Mmean(tempIndex) = mean(Magnetizations);
    MgSuspVar(tempIndex) = suspIsing(Magnetizations,kT(tempIndex));
end


figure(1);

subplot(2,2,1);
plot(kT, Emean, '.');
title('Mean Energy Per Spin vs Temperature');
xlabel('kT'); 
ylabel('<E>');

subplot(2,2,2);
plot(kT, Mmean, '.');
title('Mean Magnetization Per Spin vs Temperature');
xlabel('kT');
ylabel('<M>');

subplot(2,2,3);
plot(kT, CvVar, '.');
title('Specific Heat Per Spin vs Temperature');
xlabel('kT');
ylabel('Cv');

subplot(2,2,4);
plot(kT, MgSuspVar, '.');
title('Magnetic Susceptibility Per Spin vs Temperature');
xlabel('kT');
ylabel('X');

