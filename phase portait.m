% program to plot  phase portait for 1-D potentials

% All values are normalized and arbitrary

clc;
close all;
clear all;

m = 1;                  % Mass of the object

x = -5:0.10:5;        % Range of x-values
Nx = length(x);

 U = cos(x)+0.5*x;         % Potential Energy

E = -10:0.5:5;          % Range of energies to consider
NE = length(E);         
                        % Kinetic energy or v(x) = dx/dt
Vx1 = zeros(NE,Nx);     % Positive solution
Vx2 = zeros(NE,Nx);     % Negative solution


figure(1)
  plot(x,U,'LineWidth',2); grid on;
  axis([min(x) max(x) min(U) max(U)]);
  xlabel('X ','FontSize',14);
  ylabel('U(x) = Potential Energy','FontSize',14); 
  
  title(' Potential ' ,'FontSize',14);                            
  h=gca; 
  set(h,'FontSize',14);
  fh = figure(1);
  set(fh, 'Color', 'white'); 

figure(2)
  fh = figure(2);
  set(fh, 'Color', 'white');
  
for i = 1:NE
    
    D = E(i)-U;

    
     Vx1(i,:) = sqrt((2/m)*D); 
     Vx2(i,:) = -sqrt((2/m)*D); 


  plot(x,(Vx1(i,:)),'LineWidth',2);
  grid on; hold on;
  plot(x,(Vx2(i,:)),'LineWidth',2);
  axis([min(x) max(x) -max(E) max(E)]);
  xlabel('X','FontSize',14);
  ylabel('Vx = dx/dt','FontSize',14); 
  
  title(' Phase Portrait' ,'FontSize',14);                            
  h=gca; 
  set(h,'FontSize',14);
 
end



    
